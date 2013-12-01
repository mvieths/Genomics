/**
 * 
 */
package mvieths.SeqAlign;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Stack;

import mvieths.HMM.AlignedSegment;

/**
 * @author Foeclan
 * 
 */
public class SeqAlign {

	static boolean debug = false;
	// static boolean debug = true;

	static boolean permute = false;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String humanSequenceFile = null;
		String flySequenceFile = null;
		String anchorFile = null;
		ArrayList<String[]> anchors = null;
		ArrayList<Integer> scores = null;

		if (args.length < 2 || args.length > 3) {
			System.out.println("Usage:  seqalign humanSequence.fa flySequence.fa [anchorFile]");
			System.exit(1);
		} else if (args.length == 2) {
			humanSequenceFile = args[0];
			flySequenceFile = args[1];
		} else if (args.length <= 3) {
			humanSequenceFile = args[0];
			flySequenceFile = args[1];
			anchorFile = args[2];
		} else {
			System.out.println("Usage:  seqalign humanSequence.fa flySequence.fa [anchorFile]");
			System.exit(1);
		}

		// Create readers for each type of file
		String humanSequence = "";
		String flySequence = "";

		if (humanSequenceFile != null && flySequenceFile != null) {
			humanSequence = parseFASTA(humanSequenceFile);
			flySequence = parseFASTA(flySequenceFile);
		}
		if (anchorFile != null) {
			if (anchorFile.contains("-permute")) {
				permute = true;
			} else {
				anchors = parseAnchors(anchorFile);
			}
		}

		if (!permute) {
			if (anchors == null) {
				AlignedSegment segment = align(humanSequence, flySequence);
				System.out.println("Global Score is " + segment.getScore());
				System.out.println("Human: " + segment.getFirstAlignment());
				System.out.println("Fly  : " + segment.getSecondAlignment());
			} else {
				// Create subsequences broken up by the anchored alignments, then align them
				ArrayList<String> humanSequences = new ArrayList<String>();
				ArrayList<String> flySequences = new ArrayList<String>();
				ArrayList<Integer> humanMarkers = new ArrayList<Integer>();
				ArrayList<Integer> flyMarkers = new ArrayList<Integer>();

				// Build up a list of anchor boundaries
				for (String[] match : anchors) {
					humanMarkers.add(new Integer(match[0]));
					humanMarkers.add(new Integer(match[1]));
					flyMarkers.add(new Integer(match[2]));
					flyMarkers.add(new Integer(match[3]));
				}

				// Build up a list of substrings based on the anchor boundaries
				int beginIndex = 0;
				for (Integer marker : humanMarkers) {
					String segment = humanSequence.substring(beginIndex, marker);
					beginIndex = marker;
					humanSequences.add(segment);
				}
				beginIndex = 0;
				for (Integer marker : flyMarkers) {
					String segment = flySequence.substring(beginIndex, marker);
					beginIndex = marker;
					flySequences.add(segment);
				}

				ArrayList<AlignedSegment> segments = new ArrayList<AlignedSegment>();

				// Align all of the subsequences separately
				for (int i = 0; i < humanSequences.size(); i++) {
					String humanSeq = humanSequences.get(i);
					String flySeq = flySequences.get(i);
					AlignedSegment segment;
					segment = align(humanSeq, flySeq);
					segments.add(segment);
				}

				// Add up all the segments and display them
				int globalScore = 0;
				String assembledHumanSequence = "";
				String assembledFlySequence = "";

				for (AlignedSegment segment : segments) {
					globalScore += segment.getScore();
					// System.out.println("Human: " + segment.getFirstAlignment());
					// System.out.println("Fly  : " + segment.getSecondAlignment());
					assembledHumanSequence += segment.getFirstAlignment();
					assembledFlySequence += segment.getSecondAlignment();
				}

				System.out.println("Global Score is " + globalScore);
				System.out.println("Human: " + assembledHumanSequence);
				System.out.println("Fly  : " + assembledFlySequence);
			}
		} else {
			// Randomly alter one of the sequences and align them
			// Store the scores
			scores = new ArrayList<Integer>();
			for (int i = 0; i < 10000; i++) {
				String newHumanSeq = mutateSequence(humanSequence);
				String newFlySeq = mutateSequence(flySequence);
				AlignedSegment segment = align(newHumanSeq, newFlySeq);
				scores.add(segment.getScore());
			}

			HashMap<Integer, Integer> distribution = new HashMap<Integer, Integer>();
			for (Integer score : scores) {
				if (!distribution.containsKey(score)) {
					distribution.put(score, 1);
				} else {
					distribution.put(score, distribution.get(score) + 1);
				}
			}

			// Write the scores out to a file
			FileWriter fileWriter;
			try {
				fileWriter = new FileWriter("d:\\Users\\Foeclan\\Dropbox\\Genomics\\Homework 1\\sequences\\permutations.txt");
				BufferedWriter writer = new BufferedWriter(fileWriter);
				for (Integer score : distribution.keySet()) {
					writer.write(score + "\t" + distribution.get(score) + "\n");
				}
				writer.close();
			} catch (IOException e) {
				e.printStackTrace();
			}

		}
	}

	/**
	 * Align two sequences and print out their gap score and alignment
	 * 
	 * @param sequence1
	 * @param sequence2
	 * @return AlignedSegment
	 */
	public static AlignedSegment align(String sequence1, String sequence2) {
		char[] seq1 = sequence1.toUpperCase().toCharArray();
		char[] seq2 = sequence2.toUpperCase().toCharArray();

		// Matrix for scoring the path
		int[][] seqMatrix = new int[seq1.length + 1][seq2.length + 1];
		// Matrix for storing the path
		char[][] ptrMatrix = new char[seq1.length + 1][seq2.length + 1];

		// Initialize the scoring matrix
		for (int[] array : seqMatrix) {
			array[0] = 0;
		}
		Arrays.fill(seqMatrix[0], 0);

		// Scoring from the assignment
		int gapPenalty = -2;
		int mismatchPenalty = -3;
		int matchBonus = 1;

		// Populate the matrices, starting at 1,1
		for (int i = 1; i < seqMatrix.length; i++) {
			for (int j = 1; j < seqMatrix[0].length; j++) {
				int maxScore = seqMatrix[i - 1][j - 1];
				char direction = 'd';
				if (seq1[i - 1] == seq2[j - 1]) {
					// If we match, apply the bonus
					maxScore += matchBonus;
				} else {
					// If we mismatch, apply the penalty
					maxScore += mismatchPenalty;
				}

				// Apply the gap penalty to the value to the left
				int left = seqMatrix[i - 1][j] + gapPenalty;
				// Apply the gap penalty to the value above
				int up = seqMatrix[i][j - 1] + gapPenalty;

				// If the value to the left is higher than the diagonal, that's our score
				if (left > maxScore) {
					maxScore = left;
					direction = 'l';
				}

				// If the value above is higher than the diagonal, that's our score
				if (up > maxScore) {
					maxScore = up;
					direction = 'u';
				}

				seqMatrix[i][j] = maxScore;
				ptrMatrix[i][j] = direction;
			}
		}

		/*
		 * Print the matrices and information about them for debugging purposes
		 */
		if (debug) {
			System.out.println("Aligning\n" + sequence1 + "\n" + sequence2);
			System.out.println(seqMatrix.length + " x " + seqMatrix[0].length);

			int window = 10;

			System.out.print("          ");
			printArray(seq2, window - 1);

			for (int i = 0; i < window; i++) {
				if (i > 0) {
					System.out.printf("[%3c]", seq1[i - 1]);
				} else {
					System.out.print("     ");
				}
				printArray(seqMatrix[i], window);
			}

			System.out.print("          ");
			printArray(seq2, window - 1);

			for (int i = 0; i < window; i++) {
				if (i > 0) {
					System.out.printf("[%3c]", seq1[i - 1]);
				} else {
					System.out.print("     ");
				}
				printArray(ptrMatrix[i], window);
			}
		}

		/*
		 * Now that we've scored the matrix, we need to find a path back through it.
		 */
		Stack<Character> seqOne = new Stack<Character>();
		Stack<Character> seqTwo = new Stack<Character>();

		boolean done = false;
		int i = ptrMatrix.length - 1;
		int j = ptrMatrix[0].length - 1;
		int score = 0;
		while (!done) {
			if (ptrMatrix[i][j] == 'u') {
				// 'Up' means a gap in sequence 1
				seqOne.push('-');
				seqTwo.push(seq2[j - 1]);
				j--;
			} else if (ptrMatrix[i][j] == 'l') {
				// 'Left' means a gap in sequence 2
				seqOne.push(seq1[i - 1]);
				seqTwo.push('-');
				i--;
			} else {
				// Diagonal means a match (or at least a non-decreasing path, since mismatches are possible just penalized).
				seqOne.push(seq1[i - 1]);
				seqTwo.push(seq2[j - 1]);
				i--;
				j--;
			}
			score += seqMatrix[i][j];

			// Stop if we hit an edge
			if (j == 0 || i == 0) {
				done = true;
			}
		}

		// Print out the alignment
		String seqOneString = "";
		String seqTwoString = "";

		while (!seqOne.isEmpty()) {
			seqOneString += "" + seqOne.pop();
		}
		while (!seqTwo.isEmpty()) {
			seqTwoString += "" + seqTwo.pop();
		}

		return new AlignedSegment(score, seqOneString, seqTwoString);
	}

	/**
	 * Parse a FASTA-formatted file. Assumption is that comments start with > and all other lines are part of the sequence
	 * 
	 * @param filename
	 * @return Sequence read from the file
	 */
	public static String parseFASTA(String filename) {
		String sequence = "";

		try {
			BufferedReader reader = new BufferedReader(new FileReader(filename));

			// Read the sequence from the file
			String line;
			while ((line = reader.readLine()) != null) {
				// Ignore comments
				if (!line.startsWith(">")) {
					sequence += line;
				}
			}
		} catch (IOException e) {
			System.out.println("Failed to read file " + filename);
			e.printStackTrace();
		}
		return sequence;
	}

	/**
	 * Parse a match file
	 * 
	 * Format is (tab delimited):
	 * 
	 * humanStart humanEnd flyStart flyEnd
	 * 
	 * Where the start/end represents the borders of an anchor (pre-aligned segment)
	 * 
	 * The
	 * 
	 * @param filename
	 */
	private static ArrayList<String[]> parseAnchors(String filename) {
		ArrayList<String[]> anchors = new ArrayList<String[]>();

		try {
			BufferedReader reader = new BufferedReader(new FileReader(filename));

			String line;
			while ((line = reader.readLine()) != null) {
				// First two columns are start/finish for human, second two are start/finish for fruit fly
				anchors.add(line.split("\t"));
			}

		} catch (IOException e) {
			System.out.println("Failed to read file " + filename);
			e.printStackTrace();
		}

		return anchors;
	}

	private static String mutateSequence(String sequence) {
		String aminoAcids = "ARNDCQEGHILKMFPSTWYVBZX";
		char[] aaArray = aminoAcids.toCharArray();
		char[] seqArray = sequence.toCharArray();

		int position = (int) (Math.random() * seqArray.length);
		int newAA = (int) (Math.random() * aaArray.length);

		while (aaArray[newAA] == seqArray[position]) {
			newAA = (int) (Math.random() * aaArray.length);
		}

		seqArray[position] = aaArray[newAA];

		return new String(seqArray);
	}

	/**
	 * Print the first 'max' elements of the array
	 * 
	 * @param array
	 * @param max
	 */
	private static void printArray(int[] array, int max) {
		if (array.length < max) {
			max = array.length;
		}
		for (int i = 0; i < max; i++) {
			System.out.printf("[%3d]", array[i]);
		}
		System.out.println("");
	}

	/**
	 * Print the first 'max' elements of the array
	 * 
	 * @param array
	 * @param max
	 */
	private static void printArray(char[] array, int max) {
		if (array.length < max) {
			max = array.length;
		}
		for (int i = 0; i < max; i++) {
			System.out.printf("[%3c]", array[i]);
		}
		System.out.println("");
	}
}
