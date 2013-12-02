/**
 * 
 */
package mvieths.PhyloTree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;

import mvieths.HMM.AlignedSegment;
import mvieths.SeqAlign.SeqAlign;

/**
 * @author Foeclan
 * 
 */
public class DistanceMatrix {
	public static HashMap<String, String> sequences;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		BufferedReader reader = null;
		File parentDir = null;
		AlignedSegment[][] scoringMatrix = null;
		String inputFile = "";
		String outputFile1 = "";
		String outputFile2 = "";

		try {
		} catch (ArrayIndexOutOfBoundsException e) {
			printUsage();
		}

		try {
			if (args.length != 3) {
				printUsage();
				System.exit(1);
			} else {
				inputFile = args[0];
				outputFile1 = args[1];
				outputFile2 = args[2];

				File inFile = new File(inputFile);
				parentDir = inFile.getParentFile();
				reader = new BufferedReader(new FileReader(args[0]));
			}

			if (reader != null) {
				sequences = new HashMap<String, String>();
				String filename;
				while ((filename = reader.readLine()) != null) {
					String filepath = parentDir.getPath() + "\\" + filename;
					// System.out.println("filepath is " + filepath);
					sequences.put(filename, SeqAlign.parseFASTA(filepath));
				}

				int i = 0;
				int j = 0;
				scoringMatrix = new AlignedSegment[sequences.values().size()][sequences.values().size()];
				for (String sequence1 : sequences.values()) {
					for (String sequence2 : sequences.values()) {
						scoringMatrix[i][j] = SeqAlign.align(sequence1, sequence2);
						j++;
					}
					j = 0;
					i++;
				}

			} else {
				printUsage();
				System.exit(1);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		FileWriter fileWriter1, fileWriter2;
		try {
			String outfile1 = parentDir.getPath() + "\\" + outputFile1;
			String outfile2 = parentDir.getPath() + "\\" + outputFile2;
			fileWriter1 = new FileWriter(outfile1);
			fileWriter2 = new FileWriter(outfile2);
			BufferedWriter writer1 = new BufferedWriter(fileWriter1);
			BufferedWriter writer2 = new BufferedWriter(fileWriter2);
			for (int i = 0; i < scoringMatrix.length; i++) {
				StringBuilder d1Line = new StringBuilder();
				StringBuilder d2Line = new StringBuilder();

				for (int j = 0; j < scoringMatrix[0].length; j++) {
					double distance = 1 - scoringMatrix[i][j].getAlignmentPercent();
					d1Line.append(distance + " ");
					d2Line.append(jukesCantor(distance) + " ");
					// System.out.printf("[%.5f]", p);
				}
				fileWriter1.append(d1Line + "\r\n");
				fileWriter2.append(d2Line + "\r\n");
				// System.out.println();
			}
			writer1.close();
			writer2.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	/**
	 * Perform Jukes-Cantor correction on the provided distance
	 * 
	 * @param distance
	 * @return corrected distance
	 */
	public static double jukesCantor(double distance) {
		// -(3/4) * ln(1-(4/3)*distance)
		double correctedDistance = (-3.0 / 4.0) * Math.log(1 - ((4.0 / 3.0) * distance));

		return correctedDistance;
	}

	private static void printUsage() {
		System.out.println("Usage:\nDistanceMatrix filelist\n\nfilelist\tA text file containing a list of fasta-formatted sequence files");
	}

}
