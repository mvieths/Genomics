/**
 * 
 */
package mvieths.HMM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * @author Foeclan
 * 
 */
public class Viterbi {

	double p = 0.98;
	double q = 0.999;
	/*
	 * A + indicates a probability within a CpG island, - indicates out of CpG island
	 * 
	 * e.g. A+ to A- represents a transition from an A inside a CpG island to an A outside of a CpG island
	 */
	double[][] stateMatrix = {
			// A+ C+ G+ T+ A- C- G- T-
	{ 0.180 * p, 0.274 * p, 0.426 * p, 0.120 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4 }, // A+
	{ 0.171 * p, 0.368 * p, 0.274 * p, 0.188 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4 }, // C+
	{ 0.161 * p, 0.339 * p, 0.375 * p, 0.125 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4 }, // G+
	{ 0.079 * p, 0.355 * p, 0.384 * p, 0.182 * p, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4, (1 - p) / 4 }, // T+
	{ (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.300 * q, 0.205 * q, 0.285 * q, 0.210 * q }, // A-
	{ (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.322 * q, 0.298 * q, 0.078 * q, 0.302 * q }, // C-
	{ (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.248 * q, 0.246 * q, 0.298 * q, 0.208 * q }, // G-
	{ (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, (1 - q) / 4, 0.177 * q, 0.239 * q, 0.292 * q, 0.292 * q } // T-
	};

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String dataDirectory = "C:\\Users\\Foeclan\\Dropbox\\Genomics\\Homework 2\\";
		String toySample1 = "seq1.out";
		String toySample2 = "seq2.out";
		String chr21 = "Chr21.txt";
		String contig = "HMC21_NT_011515";
		double initialProbability = 1 / 8; // Equal chance of any state

		String toySequence1 = parseFASTA(dataDirectory + toySample1);
	}

	/**
	 * Parse a FASTA-formatted file. Assumption is that comments start with > and all other lines are part of the sequence
	 * 
	 * @param filename
	 * @return Sequence read from the file
	 */
	private static String parseFASTA(String filename) {
		String sequence = "";

		try {
			BufferedReader reader = new BufferedReader(new FileReader(filename));

			// Read the sequence from the file
			String line;
			while ((line = reader.readLine()) != null) {
				// Ignore comments
				if (!line.startsWith(">")) {
					line.replaceAll("\n", "");
					line.replaceAll("\r", "");
					sequence += line;
				}
			}
		} catch (IOException e) {
			System.out.println("Failed to read file " + filename);
			e.printStackTrace();
		}
		return sequence;
	}

}
