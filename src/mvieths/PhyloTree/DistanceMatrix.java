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
		try {
			if (args.length == 0) {
				printUsage();
			} else {
				File inFile = new File(args[0]);
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
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		FileWriter fileWriter;
		try {
			String outfile1 = parentDir.getPath() + "\\" + "d1.txt";
			fileWriter = new FileWriter(outfile1);
			BufferedWriter writer = new BufferedWriter(fileWriter);
			for (int i = 0; i < scoringMatrix.length; i++) {
				StringBuilder line = new StringBuilder();
				for (int j = 0; j < scoringMatrix[0].length; j++) {
					double p = 1 - scoringMatrix[i][j].getAlignmentPercent();
					line.append(p + " ");
					// System.out.printf("[%.5f]", p);
				}
				fileWriter.append(line + "\r\n");
				// System.out.println();
			}
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	private static void printUsage() {
		System.out.println("Usage:\nDistanceMatrix filelist\n\nfilelist\tA text file containing a list of fasta-formatted sequence files");
	}

}
