/**
 * 
 */
package mvieths.PhyloTree;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

/**
 * @author Foeclan
 * 
 */
public class NeighborJoin {
	static ArrayList<ArrayList<Double>> matrix1; // Uncorrected distance matrix
	static ArrayList<ArrayList<Double>> matrix2; // Jukes-Cantor-corrected distance matrix

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			String inputFileName1 = "";
			String inputFileName2 = "";
			String outputFileName1 = "D1-tree.txt";
			String outputFileName2 = "D2-tree.txt";

			if (args.length != 2) {
				printUsage();
				System.exit(1);
			} else {
				inputFileName1 = args[0];
				inputFileName2 = args[1];

				File inFile1 = new File(inputFileName1);
				File inFile2 = new File(inputFileName2);
				File parentDir = inFile1.getParentFile();
				BufferedReader reader1 = new BufferedReader(new FileReader(inFile1));
				BufferedReader reader2 = new BufferedReader(new FileReader(inFile2));

				matrix1 = new ArrayList<ArrayList<Double>>();
				matrix2 = new ArrayList<ArrayList<Double>>();

				String line;
				while ((line = reader1.readLine()) != null) {
					ArrayList<Double> doubleValues = new ArrayList<Double>();
					String[] values = line.split(" ");
					for (String value : values) {
						doubleValues.add(Double.parseDouble(value));
					}

					matrix1.add(doubleValues);
				}

				while ((line = reader2.readLine()) != null) {
					ArrayList<Double> doubleValues = new ArrayList<Double>();
					String[] values = line.split(" ");
					for (String value : values) {
						doubleValues.add(Double.parseDouble(value));
					}

					matrix2.add(doubleValues);
				}

				System.out.println("Matrix 1");
				// printMatrix(matrix1);

				System.out.println("\nMatrix2");
				// printMatrix(matrix2);

			}
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static void printUsage() {
		System.out.println("Usage:  NeighborJoin matrix1 matrix2");
	}

	public static void printMatrix(ArrayList<ArrayList<Double>> matrix) {
		for (ArrayList<Double> row : matrix) {
			for (Double value : row) {
				System.out.printf("[%.5f] ", value);
			}
			System.out.println();
		}
	}

	public static void buildTree(ArrayList<ArrayList<Double>> matrix) {
		for (ArrayList<Double> row : matrix) {
			for (Double value : row) {
				System.out.printf("[%.5f] ", value);
			}
			System.out.println();
		}
	}
}
