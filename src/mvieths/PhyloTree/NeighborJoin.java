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
	static PhyloNode rootNode;

	static String[] files = { "acinetobacter_sp.fasta", "bartonella_henselae.fasta", "bartonella_quintana.fasta", "bdellovibrio_bacteriovorus.fasta", "bordetella_parapertussis.fasta", "bordetella_pertussis.fasta", "bradyrhizobium_japonicum.fasta", "brucella_melitensis.fasta", "buchnera_aphidicola.fasta", "campylobacter_jejuni.fasta", "caulobacter_crescentus_cb15.fasta", "chromobacterium_violaceum.fasta", "coxiella_burnetii.fasta", "desulfovibrio_vulgaris.fasta", "escherichia_coli.fasta", "geobacter_sulfurreducens_pca.fasta", "haemophilus_influenzae.fasta", "helicobacter_hepaticus.fasta", "helicobacter_pylori.fasta", "neisseria_meningitidis.fasta", "pasteurella_multocida.fasta", "photobacterium_profundum.fasta", "photorhabdus_luminescens.fasta", "pseudomonas_aeruginosa.fasta", "pseudomonas_putida.fasta", "pseudomonas_syringae.fasta", "ralstonia_solanacearum.fasta", "rickettsia_conorii.fasta", "rickettsia_prowazekii.fasta", "salmonella_typhimurium.fasta", "shewanella_oneidensis_mr-1.fasta", "shigella_flexneri.fasta", "vibrio_cholerae.fasta", "wigglesworthia_glossinidia.fasta", "wolbachia_pipientis.fasta", "xanthomonas_axonopodis_pv._citri.fasta", "xylella_fastidiosa_9a5c.fasta", "yersinia_pestis.fasta" };

	/**
	 * @param args
	 */
	public static void main(String[] args) {
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
			initMatrices(inputFileName1, inputFileName2);
			rootNode = initTree(matrix1);

			System.out.println();
		}

	}

	public static void initMatrices(String inputFileName1, String inputFileName2) {
		try {
			File inFile1 = new File(inputFileName1);
			File inFile2 = new File(inputFileName2);
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

			// System.out.println("Matrix 1");
			// printMatrix(matrix1);

			// System.out.println("\nMatrix 2");
			// printMatrix(matrix2);
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	public static PhyloNode initTree(ArrayList<ArrayList<Double>> matrix) {
		PhyloNode newRoot = new PhyloNode();
		for (int i = 0; i < matrix.size(); i++) {
			PhyloNode node = new PhyloNode(i, files[i]);
			node.setDistances(matrix.get(i));
			newRoot.add(node);
		}
		return newRoot;
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

class PhyloNode {
	ArrayList<PhyloNode> children;
	int nodeNumber;
	int nextNodeNumber = 0;
	String nodeName;
	ArrayList<Double> distances;
	final static double BIG_DISTANCE = 10000000.0;

	public PhyloNode() {
		children = new ArrayList<PhyloNode>();
		nodeNumber = -1;
		nodeName = "";
	}

	public PhyloNode(int number, String name) {
		children = new ArrayList<PhyloNode>();
		nodeNumber = number;
		nodeName = name;
	}

	public void setNumber(int number) {
		nodeNumber = number;
	}

	public void setName(String name) {
		nodeName = name;
	}

	public void add(PhyloNode node) {
		children.add(node);
	}

	public ArrayList<PhyloNode> getNodes() {
		return children;
	}

	public void setDistances(ArrayList<Double> distances) {
		ArrayList<Double> d = new ArrayList<Double>();
		for (int i = 0; i < distances.size(); i++) {
			d.add(distances.get(i));
		}
		this.distances = d;
	}

	public ArrayList<Double> getDistances() {
		return distances;
	}

	public Double getDistance(int nodeNumber) {
		if (distances.size() < nodeNumber) {
			return 10000.0;
		} else {
			return distances.get(nodeNumber);
		}
	}

	public double getS() {
		double total = 0;
		for (Double distance : distances) {
			total += distance;
		}
		total /= (distances.size() - 2);
		return total;
	}

	// Take the child nodes with the two shortest distances and create a subnode containing them
	// If there are not two short nodes left, it won't create a subnode
	public void mergeShortest() {
		// Examine all the children
		PhyloNode short1 = null;
		PhyloNode short2 = null;
		double shortest = BIG_DISTANCE;
		double nextShortest = BIG_DISTANCE;
		for (PhyloNode child : children) {
			double nodeDistance = child.getS();
			if (nodeDistance <= shortest) {
				nextShortest = shortest;
				shortest = nodeDistance;
				short2 = short1;
				short1 = child;
			}
		}

		// If nextShortest is unchanged, we didn't find 2 short distances, so don't create a subnode
		if (nextShortest != BIG_DISTANCE) {
			// Create a new node containing the two shortest
			PhyloNode newNode = new PhyloNode();
			newNode.add(short1);
			newNode.add(short2);
			children.remove(short1);
			children.remove(short2);

			newNode.setName("U" + nextNodeNumber++);
		}
	}
}
