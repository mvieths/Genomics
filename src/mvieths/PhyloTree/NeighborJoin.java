/**
 * 
 */
package mvieths.PhyloTree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

/**
 * @author Michael Vieths
 * 
 */
public class NeighborJoin {
    static HashMap<String, HashMap<String, Double>> matrix1;                                                         // Uncorrected distance matrix
    static HashMap<String, HashMap<String, Double>> matrix2;                                                         // Jukes-Cantor-corrected distance matrix
    static PhyloRoot                                rootNode1;
    static PhyloRoot                                rootNode2;

    // This is the list of files for the DistanceMatrix portion of the assignment. This should probably be an additional argument to the program, but the assignment was explicit about only taking the two input files    
    static String[]                                 files =
                                                          { "acinetobacter_sp.fasta", "bartonella_henselae.fasta",
                                                          "bartonella_quintana.fasta",
                                                          "bdellovibrio_bacteriovorus.fasta",
                                                          "bordetella_parapertussis.fasta",
                                                          "bordetella_pertussis.fasta",
                                                          "bradyrhizobium_japonicum.fasta",
                                                          "brucella_melitensis.fasta",
                                                          "buchnera_aphidicola.fasta", "campylobacter_jejuni.fasta",
                                                          "caulobacter_crescentus_cb15.fasta",
                                                          "chromobacterium_violaceum.fasta", "coxiella_burnetii.fasta",
                                                          "desulfovibrio_vulgaris.fasta",
                                                          "escherichia_coli.fasta",
                                                          "geobacter_sulfurreducens_pca.fasta",
                                                          "haemophilus_influenzae.fasta",
                                                          "helicobacter_hepaticus.fasta", "helicobacter_pylori.fasta",
                                                          "neisseria_meningitidis.fasta",
                                                          "pasteurella_multocida.fasta",
                                                          "photobacterium_profundum.fasta",
                                                          "photorhabdus_luminescens.fasta",
                                                          "pseudomonas_aeruginosa.fasta", "pseudomonas_putida.fasta",
                                                          "pseudomonas_syringae.fasta",
                                                          "ralstonia_solanacearum.fasta", "rickettsia_conorii.fasta",
                                                          "rickettsia_prowazekii.fasta",
                                                          "salmonella_typhimurium.fasta",
                                                          "shewanella_oneidensis_mr-1.fasta",
                                                          "shigella_flexneri.fasta",
                                                          "vibrio_cholerae.fasta", "wigglesworthia_glossinidia.fasta",
                                                          "wolbachia_pipientis.fasta",
                                                          "xanthomonas_axonopodis_pv._citri.fasta",
                                                          "xylella_fastidiosa_9a5c.fasta", "yersinia_pestis.fasta" };

    /**
     * @param args
     */
    public static void main(String[] args) {
        String inputFileName1 = "";
        String outputFileName1 = "";

        if (args.length != 2) {
            printUsage();
            System.exit(1);
        }
        else {
            inputFileName1 = args[0];
            outputFileName1 = args[1];
            // Read in both matrices
            matrix1 = initMatrix(inputFileName1);

            // Create a root node for each matrix
            rootNode1 = new PhyloRoot(matrix1, files);
            // Build the tree for this matrix
            rootNode1.buildTree();
            // Get the Newick-format tree
            String newick1 = rootNode1.getNewickTree();

            // Write them out to the specified files
            try {
                FileWriter fw1 = new FileWriter(outputFileName1);

                BufferedWriter writer1 = new BufferedWriter(fw1);

                writer1.write(newick1);
                writer1.close();
            }
            catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

    }

    /**
     * Read the provided matrix file and store it in a hashmap for easier access
     * The HashMap format has the keys set to the labels of the current matrix, with the values being a further HashMap,
     *  which also has keys set to the labels of the current matrix, and values equal to what's read from the file (or, later,
     *  to the appropriately modified distance values for subnodes).  This allows easy addition and removal of nodes from the matrix
     *  without having to worry about creating new matrices when the indices change
     *  
     * @param inputFileName
     * @return
     */
    public static HashMap<String, HashMap<String, Double>> initMatrix(String inputFileName) {
        HashMap<String, HashMap<String, Double>> matrix = new HashMap<String, HashMap<String, Double>>();
        try {
            File inFile1 = new File(inputFileName);
            BufferedReader reader1 = new BufferedReader(new FileReader(inFile1));

            String line;
            int i = 0;
            while ((line = reader1.readLine()) != null) {
                String[] values = line.split(" ");
                HashMap<String, Double> row = new HashMap<String, Double>();
                for (int j = 0; j < files.length; j++) {
                    row.put(files[j], Double.parseDouble(values[j]));
                }

                matrix.put(files[i], row);
                i++;
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return matrix;
    }

    /**
     * Print a usage message if they don't provide the correct arguments
     */
    public static void printUsage() {
        System.out.println("Usage:  NeighborJoin matrix outputfile");
    }

    /**
     * Print out the matrix for debugging purposes
     * @param matrix
     */
    public static void printMatrix(HashMap<String, HashMap<String, Double>> matrix) {
        for (String file : files) {
            System.out.print(file + ":");
            for (String file2 : files) {
                System.out.printf(" [%.5f]", matrix.get(file).get(file2));
            }

            System.out.println();
        }
    }
}

/**
 * This will serve as the root node of the phylogenetic tree
 * 
 * @author Michael Vieths
 *
 */

class PhyloRoot {
    ArrayList<PhyloNode>                     children;                   // Child nodes of the root, will generally decrease over time
    HashMap<String, HashMap<String, Double>> matrix;                     // Initialized as the matrix read from the file, will change over time as nodes are added and removed
    final static double                      BIG_DISTANCE   = 10000000.0; // Pick a huge value for initialization
    int                                      nextNodeNumber = 0;         // Keep track of which Unified node we're on (U0..UN)
    ArrayList<String>                        labels;                     // Initialized as the list of files from DistanceMatrix, will change over time

    /**
     * Initialize the root node with the matrix and labels
     * Add initial child nodes
     * 
     * @param initialMatrix
     * @param files
     */
    public PhyloRoot(HashMap<String, HashMap<String, Double>> initialMatrix, String[] files) {
        children = new ArrayList<PhyloNode>();
        this.matrix = initialMatrix;
        this.labels = new ArrayList<String>(Arrays.asList(files));
        // Initialize the child nodes with their name and distances
        for (int i = 0; i < labels.size(); i++) {
            PhyloNode node = new PhyloNode(labels.get(i), initialMatrix.get(labels.get(i)));
            children.add(node);
        }
    }

    /**
     * Just keep merging until there aren't enough to merge anymore
     */
    public void buildTree() {
        while (matrix.keySet().size() > 2) {
            mergeShortest();
        }
    }

    /**
     * Print out the tree for debugging purposes
     */
    public void printTree() {
        System.out.println(getNewickTree());
    }

    /**
     * Returns the phylogenetic tree in Newick format
     */
    public String getNewickTree() {
        StringBuilder builder = new StringBuilder();
        builder.append("(");
        for (int i = 0; i < children.size(); i++) {
            builder.append(children.get(i).getNewick());
            if (i < children.size() - 1) {
                builder.append(",");
            }
        }
        builder.append(");"); // Newick format needs to end in a semicolon
        return builder.toString();
    }

    /**
     * Take the child nodes with the two shortest M-distances and create a subnode containing them
     * If there are not two short nodes left, it won't create a subnode
     * @return The merged node
     */
    public void mergeShortest() {
        double lowest = BIG_DISTANCE;
        PhyloNode lowNode1 = null;
        PhyloNode lowNode2 = null;

        // Iterate through all of the matrix pairs looking for that with the lowest M
        for (PhyloNode child : children) {
            for (PhyloNode child2 : children) {
                if (!child.equals(child2)) {
                    double m = calcM(child, child2);
                    if (m <= lowest) {
                        lowest = m;
                        lowNode1 = child;
                        lowNode2 = child2;
                    }
                }
            }
        }

        // Calculate the distances to a new node merged from these two
        double lowNodeJunction1 = calcS(lowNode1, lowNode2);

        // Name the new Union node with a unique number
        String newNodeName = "U" + nextNodeNumber++;

        // Recalculate all of the distances to this node
        HashMap<String, Double> distancesToU = new HashMap<String, Double>();
        for (String label : labels) {
            double distance;
            if (label.equals(lowNode1.getName())) {
                // If we're looking at ourself, the distances is always 0.0
                distance = 0.0;
            }
            else {
                // We know what the difference between U# and lowNode1 is, so get the distance from the target node to lowNode1 and subtract the difference
                distance = matrix.get(label).get(lowNode1.getName()) - lowNodeJunction1;
            }
            distancesToU.put(label, distance);
        }

        // Create the new Union node and perform general matrix and tree cleanup (add Union node, remove merged nodes)
        PhyloNode newNode = new PhyloNode(newNodeName, distancesToU);
        newNode.add(lowNode1);
        newNode.add(lowNode2);
        addToMatrix(newNode);
        removeFromMatrix(lowNode1);
        removeFromMatrix(lowNode2);
    }

    /**
     * Calculates the M distance, which is matrix distance minus S for both nodes
     */
    double calcM(PhyloNode i, PhyloNode j) {
        String iName = i.getName();
        String jName = j.getName();
        double distance = matrix.get(iName).get(jName);
        double s1 = i.calcS();
        double s2 = j.calcS();
        return distance - s1 - s2;
    }

    /**
     * Calculates the distance from node i and a Union node of i and j
     */
    double calcS(PhyloNode i, PhyloNode j) {
        double distance = matrix.get(i.getName()).get(j.getName());
        return ((distance / 2) + ((i.calcS() - j.calcS()) / 2));
    }

    /**
     * Remove the given node from all the children's distance matrices and ours, as well as cleaning up the labels and root's children
     */
    public void removeFromMatrix(PhyloNode node) {
        // First from ours
        matrix.remove(node.getName());
        for (HashMap<String, Double> value : matrix.values()) {
            value.remove(node.getName());
        }

        // Now from the children
        for (PhyloNode child : children) {
            child.getDistances().remove(node.getName());
        }
        children.remove(node);

        // Now from the labels
        labels.remove(labels.indexOf(node.getName()));
    }

    /**
     * Add the new node to our labels, children, distance matrix and each child node
     */
    public void addToMatrix(PhyloNode node) {
        // First add to ours
        matrix.put(node.getName(), node.getDistances());

        // Now to the children
        for (PhyloNode child : children) {
            child.getDistances().put(node.getName(), node.getDistances().get(child.getName()));
        }
        children.add(node);

        // Now to the labels
        labels.add(node.getName());
    }
}

/**
 * A node in a Phylogenetic tree
 *
 * @author Michael Vieths
 */
class PhyloNode {
    ArrayList<PhyloNode>    children; // Child nodes (2 at most)
    String                  nodeName; // Name of this node (U# if a Union node, name of a file otherwise)
    HashMap<String, Double> distances; // Map of the distances to other nodes from this one

    /**
     * Initialize the children, name, and distances
     */
    public PhyloNode(String name, HashMap<String, Double> distances) {
        children = new ArrayList<PhyloNode>();
        nodeName = name;
        setDistances(distances);
    }

    public void setName(String name) {
        nodeName = name;
    }

    public String getName() {
        return nodeName;
    }

    /**
     * Add the given node as a child to this node
     * @param node
     */
    public void add(PhyloNode node) {
        children.add(node);
    }

    /**
     * Return the child nodes of this node
     * @return
     */
    public ArrayList<PhyloNode> getNodes() {
        return children;
    }

    /**
     * Directly set the distances from this node to others
     * @param distances
     */
    public void setDistances(HashMap<String, Double> distances) {
        this.distances = distances;
    }

    /**
     * Get the distances from this node to others
     * @return
     */
    public HashMap<String, Double> getDistances() {
        return distances;
    }

    /**
     * Returns the sum of all distances for this node divided by the number of distances - 2
     * @return
     */
    public double calcS() {
        double total = 0;

        for (Double distance : distances.values()) {
            total += distance;
        }
        total /= (distances.values().size() - 2);

        return total;
    }

    /**
     * Recurse through the tree, returning Newick-format text for parsing in another tool
     * 
     * @return
     */
    public String getNewick() {
        StringBuilder builder = new StringBuilder();
        if (children.size() > 0) {
            builder.append("(");
            builder.append(children.get(0).getNewick());
            if (children.size() == 2) {
                builder.append(",");
                builder.append(children.get(1).getNewick());
            }
            builder.append(")");
        }
        else {
            builder.append(getName());
        }
        return builder.toString();
    }
}
