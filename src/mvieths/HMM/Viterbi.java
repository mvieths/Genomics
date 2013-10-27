/**
 * 
 */
package mvieths.HMM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * @author Foeclan
 * 
 */
public class Viterbi {

    final static double p                  = 0.98;
    final static double q                  = 0.999;

    /*
     * A + indicates a probability within a CpG island, - indicates out of CpG island
     * 
     * e.g. A+ to A- represents a transition from an A inside a CpG island to an A outside of a CpG island
     */
    static String[]     states             =
                                           { "A+", "C+", "G+", "T+", "A-", "C-", "G-", "T-" };

    static double[]     initialProbability =
                                           { 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0,
                                           1.0 / 8.0,
                                           1.0 / 8.0 };

    static double[][]   transitionMatrix   =
                                           {
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

    // Emission probabilities are 1 when it's the same base, 0 otherwise
    static double[][]   emissionMatrix     =
                                           {
                                           { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
                                           { 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
                                           { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
                                           { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },
                                           { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 },
                                           { 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 },
                                           { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 },
                                           { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 }
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

        String toySequence1 = parseFASTA(dataDirectory + toySample1);
        System.out.println("sequence is " + toySequence1.length() + " long");
        forwardViterbi(toySequence1);
    }

    /**
     * 
     */
    // forwardViterbi(observations, states, start_probability, transition_probability, emission_probability)
    private static void forwardViterbi(String sequence) {
        // Keep track of the path and probability for each state
        // Initialize based on initialProbability matrix
        // This will serve the Ptr function in the algorithm
        ViterbiNode[] history = new ViterbiNode[states.length];
        for (int i = 0; i < history.length; i++) {
            history[i] = new ViterbiNode(initialProbability[i], i);
            history[i].setArgMax(initialProbability[i]);
        }

        // for each observation (nucleotide), calculate its probability and store it in the appropriate ViterbiNode
        for (int nucleotide = 0; nucleotide < sequence.length(); nucleotide++) {
            char observation = sequence.charAt(nucleotide);

            double max = 0.0;

            // Allocate a new array of nodes for the current probabilities
            ViterbiNode[] current = new ViterbiNode[states.length];
            // Generate probabilities for all of the current states
            for (int j = 0; j < history.length; j++) {
                // Iterate through the previous states
                for (int k = 0; k < history.length; k++) {
                    // Get the previous probability
                    double prevProb = history[k].getTotalProbability();
                    ArrayList<Integer> pastStates = history[k].getPath();
                    // Use log to avoid underruns
                    double probability = Math.log(transitionMatrix[k][j] * emissionMatrix[k][nucleotide]);

                }
            }
        }

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
        }
        catch (IOException e) {
            System.out.println("Failed to read file " + filename);
            e.printStackTrace();
        }
        return sequence;
    }
}

class ViterbiNode {
    double             totalProbability;
    double             argMax;
    ArrayList<Integer> path;

    public ViterbiNode(double probability, int state) {
        setTotalProbability(probability);
        path = new ArrayList<Integer>();
        path.add(state);
    }

    /**
     * @return the probability
     */
    public double getTotalProbability() {
        return totalProbability;
    }

    /**
     * @param probability the probability to set
     */
    public void setTotalProbability(double probability) {
        // Since the probabilities will get very small, store the log
        this.totalProbability = Math.log(probability);
    }

    /**
     * @return the argMax
     */
    public double getArgMax() {
        return argMax;
    }

    /**
     * @param argMax the argMax to set
     */
    public void setArgMax(double argMax) {
        this.argMax = argMax;
    }

    /**
     * @return the path
     */
    public ArrayList<Integer> getPath() {
        return path;
    }

    /**
     * @param path the path to set
     */
    public void setPath(ArrayList<Integer> path) {
        this.path = path;
    }

    public void addToPath(int state) {
        path.add(state);
    }
}