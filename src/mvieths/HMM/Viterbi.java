/**
 * 
 */
package mvieths.HMM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * @author Foeclan
 * 
 */
public class Viterbi {

    //    final static double             p                  = 0.98;
    //    final static double             q                  = 0.999;
    final static double             p                = 0.9999;
    final static double             q                = 0.95;

    /*
     * A + indicates a probability within a CpG island, - indicates out of CpG island
     * 
     * e.g. A+ to A- represents a transition from an A inside a CpG island to an A outside of a CpG island
     */
    static String[]                 states           =
                                                     { "A+", "C+", "G+", "T+", "A-", "C-", "G-", "T-" };

    //    static double[]                 initialProbability =
    //                                                       { 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0,
    //                                                       1.0 / 8.0,
    //                                                       1.0 / 8.0,
    //                                                       1.0 / 8.0 };

    static double[][]               transitionMatrix =
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

    static HashMap<String, Integer> stateMap;

    /**
     * @param args
     */
    public static void main(String[] args) {
        // Create a reverse-mapping of the string state to the corresponding integer
        stateMap = new HashMap<String, Integer>();
        for (int i = 0; i < states.length; i++) {
            stateMap.put(states[i], i);
        }

        String dataDirectory = "C:\\Users\\Foeclan\\Dropbox\\Genomics\\Homework 2\\";
        String toySample1 = "seq1.out";
        String toySample2 = "seq2.out";
        String chr21 = "Chr21.txt";
        String contig = "HMC21_NT_011515.fasta";

        String toySequence1 = parseFASTA(dataDirectory + toySample1);
        //        System.out.println(toySequence1);
        System.out.println("sequence is " + toySequence1.length() + " long");
        viterbinate(toySequence1);
        //forwardViterbi(toySequence1);

        //        String chromosome21 = parseFASTA(dataDirectory + contig);
        //        forwardViterbi(chromosome21);
    }

    public static void viterbinate(String sequence) {
        char[] seqArray = sequence.toCharArray();

        // Set up the initial probability
        ViterbiNode[] initial = new ViterbiNode[states.length];
        for (int i = 0; i < states.length; i++) {
            initial[i] = new ViterbiNode();
            initial[i].total = Math.log(1.0);
        }

        ArrayList<ViterbiNode[]> history = new ArrayList<ViterbiNode[]>();
        history.add(initial);

        for (int position = 0; position < seqArray.length; position++) {
            char nucleotide = seqArray[position];
            for (int prevState = 0; prevState < states.length; prevState++) {
                double maxValue = 0.0;
                int maxState = -1;
                for (int curState = 0; curState < states.length; curState++) {
                    // Calculate the emission probability
                    double eProb = getEmissionProbability(curState, nucleotide);
                    double tProb = transitionMatrix[curState][prevState];
                    //System.out.printf("eProb [%10f] tProb [%10f]\n", eProb, tProb);
                    double probability = eProb * tProb;
                    if (probability > maxValue) {
                        maxValue = probability;
                        maxState = curState;
                    }
                }
                initial[prevState].total += Math.log(maxValue);
                initial[prevState].path.add(maxState);

            }

        }

    }

    public static double getEmissionProbability(int state, char nucleotide) {
        double eProb;
        if (states[state].contains("" + nucleotide)) {
            eProb = 1.0;
        }
        else {
            eProb = 0.0;
        }
        return eProb;
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
                    //                    line.replaceAll("\n", "");
                    //                    line.replaceAll("\r", "");
                    sequence += line;
                }
            }
            sequence = sequence.replaceAll("\\s+", "");
        }
        catch (IOException e) {
            System.out.println("Failed to read file " + filename);
            e.printStackTrace();
        }
        return sequence;
    }

}

class ViterbiNode {
    // Current probability of this path
    double total;
    // Current state
    int    state;

    public ViterbiNode(double total, int state) {
        setTotal(total);
        setState(state);
    }

    /**
     * @return the total
     */
    public double getTotal() {
        return total;
    }

    /**
     * @param total the total to set
     */
    public void setTotal(double total) {
        this.total = total;
    }

    /**
     * @return the state
     */
    public int getState() {
        return state;
    }

    /**
     * @param state the state to set
     */
    public void setState(int state) {
        this.state = state;
    }

}
