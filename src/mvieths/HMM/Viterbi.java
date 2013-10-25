/**
 * 
 */
package mvieths.HMM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 * @author Foeclan
 * 
 */
public class Viterbi {

    static double     p                  = 0.98;
    static double     q                  = 0.999;
    /*
     * A + indicates a probability within a CpG island, - indicates out of CpG island
     * 
     * e.g. A+ to A- represents a transition from an A inside a CpG island to an A outside of a CpG island
     */
    String[]          states             =
                                         { "A+", "C+", "G+", "T+", "A-", "C-", "G-", "T-" };
    String[]          observations       =
                                         { "-", "+" };                                              // In or out of a CpG island
    double[]          initialProbability =
                                         { 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8, 1 / 8 };

    static double[][] transitionMatrix   =
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
        findIslands(toySequence1);
    }

    /**
     * Parse through a provided String looking for CpG islands 
     * Note the start and end positions of the island 
     * Note any genes which occur within 500 bp of the end position of the island
     */
    public static void findIslands(String sequence) {
        Probability previous;
        boolean inIsland = false;
        char previousBase;
        char nextBase;
        double previousProbability = 1.0 / 8.0;
        System.out.println("Probability is " + previousProbability);

        HashMap<Integer, Boolean> islandMap = new HashMap<Integer, Boolean>();

        for (int i = 0; i < sequence.length() - 1; i++) {
            previousBase = sequence.charAt(i);
            nextBase = sequence.charAt(i + 1);
            previous = getMaxProbability(previousBase, nextBase, inIsland);
            if (previous.isInIsland() != inIsland) {
                // We entered or left an island
                inIsland = previous.isInIsland();
                islandMap.put(i, inIsland);
                if (inIsland) {
                    System.out.println("Entered an island at " + i);
                }
                else {
                    System.out.println("Left an island at " + i);
                }
            }
            //            System.out.println("prevProbability is " + previousProbability);
            //            System.out.println("Probability is " + previous.getProbability());

            previousProbability *= previous.getProbability();
        }

        System.out.println("Probability is " + previousProbability);
    }

    public static Probability getMaxProbability(char previous, char next, boolean inIsland) {
        Probability probability;
        double mostLikely = 0.0;
        boolean toIsland = false;

        /*
         * Need 4 probabilities:
         * 
         * Previous base is in island, next base is in island
         * 
         * Previous base is in island, next base is not
         * 
         * Previous base is not in an island, next base is
         * 
         * Previous base is not in an island, nor is the next
         */
        double curProbability;
        curProbability = getProbability(previous, inIsland, next, true);
        if (curProbability > mostLikely) {
            mostLikely = curProbability;
            toIsland = true;
        }

        curProbability = getProbability(previous, inIsland, next, false);
        if (curProbability > mostLikely) {
            mostLikely = curProbability;
            toIsland = false;
        }

        probability = new Probability(mostLikely, toIsland);

        return probability;
    }

    public static double getProbability(char previous, boolean prevInIsland, char next, boolean nextInIsland) {
        double probability = 0.0;
        int previousIndex = 0;
        int nextIndex = 0;
        int notInIslandMod = 4; // Add 4 to the index if it's not in an island

        switch (previous) {
            case 'A':
                previousIndex = 0;
                break;
            case 'C':
                previousIndex = 1;
                break;
            case 'G':
                previousIndex = 2;
                break;
            case 'T':
                previousIndex = 3;
                break;
        }
        if (!prevInIsland) {
            previousIndex += notInIslandMod;
        }

        switch (next) {
            case 'A':
                nextIndex = 0;
                break;
            case 'C':
                nextIndex = 1;
                break;
            case 'G':
                nextIndex = 2;
                break;
            case 'T':
                nextIndex = 3;
                break;
        }
        if (!nextInIsland) {
            nextIndex += notInIslandMod;
        }

        probability = transitionMatrix[previousIndex][nextIndex];
        //        System.out.println("Probability of going from " + previous + "(" + prevInIsland + ") to " + next + "("
        //                + nextInIsland + ") is " + probability);
        return probability;

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

class Probability {
    private double  probability;
    private boolean inIsland;

    public Probability(double prob, boolean inIsland) {
        setProbability(prob);
        setInIsland(inIsland);
    }

    /**
     * @return the probability
     */
    public double getProbability() {
        return probability;
    }

    /**
     * @param probability
     *            the probability to set
     */
    public void setProbability(double probability) {
        this.probability = probability;
    }

    /**
    * @return the inIsland
    */
    public boolean isInIsland() {
        return inIsland;
    }

    /**
     * @param inIsland
     *            the inIsland to set
     */
    public void setInIsland(boolean inIsland) {
        this.inIsland = inIsland;
    }
}