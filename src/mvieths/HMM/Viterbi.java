/**
 * 
 */
package mvieths.HMM;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Stack;

/**
 * @author MVieths viet0013@umn.edu
 * 
 */
public class Viterbi {

    // Starting values for p and q
    // final static double p = 0.98;
    // final static double q = 0.999;
    // Values I'm playing with
    final static double                p                  = 0.99999999;
    final static double                q                  = 0.94;

    static char[]                      observations       =
                                                          { 'A', 'C', 'G', 'T', 'N' };

    /*
     * A + indicates a probability within a CpG island, - indicates out of CpG island
     * 
     * e.g. A+ to A- represents a transition from an A inside a CpG island to an A outside of a CpG island
     */
    static String[]                    states             =
                                                          { "A+", "C+", "G+", "T+", "A-", "C-", "G-", "T-" };

    static double[][]                  transitionMatrix   =
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
    static double[][]                  emissionMatrix     =
                                                          {
                                                          // A+ C+ G+ T+ A- C- G- T-
            { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 }, // A
            { 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 }, // C
            { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, // G
            { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 }, // T
            { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 }   // N
                                                          };

    // Reverse mapping from string to integer for looking up a state
    static HashMap<String, Integer>    stateMap;
    // Reverse mapping from character to integer for looking up an observation
    static HashMap<Character, Integer> obsMap;

    static double[]                    initialProbability =
                                                          { 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0,
                                                          1.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0 };

    /**
     * @param args
     */
    public static void main(String[] args) {
        // Create a reverse-mapping of the string state to the corresponding integer
        stateMap = new HashMap<String, Integer>();
        for (int i = 0; i < states.length; i++) {
            stateMap.put(states[i], i);
        }

        // Create a reverse-mapping of the character observation to the corresponding integer
        obsMap = new HashMap<Character, Integer>();
        for (int i = 0; i < observations.length; i++) {
            obsMap.put(observations[i], i);
        }

        String chr = "Chr21.txt";
        String filename = args[0];

        String sequence = parseFASTA(filename);
        System.out.println("sequence is " + sequence.length() + " bases long");
        ArrayList<CpGIsland> islands = viterbinate(sequence);

        File inFile = new File(filename);
        File parentDir = inFile.getParentFile();
        File chrFile = new File(parentDir, chr);
        HashMap<Integer, String> forwardAnnotations = readAnnotations(chrFile);

        annotateIslands(islands, forwardAnnotations);
    }

    public static ArrayList<CpGIsland> viterbinate(String sequence) {
        double[][] probabilityTable = new double[states.length][sequence.length()];
        int[][] stateTable = new int[states.length][sequence.length()];

        // Set the initial probabilities and states
        char firstChar = sequence.charAt(0);
        for (int i = 0; i < states.length; i++) {
            // Probability is 0 if emission probability is 0 (to avoid Math.log issues)
            if (emissionMatrix[obsMap.get(firstChar)][i] != 0.0) {
                probabilityTable[i][0] = Math.log(initialProbability[i] * emissionMatrix[obsMap.get(firstChar)][i]);
            }
            else {
                probabilityTable[i][0] = 0;
            }
            stateTable[i][0] = 0;
        }

        // Continue populating the tables at sequence.charAt(1)
        for (int seqPos = 1; seqPos < sequence.length(); seqPos++) {
            char thisChar = sequence.charAt(seqPos);

            // Iterate through the states
            for (int state = 0; state < states.length; state++) {
                // Find the maximum value in the previous entry in the table
                // Since we're using Math.log() to avoid underruns and will likely have negative probabilities, use negative infinity to start
                double maxProbability = Double.NEGATIVE_INFINITY;
                int maxState = -1;
                for (int k = 0; k < states.length; k++) {
                    double prob = Double.NEGATIVE_INFINITY;
                    // Math.log(0) comes out as negative infinity, so only get the log and sum the probability if it's non-zero
                    if (thisChar == 'N') {
                        // 'N' stands for any, so provide an equal probability of any state
                        prob = probabilityTable[k][seqPos - 1] + Math.log(1.0 / 8.0)
                                + Math.log(emissionMatrix[obsMap.get(thisChar)][k]);
                    }
                    else if (emissionMatrix[obsMap.get(thisChar)][k] != 0.0) {
                        prob = probabilityTable[k][seqPos - 1] + Math.log(transitionMatrix[k][state])
                                + Math.log(emissionMatrix[obsMap.get(thisChar)][k]);
                    }

                    if (prob > maxProbability) {
                        maxProbability = prob;
                        maxState = k;
                    }
                }

                probabilityTable[state][seqPos] = maxProbability;
                stateTable[state][seqPos] = maxState;
            }
        }

        // Now that we've populated the tables, backtrace through them
        // Find the greatest probability and its state in the last entry
        Stack<Integer> backtrace = new Stack<Integer>();

        // Step backward through the table, picking out the maximum state at each position
        for (int i = sequence.length() - 1; i >= 0; i--) {
            double maxProb = Double.NEGATIVE_INFINITY;
            int maxState = -1;
            for (int j = 0; j < states.length; j++) {
                double prob = probabilityTable[j][i];
                if (prob > maxProb) {
                    maxProb = prob;
                    maxState = j;
                }
            }

            // Now that we have the max state position, push the corresponding state from the other table onto our stack
            backtrace.push(stateTable[maxState][i]);
        }

        ArrayList<CpGIsland> islands = new ArrayList<CpGIsland>();
        CpGIsland island = new CpGIsland();

        // Initialize whether we're in an island based on the state at the top of the stack
        // States < 4 are in a CpG island, >= 4 are not
        int state = backtrace.pop();
        boolean inIsland = (state < 4);
        if (inIsland) {
            island.setStart(0);
        }
        int i = 1;
        // work through the stack populating CpGIslands as we find them
        while (!backtrace.isEmpty()) {
            state = backtrace.pop();

            // Make note any time we change state (in/out of CpG island)
            if (inIsland && state > 3) {
                inIsland = false;
                island.setEnd(i);
                islands.add(island);
                island = new CpGIsland();
            }
            else if (!inIsland && state < 4) {
                inIsland = true;
                island.setStart(i);
            }
            i++;
        }

        return islands;
    }

    /**
     * Line up the islands we found with our annotation file to determine if any of them fall within 500 bp of a gene.
     * Ignore any <200bp in length.
     * 
     * @param islands
     * @param forwardAnnotations
     */
    private static void annotateIslands(ArrayList<CpGIsland> islands, HashMap<Integer, String> forwardAnnotations) {
        int window = 500;
        int i = 0;
        for (CpGIsland next : islands) {
            if (next.getLength() > 200) {
                int forwardWindowEnd = next.getEnd() + window;
                for (Integer key : forwardAnnotations.keySet()) {
                    if ((key > next.getStart()) && (key < forwardWindowEnd)) {
                        // This gene is within 500bp of the CpG window (start or end)
                        System.out.printf("CpG island %10d bp (%d-%d) %s\n", next.getLength(), next.getStart(),
                                next.getEnd(), forwardAnnotations.get(key));
                        i++;
                    }
                }
            }
        }
        System.out.println("Total CpG islands found: " + islands.size() + " ; " + i + " out of " + islands.size()
                + " are followed by a coding region");
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
                    sequence += line.replaceAll("\\s+", "");
                }
            }
        }
        catch (IOException e) {
            System.out.println("Failed to read file " + filename);
            e.printStackTrace();
        }

        return sequence;
    }

    /**
     * Parse the annotation file and return any from the forward strand.
     * 
     * @param chrFile
     * @return
     */
    private static HashMap<Integer, String> readAnnotations(File chrFile) {
        HashMap<Integer, String> annotations = new HashMap<Integer, String>();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(chrFile));

            // Read the annotations from the file
            String line;
            while ((line = reader.readLine()) != null) {
                String[] columns = line.split("\\s+");
                // Only look at forward strand
                if (columns[2].equals("+")) {
                    int startPosition = Integer.parseInt(columns[0]);
                    annotations.put(startPosition, columns[3]);
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return annotations;
    }
}

/**
 * Store the start and end of a CpG island
 * @author MVieths
 *
 */
class CpGIsland {
    // Since this particular dataset starts at 43507093, include that as an offset for our coordinates
    private final int offset = 43507093;
    private int       start;
    private int       end;

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start
     *            the start to set
     */
    public void setStart(int start) {
        this.start = start + offset;
    }

    /**
     * @return the end
     */
    public int getEnd() {
        return end;
    }

    /**
     * @param end
     *            the end to set
     */
    public void setEnd(int end) {
        this.end = end + offset;
    }

    /**
     * Return the length in bases of this island
     * @return
     */
    public int getLength() {
        return end - start;
    }
}