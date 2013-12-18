package mvieths.HMM.Util;

import java.util.HashMap;

/**
 * Stores a dataset, which consists of an amino acid sequence and its corresponding states
 * 
 * @author Michael Vieths
 * 
 */
public class HMMDataSet
{
    static String                      sequenceData;
    static String                      stateData;
    static String[]                    states;
    static char[]                      observations;
    static double[][]                  transitionMatrix;
    static double[][]                  emissionMatrix;
    static double[]                    initialProbability;
    // Reverse mapping from string to integer for looking up a state
    static HashMap<String, Integer>    stateMap;
    // Reverse mapping from character to integer for looking up an observation
    static HashMap<Character, Integer> obsMap;

    public HMMDataSet(String sequence, String stateData, char[] observations, String[] states,
            double[] initialProbability, double[][] transitionMatrix, double[][] emissionMatrix) throws Exception
    {
        if (sequence.length() != stateData.length())
        {
            throw new Exception("Sequence and state lengths do not match");
        }
        else
        {
            setSequenceData(sequence);
            setStateData(stateData);
            setObservations(observations);
            setStates(states);
            setInitialProbability(initialProbability);
            setTransitionMatrix(transitionMatrix);
            setEmissionMatrix(emissionMatrix);
        }
    }

    public static double[][] forward() {
        double[][] probabilityTable = new double[states.length][sequenceData.length()];
        int[][] stateTable = new int[states.length][sequenceData.length()];

        // Set the initial probabilities and states
        char firstChar = sequenceData.charAt(0);
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

        // Continue populating the tables at sequenceData.charAt(1)
        for (int seqPos = 1; seqPos < sequenceData.length(); seqPos++) {
            char thisChar = sequenceData.charAt(seqPos);

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
        return probabilityTable;
    }

    public static double[][] backward() {
        double[][] resultingProbabilities = new double[states.length][sequenceData.length()];

        return resultingProbabilities;
    }

    /**
     * @return the states
     */
    public static String[] getStates() {
        return states;
    }

    /**
     * @param states the states to set
     */
    public static void setStates(String[] states) {
        HMMDataSet.states = states;
        setStateMap();
    }

    /**
     * @return the observations
     */
    public static char[] getObservations() {
        return observations;
    }

    /**
     * @param observations the observations to set
     */
    public static void setObservations(char[] observations) {
        HMMDataSet.observations = observations;
        setObsMap();
    }

    /**
     * @return the transitionMatrix
     */
    public static double[][] getTransitionMatrix() {
        return transitionMatrix;
    }

    /**
     * @param transitionMatrix the transitionMatrix to set
     */
    public static void setTransitionMatrix(double[][] transitionMatrix) {
        HMMDataSet.transitionMatrix = transitionMatrix;
    }

    /**
     * @return the emissionMatrix
     */
    public static double[][] getEmissionMatrix() {
        return emissionMatrix;
    }

    /**
     * @param emissionMatrix the emissionMatrix to set
     */
    public static void setEmissionMatrix(double[][] emissionMatrix) {
        HMMDataSet.emissionMatrix = emissionMatrix;
    }

    /**
     * @return the initialProbability
     */
    public static double[] getInitialProbability() {
        return initialProbability;
    }

    /**
     * @param initialProbability the initialProbability to set
     */
    public static void setInitialProbability(double[] initialProbability) {
        HMMDataSet.initialProbability = initialProbability;
    }

    /**
     * @return the stateMap
     */
    public static HashMap<String, Integer> getStateMap() {
        return stateMap;
    }

    /**
     * @param stateMap the stateMap to set
     */
    public static void setStateMap() {
        // Create a reverse-mapping of the string state to the corresponding integer
        stateMap = new HashMap<String, Integer>();
        for (int i = 0; i < states.length; i++) {
            stateMap.put(states[i], i);
        }
    }

    /**
     * @return the obsMap
     */
    public static HashMap<Character, Integer> getObsMap() {
        return obsMap;
    }

    /**
     * @param obsMap the obsMap to set
     */
    public static void setObsMap() {
        // Create a reverse-mapping of the character observation to the corresponding integer
        obsMap = new HashMap<Character, Integer>();
        for (int i = 0; i < observations.length; i++) {
            obsMap.put(observations[i], i);
        }
    }

    public String getSequenceData()
    {
        return sequenceData;
    }

    public void setSequenceData(String sequence)
    {
        sequenceData = sequence;
    }

    public String getStateData()
    {
        return stateData;
    }

    public void setStateData(String states)
    {
        stateData = states;
    }
}