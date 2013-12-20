package mvieths.HMM.Util;

import java.util.HashMap;

/**
 * Stores a dataset, which consists of an amino acid sequence and its corresponding states.  
 * Provides methods for performing basic Hidden Markov Model algorithms on the dataset.
 * 
 * @author Michael Vieths
 * 
 */
public class HMMDataSet
{
    // The possible states for this Markov model
    static char[]                      states;

    // The possible observations for this Markov model
    static char[]                      observations;

    // The sequence we've been provided.  Must consist entirely of values from 'observations'
    static String                      sequenceData;

    // The state data we've been provided.  Must consist entirely of values from 'states'.  Must be of equal length to sequenceData.
    static String                      stateData;

    // The transition matrix contains the probability of changing from one state to another
    static double[][]                  transitionMatrix;

    // The emission matrix contains the probability of changing from one observation to another given a particular state
    static double[][]                  emissionMatrix;

    // The initial probability is the probability used for the first position in the sequence
    static double[]                    initialProbability;

    // Reverse mapping from string to integer for looking up a state
    static HashMap<Character, Integer> stateMap;

    // Reverse mapping from character to integer for looking up an observation
    static HashMap<Character, Integer> obsMap;

    public HMMDataSet(String sequence, String stateData, char[] observations, char[] states,
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

    /**
     * Performs the forward algorithm on the current dataset
     * 
     * @return resulting probability matrix
     */
    public static double[][] forward() {
        double[][] probabilityTable = new double[states.length][sequenceData.length()];

        // Set the initial probabilities starting at position 0 in the sequenceData
        for (int i = 0; i < states.length; i++) {
            probabilityTable[i][0] = initialProbability[i] * emissionMatrix[i][obsMap.get(sequenceData.charAt(0))];
        }

        // Continue populating the tables starting at sequenceData.charAt(1)
        for (int seqPos = 1; seqPos < sequenceData.length(); seqPos++) {
            // Calculate state values for each observation in the sequence
            for (int i = 0; i < states.length; i++) {
                probabilityTable[i][seqPos] = 0;
                // Add up the previous probabilities
                for (int j = 0; j < states.length; j++) {
                    probabilityTable[i][seqPos] += (probabilityTable[i][seqPos - 1] * transitionMatrix[i][j]);
                }
                // Multiply by the emission probability
                probabilityTable[i][seqPos] *= emissionMatrix[i][obsMap.get(sequenceData.charAt(seqPos))];
            }
        }
        return probabilityTable;
    }

    /**
     * Performs the backward algorithm on the current dataset
     * 
     * @return resulting probability matrix
     */
    public static double[][] backward() {
        double[][] probabilityTable = new double[states.length][sequenceData.length()];

        int lastPos = sequenceData.length() - 1;
        char lastChar = sequenceData.charAt(lastPos);
        // Set the initial probabilities starting at position 0 in the sequenceData
        for (int i = 0; i < states.length; i++) {
            probabilityTable[i][lastPos] = initialProbability[i] * emissionMatrix[i][obsMap.get(lastChar)];
        }

        // Continue populating the tables starting at sequenceData.charAt(lastPos-1)
        for (int seqPos = lastPos - 1; seqPos >= 0; seqPos--) {
            // Calculate state values for each observation in the sequence
            for (int i = 0; i < states.length; i++) {
                probabilityTable[i][seqPos] = 0;
                // Add up the previous probabilities
                for (int j = 0; j < states.length; j++) {
                    probabilityTable[i][seqPos] += (probabilityTable[i][seqPos + 1] * transitionMatrix[i][j] * emissionMatrix[i][obsMap
                            .get(sequenceData.charAt(seqPos) + 1)]);
                }
            }
        }
        return probabilityTable;
    }

    /**
     * @return the states
     */
    public static char[] getStates() {
        return states;
    }

    /**
     * @param states the states to set
     */
    public static void setStates(char[] states) {
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
    public static HashMap<Character, Integer> getStateMap() {
        return stateMap;
    }

    /**
     * @param stateMap the stateMap to set
     */
    public static void setStateMap() {
        // Create a reverse-mapping of the string state to the corresponding integer
        stateMap = new HashMap<Character, Integer>();
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