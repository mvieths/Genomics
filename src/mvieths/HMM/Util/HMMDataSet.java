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
    char[]                      states;

    // The possible observations for this Markov model
    char[]                      observations;

    // The sequence we've been provided.  Must consist entirely of values from 'observations'
    String                      sequenceData;

    // The state data we've been provided.  Must consist entirely of values from 'states'.  Must be of equal length to sequenceData.
    String                      stateData;

    // The transition matrix contains the probability of changing from one state to another
    double[][]                  transitionMatrix;

    // The emission matrix contains the probability of changing from one observation to another given a particular state
    double[][]                  emissionMatrix;

    // The initial probability is the probability used for the first position in the sequence
    double[]                    initialProbability;

    // State table for Viterbi
    int[][]                     stateTable;

    // Reverse mapping from string to integer for looking up a state
    HashMap<Character, Integer> stateMap;

    // Reverse mapping from character to integer for looking up an observation
    HashMap<Character, Integer> obsMap;

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
    public double[][] forward() {
        double[][] probabilityTable = new double[states.length][sequenceData.length()];

        // Set the initial probabilities starting at position 0 in the sequenceData
        for (int i = 0; i < states.length; i++) {
            if (initialProbability[i] == 0 || emissionMatrix[i][obsMap.get(sequenceData.charAt(0))] == 0)
            {
                probabilityTable[i][0] = 0.0;
            }
            else {
                probabilityTable[i][0] = initialProbability[i] * emissionMatrix[i][obsMap.get(sequenceData.charAt(0))];
            }
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
    public double[][] backward() {
        double[][] probabilityTable = new double[states.length][sequenceData.length()];

        int lastPos = sequenceData.length() - 1;
        char lastChar = sequenceData.charAt(lastPos);
        // Set the initial probabilities starting at position 0 in the sequenceData
        for (int i = 0; i < states.length; i++) {
            if (initialProbability[i] == 0 || emissionMatrix[i][obsMap.get(lastChar)] == 0)
            {
                probabilityTable[i][lastPos] = 0.0;
            }
            else {
                probabilityTable[i][lastPos] = initialProbability[i] * emissionMatrix[i][obsMap.get(lastChar)];
            }
        }

        // Continue populating the tables starting at sequenceData.charAt(lastPos-1)
        for (int seqPos = lastPos - 1; seqPos >= 0; seqPos--) {
            // Calculate state values for each observation in the sequence
            for (int i = 0; i < states.length; i++) {
                probabilityTable[i][seqPos] = 0;
                // Add up the previous probabilities
                for (int j = 0; j < states.length; j++) {
                    probabilityTable[i][seqPos] += (probabilityTable[i][seqPos + 1] * transitionMatrix[i][j] * emissionMatrix[i][obsMap
                            .get(sequenceData.charAt(seqPos + 1))]);
                }
            }
        }
        return probabilityTable;
    }

    public double[][] viterbi() {
        double[][] probabilityTable = new double[states.length][sequenceData.length()];
        stateTable = new int[states.length][sequenceData.length()];

        // Set the initial probabilities and states
        char firstChar = sequenceData.charAt(0);
        for (int i = 0; i < states.length; i++) {
            // Probability is 0 if emission probability is 0 (to avoid Math.log issues)
            if (emissionMatrix[i][obsMap.get(firstChar)] != 0.0) {
                probabilityTable[i][0] = Math.log(initialProbability[i] * emissionMatrix[i][obsMap.get(firstChar)]);
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
                    if (emissionMatrix[k][obsMap.get(thisChar)] != 0.0) {
                        prob = probabilityTable[k][seqPos - 1] + Math.log(transitionMatrix[k][state])
                                + Math.log(emissionMatrix[k][obsMap.get(thisChar)]);
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

    /**
     * @return the states
     */
    public char[] getStates() {
        return states;
    }

    /**
     * @param states the states to set
     */
    public void setStates(char[] states) {
        this.states = states;
        setStateMap();
    }

    /**
     * @return the observations
     */
    public char[] getObservations() {
        return observations;
    }

    /**
     * @param observations the observations to set
     */
    public void setObservations(char[] observations) {
        this.observations = observations;
        setObsMap();
    }

    /**
     * @return the transitionMatrix
     */
    public double[][] getTransitionMatrix() {
        return transitionMatrix;
    }

    /**
     * @param transitionMatrix the transitionMatrix to set
     */
    public void setTransitionMatrix(double[][] transitionMatrix) {
        this.transitionMatrix = transitionMatrix;
    }

    /**
     * @return the emissionMatrix
     */
    public double[][] getEmissionMatrix() {
        return emissionMatrix;
    }

    /**
     * @param emissionMatrix the emissionMatrix to set
     */
    public void setEmissionMatrix(double[][] emissionMatrix) {
        this.emissionMatrix = emissionMatrix;
    }

    /**
     * @return the initialProbability
     */
    public double[] getInitialProbability() {
        return initialProbability;
    }

    /**
     * @param initialProbability the initialProbability to set
     */
    public void setInitialProbability(double[] initialProbability) {
        this.initialProbability = initialProbability;
    }

    /**
     * @return the stateMap
     */
    public HashMap<Character, Integer> getStateMap() {
        return stateMap;
    }

    /**
     * @param stateMap the stateMap to set
     */
    public void setStateMap() {
        // Create a reverse-mapping of the string state to the corresponding integer
        stateMap = new HashMap<Character, Integer>();
        for (int i = 0; i < states.length; i++) {
            stateMap.put(states[i], i);
        }
    }

    /**
     * @return the obsMap
     */
    public HashMap<Character, Integer> getObsMap() {
        return obsMap;
    }

    /**
     * @param obsMap the obsMap to set
     */
    public void setObsMap() {
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