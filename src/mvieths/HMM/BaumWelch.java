/**
 * 
 */
package mvieths.HMM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Pattern;

import mvieths.HMM.Util.HMMDataSet;

/**
 * @author Michael Vieths
 * 
 */
public class BaumWelch
{
    // HMMDataSet(String sequence, String stateData, char[] observations, String[] states,
    //     double[] initialProbability, double[][] transitionMatrix, double[][] emissionMatrix)

    // Observations - An amino acid abbreviation (ARNDCQEGHILKMFPSTWYVBZX)
    char[]                      observations       = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B',
                    'Z', 'X'                      };

    HashMap<Character, Integer> obsMap;

    // States - beta sheet (e), helix (h), or loop (_)
    char[]                      states             = { 'e', 'h', '_' };

    double                      evenState          = 1.0 / 3.0;
    double                      evenObs            = 1.0 / 23.0;

    double[]                    initialProbability = { evenState, evenState, evenState };

    // Probability of the state changing (in this case double[3][3])
    //    double[][]                  transitionMatrix   =
    //                                                   {
    //                                                   { evenState, evenState, evenState },
    //                                                   { evenState, evenState, evenState },
    //                                                   { evenState, evenState, evenState }
    //                                                   };
    double[][]                  transitionMatrix   = { { 0.5, 0.3, 0.2 }, { 0.3, 0.2, 0.5 }, { 0.2, 0.5, 0.3 } };

    // Probability of the observation changing
    double[][]                  emissionMatrix     = {
                    { evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs,
                    evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs },
                    { evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs,
                    evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs },
                    { evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs,
                    evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs } };
    double[]                    currentProbability;
    double[][]                  currentTransitionMatrix;
    double[][]                  currentEmissionMatrix;

    public BaumWelch(String trainingFile, String testFile)
    {
        /**
         * @param args
         */
        currentProbability = Arrays.copyOf(initialProbability, initialProbability.length);
        currentTransitionMatrix = Arrays.copyOf(transitionMatrix, transitionMatrix.length);
        currentEmissionMatrix = Arrays.copyOf(emissionMatrix, emissionMatrix.length);
        ArrayList<HMMDataSet> trainingData;
        ArrayList<HMMDataSet> testData;

        trainingData = readSequences(trainingFile);
        testData = readSequences(testFile);
        System.out.println("There are " + trainingData.size() + " sequences in the training data");
        System.out.println("There are " + testData.size() + " sequences in the test data");

        HMMDataSet viterbiData = trainingData.get(0);
        viterbiData.setInitialProbability(initialProbability);
        viterbiData.setTransitionMatrix(transitionMatrix);
        viterbiData.setEmissionMatrix(emissionMatrix);

        double[][] viterbiProb = viterbiData.viterbi();
        char[] stateTable = buildStateTable(viterbiProb);
        System.out.println("Viterbi data:");
        System.out.println(viterbiData.getSequenceData());
        System.out.println(new String(stateTable));

        //        System.out.println("Training data");
        //        printDataset(trainingData);
        //        System.out.println("Test data");
        //        printDataset(testData);
        trainData(trainingData);
        testData(testData);
    }

    public void printDataset(ArrayList<HMMDataSet> dataset)
    {
        for (HMMDataSet data : dataset)
        {
            System.out.println(data.getSequenceData());
            System.out.println(data.getStateData());
        }
    }

    public void trainData(ArrayList<HMMDataSet> trainingData)
    {
        System.out.println("Before training:");
        print();
        double[] probability = Arrays.copyOf(initialProbability, initialProbability.length);
        double[][] tMatrix = Arrays.copyOf(transitionMatrix, transitionMatrix.length);
        double[][] eMatrix = Arrays.copyOf(emissionMatrix, emissionMatrix.length);

        int counter = 0; // Disables printing the first few results
        for (HMMDataSet data : trainingData)
        {
            // Reset the matrices and initial probabilities for the dataset
            double[][] forward = data.forward();
            double[][] backward = data.backward();
            obsMap = data.getObsMap();

            // Figure out new initial probabilities
            for (int i = 0; i < states.length; i++)
            {
                probability[i] = calcProbAtT(0, i, forward, backward);
            }

            // Recalculate the transition matrix
            for (int i = 0; i < states.length; i++)
            {
                for (int j = 0; j < states.length; j++)
                {
                    double xi = 0.0;
                    double gamma = 0.0;

                    for (int t = 0; t < data.getSequenceData().length(); t++)
                    {
                        xi += calcStateProb(t, i, j, data.getSequenceData(), forward, backward);
                        gamma += calcProbAtT(t, i, forward, backward);
                    }

                    double value = 0.0;
                    if (gamma > 0.0)
                    {
                        value = xi / gamma;
                    }

                    tMatrix[i][j] = value;
                }
            }

            // Recalculate the emission matrix
            for (int i = 0; i < states.length; i++)
            {
                for (int j = 0; j < observations.length; j++)
                {
                    double transitionsTo = 0.0;
                    double transitionsFrom = 0.0;

                    String sequence = data.getSequenceData();

                    for (int t = 0; t < sequence.length(); t++)
                    {
                        double gamma = calcProbAtT(t, i, forward, backward);
                        double multiplier = 0;
                        if (obsMap.get(sequence.charAt(t)) == j)
                        {
                            multiplier = 1.0;
                        }
                        transitionsTo += gamma * multiplier;
                        transitionsFrom += gamma;
                    }

                    double value = 0.0;
                    if (transitionsFrom > 0.0)
                    {
                        value = transitionsTo / transitionsFrom;
                    }

                    eMatrix[i][j] = value;
                }
            }

            currentProbability = probability;
            currentTransitionMatrix = tMatrix;
            currentEmissionMatrix = eMatrix;

            if (counter < 5)
            {
                System.out.println("Training values (" + counter + ") :");
                print();
            }
            counter++;
        }

        System.out.println("Final matrices:");
        print();
    }

    /**
     * Uses the currentProbability, currentTransitionMatrix, and currentEmissionMatrix to generate state maps for the test data
     */
    public void testData(ArrayList<HMMDataSet> testData)
    {
        int counter = 0;
        for (HMMDataSet data : testData)
        {
            // Reset the matrices and initial probabilities for the dataset
            data.setInitialProbability(currentProbability);
            data.setTransitionMatrix(currentTransitionMatrix);
            data.setEmissionMatrix(currentEmissionMatrix);

            double[][] forward = data.forward();
            double[][] backward = data.backward();
            char[] forwardStates = buildStateTable(forward);
            char[] backwardStates = buildStateTable(backward);

            System.out.println("Seq[" + counter + "]: " + data.getSequenceData());
            System.out.println("Fwd[" + counter + "]: " + new String(forwardStates));
            System.out.println("Bwd[" + counter + "]: " + new String(backwardStates));
            counter++;
        }
    }

    /**
     * Calculates the probability of state i at time t a_i(t)*b_i(t) / (sum of the products of a and b) aka gamma (from wikipedia article)
     * 
     * @param t
     *            The position in the sequence
     * @param i
     *            stuff
     * @param forward
     *            Forward probabilities
     * @param backward
     *            Backward probabilities
     * @return Probability of being in state i at time t
     */
    public double calcProbAtT(int t, int i, double[][] forward, double[][] backward)
    {
        // Calculate the product of the forward and backward probabilities at time t
        double ab = forward[i][t] * backward[i][t];
        double sum = 0.0;

        // Calculate the sum of the products of the forward and backward algorithms
        for (int j = 0; j < forward.length; j++)
        {
            sum += forward[j][t] * backward[j][t];
        }
        double retVal = 0.0;
        if (sum > 0.0)
        {
            retVal = ab / sum;
        }
        return retVal;
    }

    /**
     * Calculates the probability of being in state i and then state j aka xi (from wikipedia article)
     * 
     * @param t
     *            Starting position in the sequence (will also consider t+1)
     * @param i
     *            State at time t
     * @param j
     *            State at time t+1
     * @param sequence
     *            The training sequence we're considering
     * @param forward
     *            Forward probabilities
     * @param backward
     *            Backward probabilities
     * @return Probability of being in state i and then state j
     */
    public double calcStateProb(int t, int i, int j, String sequence, double[][] forward, double[][] backward)
    {
        // Find the product of the forward algorithm at t, backward at t+1, and the emission and transmission probabilities
        double ab = forward[i][t] * currentTransitionMatrix[i][j];
        // Mind the end of the array (don't try to get something that doesn't exist)
        if (t != (sequence.length() - 1))
        {
            ab *= (currentEmissionMatrix[j][obsMap.get(sequence.charAt(t + 1))] * backward[j][t + 1]);
        }

        double sum = 0.0;
        // Calculate the sum of the products of the forward and backward algorithms
        for (int k = 0; k < forward.length; k++)
        {
            sum += forward[k][t] * backward[k][t];
        }
        double retVal = 0.0;
        if (sum > 0.0)
        {
            retVal = ab / sum;
        }
        return retVal;
    }

    public char[] buildStateTable(double[][] results)
    {
        char[] stateTable = new char[results[0].length];

        for (int i = 0; i < results[0].length; i++)
        {
            double maxProb = 0.0;
            int maxIndex = 0;
            for (int j = 0; j < results.length; j++)
            {
                double currentProb = Math.pow(Math.E, results[j][i]);
                if (currentProb > maxProb)
                {
                    maxProb = currentProb;
                    maxIndex = j;
                }
            }
            stateTable[i] = states[maxIndex];
        }

        return stateTable;
    }

    /**
     * Read in the provided file, producing a list of all the datasets
     * 
     * @param filename
     *            for a training or test dataset list
     * @return An ArrayList of the data sets
     */
    private ArrayList<HMMDataSet> readSequences(String filename)
    {
        ArrayList<HMMDataSet> dataSets = new ArrayList<HMMDataSet>();
        StringBuilder sequenceData = new StringBuilder();
        StringBuilder stateData = new StringBuilder();

        try
        {
            BufferedReader reader = new BufferedReader(new FileReader(filename));

            // Read the sequence from the file
            String line;
            /*
             * Read through the provided file, storing anything between <> with a sequence and its corresponding states
             */
            while ((line = reader.readLine()) != null)
            {
                String[] columns = line.split(" ");
                if (Pattern.matches("[ARNDCQEGHILKMFPSTWYVBZX]", columns[0]))
                {
                    sequenceData.append(columns[0]);
                    stateData.append(columns[1]);
                }
                else if (columns[0].equals("<>"))
                {
                    if (sequenceData.length() > 0)
                    {
                        // Store anything between <> that's in the right format
                        dataSets.add(new HMMDataSet(sequenceData.toString(), stateData.toString(), observations, states, initialProbability, transitionMatrix,
                                                    emissionMatrix));
                        sequenceData = new StringBuilder();
                        stateData = new StringBuilder();
                    }
                }
            }
        }
        catch (Exception e)
        {
            System.out.println(e.getMessage());
        }

        return dataSets;
    }

    public void print()
    {
        // Print the initial probabilities
        System.out.println("=== Initial Probabilities ===");
        for (int i = 0; i < currentProbability.length; i++)
        {
            //            System.out.printf("[%.4f] ", currentProbability[i]);
            System.out.printf("[%.4f] ", Math.log(currentProbability[i]));
        }
        System.out.println();

        // Print the transition matrix
        System.out.println("=== Transition Matrix ===");
        for (int i = 0; i < currentTransitionMatrix.length; i++)
        {
            for (int j = 0; j < currentTransitionMatrix[0].length; j++)
            {
                //                System.out.printf("[%.4f] ", currentTransitionMatrix[i][j]);
                System.out.printf("[%.4f] ", Math.log(currentTransitionMatrix[i][j]));
            }
            System.out.println();
        }

        // Print the emission matrix
        System.out.println("=== Emission Matrix ===");
        for (int i = 0; i < currentEmissionMatrix.length; i++)
        {
            for (int j = 0; j < currentEmissionMatrix[0].length; j++)
            {
                //                System.out.printf("[%.4f] ", currentEmissionMatrix[i][j]);
                System.out.printf("[%.4f] ", Math.log(currentEmissionMatrix[i][j]));
            }
            System.out.println();
        }
    }

}