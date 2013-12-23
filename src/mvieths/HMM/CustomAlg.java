package mvieths.HMM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Pattern;

import mvieths.HMM.Util.HMMDataSet;

public class CustomAlg
{
    // Observations - An amino acid abbreviation (ARNDCQEGHILKMFPSTWYVBZX)
    char[]     observations = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X' };

    double[][] transitionMatrix;
    double[][] emissionMatrix;
    double[][] initialProbabilities;

    // States - beta sheet (e), helix (h), or loop (_)
    char[]     states       = { 'e', 'h', '_' };

    public CustomAlg(String trainingFile, String testFile)
    {
        ArrayList<HMMDataSet> trainingData;
        ArrayList<HMMDataSet> testData;

        trainingData = readSequences(trainingFile);
        testData = readSequences(testFile);
        System.out.println("There are " + trainingData.size() + " sequences in the training data");
        System.out.println("There are " + testData.size() + " sequences in the test data");

        trainData(trainingData);
        testData(testData);
        printData(testData);
    }

    public void trainData(ArrayList<HMMDataSet> trainingData)
    {
        double[][] transition = new double[states.length][states.length]; // Count how often we move from state to state
        // Initialize these with 0s
        for (int i = 0; i < transition.length; i++)
        {
            Arrays.fill(transition[i], 0.0);
        }

        double[][] emission = new double[observations.length][observations.length]; // Count how often we move from observation to observation
        // Initialize these with 0s
        for (int i = 0; i < emission.length; i++)
        {
            Arrays.fill(emission[i], 0.0);
        }

        // Initialize these with 0s
        double[][] initial = new double[states.length][observations.length];
        for (int i = 0; i < initial.length; i++)
        {
            Arrays.fill(initial[i], 0.0);
        }

        // Iterate through the training dataset, adding up the probabilities
        for (HMMDataSet data : trainingData)
        {
            data.generateProbabilities();
            double[][] tempTransition = data.getTransitionMatrix();
            double[][] tempEmission = data.getEmissionMatrix();
            HashMap<Character, Integer> obsMap = data.getObsMap();
            HashMap<Character, Integer> stateMap = data.getStateMap();

            // Add up the transition probabilities
            for (int i = 0; i < states.length; i++)
            {
                for (int j = 0; j < states.length; j++)
                {
                    transition[i][j] += tempTransition[i][j];
                }
            }

            // Add up the emission probabilities
            for (int i = 0; i < observations.length; i++)
            {
                for (int j = 0; j < observations.length; j++)
                {
                    emission[i][j] += tempEmission[i][j];
                }
            }

            // Add up the starting states
            char state = data.getStateData().charAt(0);
            char obs = data.getSequenceData().charAt(0);
            initial[stateMap.get(state)][obsMap.get(obs)]++;
        }

        // Divide each transition by the total number of datasets
        for (int i = 0; i < states.length; i++)
        {
            for (int j = 0; j < states.length; j++)
            {
                transition[i][j] /= trainingData.size();
            }
        }

        // Divide each emission by the total number of datasets
        for (int i = 0; i < observations.length; i++)
        {
            for (int j = 0; j < observations.length; j++)
            {
                emission[i][j] /= trainingData.size();
            }
        }

        transitionMatrix = transition;
        emissionMatrix = emission;
        initialProbabilities = initial;
    }

    /**
     * Utilize the matrices created by the training data to predict the states of the test data
     * 
     * @param testData
     */
    private void testData(ArrayList<HMMDataSet> testData)
    {
        for (HMMDataSet data : testData)
        {
            String sequence = data.getSequenceData();
            StringBuilder stateData = new StringBuilder();
            HashMap<Character, Integer> obsMap = data.getObsMap();
            HashMap<Character, Integer> stateMap = data.getStateMap();

            // Determine a starting state
            char firstChar = sequence.charAt(0);
            double largestVal = 0.0;
            int largestIndex = 0;
            for (int i = 0; i < states.length; i++)
            {
                if (initialProbabilities[i][obsMap.get(firstChar)] > largestVal)
                {
                    largestVal = initialProbabilities[i][obsMap.get(firstChar)];
                    largestIndex = i;
                }
            }
            char firstState = states[largestIndex];
            stateData.append(firstState);

            // Continue through the sequence
            for (int i = 1; i < sequence.length(); i++)
            {
                char prevChar = sequence.charAt(i - 1);
                char thisChar = sequence.charAt(i);
                char prevState = stateData.toString().charAt(i - 1);

                double maxProb = 0.0;
                int maxIndex = 0;
                // Figure out the most probable next state
                for (int j = 0; j < states.length; j++)
                {
                    double curProb = transitionMatrix[stateMap.get(prevState)][j] * emissionMatrix[obsMap.get(prevChar)][obsMap.get(thisChar)];
                    if (curProb > maxProb)
                    {
                        maxProb = curProb;
                        maxIndex = j;
                    }
                }
                stateData.append(states[maxIndex]);
            }
            data.setStateData(stateData.toString());
        }
    }

    private void printData(ArrayList<HMMDataSet> testData)
    {
        int counter = 0;
        for (HMMDataSet data : testData)
        {
            System.out.println("Dataset " + counter + ":");
            System.out.println(data.getSequenceData());
            System.out.println(data.getStateData());
            counter++;
        }
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
                        dataSets.add(new HMMDataSet(sequenceData.toString(), stateData.toString(), observations, states));
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

}
