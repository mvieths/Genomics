/**
 * 
 */
package mvieths.HMM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
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
    static ArrayList<HMMDataSet> trainingData;
    static ArrayList<HMMDataSet> testData;

    // Observations - An amino acid abbreviation (ARNDCQEGHILKMFPSTWYVBZX)
    static char[]                observations       =
                                                    { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
                                                    'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X' };

    // States - beta sheet (e), helix (h), or loop (_)
    static char[]                states             =
                                                    { 'e', 'h', '_' };

    static double                evenState          = 1.0 / 3.0;
    static double                evenObs            = 1.0 / 20.0;

    static double[]              initialProbability =
                                                    { evenState, evenState, evenState };
    static double[]              currentProbability;

    // Probability of the state changing (in this case double[3][3])
    static double[][]            transitionMatrix   =
                                                    {
                                                    { evenState, evenState, evenState },
                                                    { evenState, evenState, evenState },
                                                    { evenState, evenState, evenState }
                                                    };
    static double[][]            currentTransitionMatrix;

    // Probability of the observation changing (in this case double[3][20])
    static double[][]            emissionMatrix     =
                                                    {
            { evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs,
            evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs },
            { evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs,
            evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs },
            { evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs,
            evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs, evenObs }
                                                    };
    static double[][]            currentEmissionMatrix;

    /**
     * @param args
     */
    public static void main(String[] args)
    {
        currentProbability = Arrays.copyOf(initialProbability, initialProbability.length);
        currentTransitionMatrix = Arrays.copyOf(transitionMatrix, transitionMatrix.length);
        currentEmissionMatrix = Arrays.copyOf(emissionMatrix, emissionMatrix.length);

        trainingData = readSequences(args[0]);
        testData = readSequences(args[1]);

        System.out.println("There are " + trainingData.size() + " sequences in the training data");
        System.out.println("There are " + testData.size() + " sequences in the test data");
    }

    public void trainData() {
        for (HMMDataSet data : trainingData) {
            double[][] forward = data.forward();
            double[][] backward = data.backward();

            // Figure out new initial probabilities
            for (int i = 0; i < states.length; i++) {

            }

            // Recalculate the transition matrix
            for (int i = 0; i < states.length; i++) {
                for (int j = 0; j < states.length; j++) {

                }
            }

            // Recalculate the emission matrix
            for (int i = 0; i < states.length; i++) {
                for (int j = 0; j < observations.length; j++) {

                }
            }
        }
    }

    /**
     * Read in the provided file, producing a list of all the datasets
     * 
     * @param filename
     *            for a training or test dataset list
     * @return An ArrayList of the data sets
     */
    static ArrayList<HMMDataSet> readSequences(String filename)
    {
        ArrayList<HMMDataSet> dataSets = new ArrayList<HMMDataSet>();

        try
        {
            BufferedReader reader = new BufferedReader(new FileReader(filename));

            // Read the sequence from the file
            String line;
            StringBuilder sequenceData = new StringBuilder();
            StringBuilder stateData = new StringBuilder();
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
                        // HMMDataSet(String sequence, String stateData, char[] observations, String[] states,
                        //     double[] initialProbability, double[][] transitionMatrix, double[][] emissionMatrix)
                        HMMDataSet data = new HMMDataSet(sequenceData.toString(), stateData.toString(), observations,
                                states, currentProbability, currentTransitionMatrix, currentEmissionMatrix);
                        dataSets.add(data);
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