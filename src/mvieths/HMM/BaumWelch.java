/**
 * 
 */
package mvieths.HMM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.regex.Pattern;

import mvieths.HMM.Util.HMMDataSet;

/**
 * @author Michael Vieths
 * 
 */
public class BaumWelch
{
    static ArrayList<HMMDataSet> trainingData;
    static ArrayList<HMMDataSet> testData;

    /**
     * @param args
     */
    public static void main(String[] args)
    {
        trainingData = readSequences(args[0]);
        testData = readSequences(args[1]);

        System.out.println("There are " + trainingData.size() + " sequences in the training data");
        System.out.println("There are " + testData.size() + " sequences in the test data");
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
            StringBuilder sequence = new StringBuilder();
            StringBuilder states = new StringBuilder();
            /*
             * Read through the provided file, storing anything between <> with a sequence and its corresponding states
             */
            while ((line = reader.readLine()) != null)
            {
                String[] columns = line.split(" ");
                if (Pattern.matches("[ARNDCQEGHILKMFPSTWYVBZX]", columns[0]))
                {
                    sequence.append(columns[0]);
                    states.append(columns[1]);
                }
                else if (columns[0].equals("<>"))
                {
                    if (sequence.length() > 0)
                    {
                        // Store anything between <> that's in the right format
                        HMMDataSet data = new HMMDataSet(sequence.toString(), states.toString());
                        dataSets.add(data);
                        sequence = new StringBuilder();
                        states = new StringBuilder();
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