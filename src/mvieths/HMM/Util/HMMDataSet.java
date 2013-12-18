package mvieths.HMM.Util;

/**
 * Stores a dataset, which consists of an amino acid sequence and its corresponding states
 * 
 * @author Michael Vieths
 * 
 */
public class HMMDataSet
{
    String   sequenceData;
    String   stateData;
    String[] states;

    public HMMDataSet(String sequence, String states) throws Exception
    {
        if (sequence.length() != states.length())
        {
            throw new Exception("Sequence and state lengths do not match");
        }
        else
        {
            setSequenceData(sequence);
            setStateData(states);
        }
    }

    public String getSequenceData()
    {
        return sequenceData;
    }

    public void setSequenceData(String sequence)
    {
        this.sequenceData = sequence;
    }

    public String getStateData()
    {
        return stateData;
    }

    public void setStateData(String states)
    {
        this.stateData = states;
    }
}