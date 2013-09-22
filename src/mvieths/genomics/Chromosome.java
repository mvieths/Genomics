package mvieths.genomics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Chromosome {
    /**
     * This class will take a FASTA chromosome file and provide helper methods for analyzing it
     */
    String                       seqFileName;
    String                       geneFileName;
    long                         gContent;             // Running count of guanine nucleotides
    long                         cContent;             // Running count of cytosine nucleotides
    long                         aContent;             // Running count of adenine nucleotides
    long                         tContent;             // Running count of thymine nucleotides
    long                         chromosomeLength;     // Length of G+C+A+T
    long                         totalChromosomeLength; // Length of G+C+A+T+N
    static ArrayList<GCATWindow> windows;              // Windows of set length in the chromosome and their g/c/a/t
    static ArrayList<GeneWindow> genes;                // Locations of genes in this chromosome

    public Chromosome(String seqFile, String geneFile) {
        gContent = 0;
        cContent = 0;
        aContent = 0;
        tContent = 0;
        chromosomeLength = 0;
        totalChromosomeLength = 0;

        try {
            seqFileName = seqFile;
            geneFileName = geneFile;
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public Chromosome(String seqFile) {
        this(seqFile, "");
    }

    /**
     * Calculates the G, C, A and T content of the chromosome
     */
    void calcContent() {
        System.out.println("Calculating nucleotide content...");
        String line;
        try {
            BufferedReader reader = new BufferedReader(new FileReader(seqFileName));

            while ((line = reader.readLine()) != null) {
                countLine(line);
            }
            // Add up all of the nucleotides we count
            chromosomeLength = gContent + cContent + aContent + tContent;
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Calculates the G, C, A and T content in increments of the specified window size
     */
    void calcContentByWindow(long windowSize) {
        System.out.println("Building nucleotide content list...");
        try {
            windows = new ArrayList<GCATWindow>();
            FileReader fileReader = new FileReader(seqFileName);
            BufferedReader lineReader = new BufferedReader(fileReader);
            lineReader.readLine(); // Get rid of the first line, which is a comment
            GCATWindow window;
            long counter = 0;
            int c;
            long gValue = 0;
            long cValue = 0;
            long aValue = 0;
            long tValue = 0;

            while ((c = fileReader.read()) != -1) {
                counter++;
                switch (c) {
                    case 'g':
                    case 'G':
                        gValue++;
                        break;
                    case 'c':
                    case 'C':
                        cValue++;
                        break;
                    case 'a':
                    case 'A':
                        aValue++;
                        break;
                    case 't':
                    case 'T':
                        tValue++;
                        break;
                }

                if (counter % windowSize == 0) {
                    window = new GCATWindow(counter, gValue, cValue, aValue, tValue);
                    windows.add(window);

                    gValue = 0;
                    cValue = 0;
                    aValue = 0;
                    tValue = 0;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Processes the geneFile for genes and stores their name, starting and ending positions in an ArrayList of
     * GeneWindow objects
     */
    void calcGenes() {
        System.out.println("Building gene list...");
        FileReader fileReader;
        try {
            if (geneFileName.equals("")) {
                System.out.println("No gene file specified");
                return;
            } else {
                fileReader = new FileReader(geneFileName);
                BufferedReader lineReader = new BufferedReader(fileReader);
                String line;
                genes = new ArrayList<GeneWindow>();
                GeneWindow gene;
                while ((line = lineReader.readLine()) != null) {
                    String[] tokens = line.split("\t");
                    if (tokens[5].equals("'+'")) {
                        gene = new GeneWindow(Long.parseLong(tokens[2]), Long.parseLong(tokens[3]), tokens[0]);
                    } else {
                        // Change this to take into account counting from the end of the chromosome
                        long start = totalChromosomeLength - Long.parseLong(tokens[3]);
                        long end = totalChromosomeLength - Long.parseLong(tokens[2]);
                        gene = new GeneWindow(start, end, tokens[0]);
                    }
                    genes.add(gene);
                }
            }
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    /**
     * Counts the G/C/A/T content in a given line of text. Skips anything that's not a G/C/A/T, and comment lines (those
     * starting with a >). Case insensitive. Stores the content in the global variables gContent, cContent, aContent,
     * and tContent, with the total stored in chromosomeLength
     * 
     * @param line
     *            A line of text to analyze
     */
    private void countLine(String line) {
        if (line.startsWith(">")) {
            // Comment, ignore it
        } else {
            totalChromosomeLength += line.length();
            for (int i = 0; i < line.length(); i++) {
                char c = line.charAt(i);
                switch (c) {
                    case 'g':
                    case 'G':
                        gContent++;
                        break;
                    case 'c':
                    case 'C':
                        cContent++;
                        break;
                    case 'a':
                    case 'A':
                        aContent++;
                        break;
                    case 't':
                    case 'T':
                        tContent++;
                        break;
                }
            }
        }
    }

    /**
     * Print the percentage of G+C nucleotides in this chromosome
     */
    void printGCPercent() {
        double gcPercent = (((double) gContent + (double) cContent) / chromosomeLength) * 100;

        System.out.println("G = " + gContent + " C = " + cContent + " Total Length = " + chromosomeLength
                + "\nGC Percent = " + gcPercent);
    }

    /**
     * Print the percentage of G/C/A/T nucleotides in this chromosome
     */
    void printAllPercents() {
        double gPercent = ((double) gContent / chromosomeLength) * 100;
        double cPercent = ((double) cContent / chromosomeLength) * 100;
        double aPercent = ((double) aContent / chromosomeLength) * 100;
        double tPercent = ((double) tContent / chromosomeLength) * 100;
        double totalPercent = gPercent + cPercent + aPercent + tPercent;
        System.out.println("gPercent = " + gPercent + "\ncPercent = " + cPercent + "\naPercent = " + aPercent
                + "\ntPercent = " + tPercent + "\nTotal Percent = " + totalPercent);

    }

    /**
     * Provides the number of G+C+A+T nucleotides. Note that 'N' or other values are not included in this value
     * 
     * @return Number of G+C+A+T nucleotides
     */
    public long getChromosomeLength() {
        return chromosomeLength;
    }

    /**
     * Provides the list of G/C/A/T windows in this chromosome
     * 
     * @return An ArrayList of the G/C/A/T statistics for this chromosome
     */
    public ArrayList<GCATWindow> getWindows() {
        return windows;
    }

    /**
     * Provides the list of genes in this chromosome and their positions
     * 
     * @return Provides an ArrayList of GeneWindow objects in this chromosome and their positions
     */
    public ArrayList<GeneWindow> getGenes() {
        return genes;
    }
}

/**
 * 
 * @author Foeclan
 * 
 */
class GCATWindow {
    private long       gContent;
    private long       cContent;
    private long       aContent;
    private long       tContent;
    private final long total;
    private final long windowNumber;

    public GCATWindow(long number, long g, long c, long a, long t) {
        windowNumber = number;
        setgContent(g);
        setcContent(c);
        setaContent(a);
        settContent(t);
        total = g + c + a + t;
    }

    /**
     * @param gContent
     *            the gContent to set
     */
    public void setgContent(long gContent) {
        this.gContent = gContent;
    }

    /**
     * @return the gContent
     */
    public long getgContent() {
        return gContent;
    }

    public double getgPercent() {
        return ((double) gContent / total) * 100;
    }

    /**
     * @param cContent
     *            the cContent to set
     */
    public void setcContent(long cContent) {
        this.cContent = cContent;
    }

    /**
     * @return the cContent
     */
    public long getcContent() {
        return cContent;
    }

    public double getcPercent() {
        return ((double) cContent / total) * 100;
    }

    /**
     * @param aContent
     *            the aContent to set
     */
    public void setaContent(long aContent) {
        this.aContent = aContent;
    }

    /**
     * @return the aContent
     */
    public long getaContent() {
        return aContent;
    }

    public double getaPercent() {
        return ((double) aContent / total) * 100;
    }

    /**
     * @param tContent
     *            the tContent to set
     */
    public void settContent(long tContent) {
        this.tContent = tContent;
    }

    /**
     * @return the tContent
     */
    public long gettContent() {
        return tContent;
    }

    public long gettPercent() {
        double percent = ((double) tContent / total) * 100;
        if (Double.isNaN(percent)) {
            percent = 0;
        }
        if (percent > 100.0) {
            System.out.println("Invalid percentage " + percent);
        }
        return (long) percent;
    }

    public long getGCContent() {
        return gContent + cContent;
    }

    public double getGCPercent() {
        long gcContent = gContent + cContent;
        double gcPercent = ((double) gcContent / total) * 100;

        return gcPercent;
    }

    /**
     * @return the windowNumber
     */
    public String getWindowNumberAsString() {
        return "" + windowNumber;
    }

    public long getWindowNumber() {
        return windowNumber;
    }
}

class GeneWindow {
    private long   start;
    private long   end;
    private String name;

    public GeneWindow(long start, long end, String name) {
        setStart(start);
        setEnd(end);
        setName(name);
    }

    /**
     * @param start
     *            the start to set
     */
    public void setStart(long start) {
        this.start = start;
    }

    /**
     * @return the start
     */
    public long getStart() {
        return start;
    }

    /**
     * @param end
     *            the end to set
     */
    public void setEnd(long end) {
        this.end = end;
    }

    /**
     * @return the end
     */
    public long getEnd() {
        return end;
    }

    /**
     * @param name
     *            the name to set
     */
    public void setName(String name) {
        this.name = name;
    }

    /**
     * @return the name
     */
    public String getName() {
        return name;
    }
}