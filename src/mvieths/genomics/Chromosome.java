package mvieths.genomics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Chromosome {
    /**
     * This class will take a FASTA chromosome file and provide helper methods for analyzing it
     */
    String                       fileName;
    long                         gContent;             // Running count of guanine nucleotides
    long                         cContent;             // Running count of cytosine nucleotides
    long                         aContent;             // Running count of adenine nucleotides
    long                         tContent;             // Running count of thymine nucleotides
    long                         chromosomeLength;     // Length of G+C+A+T
    long                         totalChromosomeLength; // Length of G+C+A+T+N
    static ArrayList<GCATWindow> windows;              // Windows of set length in the chromosome and their g/c/a/t

    public Chromosome(String filename) {
        gContent = 0;
        cContent = 0;
        aContent = 0;
        tContent = 0;
        chromosomeLength = 0;
        totalChromosomeLength = 0;

        try {
            fileName = filename;
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * @param args
     */
    // public static void main(String[] args) {
    // Chromosome chr = new Chromosome(args[0]);
    // chr.calcContent();
    // chr.printGCPercent();
    // chr.printAllPercents();
    // chr.calcContentByWindow(200);
    // System.out.println("There are " + windows.size() + " windows");
    // }

    /**
     * Calculates the G, C, A and T content of the chromosome
     */
    void calcContent() {
        String line;
        try {
            BufferedReader reader = new BufferedReader(new FileReader(fileName));

            while ((line = reader.readLine()) != null) {
                countLine(line);
            }
            // Add up all of the nucleotides we count
            chromosomeLength = gContent + cContent + aContent + tContent;
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public long getChromosomeLength() {
        return chromosomeLength;
    }

    /**
     * Calculates the G, C, A and T content in increments of the specified window size
     */
    void calcContentByWindow(long windowSize) {
        try {
            windows = new ArrayList<GCATWindow>();
            FileReader fileReader = new FileReader(fileName);
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

    void printGCPercent() {
        double gcPercent = (((double) gContent + (double) cContent) / chromosomeLength) * 100;

        System.out.println("G = " + gContent + " C = " + cContent + " Total Length = " + chromosomeLength
                + "\nGC Percent = " + gcPercent);
    }

    void printAllPercents() {
        double gPercent = ((double) gContent / chromosomeLength) * 100;
        double cPercent = ((double) cContent / chromosomeLength) * 100;
        double aPercent = ((double) aContent / chromosomeLength) * 100;
        double tPercent = ((double) tContent / chromosomeLength) * 100;
        double totalPercent = gPercent + cPercent + aPercent + tPercent;
        System.out.println("gPercent = " + gPercent + "\ncPercent = " + cPercent + "\naPercent = " + aPercent
                + "\ntPercent = " + tPercent + "\nTotal Percent = " + totalPercent);

    }

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

    public ArrayList<GCATWindow> getWindows() {
        return windows;
    }
}

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

    /**
     * @return the windowNumber
     */
    public String getWindowNumber() {
        return "" + windowNumber;
    }
}