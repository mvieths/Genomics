package mvieths.HMM;

/**
 * Stores a score for the aligned segment and the string version of both alignments
 * 
 * @author Foeclan
 * 
 */
public class AlignedSegment {
	private int score;
	private String firstAlignment;
	private String secondAlignment;

	public AlignedSegment(int score, String alignment1, String alignment2) {
		setScore(score);
		setFirstAlignment(alignment1);
		setSecondAlignment(alignment2);
	}

	public int getScore() {
		return score;
	}

	public void setScore(int score) {
		this.score = score;
	}

	public String getFirstAlignment() {
		return firstAlignment;
	}

	public void setFirstAlignment(String alignment) {
		this.firstAlignment = alignment;
	}

	public String getSecondAlignment() {
		return secondAlignment;
	}

	public void setSecondAlignment(String secondAlignment) {
		this.secondAlignment = secondAlignment;
	}

	public double getAlignmentPercent() {
		double perfect = 0;
		for (int i = 0; i < firstAlignment.length(); i++) {
			if (firstAlignment.charAt(i) == secondAlignment.charAt(i)) {
				perfect++;
			}
		}

		double percent = perfect / firstAlignment.length();
		return percent;
	}
}
