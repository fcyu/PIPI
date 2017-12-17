package proteomics.Types;

public class ResultEntry implements Comparable<ResultEntry> {

    public final double score;
    public final String peptide;
    public final char leftFlank;
    public final char rightFlank;
    private final int hashCode;
    private final boolean isDecoy;
    private final boolean containPtm;

    public ResultEntry(double score, String peptide, char leftFlank, char rightFlank, boolean isDecoy, boolean containPtm) {
        this.score = score;
        this.peptide = peptide;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
        String toString = peptide + "-" + score;
        hashCode = toString.hashCode();
        this.isDecoy = isDecoy;
        this.containPtm = containPtm;
    }

    public int compareTo(ResultEntry other) {
        if (score > other.score) {
            return 1;
        } else if (score < other.score) {
            return -1;
        } else {
            return 0;
        }
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public boolean containPTM() {
        return containPtm;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof ResultEntry) {
            ResultEntry temp = (ResultEntry) other;
            return ((temp.peptide.contentEquals(peptide)) && (temp.score == score));
        } else {
            return false;
        }
    }
}