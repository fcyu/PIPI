package proteomics.Types;

public class ResultEntry implements Comparable<ResultEntry> {

    public final double score;
    public final String peptide;
    private final int hashCode;
    private final boolean isDecoy;

    public ResultEntry(double score, String peptide, boolean isDecoy) {
        this.score = score;
        this.peptide = peptide;
        String toString = peptide + "-" + score;
        hashCode = toString.hashCode();
        this.isDecoy = isDecoy;
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