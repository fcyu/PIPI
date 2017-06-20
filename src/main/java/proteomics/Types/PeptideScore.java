package proteomics.Types;


public class PeptideScore implements Comparable<PeptideScore> {

    public final double score;
    public final Peptide peptide;
    private final String toString;
    private final int hashCode;

    public PeptideScore(double score, Peptide peptide) {
        this.score = score;
        this.peptide = peptide;
        toString = peptide.toString() + "-" + score;
        hashCode = toString.hashCode();
    }

    public String toString() {
        return toString;
    }

    public int hashCode() {
        return hashCode;
    }

    public int compareTo(PeptideScore other) {
        if (score > other.score) {
            return 1;
        } else if (score < other.score) {
            return -1;
        } else {
            if (peptide.getVarPTMNum() < other.peptide.getVarPTMNum()) {
                return 1;
            } else if (peptide.getVarPTMNum() > other.peptide.getVarPTMNum()) {
                return -1;
            } else {
                return 0;
            }
        }
    }

    public boolean equals(Object other) {
        if (other instanceof PeptideScore) {
            PeptideScore temp = (PeptideScore) other;
            return peptide.equals(temp.peptide) && (score == temp.score);
        } else {
            return false;
        }
    }
}
