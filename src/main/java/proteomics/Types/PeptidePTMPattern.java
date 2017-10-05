package proteomics.Types;

public class PeptidePTMPattern {

    public final String ptmFreeSequence;

    private Entry topEntry = null; // May be we will maintain a list to hold top N candidates in the future.

    public PeptidePTMPattern(String ptmFreeSequence) {
        this.ptmFreeSequence = ptmFreeSequence;
    }

    public void update(Peptide peptide, int K1) {
        Entry entry = new Entry(K1, peptide);
        if (topEntry == null) {
            topEntry = entry;
        } else if (topEntry.compareTo(entry) < 0) {
            topEntry = entry;
        }
    }

    public Entry getTopEntry() {
        return topEntry;
    }

    public class Entry implements Comparable<Entry> {

        public final int matchedPeakNum;
        public final Peptide peptide ;

        private final int hashCode;

        Entry(int matchedPeakNum, Peptide peptide) {
            this.matchedPeakNum = matchedPeakNum;
            this.peptide = peptide;

            hashCode = peptide.hashCode();
        }

        public int hashCode() {
            return hashCode;
        }

        public int compareTo(Entry other) {
            if (peptide.getScore() > other.peptide.getScore()) {
                return 1;
            } else if (peptide.getScore() < other.peptide.getScore()) {
                return -1;
            } else {
                if (this.matchedPeakNum > other.matchedPeakNum) {
                    return 1;
                } else if (this.matchedPeakNum < other.matchedPeakNum) {
                    return -1;
                } else {
                    return 0;
                }
            }
        }

        public boolean equals(Object other) {
            if (other instanceof Entry) {
                return hashCode == ((Entry) other).hashCode;
            } else {
                return false;
            }
        }
    }
}
