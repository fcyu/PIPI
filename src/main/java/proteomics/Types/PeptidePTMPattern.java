package proteomics.Types;

import proteomics.PIPI;

import java.util.TreeMap;

public class PeptidePTMPattern {

    public final String ptmFreeSequence;

    private final static int topNum = 3;

    private double pValueLowerBound;
    private TreeMap<Integer, Entry> candidateNumEntryMap = new TreeMap<>();
    private Entry topEntry = null;

    private int candidateNum = 0;
    private boolean stopped = false;

    public PeptidePTMPattern(String ptmFreeSequence, double pValueLowerBound) {
        this.ptmFreeSequence = ptmFreeSequence;
        this.pValueLowerBound = pValueLowerBound;
    }

    public void update(double pValue, Peptide peptide, int K1) {
        Entry entry = new Entry(candidateNum, pValueLowerBound, pValue, K1, peptide);
        if (topEntry == null) {
            topEntry = entry;
        } else if (topEntry.compareTo(entry) < 0) {
            topEntry = entry;
        }

        if (PIPI.DEBUG) {
            candidateNumEntryMap.put(candidateNum, entry);
        }
    }

    public void addCandidateNum(int candidateNum) {
        this.candidateNum += candidateNum;
    }

    public boolean reachLowerBound() {
        return topEntry.eValue <= pValueLowerBound * candidateNum;
    }

    public void setPValueLowerBound(double pValueLowerBound) {
        this.pValueLowerBound = pValueLowerBound;
    }

    public double getPValueLowerBound() {
        return pValueLowerBound;
    }

    public TreeMap<Integer, Entry> getCandidateNumEntryMap() {
        return candidateNumEntryMap;
    }

    public Entry getTopEntry() {
        return topEntry;
    }

    public void setStopped() {
        stopped = true;
    }

    public boolean isStopped() {
        return stopped;
    }

    public int getCandidateNum() {
        return candidateNum;
    }

    public class Entry implements Comparable<Entry> {

        public final int candidateNum;
        public final double pValueLowerBound;
        public final double pValue;
        public final int matchedPeakNum;
        public final Peptide peptide ;
        public final double eValue;
        public final double eValueLowerBound;

        private final int hashCode;

        Entry(int candidateNum, double pValueLowerBound, double pValue, int matchedPeakNum, Peptide peptide) {
            this.candidateNum = candidateNum;
            this.pValueLowerBound = pValueLowerBound;
            this.pValue = pValue;
            this.matchedPeakNum = matchedPeakNum;
            this.peptide = peptide;
            eValue = candidateNum * pValue;
            eValueLowerBound = candidateNum * pValueLowerBound;

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
                    if (this.candidateNum < other.candidateNum) {
                        return 1;
                    } else if (this.candidateNum > other.candidateNum) {
                        return -1;
                    } else {
                        return 0;
                    }
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
