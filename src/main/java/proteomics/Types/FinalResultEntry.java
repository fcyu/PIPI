package proteomics.Types;


import java.util.*;

public class FinalResultEntry {

    private final int scanNum;
    private final int charge;
    private final float precursorMz;
    private final String mgtTitle;
    private final String labeling; // N14, N15;

    private TreeSet<Peptide> peptideSet = new TreeSet<>(Collections.reverseOrder());
    private Map<String, TreeSet<Peptide>> ptmPatterns;

    public FinalResultEntry(int scanNum, int charge, float precursorMz, String mgtTitle, String labeling) {
        this.scanNum = scanNum;
        this.charge = charge;
        this.precursorMz = precursorMz;
        this.mgtTitle = mgtTitle;
        this.labeling = labeling;
    }

    public String getMgtTitle() {
        return mgtTitle;
    }

    public boolean hasHit() {
        return !peptideSet.isEmpty();
    }

    public int getScanNum() {
        return scanNum;
    }

    public int getCharge() {
        return charge;
    }

    public float getPrecursorMz() {
        return precursorMz;
    }

    public String getLabeling() {
        return labeling;
    }

    public void addScore(Peptide peptide) {
        if (peptideSet.size() < 5) {
            peptideSet.add(peptide);
        } else if (peptide.getScore() > peptideSet.last().getScore()) {
            peptideSet.pollLast();
            peptideSet.add(peptide);
        }
    }

    public TreeSet<Peptide> getPeptideSet() {
        return peptideSet;
    }

    public void setPtmPatterns(Map<String, TreeSet<Peptide>> ptmPatterns) {
        this.ptmPatterns = ptmPatterns;
    }

    public Map<String, TreeSet<Peptide>> getPtmPatterns() {
        return ptmPatterns;
    }
}
