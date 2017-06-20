package proteomics.Types;


import java.util.*;

public class FinalResultEntry {

    private static final int scoreNum = 5;

    private final int scanNum;
    private final int charge;
    private final float precursorMz;
    private LinkedList<Double> scoreList = new LinkedList<>();
    private Peptide peptide;
    private double normalizedCrossXcorr;
    private int globalSearchRank;
    private double ionFrac;
    private double matchedHighestIntensityFrac;
    private double ptmDeltasScore;
    private LinkedList<PeptideScore> ptmPatterns;

    private float qValue = -1;

    public FinalResultEntry(int scanNum, int charge, float precursorMz) {
        this.scanNum = scanNum;
        this.charge = charge;
        this.precursorMz = precursorMz;
    }

    public double getScore() {
        return scoreList.peekFirst();
    }

    public boolean noScore() {
        return scoreList.isEmpty();
    }

    public Peptide getPeptide() {
        return peptide;
    }
    
    public boolean isDecoy() {
        return peptide.isDecoy();
    }

    public double getDeltaC() {
        if (scoreList.size() < 2) {
            return 1;
        } else {
            return (scoreList.peekFirst() - scoreList.get(1)) / scoreList.peekFirst();
        }
    }

    public double getDeltaLC() {
        if (scoreList.size() < 2) {
            return 1;
        } else {
            return (scoreList.peekFirst() - scoreList.peekLast()) / scoreList.peekFirst();
        }
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

    public void addScore(double score) {
        if (scoreList.size() < scoreNum) {
            scoreList.add(score);
            scoreList.sort(Collections.reverseOrder());
        } else if (score > scoreList.peekLast()) {
            scoreList.pollLast();
            scoreList.add(score);
            scoreList.sort(Collections.reverseOrder());
        }
    }

    public void setPeptide(Peptide peptide) {
        this.peptide = peptide;
    }

    public void setQValue(float qValue) {
        this.qValue = qValue;
    }

    public float getQValue() {
        return qValue;
    }

    public void setNormalizedCrossXcorr(double normalizedCrossXcorr) {
        this.normalizedCrossXcorr = normalizedCrossXcorr;
    }

    public void setGlobalSearchRank(int rank) {
        globalSearchRank = rank;
    }

    public double getNormalizedCrossXcorr() {
        return normalizedCrossXcorr;
    }

    public int getGlobalSearchRank() {
        return globalSearchRank;
    }

    public void setIonFrac(double ionFrac) {
        this.ionFrac = ionFrac;
    }

    public void setMatchedHighestIntensityFrac(double matchedHighestIntensityFrac) {
        this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
    }

    public double getIonFrac() {
        return ionFrac;
    }

    public double getMatchedHighestIntensityFrac() {
        return matchedHighestIntensityFrac;
    }

    public void setPtmDeltasScore(double ptmDeltasScore) {
        this.ptmDeltasScore = ptmDeltasScore;
    }

    public double getPtmDeltasScore() {
        return ptmDeltasScore;
    }

    public void setPtmPatterns(LinkedList<PeptideScore> ptmPatterns) {
        this.ptmPatterns = ptmPatterns;
    }

    public LinkedList<PeptideScore> getPtmPatterns() {
        return ptmPatterns;
    }
}
