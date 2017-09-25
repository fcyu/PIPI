package proteomics.Types;


import java.util.*;

public class FinalResultEntry {

    private static final int scoreNum = 5;

    private final int scanNum;
    private final int charge;
    private final float precursorMz;
    private final String mgtTitle;
    private TreeSet<PeptideScore> peptideScoreList = new TreeSet<>(Collections.reverseOrder());
    private double ionFrac;
    private double matchedHighestIntensityFrac;
    private double ptmDeltasScore;
    private TreeSet<PeptideScore> ptmPatterns;
    private int explainedAaNum;

    private float qValue = -1;

    public FinalResultEntry(int scanNum, int charge, float precursorMz, String mgtTitle) {
        this.scanNum = scanNum;
        this.charge = charge;
        this.precursorMz = precursorMz;
        this.mgtTitle = mgtTitle;
    }

    public String getMgtTitle() {
        return mgtTitle;
    }

    public double getScore() {
        return peptideScoreList.first().score;
    }

    public boolean hasHit() {
        return !peptideScoreList.isEmpty();
    }

    public Peptide getPeptide() {
        return peptideScoreList.first().peptide;
    }
    
    public boolean isDecoy() {
        return peptideScoreList.first().peptide.isDecoy();
    }

    public double getDeltaC() {
        if (peptideScoreList.size() < 2) {
            return 1;
        } else {
            Iterator<PeptideScore> it = peptideScoreList.iterator();
            return (it.next().score - it.next().score) / peptideScoreList.first().score;
        }
    }

    public double getDeltaLC() {
        if (peptideScoreList.size() < 2) {
            return 1;
        } else {
            return (peptideScoreList.first().score - peptideScoreList.last().score) / peptideScoreList.first().score;
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

    public void addScore(PeptideScore peptideScore) {
        if (peptideScoreList.size() < scoreNum) {
            peptideScoreList.add(peptideScore);
        } else if (peptideScore.score > peptideScoreList.last().score) {
            peptideScoreList.pollLast();
            peptideScoreList.add(peptideScore);
        }
    }

    public void setQValue(float qValue) {
        this.qValue = qValue;
    }

    public float getQValue() {
        return qValue;
    }

    public double getNormalizedCrossCorrelationCoefficient() {
        return peptideScoreList.first().peptide.getNormalizedCrossCorr();
    }

    public int getGlobalSearchRank() {
        return peptideScoreList.first().peptide.getGlobalRank();
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

    public void setPtmPatterns(TreeSet<PeptideScore> ptmPatterns) {
        this.ptmPatterns = ptmPatterns;
    }

    public TreeSet<PeptideScore> getPtmPatterns() {
        return ptmPatterns;
    }

    public void setExplainedAaNum(int explainedAaNum) {
        this.explainedAaNum = explainedAaNum;
    }

    public int getExplainedAaNum() {
        return explainedAaNum;
    }
}
