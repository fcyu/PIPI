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

    private float qValue = -1;

    public FinalResultEntry(int scanNum, int charge, float precursorMz) {
        this.scanNum = scanNum;
        this.charge = charge;
        this.precursorMz = precursorMz;
    }

    public double getScore() {
        return scoreList.get(0);
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
            return (scoreList.get(0) - scoreList.get(1)) / scoreList.get(0);
        }
    }

    public double getDeltaLC() {
        if (scoreList.size() < 3) {
            return 1;
        } else {
            return (scoreList.get(0) - scoreList.get(scoreList.size() - 1)) / scoreList.get(0);
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
}
