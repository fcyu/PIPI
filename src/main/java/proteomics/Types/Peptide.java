package proteomics.Types;

import ProteomicsLibrary.Types.Coordinate;
import proteomics.Segment.InferSegment;
import ProteomicsLibrary.MassTool;

import java.util.Locale;
import java.util.Map;

public class Peptide implements Comparable<Peptide> {

    private final String ptmFreePeptide;
    private final boolean isDecoy;
    private final String normalizedPeptideString;
    private final MassTool massTool;
    private final int maxMs2Charge;
    private final int globalRank;
    private final double normalizedCrossCorrelationCoefficient;

    private int hashCode;

    // these fields need to be changed every time PTM changed.
    private PositionDeltaMassMap varPTMMap = null;
    private double theoMass = -1;
    private double[][] ionMatrix = null;
    private double[] chargeOneBIonArray = null;
    private String varPtmContainingSeq = null;

    private String ptmContainingSeq = null;

    // score part
    private double score = -1;
    private int matchedPeakNum = -1;
    private double ionFrac = -1;
    private double matchedHighestIntensityFrac = -1;
    private double explainedAaFrac = -1;
    private double qValue = -1;
    private String aScore = "-";

    public Peptide(String ptmFreePeptide, boolean isDecoy, MassTool massTool, int maxMs2Charge, double normalizedCrossCorrelationCoefficient, int globalRank) {
        this.ptmFreePeptide = ptmFreePeptide;
        this.isDecoy = isDecoy;
        this.normalizedPeptideString = InferSegment.normalizeSequence(ptmFreePeptide);
        this.normalizedCrossCorrelationCoefficient = normalizedCrossCorrelationCoefficient;
        this.massTool = massTool;
        this.maxMs2Charge = maxMs2Charge;
        this.globalRank = globalRank;

        hashCode = ptmFreePeptide.hashCode();
    }

    public int getGlobalRank() {
        return globalRank;
    }

    public double[][] getIonMatrix() {
        if (ionMatrix == null) {
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massTool.buildIonArray(varPtmContainingSeq, maxMs2Charge);
            theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
            chargeOneBIonArray = ionMatrix[0];
        }
        return ionMatrix;
    }

    public String getNormalizedPeptideString() {
        return normalizedPeptideString;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public double getTheoMass() {
        if (theoMass < 0) {
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massTool.buildIonArray(varPtmContainingSeq, maxMs2Charge);
            theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
            chargeOneBIonArray = ionMatrix[0];
        }
        return theoMass;
    }

    public double[] getChargeOneBIonArray() {
        if (chargeOneBIonArray == null) {
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massTool.buildIonArray(varPtmContainingSeq, maxMs2Charge);
            theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
            chargeOneBIonArray = ionMatrix[0];
        }
        return chargeOneBIonArray;
    }

    public boolean equals(Object other) {
        if (!(other instanceof Peptide)) {
            return false;
        }

        Peptide otherPeptide = (Peptide) other;
        return this.hashCode == otherPeptide.hashCode;
    }

    public Peptide clone() throws CloneNotSupportedException {
        super.clone();
        Peptide other = new Peptide(ptmFreePeptide, isDecoy, massTool, maxMs2Charge, normalizedCrossCorrelationCoefficient, globalRank);
        if (varPTMMap != null) {
            other.setVarPTM(varPTMMap.clone());
            other.setScore(score);
            other.setMatchedHighestIntensityFrac(matchedHighestIntensityFrac);
            other.setExplainedAaFrac(explainedAaFrac);
            other.setIonFrac(ionFrac);
            other.setaScore(aScore);
            other.setQValue(qValue);
        }

        return other;
    }

    public int hashCode() {
        return hashCode;
    }

    public int length() {
        return ptmFreePeptide.length();
    }

    public void setVarPTM(PositionDeltaMassMap ptmMap) {
        this.varPTMMap = ptmMap;
        if (ptmMap != null) {
            // reset these fields to make them being regenerated again.
            theoMass = -1;
            ionMatrix = null;
            chargeOneBIonArray = null;
            varPtmContainingSeq = null;
            ptmContainingSeq = null;

            String toString = ptmFreePeptide + "." + ptmMap.toString();
            hashCode = toString.hashCode();
        }
    }

    public String toString() {
        return ptmFreePeptide + "." + varPTMMap.toString();
    }

    public boolean hasVarPTM() {
        return varPTMMap != null;
    }

    private int getVarPTMNum() {
        if (hasVarPTM()) {
            return varPTMMap.size();
        } else {
            return 0;
        }
    }

    public String getPTMFreePeptide() {
        return ptmFreePeptide;
    }

    public PositionDeltaMassMap getVarPTMs() {
        return varPTMMap;
    }

    private String getVarPtmContainingSeq() {
        if (varPtmContainingSeq == null) {
            if (varPTMMap != null) {
                StringBuilder sb = new StringBuilder(ptmFreePeptide.length() * 5);
                int tempIdx = varPTMMap.firstKey().y;
                if (tempIdx > 1) {
                    sb.append(ptmFreePeptide.substring(0, tempIdx - 1));
                }
                int i = tempIdx - 1;
                tempIdx = varPTMMap.lastKey().y;
                while (i < ptmFreePeptide.length()) {
                    boolean hasMod = false;
                    if (tempIdx > i) {
                        for (Coordinate co : varPTMMap.keySet()) {
                            if (co.y - 1 == i) {
                                sb.append(String.format(Locale.US, "%c(%.3f)", ptmFreePeptide.charAt(i), varPTMMap.get(co)));
                                hasMod = true;
                                ++i;
                                break;
                            }
                        }
                        if (!hasMod) {
                            sb.append(ptmFreePeptide.charAt(i));
                            ++i;
                        }
                    } else {
                        break;
                    }
                }
                if (tempIdx < ptmFreePeptide.length()) {
                    sb.append(ptmFreePeptide.substring(tempIdx));
                }
                varPtmContainingSeq = sb.toString();
            } else {
                varPtmContainingSeq = ptmFreePeptide;
            }
        }

        return varPtmContainingSeq;
    }

    public String getPtmContainingSeq(Map<Character, Double> fixModMap) { // caution: containing fix modification. Calculating ion masses based on it is incorrect.
        if (ptmContainingSeq == null) {
            ptmContainingSeq = getVarPtmContainingSeq();
            for (char aa : fixModMap.keySet()) {
                if (Math.abs(fixModMap.get(aa)) > 0.01) {
                    ptmContainingSeq = ptmContainingSeq.replaceAll(String.valueOf(aa), String.format(Locale.US, "%c(%.3f)", aa, fixModMap.get(aa)));
                }
            }
        }

        return ptmContainingSeq;
    }

    public double getNormalizedCrossCorr() {
        return normalizedCrossCorrelationCoefficient;
    }

    public void setScore(double score) {
        this.score = score;
    }

    public void setMatchedPeakNum(int matchedPeakNum) {
        this.matchedPeakNum = matchedPeakNum;
    }

    public void setIonFrac(double ionFrac) {
        this.ionFrac = ionFrac;
    }

    public void setMatchedHighestIntensityFrac(double matchedHighestIntensityFrac) {
        this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
    }

    public void setExplainedAaFrac(double explainedAaFrac) {
        this.explainedAaFrac = explainedAaFrac;
    }

    public void setaScore(String aScore) {
        this.aScore = aScore;
    }

    private void setQValue(double qValue) {
        this.qValue = qValue;
    }

    public double getScore() {
        return score;
    }

    public int getMatchedPeakNum() {
        return matchedPeakNum;
    }

    public double getIonFrac() {
        return ionFrac;
    }

    public double getMatchedHighestIntensityFrac() {
        return matchedHighestIntensityFrac;
    }

    public double getExplainedAaFrac() {
        return explainedAaFrac;
    }

    public String getaScore() {
        return aScore;
    }

    public double getQValue() {
        return qValue;
    }

    public int compareTo(Peptide peptide) {
        if (score > peptide.getScore()) {
            return 1;
        } else if (score < peptide.getScore()) {
            return -1;
        } else {
            if (matchedPeakNum > peptide.getMatchedPeakNum()) {
                return 1;
            } else if (matchedPeakNum < peptide.getMatchedPeakNum()) {
                return -1;
            } else {
                if (explainedAaFrac > peptide.getExplainedAaFrac()) {
                    return 1;
                } else if (explainedAaFrac < peptide.getExplainedAaFrac()) {
                    return -1;
                } else {
                    if (getVarPTMNum() < peptide.getVarPTMNum()) {
                        return 1;
                    } else if (getVarPTMNum() > peptide.getVarPTMNum()) {
                        return -1;
                    } else if (normalizedCrossCorrelationCoefficient > peptide.getNormalizedCrossCorr()) {
                        return 1;
                    } else if (normalizedCrossCorrelationCoefficient < peptide.getNormalizedCrossCorr()) {
                        return -1;
                    } else {
                        if (!isDecoy && peptide.isDecoy()) {
                            return 1;
                        } else if (isDecoy && !peptide.isDecoy()) {
                            return -1;
                        } else{
                            return 0;
                        }
                    }
                }
            }
        }
    }
}
