package proteomics.Types;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Segment.InferenceSegment;
import proteomics.TheoSeq.MassTool;

import java.util.Map;

public class Peptide implements Comparable<Peptide> {

    private static Logger logger = LoggerFactory.getLogger(Peptide.class);

    private final String ptmFreeSeq;
    private final boolean isDecoy;
    private final String normalizedPeptideString;
    private final MassTool massToolObj;
    private final int maxMs2Charge;
    private final String leftFlank;
    private final String rightFlank;
    private final int globalRank;
    private final double normalizedCrossXcorr;
    private PositionDeltaMassMap varPTMMap = null;

    private String toString;
    private int hashCode;

    // these fields need to be changed every time PTM changed.
    private float precursorMass = -1;
    private float[][] ionMatrix = null;
    private float[] chargeOneBIonArray = null;
    private String varPtmContainingSeq = null;
    private String ptmContainingSeq = null;

    public Peptide(String ptmFreeSeq, boolean isDecoy, MassTool massToolObj, int maxMs2Charge, double normalizedCrossXcorr, String leftFlank, String rightFlank, int globalRank) {
        this.ptmFreeSeq = ptmFreeSeq;
        this.isDecoy = isDecoy;
        this.normalizedPeptideString = InferenceSegment.normalizeSequence(ptmFreeSeq);
        this.normalizedCrossXcorr = normalizedCrossXcorr;
        this.massToolObj = massToolObj;
        this.maxMs2Charge = maxMs2Charge;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
        this.globalRank = globalRank;

        toString = leftFlank + "." + ptmFreeSeq + "." + rightFlank;
        hashCode = toString.hashCode();
    }

    public int getGlobalRank() {
        return globalRank;
    }

    public float[][] getIonMatrix() {
        if (ionMatrix == null) {
            ionMatrix = massToolObj.buildIonArray(getVarPtmContainingSeq(), maxMs2Charge);
        }
        return ionMatrix;
    }

    public String getNormalizedPeptideString() {
        return normalizedPeptideString;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public float getPrecursorMass() {
        if (precursorMass < 0) {
            precursorMass = massToolObj.calResidueMass(getVarPtmContainingSeq()) + MassTool.H2O;
        }
        return precursorMass;
    }

    public float[] getChargeOneBIonArray() {
        if (chargeOneBIonArray == null) {
            float[][] temp = massToolObj.buildIonArray(getVarPtmContainingSeq(), 1);
            chargeOneBIonArray = temp[0];
        }
        return chargeOneBIonArray;
    }

    public String toString() {
        return toString;
    }

    public boolean equals(Object other) {
        if (!(other instanceof Peptide)) {
            return false;
        }

        Peptide otherPeptide = (Peptide) other;
        return this.hashCode == otherPeptide.hashCode;
    }

    public Peptide clone() {
        Peptide other = null;
        try {
            other = new Peptide(ptmFreeSeq, isDecoy, massToolObj, maxMs2Charge, normalizedCrossXcorr, leftFlank, rightFlank, globalRank);
            if (varPTMMap != null) {
                other.setVarPTM(varPTMMap.clone());
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }

        return other;
    }

    public int hashCode() {
        return hashCode;
    }

    public int length() {
        return ptmFreeSeq.length();
    }

    public String getLeftFlank() {
        return leftFlank;
    }

    public String getRightFlank() {
        return rightFlank;
    }

    public void setVarPTM(PositionDeltaMassMap ptmMap) {
        this.varPTMMap = ptmMap;
        if (ptmMap != null) {
            // reset these fields to make them being regenerated again.
            precursorMass = -1;
            ionMatrix = null;
            chargeOneBIonArray = null;
            varPtmContainingSeq = null;
            ptmContainingSeq = null;

            toString = leftFlank + "." + ptmFreeSeq + "." + rightFlank + "." + ptmMap.toString();
            hashCode = toString.hashCode();
        }
    }

    public boolean hasVarPTM() {
        return varPTMMap != null;
    }

    public int getVarPTMNum() {
        if (hasVarPTM()) {
            return varPTMMap.size();
        } else {
            return 0;
        }
    }

    public String getPTMFreeSeq() {
        return ptmFreeSeq;
    }

    public int getUnexplainedAaNum() {
        if (varPTMMap == null) {
            return 0;
        } else {
            int total = 0;
            for (Coordinate co : varPTMMap.keySet()) {
                if (co.y - co.x > 1) {
                    total += co.y - co.x;
                }
            }
            return total;
        }
    }

    public PositionDeltaMassMap getVarPTMs() {
        return varPTMMap;
    }

    public String getVarPtmContainingSeq() {
        if (varPtmContainingSeq == null) {
            if (varPTMMap != null) {
                StringBuilder sb = new StringBuilder(ptmFreeSeq.length() * 5);
                int tempIdx = varPTMMap.firstKey().y;
                if (tempIdx > 1) {
                    sb.append(ptmFreeSeq.substring(0, tempIdx - 1));
                }
                int i = tempIdx - 1;
                tempIdx = varPTMMap.lastKey().y;
                while (i < ptmFreeSeq.length()) {
                    boolean hasMod = false;
                    if (tempIdx > i) {
                        for (Coordinate co : varPTMMap.keySet()) {
                            if (co.y - 1 == i) {
                                sb.append(String.format("%c(%.1f)", ptmFreeSeq.charAt(i), varPTMMap.get(co)));
                                hasMod = true;
                                ++i;
                                break;
                            }
                        }
                        if (!hasMod) {
                            sb.append(ptmFreeSeq.charAt(i));
                            ++i;
                        }
                    } else {
                        break;
                    }
                }
                if (tempIdx < ptmFreeSeq.length()) {
                    sb.append(ptmFreeSeq.substring(tempIdx));
                }
                varPtmContainingSeq = sb.toString();
            } else {
                varPtmContainingSeq = ptmFreeSeq;
            }
        }

        return varPtmContainingSeq;
    }

    public String getPtmContainingSeq(Map<Character, Float> fixModMap) { // caution: containing fix modification. Calculating ion masses based on it is incorrect.
        if (ptmContainingSeq == null) {
            ptmContainingSeq = getVarPtmContainingSeq();
            for (char aa : fixModMap.keySet()) {
                if (Math.abs(fixModMap.get(aa)) > 0.01) {
                    ptmContainingSeq = ptmContainingSeq.replaceAll(String.valueOf(aa), String.format("%c(%.1f)", aa, fixModMap.get(aa)));
                }
            }
        }

        return ptmContainingSeq;
    }

    public double getNormalizedCrossCorr() {
        return normalizedCrossXcorr;
    }

    public int compareTo(Peptide peptide) {
        if (normalizedCrossXcorr > peptide.getNormalizedCrossCorr()) {
            return 1;
        } else if (normalizedCrossXcorr < peptide.getNormalizedCrossCorr()) {
            return -1;
        } else {
            return 0;
        }
    }
}
