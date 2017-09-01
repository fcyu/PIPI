package proteomics.Types;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Segment.InferenceSegment;
import proteomics.TheoSeq.MassTool;

import java.util.Locale;
import java.util.Map;

public class Peptide implements Comparable<Peptide> {

    private static Logger logger = LoggerFactory.getLogger(Peptide.class);

    private final String ptmFreeSeq;
    private final boolean isDecoy;
    private final String normalizedPeptideString;
    private final MassTool massToolObj;
    private final int maxMs2Charge;
    private final char leftFlank;
    private final char rightFlank;
    private final int globalRank;
    private final double normalizedCrossCorrelationCoefficient;

    private String toString;
    private int hashCode;

    private int unknownPtmNum = 0;

    // these fields need to be changed every time PTM changed.
    private PositionDeltaMassMap varPTMMap = null;
    private float precursorMass = -1;
    private float[][] ionMatrix = null;
    private float[] chargeOneBIonArray = null;
    private String varPtmContainingSeq = null;

    private String ptmContainingSeq = null;

    public Peptide(String ptmFreeSeq, boolean isDecoy, MassTool massToolObj, int maxMs2Charge, double normalizedCrossCorrelationCoefficient, char leftFlank, char rightFlank, int globalRank) {
        this.ptmFreeSeq = ptmFreeSeq;
        this.isDecoy = isDecoy;
        this.normalizedPeptideString = InferenceSegment.normalizeSequence(ptmFreeSeq);
        this.normalizedCrossCorrelationCoefficient = normalizedCrossCorrelationCoefficient;
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
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massToolObj.buildIonArray(varPtmContainingSeq, maxMs2Charge);
            precursorMass = massToolObj.calResidueMass(varPtmContainingSeq);
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

    public float getPrecursorMass() {
        if (precursorMass < 0) {
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massToolObj.buildIonArray(varPtmContainingSeq, maxMs2Charge);
            precursorMass = massToolObj.calResidueMass(varPtmContainingSeq);
            chargeOneBIonArray = ionMatrix[0];
        }
        return precursorMass;
    }

    public float[] getChargeOneBIonArray() {
        if (chargeOneBIonArray == null) {
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massToolObj.buildIonArray(varPtmContainingSeq, maxMs2Charge);
            precursorMass = massToolObj.calResidueMass(varPtmContainingSeq);
            chargeOneBIonArray = ionMatrix[0];
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
            other = new Peptide(ptmFreeSeq, isDecoy, massToolObj, maxMs2Charge, normalizedCrossCorrelationCoefficient, leftFlank, rightFlank, globalRank);
            if (varPTMMap != null) {
                other.setVarPTM(varPTMMap.clone());
                other.setUnknownPtmNum(unknownPtmNum);
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

    public char getLeftFlank() {
        return leftFlank;
    }

    public char getRightFlank() {
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

    public void setUnknownPtmNum(int unknownPtmNum) {
        this.unknownPtmNum = unknownPtmNum;
    }

    public int getUnknownPtmNum() {
        return unknownPtmNum;
    }

    public String getPTMFreeSeq() {
        return ptmFreeSeq;
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
                                sb.append(String.format(Locale.US, "%c(%.1f)", ptmFreeSeq.charAt(i), varPTMMap.get(co)));
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
                    ptmContainingSeq = ptmContainingSeq.replaceAll(String.valueOf(aa), String.format(Locale.US, "%c(%.1f)", aa, fixModMap.get(aa)));
                }
            }
        }

        return ptmContainingSeq;
    }

    public double getNormalizedCrossCorr() {
        return normalizedCrossCorrelationCoefficient;
    }

    public int compareTo(Peptide peptide) {
        if (normalizedCrossCorrelationCoefficient > peptide.getNormalizedCrossCorr()) {
            return 1;
        } else if (normalizedCrossCorrelationCoefficient < peptide.getNormalizedCrossCorr()) {
            return -1;
        } else {
            return 0;
        }
    }
}
