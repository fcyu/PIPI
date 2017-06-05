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
    private float precursorMass = 0;
    private final float[][] ionMatrix;
    private final float[] chargeOneBIonArray;
    private PositionDeltaMassMap varPTMMap = null;
    private final MassTool massToolObj;
    private final int maxMs2Charge;
    private final String leftFlank;
    private final String rightFlank;
    private final int globalRank;
    private String toString;
    private int hashCode;

    private final double normalizedCrossXcorr;

    public Peptide(String ptmFreeSeq, boolean isDecoy, MassTool massToolObj, int maxMs2Charge, double normalizedCrossXcorr, String leftFlank, String rightFlank, int globalRank) {
        this.ptmFreeSeq = ptmFreeSeq;
        this.isDecoy = isDecoy;
        this.normalizedPeptideString = InferenceSegment.normalizeSequence(ptmFreeSeq);
        this.normalizedCrossXcorr = normalizedCrossXcorr;
        precursorMass = massToolObj.calResidueMass(ptmFreeSeq) + massToolObj.returnMassTable().get("H2O");
        ionMatrix = massToolObj.buildIonArray(ptmFreeSeq, maxMs2Charge);
        chargeOneBIonArray = ionMatrix[0];

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
        return ionMatrix;
    }

    public String getNormalizedPeptideString() {
        return normalizedPeptideString;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public float getPrecursorMass() {
        return precursorMass;
    }

    public float[] getChargeOneBIonArray() {
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
            for (Coordinate co : ptmMap.keySet()) {
                float deltaMass = ptmMap.get(co);
                for (int i = 0; i < ionMatrix.length / 2; ++i) {
                    float deltaMz = deltaMass / (i + 1);
                    // if the PTM can be pin-pointed to a single amino acid, change the ions' mz values as what they should be
                    // the PTM cannot be pin-pointed to a single amino acid, all peaks in the block except for head or tail are treated as unmodified peaks
                    if (co.x == 0) {
                        for (int j = Math.max(0, co.y - 2); j < ionMatrix[0].length; ++j) {
                            ionMatrix[i * 2][j] += deltaMz;
                        }
                        ionMatrix[i * 2 + 1][0] += deltaMz;
                    } else if (co.y > ptmFreeSeq.length() - 2) {
                        ionMatrix[i * 2][ionMatrix[0].length - 1] += deltaMz;
                        for (int j = 0; j < co.x; ++j) {
                            ionMatrix[i * 2 + 1][j] += deltaMz;
                        }
                    } else {
                        for (int j = co.y - 2; j < ionMatrix[0].length; ++j) {
                            ionMatrix[i * 2][j] += deltaMz;
                        }
                        for (int j = 0; j < co.x; ++j) {
                            ionMatrix[i * 2 + 1][j] += deltaMz;
                        }
                    }
                }
            }
            float totalDeltaMass = 0;
            for (float deltaMass : ptmMap.values()) {
                totalDeltaMass += deltaMass;
            }
            precursorMass += totalDeltaMass;

            toString += "." + varPTMMap.toString();
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

    public String getPTMContainedString(Map<String, Float> fixModMap, int decimalPoint) { // include fix modification
        if (hasVarPTM()) {
            StringBuilder sb = new StringBuilder(ptmFreeSeq.length() * 5);
            int i = 0;
            while (i < ptmFreeSeq.length()) {
                boolean ok = false;

                // add variable modification, and maybe fix modification
                for (Coordinate co : varPTMMap.keySet()) {
                    if (co.x == i) {
                        // sb.append("(");
                        while (i < co.y) {
                            sb.append(ptmFreeSeq.charAt(i));
                            ++i;
                        }
                        if (decimalPoint == 0) {
                            sb.append(String.format("(%d)", Math.round(varPTMMap.get(co) + fixModMap.get(String.valueOf(ptmFreeSeq.charAt(i - 1)))))); // add fix modification to variable modification.
                        } else if (decimalPoint == 1) {
                            sb.append(String.format("(%.1f)", varPTMMap.get(co) + fixModMap.get(String.valueOf(ptmFreeSeq.charAt(i - 1))))); // add fix modification to variable modification.
                        } else {
                            sb.append(String.format("(%.2f)", varPTMMap.get(co) + fixModMap.get(String.valueOf(ptmFreeSeq.charAt(i - 1))))); // add fix modification to variable modification.
                        }
                        ok = true;
                        break;
                    }
                }

                // add fix modification or not
                if (!ok) {
                    float deltaMass = fixModMap.get(String.valueOf(ptmFreeSeq.charAt(i)));
                    if (Math.abs(deltaMass) > 1e-6) {
                        sb.append(ptmFreeSeq.charAt(i));
                        sb.append(String.format("(%.2f)", deltaMass)); // for fix modification, the decimal point is always 2
                    } else {
                        sb.append(ptmFreeSeq.charAt(i));
                    }
                    ++i;
                }
            }
            return sb.toString();
        } else {
            StringBuilder sb = new StringBuilder(ptmFreeSeq.length() * 5);
            int i = 0;
            while (i < ptmFreeSeq.length()) {
                // add fix modification
                float deltaMass = fixModMap.get(String.valueOf(ptmFreeSeq.charAt(i)));
                if (Math.abs(deltaMass) > 1e-6) {
                    sb.append(ptmFreeSeq.charAt(i));
                    sb.append(String.format("(%.2f)", deltaMass));
                } else {
                    sb.append(ptmFreeSeq.charAt(i));
                }
                ++i;
            }
            return sb.toString();
        }
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
