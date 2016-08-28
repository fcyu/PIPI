package proteomics.Types;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Segment.InferenceSegment;
import proteomics.TheoSeq.MassTool;

import java.util.Map;

public class Peptide implements Comparable<Peptide> {

    private static Logger logger = LoggerFactory.getLogger(Peptide.class);

    private final String peptideString;
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

    private final double normalizedCrossXcorr;

    public Peptide(String peptideString, boolean isDecoy, MassTool massToolObj, int maxMs2Charge, double normalizedCrossXcorr, String leftFlank, String rightFlank, int globalRank) {
        this.peptideString = peptideString;
        this.isDecoy = isDecoy;
        this.normalizedPeptideString = InferenceSegment.normalizeSequence(false, peptideString);
        this.normalizedCrossXcorr = normalizedCrossXcorr;
        precursorMass = massToolObj.calResidueMass(peptideString) + massToolObj.returnMassTable().get("H2O");
        ionMatrix = massToolObj.buildIonArray(peptideString, maxMs2Charge);
        chargeOneBIonArray = ionMatrix[0];

        this.massToolObj = massToolObj;
        this.maxMs2Charge = maxMs2Charge;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
        this.globalRank = globalRank;
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
        if (hasVarPTM()) {
            return leftFlank + "." + peptideString + "." + rightFlank + "," + varPTMMap.toString();
        } else {
            return leftFlank + "." + peptideString + "." + rightFlank;
        }
    }

    public boolean equals(Object other) {
        if (!(other instanceof Peptide)) {
            return false;
        }

        Peptide otherPeptide = (Peptide) other;
        return this.hashCode() == otherPeptide.hashCode();
    }

    public Peptide clone() {
        Peptide other = null;
        try {
            other = new Peptide(peptideString, isDecoy, massToolObj, maxMs2Charge, normalizedCrossXcorr, leftFlank, rightFlank, globalRank);
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
        return this.toString().hashCode();
    }

    public int length() {
        return peptideString.length();
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
                    // the PTM cannot be pin-pointed to a single amino acid, all peaks in the block except for head or tail are treated as missing peaks
                    for (int j = co.x; j < co.y - 1; ++j) {
                        ionMatrix[i * 2][j] = 0;
                    }
                    for (int j = co.y - 1; j < ionMatrix[0].length; ++j) {
                        ionMatrix[i * 2][j] += deltaMz;
                    }
                    for (int j = co.x + 1; j < co.y; ++j) {
                        ionMatrix[i * 2 + 1][j] = 0;
                    }
                    for (int j = 0; j <= co.x; ++j) {
                        ionMatrix[i * 2 + 1][j] += deltaMz;
                    }
                }
            }
            float totalDeltaMass = 0;
            for (float deltaMass : ptmMap.values()) {
                totalDeltaMass += deltaMass;
            }
            precursorMass += totalDeltaMass;
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
        return peptideString;
    }

    public PositionDeltaMassMap getVarPTMs() {
        return varPTMMap;
    }

    public String getPTMContainedString(Map<String, Float> fixModMap) { // include fix modification
        if (hasVarPTM()) {
            StringBuilder sb = new StringBuilder(peptideString.length() * 5);
            int i = 0;
            while (i < peptideString.length()) {
                boolean ok = false;

                // add variable modification, and maybe fix modification
                for (Coordinate co : varPTMMap.keySet()) {
                    if (co.x == i) {
                        // sb.append("(");
                        while (i < co.y) {
                            sb.append(peptideString.charAt(i));
                            ++i;
                        }
                        sb.append(String.format("(%.2f)", varPTMMap.get(co) + fixModMap.get(String.valueOf(peptideString.charAt(i - 1))))); // add fix modification to variable modification.
                        ok = true;
                        break;
                    }
                }

                // add fix modification or not
                if (!ok) {
                    float deltaMass = fixModMap.get(String.valueOf(peptideString.charAt(i)));
                    if (Math.abs(deltaMass) > 1e-6) {
                        sb.append(peptideString.charAt(i));
                        sb.append(String.format("(%.2f)", deltaMass));
                    } else {
                        sb.append(peptideString.charAt(i));
                    }
                    ++i;
                }
            }
            return sb.toString();
        } else {
            StringBuilder sb = new StringBuilder(peptideString.length() * 5);
            int i = 0;
            while (i < peptideString.length()) {
                // add fix modification
                float deltaMass = fixModMap.get(String.valueOf(peptideString.charAt(i)));
                if (Math.abs(deltaMass) > 1e-6) {
                    sb.append(peptideString.charAt(i));
                    sb.append(String.format("(%.2f)", deltaMass));
                } else {
                    sb.append(peptideString.charAt(i));
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
