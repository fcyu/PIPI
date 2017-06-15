package proteomics.TheoSeq;

import proteomics.Types.AA;
import proteomics.Types.SparseBooleanVector;

import java.util.*;
import java.util.regex.*;

public class MassTool {

    private static final Pattern mod_aa_pattern = Pattern.compile("([A-Znc])(\\(([0-9\\.\\-]+)\\))?");

    public static final float C13_DIFF = 1.00335483f;
    public static final float H2O = 18.010564684f;
    public static final float PROTON = 1.00727646688f;
    public static final float NH3 = 17.026549106f;

    private final Map<Character, Float> massTable = new HashMap<>(30, 1);
    private final int missedCleavage;
    private float ms2Tolerance = 1.0005f;
    private float oneMinusBinOffset = 0.4f;
    private Pattern digestSitePattern;
    private final boolean cleavageFromCTerm;

    public MassTool(final int missedCleavage, Map<Character, Float> fixModMap, String cleavageSite, String protectionSite, boolean cleavageFromCTerm, float ms2Tolerance, float oneMinusBinOffset) {
        this.missedCleavage = missedCleavage;
        this.ms2Tolerance = ms2Tolerance;
        this.oneMinusBinOffset = oneMinusBinOffset;
        this.cleavageFromCTerm = cleavageFromCTerm;
        massTable.put('G', 57.021464f + fixModMap.get('G'));
        massTable.put('A', 71.037114f + fixModMap.get('A'));
        massTable.put('S', 87.032028f + fixModMap.get('S'));
        massTable.put('P', 97.052764f + fixModMap.get('P'));
        massTable.put('V', 99.068414f + fixModMap.get('V'));
        massTable.put('T', 101.047678f + fixModMap.get('I'));
        massTable.put('C', 103.009184f + fixModMap.get('C'));
        massTable.put('I', 113.084064f + fixModMap.get('I'));
        massTable.put('L', 113.084064f + fixModMap.get('L'));
        massTable.put('N', 114.042927f + fixModMap.get('N'));
        massTable.put('D', 115.026943f + fixModMap.get('D'));
        massTable.put('Q', 128.058578f + fixModMap.get('Q'));
        massTable.put('K', 128.094963f + fixModMap.get('K'));
        massTable.put('E', 129.042593f + fixModMap.get('E'));
        massTable.put('M', 131.040485f + fixModMap.get('M'));
        massTable.put('H', 137.058912f + fixModMap.get('H'));
        massTable.put('F', 147.068414f + fixModMap.get('F'));
        massTable.put('R', 156.101111f + fixModMap.get('R'));
        massTable.put('Y', 163.063329f + fixModMap.get('Y'));
        massTable.put('W', 186.079313f + fixModMap.get('W'));
        massTable.put('U', 150.953636f + fixModMap.get('U'));
        massTable.put('O', 132.08988f + fixModMap.get('O'));
        massTable.put('n', fixModMap.get('n'));
        massTable.put('c', fixModMap.get('c'));
        massTable.put('#', 113.084064f); // for I and L.
        massTable.put('$', 128.0767705f); // for Q and K.

        if (protectionSite.contentEquals("-")) {
            digestSitePattern = Pattern.compile("[" + cleavageSite + "]");
        } else if (cleavageFromCTerm) {
            digestSitePattern = Pattern.compile("[" + cleavageSite + "](?=[^" + protectionSite + "])");
        } else {
            digestSitePattern = Pattern.compile("(?<=[^" + protectionSite + "])" + "[" + cleavageSite + "]");
        }
    }

    public float calResidueMass(String seq) { // n and c are also AA.
        double total_mass = 0;
        Matcher matcher = mod_aa_pattern.matcher(seq);
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            float delta_mass = 0;
            if (matcher.group(3) != null) {
                delta_mass = Float.valueOf(matcher.group(3));
            }
            total_mass += massTable.get(aa) + delta_mass;
        }

        return (float) total_mass;
    }

    public AA[] seqToAAList(String seq) {
        Matcher matcher = mod_aa_pattern.matcher(seq);
        List<AA> temp = new LinkedList<>();
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            float delta_mass = 0;
            if (matcher.group(3) != null) {
                delta_mass = Float.valueOf(matcher.group(3));
            }
            temp.add(new AA(aa, delta_mass));
        }
        return temp.toArray(new AA[temp.size()]);
    }

    public Set<String> buildPeptideSet(String proSeq) {
        Map<Integer, List<int[]>> digestRangeMap = digestTrypsin(proSeq);
        Set<String> peptideSeqSet = new HashSet<>();

        for (int i = 0; i <= missedCleavage; ++i) {
            for (int[] digestRange1 : digestRangeMap.get(i)) {
                String subString = proSeq.substring(digestRange1[0], digestRange1[1]);
                peptideSeqSet.add("n" + subString + "c");
            }
        }
        return peptideSeqSet;
    }

    public float[][] buildIonArray(String seq, int maxCharge) {
        AA[] aaArray = seqToAAList(seq);

        float[][] peptideIonArray = new float[2 * maxCharge][seq.length() - 2];
        // traverse the sequence to get b-ion
        float bIonMass = massTable.get(aaArray[0].aa) + aaArray[0].ptmDeltaMass; // add N-term modification
        for (int i = 1; i < aaArray.length - 2; ++i) {
            bIonMass += massTable.get(aaArray[i].aa) + aaArray[i].ptmDeltaMass;
            for (int charge = 1; charge <= maxCharge; ++charge) {
                peptideIonArray[2 * (charge - 1)][i - 1]  = bIonMass / charge + 1.00727646688f;
            }
        }
        // calculate the last b-ion with C-term modification
        bIonMass +=  massTable.get(aaArray[aaArray.length - 2].aa) + aaArray[aaArray.length - 2].ptmDeltaMass + massTable.get(aaArray[aaArray.length - 1].aa) + aaArray[aaArray.length - 1].ptmDeltaMass;
        for (int charge = 1; charge <= maxCharge; ++charge) {
            peptideIonArray[2 * (charge - 1)][aaArray.length - 3] = bIonMass / charge + 1.00727646688f;
        }

        // traverse the sequence with reversed order to get y-ion
        // the whole sequence
        float yIonMass = bIonMass + H2O;
        for (int charge = 1; charge <= maxCharge; ++charge) {
            peptideIonArray[2 * (charge - 1) + 1][0] = yIonMass  / charge + 1.00727646688f;
        }
        // delete the first amino acid and N-term modification
        yIonMass -= massTable.get(aaArray[0].aa) + aaArray[0].ptmDeltaMass + massTable.get(aaArray[1].aa) + aaArray[1].ptmDeltaMass;
        for (int charge = 1; charge <= maxCharge; ++charge) {
            peptideIonArray[2 * (charge - 1) + 1][1] = yIonMass / charge + 1.00727646688f;
        }

        // rest of the sequence
        for (int i = 2; i < aaArray.length - 2; ++i) {
            yIonMass -= massTable.get(aaArray[i].aa) + aaArray[i].ptmDeltaMass;
            for (int charge = 1; charge <= maxCharge; ++charge) {
                peptideIonArray[2 * (charge - 1) + 1][i] = yIonMass / charge + 1.00727646688f;
            }
        }

        return peptideIonArray;
    }

    public SparseBooleanVector buildVector(float[][] ionMatrix, int precursorCharge) {
        int colNum = ionMatrix[0].length;
        int rowNum = Math.min(ionMatrix.length / 2, precursorCharge - 1) * 2;
        if (precursorCharge == 1) {
            rowNum = 2;
        }

        SparseBooleanVector theoIonVector = new SparseBooleanVector();
        for (int i = 0; i < rowNum; ++i) {
            for (int j = 0; j < colNum; ++j) {
                if (ionMatrix[i][j] > 1e-6) {
                    int idx = mzToBin(ionMatrix[i][j]);
                    theoIonVector.put(idx);
                }
            }
        }

        return theoIonVector;
    }

    public Map<Character, Float> returnMassTable() {
        return massTable;
    }

    public int mzToBin(float mz) {
        return (int) Math.floor(mz / (2 * ms2Tolerance) + oneMinusBinOffset);
    }

    public float binToMz(int idx) {
        return (idx - oneMinusBinOffset) * 2 * ms2Tolerance;
    }

    private Map<Integer, List<int[]>> digestTrypsin(String proSeq) {
        // Cut a protein
        List<Integer> cutPointList = new ArrayList<>(200);
        int length = proSeq.length();
        int idxStart = 0;
        Matcher matchObj = digestSitePattern.matcher(proSeq);
        cutPointList.add(0);
        while (idxStart < length) {
            if (matchObj.find()) {
                int cutPoint;
                if (cleavageFromCTerm) {
                    cutPoint = matchObj.end();
                } else {
                    cutPoint = matchObj.start();
                }
                cutPointList.add(cutPoint);
                idxStart = cutPoint;
            } else {
                cutPointList.add(length);
                break;
            }
        }

        Collections.sort(cutPointList);

        // Deal with missed cleavage
        Map<Integer, List<int[]>> digestRangeMap = new HashMap<>();
        for (int time = 0; time <= missedCleavage; ++time) {
            List<int[]> temp = new LinkedList<>();
            int leftPoint;
            int rightPoint;
            for (int i = 0; i + 1 + time < cutPointList.size(); ++i) {
                leftPoint = cutPointList.get(i);
                rightPoint = cutPointList.get(i + 1 + time);
                temp.add(new int[]{leftPoint, rightPoint});
            }
            digestRangeMap.put(time, temp);
        }

        return digestRangeMap;
    }
}
