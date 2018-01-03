package proteomics.TheoSeq;

import proteomics.Types.AA;
import proteomics.Types.SparseBooleanVector;

import java.util.*;
import java.util.regex.*;

public class MassTool {

    private static final Pattern mod_aa_pattern = Pattern.compile("([A-Znc])(\\(([0-9.\\-]+)\\))?");
    public static final float PROTON = 1.00727646688f;
    public static final float C13_DIFF = 1.00335483f;

    public final float H2O;
    public Map<String, Double> elementTable = new HashMap<>();

    private final Map<Character, Float> massTable = new HashMap<>(30, 1);
    private final int missedCleavage;
    private float ms2Tolerance;
    private float oneMinusBinOffset;
    private Pattern digestSitePattern;
    private final boolean cleavageFromCTerm;
    private final String labelling;
    private final Map<Character, Float> fixModMap;

    public MassTool(final int missedCleavage, Map<Character, Float> fixModMap, String cleavageSite, String protectionSite, boolean cleavageFromCTerm, float ms2Tolerance, float oneMinusBinOffset, String labelling) {
        this.labelling = labelling;
        this.fixModMap = fixModMap;

        elementTable.put("-", 0d);
        elementTable.put("H", 1.0078246);
        elementTable.put("He", 3.01603);
        elementTable.put("Li", 6.015121);
        elementTable.put("Be", 9.012182);
        elementTable.put("B", 10.012937);
        elementTable.put("C", 12.0000000);
        elementTable.put("N", 14.0030732);
        if (labelling.contentEquals("N15")) {
            elementTable.put("N", 15.0001088);
        }
        elementTable.put("O", 15.9949141);
        elementTable.put("F", 18.9984032);
        elementTable.put("Ne", 19.992435);
        elementTable.put("Na", 22.989767);
        elementTable.put("Mg", 23.985042);
        elementTable.put("Al", 26.981539);
        elementTable.put("Si", 27.976927);
        elementTable.put("P", 30.973762);
        elementTable.put("S", 31.972070);
        elementTable.put("Cl", 34.9688531);
        elementTable.put("Ar", 35.967545);
        elementTable.put("K", 38.963707);
        elementTable.put("Ca", 39.962591);
        elementTable.put("Sc", 44.955910);
        elementTable.put("Ti", 45.952629);
        elementTable.put("V", 49.947161);
        elementTable.put("Cr", 49.946046);
        elementTable.put("Mn", 54.938047);
        elementTable.put("Fe", 53.939612);
        elementTable.put("Co", 58.933198);
        elementTable.put("Ni", 57.935346);
        elementTable.put("Cu", 62.939598);
        elementTable.put("Zn", 63.929145);
        elementTable.put("Ga", 68.925580);
        elementTable.put("Ge", 69.924250);
        elementTable.put("As", 74.921594);
        elementTable.put("Se", 73.922475);
        elementTable.put("Br", 78.918336);
        elementTable.put("Kr", 77.914);
        elementTable.put("Rb", 84.911794);
        elementTable.put("Sr", 83.913430);
        elementTable.put("Y", 88.905849);
        elementTable.put("Zr", 89.904703);
        elementTable.put("Nb", 92.906377);
        elementTable.put("Mo", 91.906808);
        elementTable.put("Tc", 98.0);
        elementTable.put("Ru", 95.907599);
        elementTable.put("Rh", 102.905500);
        elementTable.put("Pd", 101.905634);
        elementTable.put("Ag", 106.905092);
        elementTable.put("Cd", 105.906461);
        elementTable.put("In", 112.904061);
        elementTable.put("Sn", 111.904826);
        elementTable.put("Sb", 120.903821);
        elementTable.put("Te", 119.904048);
        elementTable.put("I", 126.904473);
        elementTable.put("Xe", 123.905894);
        elementTable.put("Cs", 132.905429);
        elementTable.put("Ba", 129.906282);
        elementTable.put("La", 137.90711);
        elementTable.put("Ce", 135.907140);
        elementTable.put("Pr", 140.907647);
        elementTable.put("Nd", 141.907719);
        elementTable.put("Pm", 145.0);
        elementTable.put("Sm", 143.911998);
        elementTable.put("Eu", 150.919847);
        elementTable.put("Gd", 151.919786);
        elementTable.put("Tb", 158.925342);
        elementTable.put("Dy", 155.925277);
        elementTable.put("Ho", 164.930319);
        elementTable.put("Er", 161.928775);
        elementTable.put("Tm", 168.934212);
        elementTable.put("Yb", 167.933894);
        elementTable.put("Lu", 174.940770);
        elementTable.put("Hf", 173.940044);
        elementTable.put("Ta", 179.947462);
        elementTable.put("W", 179.946701);
        elementTable.put("Re", 184.952951);
        elementTable.put("Os", 183.952488);
        elementTable.put("Ir", 190.960584);
        elementTable.put("Pt", 189.959917);
        elementTable.put("Au", 196.966543);
        elementTable.put("Hg", 195.965807);
        elementTable.put("Tl", 202.972320);
        elementTable.put("Pb", 203.973020);
        elementTable.put("Bi", 208.980374);
        elementTable.put("Po", 209.0);
        elementTable.put("At", 210.0);
        elementTable.put("Rn", 222.0);
        elementTable.put("Fr", 223.0);
        elementTable.put("Ra", 226.025);
        // elementTable.put("Ac", 227.028); // conflict with Unimod bricks
        elementTable.put("Th", 232.038054);
        elementTable.put("Pa", 231.0359);
        elementTable.put("U", 234.040946);
        elementTable.put("Np", 237.048);
        elementTable.put("Pu", 244.0);
        elementTable.put("Am", 243.0);
        elementTable.put("Cm", 247.0);
        elementTable.put("Bk", 247.0);
        elementTable.put("Cf", 251.0);
        elementTable.put("Es", 252.0);
        elementTable.put("Fm", 257.0);
        elementTable.put("Md", 258.0);
        elementTable.put("No", 259.0);
        elementTable.put("Lr", 260.0);
        elementTable.put("13C", 13.0033554);
        elementTable.put("15N", 15.0001088);
        elementTable.put("18O", 17.9991616);
        elementTable.put("2H", 2.0141021);
        elementTable.put("dHex", elementTable.get("C") * 6 + elementTable.get("O") * 4 + elementTable.get("H") * 10);
        elementTable.put("Hep", elementTable.get("C") * 7 + elementTable.get("O") * 6 + elementTable.get("H") * 12);
        elementTable.put("Hex", elementTable.get("C") * 6 + elementTable.get("O") * 5 + elementTable.get("H") * 10);
        elementTable.put("HexA", elementTable.get("C") * 6 + elementTable.get("O") * 6 + elementTable.get("H") * 8);
        elementTable.put("HexN", elementTable.get("C") * 6 + elementTable.get("O") * 4 + elementTable.get("H") * 11 + elementTable.get("N"));
        elementTable.put("HexNAc", elementTable.get("C") * 8 + elementTable.get("O") * 5 + + elementTable.get("N") + elementTable.get("H") * 13);
        elementTable.put("Kdn", elementTable.get("C") * 9 + elementTable.get("H") * 14 + elementTable.get("O") * 8);
        elementTable.put("Kdo", elementTable.get("C") * 8 + elementTable.get("H") * 12 + elementTable.get("O") * 7);
        elementTable.put("NeuAc", elementTable.get("C") * 11 + elementTable.get("H") * 17 + elementTable.get("O") * 8 + elementTable.get("N"));
        elementTable.put("NeuGc", elementTable.get("C") * 11 + elementTable.get("H") * 17 + elementTable.get("O") * 9 + elementTable.get("N"));
        elementTable.put("Pent", elementTable.get("C") * 5 + elementTable.get("O") * 4 + elementTable.get("H") * 8);
        elementTable.put("Phos", elementTable.get("O") * 3 + elementTable.get("H") + elementTable.get("P"));
        elementTable.put("Sulf", elementTable.get("S") + elementTable.get("O") * 3);
        elementTable.put("Water", elementTable.get("H") * 2 + elementTable.get("O"));
        elementTable.put("Me", elementTable.get("C") + elementTable.get("H") * 2);
        elementTable.put("Ac", elementTable.get("C") * 2 + elementTable.get("H") * 2 + elementTable.get("O")); // Caution! This is not Actinium

        this.missedCleavage = missedCleavage;
        this.ms2Tolerance = ms2Tolerance;
        this.oneMinusBinOffset = oneMinusBinOffset;
        this.cleavageFromCTerm = cleavageFromCTerm;
        massTable.put('G', (float) (elementTable.get("C") * 2 + elementTable.get("H") * 3 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('G')));
        massTable.put('A', (float) (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('A')));
        massTable.put('S', (float) (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('S')));
        massTable.put('P', (float) (elementTable.get("C") * 5 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('P')));
        massTable.put('V', (float) (elementTable.get("C") * 5 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('V')));
        massTable.put('T', (float) (elementTable.get("C") * 4 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('I')));
        massTable.put('C', (float) (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") + elementTable.get("S") + fixModMap.get('C')));
        massTable.put('I', (float) (elementTable.get("C") * 6 + elementTable.get("H") * 11 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('I')));
        massTable.put('L', (float) (elementTable.get("C") * 6 + elementTable.get("H") * 11 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('L')));
        massTable.put('N', (float) (elementTable.get("C") * 4 + elementTable.get("H") * 6 + elementTable.get("N") * 2 + elementTable.get("O") * 2 + fixModMap.get('N')));
        massTable.put('D', (float) (elementTable.get("C") * 4 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") * 3 + fixModMap.get('D')));
        massTable.put('Q', (float) (elementTable.get("C") * 5 + elementTable.get("H") * 8 + elementTable.get("N") * 2 + elementTable.get("O") * 2 + fixModMap.get('Q')));
        massTable.put('K', (float) (elementTable.get("C") * 6 + elementTable.get("H") * 12 + elementTable.get("N") * 2 + elementTable.get("O") + fixModMap.get('K')));
        massTable.put('E', (float) (elementTable.get("C") * 5 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 3 + fixModMap.get('E')));
        massTable.put('M', (float) (elementTable.get("C") * 5 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + elementTable.get("S") + fixModMap.get('M')));
        massTable.put('H', (float) (elementTable.get("C") * 6 + elementTable.get("H") * 7 + elementTable.get("N") * 3 + elementTable.get("O") + fixModMap.get('H')));
        massTable.put('F', (float) (elementTable.get("C") * 9 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('F')));
        massTable.put('R', (float) (elementTable.get("C") * 6 + elementTable.get("H") * 12 + elementTable.get("N") * 4 + elementTable.get("O") + fixModMap.get('R')));
        massTable.put('Y', (float) (elementTable.get("C") * 9 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('Y')));
        massTable.put('W', (float) (elementTable.get("C") * 11 + elementTable.get("H") * 10 + elementTable.get("N") * 2 + elementTable.get("O") + fixModMap.get('W')));
        massTable.put('U', (float) (elementTable.get("C") * 3 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 2 + elementTable.get("Se") + fixModMap.get('U')));
        massTable.put('O', (float) (elementTable.get("C") * 12 + elementTable.get("H") * 21 + elementTable.get("N") * 3 + elementTable.get("O") * 3 + fixModMap.get('O')));
        massTable.put('n', fixModMap.get('n'));
        massTable.put('c', fixModMap.get('c'));
        massTable.put('#', massTable.get('I')); // for I and L.
        massTable.put('$', (massTable.get('Q') + massTable.get('K')) / 2); // for Q and K.
        H2O = (float) (elementTable.get("H") * 2 + elementTable.get("O"));

        if (protectionSite.contentEquals("-")) {
            digestSitePattern = Pattern.compile("[" + cleavageSite + "]");
        } else if (cleavageFromCTerm) {
            digestSitePattern = Pattern.compile("[" + cleavageSite + "](?=[^" + protectionSite + "])");
        } else {
            digestSitePattern = Pattern.compile("(?<=[^" + protectionSite + "])" + "[" + cleavageSite + "]");
        }
    }

    public float calResidueMass(String seq) { // n and c are also AA. Consider fixed modification automatically
        float total_mass = 0;
        Matcher matcher = mod_aa_pattern.matcher(seq);
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            float delta_mass = 0;
            if (matcher.group(3) != null) {
                delta_mass = Float.valueOf(matcher.group(3));
            }
            total_mass += massTable.get(aa) + delta_mass;
        }

        return total_mass;
    }

    public float calResidueMass2(String seq) { // n and c are also AA. Don't consider fixed modification
        float total_mass = 0;
        Matcher matcher = mod_aa_pattern.matcher(seq);
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            float delta_mass = 0;
            if (matcher.group(3) != null) {
                delta_mass = Float.valueOf(matcher.group(3));
            }
            if (Math.abs(fixModMap.get(aa) - delta_mass) < 0.01) {
                total_mass += massTable.get(aa);
            } else {
                total_mass += massTable.get(aa) + delta_mass;
            }
        }

        return total_mass;
    }

    public static AA[] seqToAAList(String seq) { // n and c are also AA.
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

        float[][] peptideIonArray = new float[2 * maxCharge][aaArray.length - 2];
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

        Set<Integer> tempSet = new HashSet<>();
        for (int i = 0; i < rowNum; ++i) {
            for (int j = 0; j < colNum; ++j) {
                if (ionMatrix[i][j] > 1e-6) {
                    int idx = mzToBin(ionMatrix[i][j]);
                    tempSet.add(idx);
                }
            }
        }

        return new SparseBooleanVector(tempSet);
    }

    public Map<Character, Float> returnMassTable() {
        return massTable;
    }

    public Map<String, Double> getElementTable() {
        return elementTable;
    }

    public int mzToBin(float mz) {
        return (int) Math.floor(mz / (2 * ms2Tolerance) + oneMinusBinOffset);
    }

    public float binToMz(int idx) {
        return (idx - oneMinusBinOffset) * 2 * ms2Tolerance;
    }

    public String getLabelling() {
        return labelling;
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
