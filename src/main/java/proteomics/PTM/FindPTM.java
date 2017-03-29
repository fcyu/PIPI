package proteomics.PTM;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.DynamicProgramming.TagAlignment;
import proteomics.PIPI;
import proteomics.Types.*;
import proteomics.TheoSeq.MassTool;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class FindPTM {

    private static final Logger logger = LoggerFactory.getLogger(FindPTM.class);

    private final float ms1Tolerance;
    private final int ms1ToleranceUnit;
    private final float ms2Tolerance;
    private final float minPtmMass;
    private final float maxPtmMass;
    private final Map<String, Float> massTable;
    private List<Peptide> peptideWithPtmList = new LinkedList<>();
    private float totalResidueMz;
    private final Map<String, TreeSet<Integer>> siteMass100Map;
    private final Set<String> noResultScanPeptide = new HashSet<>();

    public FindPTM(List<Peptide> peptideList, SpectrumEntry spectrumEntry, List<ThreeExpAA> exp3aaLists, MassTool massToolObj, Map<String, TreeSet<Integer>> siteMass100Map, float minPtmMass, float maxPtmMass, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance) {
        this.ms1Tolerance = ms1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.massTable = massToolObj.returnMassTable();
        this.siteMass100Map = siteMass100Map;

        totalResidueMz = spectrumEntry.precursorMass - massTable.get("H2O") + massTable.get("PROTON");
        for (Peptide peptide : peptideList) {
            PositionDeltaMassMap gaps = inferGaps(peptide, spectrumEntry, massTable.get("PROTON"), totalResidueMz, exp3aaLists);
            if ((gaps != null) && (!gaps.isEmpty())) {
                PositionDeltaMassMap ptms = gaps2Ptms(peptide, gaps);
                if (ptms != null) {
                    Peptide peptideWithPtm = peptide.clone();
                    peptideWithPtm.setVarPTM(ptms);
                    peptideWithPtmList.add(peptideWithPtm);
                } else if (PIPI.DEV) {
                    noResultScanPeptide.add(peptide.getPTMFreeSeq());
                }
            } else if (PIPI.DEV) {
                noResultScanPeptide.add(peptide.getPTMFreeSeq());
            }
        }

        if (PIPI.DEV) {
            try (BufferedWriter writer = new BufferedWriter(new FileWriter("no_result_scans" + "." + spectrumEntry.scanNum + "." + spectrumEntry.precursorCharge + ".csv"))) {
                writer.write("total_candidate,no_result_candidate\n");
                for (String seq : noResultScanPeptide) {
                    writer.write(String.format("%d,%s\n", peptideList.size(), seq));
                }
            } catch (IOException ex) {
                ex.printStackTrace();
                logger.error(ex.getMessage());
                System.exit(1);
            }
        }
    }

    public List<Peptide> getPeptidesWithPTMs() {
        return peptideWithPtmList;
    }

    private PositionDeltaMassMap inferGaps(Peptide peptide, SpectrumEntry spectrumEntry, float NTermMz, float CTermMz, List<ThreeExpAA> exp3aaListSet) {
        float expPrecursorMass = spectrumEntry.precursorMass;
        float peptideMass = peptide.getPrecursorMass();
        float deltaMass = expPrecursorMass - peptideMass;
        float[] chargeOneBMzArray = peptide.getChargeOneBIonArray();

        // Get useful ExpAALists and their theoretical locations.
        List<ThreeExpAA> usefulExpAas = new LinkedList<>(cleanLocateTags(exp3aaListSet, peptide.getChargeOneBIonArray(), peptide.getNormalizedPeptideString()));
        if (usefulExpAas.isEmpty()) {
            return null;
        }

        Collections.sort(usefulExpAas);

        // Tag alignment
        TagAlignment tagAlignmentObj = new TagAlignment(usefulExpAas.toArray(new ThreeExpAA[usefulExpAas.size()]), peptide, siteMass100Map, NTermMz, CTermMz, ms2Tolerance);
        List<ExpAA> alignedResult = tagAlignmentObj.getAlignedList();

        if ((alignedResult == null) || (alignedResult.isEmpty())) {
            return null;
        }

        // find gaps
        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap();
        float totalDeltaMass = 0;
        for (ExpAA aa : alignedResult) {
            if (aa.getMod() > 0) {
                totalDeltaMass += aa.getMod();
                positionDeltaMassMap.put(new Coordinate(aa.getTheoLocation(), aa.getTheoLocation() + 1), aa.getMod());
            }
        }

        // check N-term first
        ExpAA firstExpAa = alignedResult.get(0);
        if (Math.abs(firstExpAa.getHeadLocation() - massTable.get("PROTON")) > ms2Tolerance) {
            float expMass = firstExpAa.getHeadLocation() - massTable.get("PROTON");
            float theoMass = chargeOneBMzArray[firstExpAa.getTheoLocation() - 1] - massTable.get("PROTON");
            float massDiff = expMass - theoMass;
            if (Math.abs(massDiff) > ms2Tolerance) {
                positionDeltaMassMap.put(new Coordinate(0, firstExpAa.getTheoLocation()), massDiff);
                totalDeltaMass += massDiff;
            }
        }

        // check C-term
        ExpAA lastExpAa = alignedResult.get(alignedResult.size() - 1);
        if (Math.abs(totalResidueMz -lastExpAa.getTailLocation()) > ms2Tolerance) {
            float expMass =  totalResidueMz - lastExpAa.getTailLocation();
            float theoMass = chargeOneBMzArray[chargeOneBMzArray.length - 1] - chargeOneBMzArray[lastExpAa.getTheoLocation()];
            float massDiff = expMass - theoMass;
            if (Math.abs(massDiff) > ms2Tolerance) {
                positionDeltaMassMap.put(new Coordinate(lastExpAa.getTheoLocation() + 1, peptide.length()), massDiff);
                totalDeltaMass += massDiff;
            }
        }

        for (int i = 1; i < alignedResult.size(); ++i) {
            ExpAA formerExpAa = alignedResult.get(i - 1);
            ExpAA laterExpAa = alignedResult.get(i);
            float expMassGap = laterExpAa.getHeadLocation() - formerExpAa.getTailLocation();
            float theoMass = chargeOneBMzArray[laterExpAa.getTheoLocation() - 1] - chargeOneBMzArray[formerExpAa.getTheoLocation()];
            float massDiff = expMassGap - theoMass;
            if (Math.abs(massDiff) > 2 * ms2Tolerance) {
                // There is a gap in exp and theo spectrum.
                if (Math.abs(massDiff) > ms2Tolerance) {
                    positionDeltaMassMap.put(new Coordinate(formerExpAa.getTheoLocation() + 1, laterExpAa.getTheoLocation()), massDiff);
                    totalDeltaMass += massDiff;
                }
            }
        }

        // check the total delta mass of modification
        float tolerance = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            tolerance = peptideMass * ms1Tolerance / 1e6f;
        }
        tolerance = (float) Math.max(tolerance, 0.1);
        if (Math.abs(deltaMass - totalDeltaMass) > tolerance) {
            logger.debug("Spectrum {}, peptide {}, the precursor mass difference is not balance. There are {} .", spectrumEntry.scanNum, peptide.getPTMFreeSeq(), deltaMass - totalDeltaMass);
        }

        return positionDeltaMassMap;
    }

    private PositionDeltaMassMap gaps2Ptms(Peptide peptide, PositionDeltaMassMap positionGapMap) {
        String peptideString = peptide.getPTMFreeSeq();

        // check if there are complementary modification (|diffMass| are equal)
        // If there are, eliminate both of them. (Occam's razor)
        Coordinate[] coordinateArray = positionGapMap.keySet().toArray(new Coordinate[positionGapMap.size()]);
        Set<Coordinate> delCoordinateSet = new HashSet<>();
        for (int i = 0; i < coordinateArray.length - 1; ++i) {
            for (int j = i + 1; j < coordinateArray.length; ++j) {
                float deltaMass1 = positionGapMap.get(coordinateArray[i]);
                float deltaMass2 = positionGapMap.get(coordinateArray[j]);
                float tolerance = 2 * ms2Tolerance;
                if (Math.abs(deltaMass1 + deltaMass2) <= tolerance) {
                    delCoordinateSet.add(coordinateArray[i]);
                    delCoordinateSet.add(coordinateArray[j]);
                }
            }
        }

        // delete gaps with abs PTM mass smaller than a threshold
        for (Coordinate coordinate : positionGapMap.keySet()) {
            if (Math.abs(positionGapMap.get(coordinate)) < Math.max(ms2Tolerance, 0.1)) {
                delCoordinateSet.add(coordinate);
            }
        }

        for (Coordinate coordinate : delCoordinateSet) {
            positionGapMap.remove(coordinate);
        }

        if (!positionGapMap.isEmpty()) {
            // Infer and validate PTMs
            PositionDeltaMassMap newPositionGapMap = new PositionDeltaMassMap();
            for (Coordinate coordinate : positionGapMap.keySet()) {
                if (coordinate.y - coordinate.x == 1) {
                    newPositionGapMap.put(new Coordinate(coordinate.x, coordinate.y), positionGapMap.get(coordinate));
                } else {
                    String aaSeq = peptideString.substring(coordinate.x, coordinate.y);
                    float deltaMass = positionGapMap.get(coordinate);
                    boolean peptideN = false;
                    boolean peptideC = false;
                    boolean proteinN = false;
                    boolean proteinC = false;
                    if (coordinate.x == 0) {
                        if (peptide.getLeftFlank().contentEquals("-")) {
                            proteinN = true;
                        } else {
                            peptideN = true;
                        }
                    } else if (coordinate.y == peptideString.length()) {
                        if (peptide.getRightFlank().contentEquals("-")) {
                            proteinC = true;
                        } else {
                            peptideC = true;
                        }
                    }

                    Map<Integer, Float> tempMap = pinPointPTM(siteMass100Map, aaSeq, peptideN, peptideC, proteinN, proteinC, deltaMass, 2 * ms2Tolerance);
                    for (int idx : tempMap.keySet()) {
                        newPositionGapMap.put(new Coordinate(coordinate.x + idx, coordinate.x + idx + 1), tempMap.get(idx));
                    }
                }
            }
            return newPositionGapMap;
        } else {
            return null;
        }
    }

    private Set<ThreeExpAA> cleanLocateTags(List<ThreeExpAA> inputSet, float[] peptideBIonArray, String normalizedPeptideString) {
        Set<ThreeExpAA> outputSet = new HashSet<>();
        for (ThreeExpAA expAaList : inputSet) {
            int idx = 0;
            while (idx != -1) {
                idx = normalizedPeptideString.indexOf(expAaList.getPtmFreeAAString(), idx);
                if (idx != -1) {
                    if ((expAaList.get(0).getnTermMod() == 0) || ((idx == 1) && (expAaList.get(2).getcTermMod() == 0))) { // check N-term modification
                        if ((expAaList.get(2).getcTermMod() == 0) || ((idx == normalizedPeptideString.length() - 4) && (expAaList.get(0).getnTermMod() == 0))) { // check C-term modification
                            // check if the segment is in the tolerance boundary.
                            float temp = expAaList.getTailLocation() - peptideBIonArray[idx - 1 + expAaList.size() - 1]; // there is a "n" at the beginning
                            if ((temp <= maxPtmMass) && (temp >= minPtmMass)) {
                                ThreeExpAA usefulAaList = expAaList.clone();
                                for (int i = 0; i < usefulAaList.size(); ++i) {
                                    usefulAaList.setTheoLocation(i, idx + i);
                                }
                                outputSet.add(usefulAaList);
                            }
                        }
                    }
                    ++idx;
                }
            }
        }

    private Map<Integer, Float> pinPointPTM(Map<String, TreeSet<Integer>> siteMass1000Map, String seq, boolean peptideN, boolean peptideC, boolean proteinN, boolean proteinC, float deltaMass, float tolerance) { // todo: check this function is from TagAlignment.checkJumpConstrain
        Map<Integer, int[]> localSiteMass100Map = new HashMap<>();
        for (int i = 0; i < seq.length(); ++i) {
            String aa = seq.substring(i, i + 1);
            TreeSet<Integer> tempSet = new TreeSet<>();
            // Consider the amino acid only.
            if (siteMass1000Map.containsKey(aa)) {
                tempSet.addAll(siteMass1000Map.get(aa));
            }

            // When the amino acid is in the terminal, consider corresponding modification.
            if ((i == 0) && peptideN && siteMass1000Map.containsKey(aa + "-PEPTIDE_N")) {
                tempSet.addAll(siteMass1000Map.get(aa + "-PEPTIDE_N"));
            }
            if ((i == 0) && proteinN && siteMass1000Map.containsKey(aa + "-PROTEIN_N")) {
                tempSet.addAll(siteMass1000Map.get(aa + "-PROTEIN_N"));
            }
            if ((i == seq.length() - 1) && peptideC && siteMass1000Map.containsKey(aa + "-PEPTIDE_C")) {
                tempSet.addAll(siteMass1000Map.get(aa + "-PEPTIDE_C"));
            }
            if ((i == seq.length() - 1) && proteinC && siteMass1000Map.containsKey(aa + "-PROTEIN_C")) {
                tempSet.addAll(siteMass1000Map.get(aa + "-PROTEIN_C"));
            }

            // When in the terminal, consider the terminal modification.
            if ((i == 0) && peptideN && siteMass1000Map.containsKey("PEPTIDE_N")) {
                tempSet.addAll(siteMass1000Map.get("PEPTIDE_N"));
            }
            if ((i == 0) && proteinN && siteMass1000Map.containsKey("PROTEIN_N")) {
                tempSet.addAll(siteMass1000Map.get("PROTEIN_N"));
            }
            if ((i == seq.length() - 1) && peptideC && siteMass1000Map.containsKey("PEPTIDE_C")) {
                tempSet.addAll(siteMass1000Map.get("PEPTIDE_C"));
            }
            if ((i == seq.length() - 1) && proteinC && siteMass1000Map.containsKey("PEPTIDE_C")) {
                tempSet.addAll(siteMass1000Map.get("PEPTIDE_C"));
            }

            if (!tempSet.isEmpty()) {
                int[] tempArray = new int[tempSet.size()];
                int tempIdx = 0;
                for (int mass1000 : tempSet) {
                    tempArray[tempIdx] = mass1000;
                    ++tempIdx;
                }
                localSiteMass100Map.put(i, tempArray);
            }
        }

        int leftBoard = (int) ((deltaMass - tolerance) * 1000);
        int rightBoard = (int) ((deltaMass + tolerance) * 1000);

        int[] idxArray = new int[localSiteMass100Map.size()];
        int tempIdx = 0;
        for (int idx : localSiteMass100Map.keySet()) {
            idxArray[tempIdx] = idx;
            ++tempIdx;
        }
        Arrays.sort(idxArray);

        boolean ok = false;
        Map<Integer, Float> modifiedAaMap = new HashMap<>();

        // one PTM
        for (int idx : idxArray) {
            int[] mass100Array = localSiteMass100Map.get(idx);
            for (int mass1000 : mass100Array) {
                if ((mass1000 >= leftBoard) && (mass1000 <= rightBoard)) {
                    modifiedAaMap.put(idx, mass1000 * 0.001f);
                    ok = true;
                    break;
                }
            }
            if (ok) {
                break;
            }
        }

        // two PTM
        if ((!ok) && (idxArray.length > 1)) {
            for (int i = 0; i < idxArray.length - 1; ++i) {
                int[] mass100Array1 = localSiteMass100Map.get(idxArray[i]);
                for (int j = i + 1; j < idxArray.length; ++j) {
                    int[] mass100Array2 = localSiteMass100Map.get(idxArray[j]);
                    for (int mass100_1 : mass100Array1) {
                        for (int mass100_2 : mass100Array2) {
                            int mass100T = mass100_1 + mass100_2;
                            if ((mass100T >= leftBoard) && (mass100T <= rightBoard)) {
                                modifiedAaMap.put(idxArray[i], mass100_1 * 0.001f);
                                modifiedAaMap.put(idxArray[j], mass100_2 * 0.001f);
                                ok = true;
                                break;
                            }
                        }
                        if (ok) {
                            break;
                        }
                    }
                    if (ok) {
                        break;
                    }
                }
                if (ok) {
                    break;
                }
            }
        }

        if (!ok) {
            logger.debug("There is something wrong in FindPTM (line: 345).");
        }
        return modifiedAaMap;
    }
}
