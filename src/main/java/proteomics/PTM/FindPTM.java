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
    private List<Peptide> peptideWithPtmList = new LinkedList<>();
    private float totalResidueMz;
    private final Map<String, Float> modifiedAAMassMap;
    private final float[] pepNTermPossibleMod;
    private final float[] pepCTermPossibleMod;
    private final float[] proNTermPossibleMod;
    private final float[] proCTermPossibleMod;
    private final Set<String> noResultScanPeptide = new HashSet<>();

    public FindPTM(List<Peptide> peptideList, SpectrumEntry spectrumEntry, List<ThreeExpAA> exp3aaLists, Map<String, Float> modifiedAAMassMap, float[] pepNTermPossibleMod, float[] pepCTermPossibleMod, float[] proNTermPossibleMod, float[] proCTermPossibleMod, Set<VarModParam> varModParamSet, float minPtmMass, float maxPtmMass, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance) {
        this.ms1Tolerance = ms1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.modifiedAAMassMap = modifiedAAMassMap;
        this.pepNTermPossibleMod = pepNTermPossibleMod;
        this.pepCTermPossibleMod = pepCTermPossibleMod;
        this.proNTermPossibleMod = proNTermPossibleMod;
        this.proCTermPossibleMod = proCTermPossibleMod;

        totalResidueMz = spectrumEntry.precursorMass - MassTool.H2O + MassTool.PROTON;
        for (Peptide peptide : peptideList) {
            PositionDeltaMassMap gaps = inferGaps(peptide, spectrumEntry, MassTool.PROTON, totalResidueMz, exp3aaLists);
            if ((gaps != null) && (!gaps.isEmpty())) {
                PositionDeltaMassMap ptms = refineGaps(peptide, gaps);
                if (ptms != null) {
                    Peptide peptideWithPtm = peptide.clone();
                    peptideWithPtm.setVarPTM(ptms);

                    int unknownPtmNum = 0;
                    for (float ptmMass : ptms.values()) {
                        if (GeneratePtmCandidatesCalculateScore.isNewPtmMass(varModParamSet, ptmMass, Math.max(ms2Tolerance, 0.1f))) {
                            ++unknownPtmNum;
                        }
                    }
                    peptideWithPtm.setUnknownPtmNum(unknownPtmNum);

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
                    writer.write(String.format(Locale.US, "%d,%s\n", peptideList.size(), seq));
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
        TagAlignment tagAlignmentObj = new TagAlignment(usefulExpAas.toArray(new ThreeExpAA[usefulExpAas.size()]), peptide, NTermMz, CTermMz, ms2Tolerance);
        List<ExpAA> alignedResult = tagAlignmentObj.getAlignedList();

        if ((alignedResult == null) || (alignedResult.isEmpty())) {
            return null;
        }

        // find gaps
        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(peptide.getPTMFreeSeq().length());
        float totalDeltaMass = 0;
        for (ExpAA aa : alignedResult) {
            if (aa.getMod() > 0) {
                totalDeltaMass += aa.getMod();
                positionDeltaMassMap.put(new Coordinate(aa.getTheoLocation(), aa.getTheoLocation() + 1), aa.getMod());
            }
            if (aa.getnTermMod() > 0) {
                totalDeltaMass += aa.getnTermMod();
                positionDeltaMassMap.put(new Coordinate(0, 1), aa.getnTermMod());
            }
            if (aa.getcTermMod() > 0) {
                totalDeltaMass += aa.getcTermMod();
                positionDeltaMassMap.put(new Coordinate(aa.getTheoLocation(), aa.getTheoLocation() + 1), aa.getcTermMod());
            }
        }

        // check N-term first
        ExpAA firstExpAa = alignedResult.get(0);
        if (Math.abs(firstExpAa.getHeadLocation() - MassTool.PROTON) > ms2Tolerance) {
            float expMass = firstExpAa.getHeadLocation() - MassTool.PROTON;
            float theoMass = chargeOneBMzArray[firstExpAa.getTheoLocation() - 2] - MassTool.PROTON;
            float massDiff = expMass - theoMass;
            if (Math.abs(massDiff) > ms2Tolerance) {
                positionDeltaMassMap.put(new Coordinate(0, firstExpAa.getTheoLocation()), massDiff);
                totalDeltaMass += massDiff;
            }
        }

        // check C-term
        ExpAA lastExpAa = alignedResult.get(alignedResult.size() - 1);
        if (Math.abs(totalResidueMz - lastExpAa.getTailLocation()) > ms2Tolerance) {
            float expMass =  totalResidueMz - lastExpAa.getTailLocation();
            float theoMass = chargeOneBMzArray[chargeOneBMzArray.length - 1] - chargeOneBMzArray[lastExpAa.getTheoLocation() - 1];
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
            float theoMass = chargeOneBMzArray[laterExpAa.getTheoLocation() - 2] - chargeOneBMzArray[formerExpAa.getTheoLocation() - 1];
            float massDiff = expMassGap - theoMass;
            if (Math.abs(massDiff) > 2 * ms2Tolerance) {
                // There is a gap in exp and theo spectrum.
                positionDeltaMassMap.put(new Coordinate(formerExpAa.getTheoLocation() + 1, laterExpAa.getTheoLocation()), massDiff);
                totalDeltaMass += massDiff;
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

    private PositionDeltaMassMap refineGaps(Peptide peptide, PositionDeltaMassMap positionGapMap) {
        String peptideString = peptide.getPTMFreeSeq();

        // check if there are complementary modification (|diffMass| are equal)
        // If there are, eliminate both of them. (Occam's razor)
        Coordinate[] coordinateArray = positionGapMap.keySet().toArray(new Coordinate[positionGapMap.size()]);
        Set<Coordinate> delCoordinateSet = new HashSet<>();
        for (int i = 0; i < coordinateArray.length - 1; ++i) {
            for (int j = i + 1; j < coordinateArray.length; ++j) {
                float deltaMass1 = positionGapMap.get(coordinateArray[i]);
                float deltaMass2 = positionGapMap.get(coordinateArray[j]);
                if (Math.abs(deltaMass1 + deltaMass2) <= 2 * ms2Tolerance) {
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
            PositionDeltaMassMap newPositionGapMap = new PositionDeltaMassMap(positionGapMap.peptideLength);
            for (Coordinate coordinate : positionGapMap.keySet()) {
                if (coordinate.y - coordinate.x == 1) {
                    newPositionGapMap.put(new Coordinate(coordinate.x, coordinate.y), positionGapMap.get(coordinate));
                } else {
                    String aaSeq = peptideString.substring(coordinate.x, coordinate.y);
                    float deltaMass = positionGapMap.get(coordinate);
                    char[] aaSeqCharArray = aaSeq.toCharArray();
                    boolean ok = false;
                    // try protein-terminal modification first because it's before/important than the peptide-terminal modification in most experiments.
                    // relax the tolerance in high resolution mode
                    // try one known PTM
                    for (int i = 0; i < aaSeqCharArray.length; ++i) {
                        if (coordinate.x + i == 0) {
                            if ((peptide.getLeftFlank() == '-') && (proNTermPossibleMod != null)) {
                                for (float nTermMod : proNTermPossibleMod) {
                                    if (Math.abs(nTermMod - deltaMass) <= Math.min(1.5, 4 * ms2Tolerance)) {
                                        newPositionGapMap.put(new Coordinate(coordinate.x + i, coordinate.x + i + 1), nTermMod);
                                        ok = true;
                                        break;
                                    }
                                }
                            }
                            if (!ok && (pepNTermPossibleMod != null)) {
                                for (float nTermMod : pepNTermPossibleMod) {
                                    if (Math.abs(nTermMod - deltaMass) <= Math.min(1.5, 4 * ms2Tolerance)) {
                                        newPositionGapMap.put(new Coordinate(coordinate.x + i, coordinate.x + i + 1), nTermMod);
                                        ok = true;
                                        break;
                                    }
                                }
                            }
                        } else if (coordinate.x + i == peptideString.length() - 1) {
                            if ((peptide.getRightFlank() == '-') && (proCTermPossibleMod != null)) {
                                for (float cTermMod : proCTermPossibleMod) {
                                    if (Math.abs(cTermMod - deltaMass) <= Math.min(1.5, 4 * ms2Tolerance)) {
                                        newPositionGapMap.put(new Coordinate(coordinate.x + i, coordinate.x + i + 1), cTermMod);
                                        ok = true;
                                        break;
                                    }
                                }
                            }
                            if (!ok && (pepCTermPossibleMod != null)) {
                                for (float cTermMod : pepCTermPossibleMod) {
                                    if (Math.abs(cTermMod - deltaMass) <= Math.min(1.5, 4 * ms2Tolerance)) {
                                        newPositionGapMap.put(new Coordinate(coordinate.x + i, coordinate.x + i + 1), cTermMod);
                                        ok = true;
                                        break;
                                    }
                                }
                            }
                        } else {
                            for (String modifiedAA : modifiedAAMassMap.keySet()) {
                                if ((modifiedAA.charAt(0) == aaSeqCharArray[i]) && (Math.abs(modifiedAAMassMap.get(modifiedAA) - deltaMass) <= Math.min(1.5, 4 * ms2Tolerance))) {
                                    newPositionGapMap.put(new Coordinate(coordinate.x + i, coordinate.x + i + 1), modifiedAAMassMap.get(modifiedAA));
                                    ok = true;
                                    break;
                                }
                            }
                        }
                        if (ok) {
                            break;
                        }
                    }

                    if (!ok) {
                        // try two known PTMs
                        for (int i = 0; i < aaSeqCharArray.length - 1; ++i) {
                            for (int j = i + 1; j < aaSeqCharArray.length; ++j) {
                                if (coordinate.x + i == 0) {
                                    if ((peptide.getLeftFlank() == '-') && (proNTermPossibleMod != null)) {
                                        for (float nTermMod : proNTermPossibleMod) {
                                            for (String modifiedAA : modifiedAAMassMap.keySet()) {
                                                if (modifiedAA.charAt(0) == aaSeqCharArray[j]) {
                                                    if (Math.abs(nTermMod + modifiedAAMassMap.get(modifiedAA) - deltaMass) <= Math.min(1.5, 4 * ms2Tolerance)) {
                                                        newPositionGapMap.put(new Coordinate(coordinate.x + i, coordinate.x + i + 1), nTermMod);
                                                        newPositionGapMap.put(new Coordinate(coordinate.x + j, coordinate.x + j + 1), modifiedAAMassMap.get(modifiedAA));
                                                        ok = true;
                                                        break;
                                                    }
                                                }
                                            }
                                            if (ok) {
                                                break;
                                            }
                                        }
                                    }
                                    if (!ok && (pepNTermPossibleMod != null)) {
                                        for (float nTermMod : pepNTermPossibleMod) {
                                            for (String modifiedAA : modifiedAAMassMap.keySet()) {
                                                if (modifiedAA.charAt(0) == aaSeqCharArray[j]) {
                                                    if (Math.abs(nTermMod + modifiedAAMassMap.get(modifiedAA) - deltaMass) <= Math.min(1.5, 4 * ms2Tolerance)) {
                                                        newPositionGapMap.put(new Coordinate(coordinate.x + i, coordinate.x + i + 1), nTermMod);
                                                        newPositionGapMap.put(new Coordinate(coordinate.x + j, coordinate.x + j + 1), modifiedAAMassMap.get(modifiedAA));
                                                        ok = true;
                                                        break;
                                                    }
                                                }
                                            }
                                            if (ok) {
                                                break;
                                            }
                                        }
                                    }
                                } else if (coordinate.x + j == peptideString.length() - 1) {
                                    if ((peptide.getRightFlank() == '-') && (proCTermPossibleMod != null)) {
                                        for (float cTermMod : proCTermPossibleMod) {
                                            for (String modifiedAA : modifiedAAMassMap.keySet()) {
                                                if (modifiedAA.charAt(0) == aaSeqCharArray[i]) {
                                                    if (Math.abs(cTermMod + modifiedAAMassMap.get(modifiedAA) - deltaMass) <= Math.min(1.5, 4 * ms2Tolerance)) {
                                                        newPositionGapMap.put(new Coordinate(coordinate.x + i, coordinate.x + i + 1), modifiedAAMassMap.get(modifiedAA));
                                                        newPositionGapMap.put(new Coordinate(coordinate.x + j, coordinate.x + j + 1), cTermMod);
                                                        ok = true;
                                                        break;
                                                    }
                                                }
                                            }
                                            if (ok) {
                                                break;
                                            }
                                        }
                                    }
                                    if (!ok && (pepCTermPossibleMod != null)) {
                                        for (float cTermMod : pepCTermPossibleMod) {
                                            for (String modifiedAA : modifiedAAMassMap.keySet()) {
                                                if (modifiedAA.charAt(0) == aaSeqCharArray[i]) {
                                                    if (Math.abs(cTermMod + modifiedAAMassMap.get(modifiedAA) - deltaMass) <= Math.min(1.5, 4 * ms2Tolerance)) {
                                                        newPositionGapMap.put(new Coordinate(coordinate.x + i, coordinate.x + i + 1), modifiedAAMassMap.get(modifiedAA));
                                                        newPositionGapMap.put(new Coordinate(coordinate.x + j, coordinate.x + j + 1), cTermMod);
                                                        ok = true;
                                                        break;
                                                    }
                                                }
                                            }
                                            if (ok) {
                                                break;
                                            }
                                        }
                                    }
                                } else {
                                    for (String modifiedAA1 : modifiedAAMassMap.keySet()) {
                                        if (modifiedAA1.charAt(0) == aaSeqCharArray[i]) {
                                            for (String modifiedAA2 : modifiedAAMassMap.keySet()) {
                                                if (modifiedAA2.charAt(0) == aaSeqCharArray[j]) {
                                                    if (Math.abs(modifiedAAMassMap.get(modifiedAA1) + modifiedAAMassMap.get(modifiedAA2) - deltaMass) <= Math.min(1.5, 5 * ms2Tolerance)) {
                                                        newPositionGapMap.put(new Coordinate(coordinate.x + i, coordinate.x + i + 1), modifiedAAMassMap.get(modifiedAA1));
                                                        newPositionGapMap.put(new Coordinate(coordinate.x + j, coordinate.x + j + 1), modifiedAAMassMap.get(modifiedAA2));
                                                        ok = true;
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                        if (ok) {
                                            break;
                                        }
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
                        // we cannot pin-point the PTM
                        newPositionGapMap.put(new Coordinate(coordinate.x, coordinate.y), positionGapMap.get(coordinate));
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
        return outputSet;
    }
}
