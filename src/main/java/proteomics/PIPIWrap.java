package proteomics;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Search.CalSubscores;
import proteomics.Search.CalScore;
import proteomics.Search.Search;
import proteomics.Segment.InferenceSegment;
import proteomics.Spectrum.PreSpectrum;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;

public class PIPIWrap implements Callable<FinalResultEntry> {

    private static final Logger logger = LoggerFactory.getLogger(PIPIWrap.class);

    private final BuildIndex buildIndexObj;
    private final MassTool massToolObj;
    private final InferenceSegment inference3SegmentObj;
    private final SpectrumEntry spectrumEntry;
    private final float ms1Tolerance;
    private final int ms1ToleranceUnit;
    private final float ms2Tolerance;
    private final float minPtmMass;
    private final float maxPtmMass;
    private final int maxMs2Charge;
    private final InferPTM inferPTM;

    public Map<Integer, List<DevEntry>> scanDevMap;

    public PIPIWrap(BuildIndex buildIndexObj, MassTool massToolObj, SpectrumEntry spectrumEntry, float ms1Tolerance, int ms1ToleranceUnit, float ms2Tolerance, float minPtmMass, float maxPtmMass, int maxMs2Charge, Map<Integer, List<DevEntry>> scanDevMap) {
        this.buildIndexObj = buildIndexObj;
        this.massToolObj = massToolObj;
        this.spectrumEntry = spectrumEntry;
        this.ms1Tolerance = ms1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.maxMs2Charge = maxMs2Charge;
        inference3SegmentObj = buildIndexObj.getInference3SegmentObj();
        inferPTM = new InferPTM(massToolObj, maxMs2Charge, buildIndexObj.returnFixModMap(), inference3SegmentObj.getVarModParamSet(), minPtmMass, maxPtmMass);
        this.scanDevMap = scanDevMap;
    }

    @Override
    public FinalResultEntry call() {
        // Coding
        List<ThreeExpAA> expAaLists = inference3SegmentObj.inferSegmentLocationFromSpectrum(spectrumEntry);
        if (!expAaLists.isEmpty()) {
            SparseVector scanCode = inference3SegmentObj.generateSegmentIntensityVector(expAaLists);

            // Begin search.
            Search searchObj = new Search(buildIndexObj, spectrumEntry, scanCode, massToolObj, ms1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, maxMs2Charge);

            // prepare the spectrum
            PreSpectrum preSpectrumObj = new PreSpectrum(massToolObj);
            SparseVector expProcessedPL;
            if (PIPI.useXcorr) {
                expProcessedPL = preSpectrumObj.prepareXcorr(spectrumEntry.plMap, false);
            } else {
                expProcessedPL = preSpectrumObj.prepareDigitizedPL(spectrumEntry.plMap, false);
            }

            FinalResultEntry psm = new FinalResultEntry(spectrumEntry.scanNum, spectrumEntry.precursorCharge, spectrumEntry.precursorMz, spectrumEntry.mgfTitle);

            float precursorMass = spectrumEntry.precursorMass;
            float localMS1ToleranceL = -1 * ms1Tolerance;
            float localMS1ToleranceR = ms1Tolerance;
            if (ms1ToleranceUnit == 1) {
                localMS1ToleranceL = (precursorMass / (1 + ms1Tolerance * 1e-6f)) - precursorMass;
                localMS1ToleranceR = (precursorMass / (1 - ms1Tolerance * 1e-6f)) - precursorMass;
            }
            // TreeMap<Float, Float> expPL = preSpectrumObj.selectTopN(plMap, topN, range);
            double p = (double) spectrumEntry.plMap.size() / (double) ((int) precursorMass + 1);

            // infer PTM using the new approach
            List<DevEntry> devEntryList = new LinkedList<>();
            Map<String, TreeSet<Peptide>> modSequences = new TreeMap<>();
            for (Peptide peptide : searchObj.getPTMOnlyResult()) {
                PeptidePTMPattern peptidePTMPattern = inferPTM.tryPTM(expProcessedPL, precursorMass, peptide.getPTMFreeSeq(), peptide.isDecoy(), peptide.getNormalizedCrossCorr(), peptide.getLeftFlank(), peptide.getRightFlank(), peptide.getGlobalRank(), maxMs2Charge, p, localMS1ToleranceL, localMS1ToleranceR);
                if (peptidePTMPattern.getTopEntry() != null) {
                    psm.addScore(peptidePTMPattern.getTopEntry().peptide);
                    // record scores with different PTM patterns for calculating PTM delta score.
                    if (modSequences.containsKey(peptidePTMPattern.getTopEntry().peptide.getPTMFreeSeq())) {
                        TreeSet<Peptide> temp = modSequences.get(peptidePTMPattern.getTopEntry().peptide.getPTMFreeSeq());
                        if (temp.size() < 5) {
                            temp.add(peptidePTMPattern.getTopEntry().peptide);
                        } else if (peptidePTMPattern.getTopEntry().peptide.getScore() > temp.last().getScore()) {
                            temp.pollLast();
                            temp.add(peptidePTMPattern.getTopEntry().peptide);
                        }
                    } else {
                        TreeSet<Peptide> temp = new TreeSet<>(Collections.reverseOrder());
                        temp.add(peptidePTMPattern.getTopEntry().peptide);
                        modSequences.put(peptidePTMPattern.getTopEntry().peptide.getPTMFreeSeq(), temp);
                    }

                    if (PIPI.DEV) {
                        devEntryList.add(new DevEntry(peptidePTMPattern.getCandidateNum(), peptidePTMPattern.isStopped(), peptidePTMPattern.ptmFreeSequence, peptidePTMPattern.getTopEntry().peptide.getVarPtmContainingSeq()));
                    }
                }

                if (PIPI.DEBUG) {
                    try (BufferedWriter writer = new BufferedWriter(new FileWriter(spectrumEntry.scanNum + "_" + peptide.getPTMFreeSeq() + "_PtmInferDev.csv"))) {
                        writer.write("candidateNum,log10CandidateNum,Log10PValueLowerBound,Log10EValueLowerBound,log10PValue,Log10EValue,matchedPeakNum,ptmContainingSequence,score\n");
                        TreeMap<Integer, PeptidePTMPattern.Entry> tempMap = peptidePTMPattern.getCandidateNumEntryMap();
                        TreeSet<Peptide> tempSet = new TreeSet<>(Comparator.reverseOrder());
                        if (modSequences.containsKey(peptide.getPTMFreeSeq())) {
                            tempSet = modSequences.get(peptide.getPTMFreeSeq());
                        }
                        for (int candidateNum : tempMap.keySet()) {
                            PeptidePTMPattern.Entry entry = tempMap.get(candidateNum);
                            double score = -1;
                            for(Peptide peptide1 : tempSet) {
                                if (peptide1.equals(entry.peptide)) {
                                    score = peptide1.getScore();
                                    break;
                                }
                            }
                            writer.write(candidateNum + "," + Math.log10(candidateNum) + "," + Math.log10(entry.pValueLowerBound) + "," + Math.log10(entry.eValueLowerBound) + "," + Math.log10(entry.pValue) + "," + Math.log10(entry.eValue) + "," + entry.matchedPeakNum + "," + entry.peptide.getVarPtmContainingSeq() + "," + score + "\n");
                        }
                    } catch (IOException ex) {
                        ex.printStackTrace();
                        logger.error(ex.getMessage());
                        System.exit(1);
                    }
                }
            }

            if (PIPI.DEV) {
                scanDevMap.put(spectrumEntry.scanNum, devEntryList);
            }

            // Calculate Score for PTM free peptide
            for (Peptide peptide : searchObj.getPTMFreeResult()) {
                CalScore.calScore(peptide, expProcessedPL, psm, massToolObj, null);
            }

            if (psm.hasHit()) {
                psm.setPtmPatterns(modSequences);
                for (Peptide peptide : psm.getPeptideSet()) {
                    new CalSubscores(peptide, spectrumEntry, ms2Tolerance);
                }
                return psm;
            } else {
                return null;
            }
        } else {
            return null;
        }
    }

    public class DevEntry {

        public final int totalCheckedNum;
        public final boolean stopped;
        public final String ptmFreePeptide;
        public final String peptide;

        public boolean isFinalResult = false;

        public DevEntry(int totalCheckedNum, boolean stopped, String ptmFreePeptide, String peptide) {
            this.totalCheckedNum = totalCheckedNum;
            this.stopped = stopped;
            this.ptmFreePeptide = ptmFreePeptide;
            this.peptide = peptide;
        }
    }
}
