package proteomics.Search;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.PIPI;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Search {

    private static final Logger logger = LoggerFactory.getLogger(Search.class);
    private static final int rankNum = 10;

    private List<Peptide> ptmOnlyResult = new LinkedList<>();
    private List<Peptide> ptmFreeResult = new LinkedList<>();


    public Search(BuildIndex buildIndexObj, SpectrumEntry spectrumEntry, SparseVector scanCode, MassTool massToolObj, float ms1Tolerance, int ms1ToleranceUnit, float minPtmMass, float maxPtmMass, int maxMs2Charge) {
        PriorityQueue<ResultEntry> ptmFreeQueue = new PriorityQueue<>(rankNum * 2);
        PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum * 2);
        double scanNormSquare = scanCode.norm2square();
        float leftTol = ms1Tolerance;
        float rightTol = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            leftTol = spectrumEntry.precursorMass - (spectrumEntry.precursorMass / (1 + ms1Tolerance * 1e-6f));
            rightTol = (spectrumEntry.precursorMass / (1 - ms1Tolerance * 1e-6f)) - spectrumEntry.precursorMass;
        }
        float leftMass = Math.max(spectrumEntry.precursorMass + minPtmMass, buildIndexObj.getMinPeptideMass());
        float rightMass = Math.min(spectrumEntry.precursorMass + maxPtmMass, buildIndexObj.getMaxPeptideMass());

        if (leftMass >= rightMass) {
            return;
        }

        Map<String, Peptide0> peptide0Map = buildIndexObj.getPeptide0Map();
        TreeMap<Float, Set<String>> massPeptideMap = buildIndexObj.getMassPeptideMap();

        NavigableMap<Float, Set<String>> subMassPeptideMap = massPeptideMap.subMap(leftMass, true, rightMass, true);

        if (!subMassPeptideMap.isEmpty()) {
            for (float mass : subMassPeptideMap.keySet()) {
                for (String sequence : massPeptideMap.get(mass)) {
                    Peptide0 peptide0 = peptide0Map.get(sequence);
                    double score = 0;
                    double temp1 = Math.sqrt(peptide0.code.norm2square() * scanNormSquare);
                    if (temp1 > 1e-6) {
                        score = peptide0.code.dot(scanCode) / temp1;
                    }
                    float deltaMass = peptide0.mass - spectrumEntry.precursorMass; // caution: the order matters under ms1ToleranceUnit == 1 situation

                    if (peptide0.isTarget) {
                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum) {
                                ptmFreeQueue.add(new ResultEntry(score, peptide0.sequence, peptide0.leftFlank, peptide0.rightFlank, false, false));
                            } else {
                                if (score > ptmFreeQueue.peek().score) {
                                    ptmFreeQueue.poll();
                                    ptmFreeQueue.add(new ResultEntry(score, peptide0.sequence, peptide0.leftFlank, peptide0.rightFlank, false, false));
                                }
                            }
                        }

                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum) {
                                ptmOnlyQueue.add(new ResultEntry(score, peptide0.sequence, peptide0.leftFlank, peptide0.rightFlank, false, true));
                            } else {
                                if (score > ptmOnlyQueue.peek().score) {
                                    ptmOnlyQueue.poll();
                                    ptmOnlyQueue.add(new ResultEntry(score, peptide0.sequence, peptide0.leftFlank, peptide0.rightFlank, false, true));
                                }
                            }
                        }
                    } else {
                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum) {
                                ptmFreeQueue.add(new ResultEntry(score, peptide0.sequence, peptide0.leftFlank, peptide0.rightFlank, true, false));
                            } else {
                                if (score > ptmFreeQueue.peek().score) {
                                    ptmFreeQueue.poll();
                                    ptmFreeQueue.add(new ResultEntry(score, peptide0.sequence, peptide0.leftFlank, peptide0.rightFlank, true, false));
                                }
                            }
                        }

                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum) {
                                ptmOnlyQueue.add(new ResultEntry(score, peptide0.sequence, peptide0.leftFlank, peptide0.rightFlank, true, true));
                            } else {
                                if (score > ptmOnlyQueue.peek().score) {
                                    ptmOnlyQueue.poll();
                                    ptmOnlyQueue.add(new ResultEntry(score, peptide0.sequence, peptide0.leftFlank, peptide0.rightFlank, true, true));
                                }
                            }
                        }
                    }
                }
            }
        }

        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            ptmFreeResult = convertResult(ptmFreeQueue, massToolObj, maxMs2Charge);
            ptmOnlyResult = convertResult(ptmOnlyQueue, massToolObj, maxMs2Charge);

            if (PIPI.DEV) {
                try {
                    BufferedWriter writer = new BufferedWriter(new FileWriter("PTMFreeContainingGlobalScoreBound" + "." + spectrumEntry.scanNum + "." + spectrumEntry.precursorCharge + ".csv"));
                    writer.write("PTM_free_low,PTM_free_high,PTM_containing_low,PTM_containing_high\n");
                    if (!ptmFreeQueue.isEmpty()) {
                        double highScore = 0;
                        double lowScore = 1;
                        for (ResultEntry temp : ptmFreeQueue) {
                            if (temp.score > highScore) {
                                highScore = temp.score;
                            }
                            if (temp.score < lowScore) {
                                lowScore = temp.score;
                            }
                        }
                        writer.write(String.format(Locale.US, "%f,%f,", lowScore, highScore));
                    } else {
                        writer.write("-,-,");
                    }
                    if (!ptmOnlyQueue.isEmpty()) {
                        double highScore = 0;
                        double lowScore = 1;
                        for (ResultEntry temp : ptmOnlyQueue) {
                            if (temp.score > highScore) {
                                highScore = temp.score;
                            }
                            if (temp.score < lowScore) {
                                lowScore = temp.score;
                            }
                        }
                        writer.write(String.format(Locale.US, "%f,%f\n", lowScore, highScore));
                    } else {
                        writer.write("-,-\n");
                    }
                    writer.close();

                    writer = new BufferedWriter(new FileWriter("ptm_only_candidates" + "." + spectrumEntry.scanNum + "." + spectrumEntry.precursorCharge + ".csv"));
                    writer.write("peptide,globalRank,is_decoy\n");
                    for (Peptide peptide : ptmOnlyResult) {
                        if (peptide.isDecoy()) {
                            writer.write(peptide.getPTMFreeSeq() + "," + peptide.getGlobalRank() + ",1\n");
                        } else {
                            writer.write(peptide.getPTMFreeSeq() + "," + peptide.getGlobalRank() + ",0\n");
                        }
                    }
                    writer.close();
                } catch (IOException ex) {
                    ex.printStackTrace();
                    System.exit(1);
                }
            }
        }
    }

    private List<Peptide> convertResult(PriorityQueue<ResultEntry> inputQueue, MassTool massToolObj, int maxMs2Charge) {
        List<Peptide> peptideList = new LinkedList<>();
        int globalRank = inputQueue.size();
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massToolObj, maxMs2Charge, temp.score, temp.leftFlank, temp.rightFlank, globalRank));
            --globalRank;
        }

        return peptideList;
    }

    public List<Peptide> getPTMOnlyResult() {
        return ptmOnlyResult;
    }

    public List<Peptide> getPTMFreeResult() {
        return ptmFreeResult;
    }
}
