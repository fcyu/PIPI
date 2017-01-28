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


    public Search(BuildIndex buildIndexObj, SpectrumEntry spectrumEntry, SparseVector scanCode, Map<String, Peptide0> peptideCodeMap, TreeMap<Float, Set<String>> massPeptideMap, MassTool massToolObj, float ms1Tolerance, int ms1ToleranceUnit, float minPtmMass, float maxPtmMass, int maxMs2Charge) {
        int scanNum = spectrumEntry.scanNum;
        float precursorMass = spectrumEntry.precursorMass;

        Map<String, Set<String>> peptideProteinMap = buildIndexObj.returnPepProMap();
        Map<String, String> proteinSeqMap = buildIndexObj.returnProPepMap();

        float minPeptideMass = massPeptideMap.firstKey();
        float maxPeptideMass = massPeptideMap.lastKey();

        double scanNormSquare = scanCode.norm2square();
        float leftTol = ms1Tolerance;
        float rightTol = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            leftTol = precursorMass - (precursorMass / (1 + ms1Tolerance * 1e-6f));
            rightTol = (precursorMass / (1 - ms1Tolerance * 1e-6f)) - precursorMass;
        }
        float leftMass = Math.max(precursorMass + minPtmMass, minPeptideMass);
        float rightMass = Math.min(precursorMass + maxPtmMass, maxPeptideMass);
        NavigableMap<Float, Set<String>> subMassPeptideMap = massPeptideMap.subMap(leftMass, true, rightMass, true);

        if (!subMassPeptideMap.isEmpty()) {
            PriorityQueue<ResultEntry> ptmFreeQueue = new PriorityQueue<>(rankNum * 2);
            PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum * 2);
            for (Set<String> peptideSet: subMassPeptideMap.values()) {
                for (String peptide : peptideSet) {
                    Peptide0 peptideOObj = peptideCodeMap.get(peptide);
                    double score = 0;
                    double temp1 = Math.sqrt(peptideOObj.codeNormSquare * scanNormSquare);
                    if (temp1 > 1e-6) {
                        score = peptideOObj.code.dot(scanCode) / temp1;
                    }
                    float deltaMass = peptideOObj.peptideMass - spectrumEntry.precursorMass; // caution: the order matters under ms1ToleranceUnit == 1 situation

                    if (peptideOObj.isTarget) {
                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum) {
                                ptmFreeQueue.add(new ResultEntry(score, peptide, false, false));
                            } else {
                                if (score > ptmFreeQueue.peek().score) {
                                    ptmFreeQueue.poll();
                                    ptmFreeQueue.add(new ResultEntry(score, peptide, false, false));
                                }
                            }
                        }

                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum) {
                                ptmOnlyQueue.add(new ResultEntry(score, peptide, false, true));
                            } else {
                                if (score > ptmOnlyQueue.peek().score) {
                                    ptmOnlyQueue.poll();
                                    ptmOnlyQueue.add(new ResultEntry(score, peptide, false, true));
                                }
                            }
                        }
                    } else {
                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum) {
                                ptmFreeQueue.add(new ResultEntry(score, peptide, true, false));
                            } else {
                                if (score > ptmFreeQueue.peek().score) {
                                    ptmFreeQueue.poll();
                                    ptmFreeQueue.add(new ResultEntry(score, peptide, true, false));
                                }
                            }
                        }

                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum) {
                                ptmOnlyQueue.add(new ResultEntry(score, peptide, true, true));
                            } else {
                                if (score > ptmOnlyQueue.peek().score) {
                                    ptmOnlyQueue.poll();
                                    ptmOnlyQueue.add(new ResultEntry(score, peptide, true, true));
                                }
                            }
                        }
                    }
                }
            }

            ptmFreeResult = convertResult(ptmFreeQueue, peptideProteinMap, proteinSeqMap, massToolObj, maxMs2Charge, scanNum);
            ptmOnlyResult = convertResult(ptmOnlyQueue, peptideProteinMap, proteinSeqMap, massToolObj, maxMs2Charge, scanNum);

            if (PIPI.DEV) {
                try {
                    BufferedWriter writer = new BufferedWriter(new FileWriter("PTMFreeContainingGlobalScoreBound" + "." + scanNum + ".csv"));
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
                        writer.write(String.format("%f,%f,", lowScore, highScore));
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
                        writer.write(String.format("%f,%f\n", lowScore, highScore));
                    } else {
                        writer.write("-,-\n");
                    }
                    writer.close();

                    writer = new BufferedWriter(new FileWriter("ptm_only_candidates" + "." + scanNum + ".csv"));
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

    private List<Peptide> convertResult(PriorityQueue<ResultEntry> inputQueue, Map<String, Set<String>> peptideProteinMap, Map<String, String> proteinSeqMap, MassTool massToolObj, int maxMs2Charge, int scanNum) {
        List<Peptide> peptideList = new LinkedList<>();
        int globalRank = inputQueue.size();
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            String leftFlank = "-";
            String rightFlank = "-";
            String peptideString = temp.peptide;
            if (peptideProteinMap.containsKey(peptideString)) {
                String proSeq = proteinSeqMap.get(peptideProteinMap.get(peptideString).iterator().next());
                int startIdx = proSeq.indexOf(peptideString);
                if (startIdx == -1) {
                    logger.warn("Something wrong happened in Search.java (line: 223), scan num = {}, peptide = {}.", scanNum, peptideString);
                } else if (startIdx == 0) {
                    int tempIdx = peptideString.length();
                    if (tempIdx < proSeq.length()) {
                        rightFlank = proSeq.substring(tempIdx, tempIdx + 1);
                    } else {
                        rightFlank = "-";
                    }
                } else if (startIdx == proSeq.length() - peptideString.length()) {
                    leftFlank = proSeq.substring(startIdx - 1, startIdx);
                } else {
                    leftFlank = proSeq.substring(startIdx - 1, startIdx);
                    rightFlank = proSeq.substring(startIdx + peptideString.length(), startIdx + peptideString.length() + 1);
                }
            }
            peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massToolObj, maxMs2Charge, temp.score, leftFlank, rightFlank, globalRank));
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
