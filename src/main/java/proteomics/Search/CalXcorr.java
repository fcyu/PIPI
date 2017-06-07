package proteomics.Search;


import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Spectrum.PreSpectrum;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.*;

import java.util.*;

public class CalXcorr {

    private static final Logger logger = LoggerFactory.getLogger(CalXcorr.class);

    private FinalResultEntry psm;

    public CalXcorr(Set<Peptide> candidateSet, SpectrumEntry spectrum, MassTool massToolObj) {
        PreSpectrum preSpectrumObj = new PreSpectrum(massToolObj);

        // prepare the XCORR vector
        SparseVector expXcorrPl = preSpectrumObj.prepareXcorr(spectrum.unprocessedPlMap);

        Map<String, LinkedList<PeptideScore>> modSequences = new HashMap<>();

        // calculate Xcorr
        psm = new FinalResultEntry(spectrum.scanNum, spectrum.precursorCharge, spectrum.precursorMz);
        for (Peptide peptide : candidateSet) {
            SparseBooleanVector theoIonVector = massToolObj.buildVector(peptide.getIonMatrix(), spectrum.precursorCharge);
            double xcorr = theoIonVector.dot(expXcorrPl) * 0.25; // scaling the xcorr to original SEQUEST type.
            if (xcorr > 0) {
                if (psm.noScore() || (xcorr > psm.getScore()) || ((xcorr == psm.getScore()) && psm.isDecoy() && (!peptide.isDecoy()))) {
                    psm.setPeptide(peptide);
                    psm.setGlobalSearchRank(peptide.getGlobalRank());
                    psm.setNormalizedCrossXcorr(peptide.getNormalizedCrossCorr());
                }
                psm.addScore(xcorr);

                if (peptide.hasVarPTM()) {
                    // record scores with different PTM patterns for calculating PTM delta score.
                    if (modSequences.containsKey(peptide.getPTMFreeSeq())) {
                        LinkedList<PeptideScore> temp = modSequences.get(peptide.getPTMFreeSeq());
                        if (temp.size() < 2) {
                            temp.add(new PeptideScore(xcorr, peptide));
                            temp.sort(Collections.reverseOrder());
                        } else if (xcorr > temp.peekLast().score) {
                            temp.pollLast();
                            temp.add(new PeptideScore(xcorr, peptide));
                            temp.sort(Collections.reverseOrder());
                        }
                    } else {
                        LinkedList<PeptideScore> temp = new LinkedList<>();
                        temp.add(new PeptideScore(xcorr, peptide));
                        modSequences.put(peptide.getPTMFreeSeq(), temp);
                    }
                }
            }
        }

        if (psm.getPeptide() == null) {
            psm = null;
        } else {
            if (psm.getPeptide().hasVarPTM()) {
                LinkedList<PeptideScore> temp = modSequences.get(psm.getPeptide().getPTMFreeSeq());
                if (temp.size() == 1) {
                    psm.setPtmDeltasScore(temp.peekFirst().score);
                } else {
                    psm.setPtmDeltasScore(temp.peekFirst().score - temp.get(1).score);
                }
            } else {
                psm.setPtmDeltasScore(psm.getScore());
            }
        }
    }

    public FinalResultEntry getScoredPSM() {
        return psm;
    }
}
