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

        // calculate Xcorr
        psm = new FinalResultEntry(spectrum.scanNum, spectrum.precursorCharge, spectrum.precursorMz);
        for (Peptide peptide : candidateList) {
            SparseBooleanVector theoIonVector = massToolObj.buildVector(peptide.getIonMatrix(), spectrum.precursorCharge);
            double xcorr = theoIonVector.dot(expXcorrPl) * 0.25; // scaling the xcorr to original SEQUEST type.
            if (xcorr > 0) {
                if (psm.noScore() || (xcorr > psm.getScore()) || ((xcorr == psm.getScore()) && psm.isDecoy() && (!peptide.isDecoy()))) {
                    psm.setPeptide(peptide);
                    psm.setGlobalSearchRank(peptide.getGlobalRank());
                    psm.setNormalizedCrossXcorr(peptide.getNormalizedCrossCorr());
                }
                psm.addScore(xcorr);
            }
        }

        if (psm.getPeptide() == null) {
            psm = null;
        }
    }

    public FinalResultEntry getScoredPSM() {
        return psm;
    }
}
