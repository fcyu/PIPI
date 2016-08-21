package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.SpectrumEntry;
import uk.ac.ebi.pride.tools.jmzreader.*;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.util.*;

public class PreSpectraOld {

    private static final Logger logger = LoggerFactory.getLogger(PreSpectra.class);

    private Map<Integer, SpectrumEntry> numSpectrumMap = new HashMap<>();
    private TreeMap<Float, List<Integer>> massNumMap = new TreeMap<>();

    public PreSpectraOld(JMzReader spectraParser, Map<String, String> parameterMap, MassTool massToolObj, String ext) {
        int minMs1Charge = Integer.valueOf(parameterMap.get("min_ms1_charge"));
        int maxMs1Charge = Math.min(6, Integer.valueOf(parameterMap.get("max_ms1_charge")));
        float ms2Tolerance = Float.valueOf(parameterMap.get("ms2_tolerance"));
        float minPrecursorMass =  Float.valueOf(parameterMap.get("min_precursor_mass"));
        float maxPrecursorMass = Float.valueOf(parameterMap.get("max_precursor_mass"));
        int minPeakNum = Integer.valueOf(parameterMap.get("min_peak_num"));
        Map<String, Float> massTable = massToolObj.returnMassTable();

        PreSpectrum preSpectrumObj = new PreSpectrum(massToolObj);
        Iterator<Spectrum> spectrumIterator = spectraParser.getSpectrumIterator();
        while (spectrumIterator.hasNext()) {
            Spectrum spectrum = spectrumIterator.next();

            if (spectrum.getMsLevel() != 2) {
                continue;
            }

            float precursorMz = spectrum.getPrecursorMZ().floatValue();


            if (spectrum.getPrecursorCharge() == null) {
                logger.debug("Scan {} doesn't contain precursor charge information.", Integer.valueOf(spectrum.getId()));
                continue;
            }

            int precursorCharge = spectrum.getPrecursorCharge();
            float precursorMass = precursorMz * precursorCharge - precursorCharge * massTable.get("PROTON");

            if ((precursorCharge < minMs1Charge) || (precursorCharge > maxMs1Charge)) {
                continue;
            }

            if ((precursorMass > maxPrecursorMass) || (precursorMass < minPrecursorMass)) {
                continue;
            }

            Map<Double, Double> rawMzIntensityMap = spectrum.getPeakList();
            if (rawMzIntensityMap.size() < minPeakNum) {
                continue;
            }

            TreeMap<Float, Float> plMap = preSpectrumObj.preSpectrum(rawMzIntensityMap, precursorMass, precursorCharge, ms2Tolerance);
            if (plMap.size() <= minPeakNum) {
                continue;
            }

            int scanNum = Integer.valueOf(spectrum.getId());
            try {
                if (ext.contentEquals("mgf")) {
                    String title = ((Ms2Query) spectrum).getTitle();
                    String[] temp = title.split("\\.");
                    scanNum = Integer.valueOf(temp[temp.length - 2]);
                }
            } catch (Exception ex) {}

            SpectrumEntry spectrumEntry = new SpectrumEntry(scanNum, precursorMz, precursorMass, precursorCharge, plMap, null);

            if (massNumMap.containsKey(precursorMass)) {
                List<Integer> spectrumList = massNumMap.get(precursorMass);
                spectrumList.add(scanNum);
                massNumMap.put(precursorMass, spectrumList);
            } else {
                List<Integer> scanNumList = new LinkedList<>();
                scanNumList.add(scanNum);
                massNumMap.put(precursorMass, scanNumList);
            }

            numSpectrumMap.put(scanNum, spectrumEntry);
        }
    }

    public Map<Integer, SpectrumEntry> returnNumSpectrumMap() {
        return numSpectrumMap;
    }

    public TreeMap<Float, List<Integer>> returnMassNumMap() {
        return massNumMap;
    }
}
