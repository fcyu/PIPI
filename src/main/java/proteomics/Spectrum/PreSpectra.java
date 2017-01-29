package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.TheoSeq.MassTool;
import proteomics.Types.SpectrumEntry;
import uk.ac.ebi.pride.tools.jmzreader.*;
import uk.ac.ebi.pride.tools.jmzreader.model.*;
import uk.ac.ebi.pride.tools.mgf_parser.model.Ms2Query;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;

public class PreSpectra {

    private static final Logger logger = LoggerFactory.getLogger(PreSpectra.class);

    private Map<Integer, SpectrumEntry> numSpectrumMap = new HashMap<>();

    public PreSpectra(JMzReader spectraParser, Map<String, String> parameterMap, MassTool massToolObj, String ext) {
        int minMs1Charge = Integer.valueOf(parameterMap.get("min_ms1_charge"));
        int maxMs1Charge = Integer.valueOf(parameterMap.get("max_ms1_charge"));
        float minPrecursorMass =  Float.valueOf(parameterMap.get("min_precursor_mass"));
        float maxPrecursorMass = Float.valueOf(parameterMap.get("max_precursor_mass"));
        int minPeakNum = Integer.valueOf(parameterMap.get("min_peak_num"));
        float ms2Tolerance = Float.valueOf(parameterMap.get("ms2_tolerance"));
        Map<String, Float> massTable = massToolObj.returnMassTable();

        PreSpectrum preSpectrumObj = new PreSpectrum(massToolObj);
        PrintStream originalStream = System.out;
        PrintStream nullStream = new PrintStream(new OutputStream() {
            @Override
            public void write(int b) throws IOException {}
        });
        System.setOut(nullStream);

        Iterator<Spectrum> spectrumIterator = spectraParser.getSpectrumIterator();
        while (spectrumIterator.hasNext()) {
            Spectrum spectrum = spectrumIterator.next();

            if (spectrum.getMsLevel() != 2) {
                continue;
            }

            Map<Double, Double> rawMzIntensityMap = spectrum.getPeakList();
            if (rawMzIntensityMap.size() < minPeakNum) {
                logger.debug("Scan {} doesn't contain enough peak number ({}). Skip.", spectrum.getId(), minPeakNum);
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

            float precursorMz = spectrum.getPrecursorMZ().floatValue();

            SpectrumEntry spectrumEntry;
            int precursorCharge;
            float precursorMass;
            if (spectrum.getPrecursorCharge() == null) {
                logger.debug("Scan {} doesn't have charge information.", spectrum.getId());
                precursorCharge = 0;
                precursorMass = 0;
                TreeMap<Float, Float> plMap = preSpectrumObj.preSpectrum(rawMzIntensityMap, precursorMass, precursorCharge, ms2Tolerance);
                if (plMap.size() <= minPeakNum) {
                    continue;
                }
                spectrumEntry = new SpectrumEntry(scanNum, precursorMz, precursorMass, precursorCharge, plMap);
            } else {
                precursorCharge = spectrum.getPrecursorCharge();
                if ((precursorCharge < minMs1Charge) || (precursorCharge > maxMs1Charge)) {
                    continue;
                }
                precursorMass = precursorMz * precursorCharge - precursorCharge * massTable.get("PROTON");
                if ((precursorMass > maxPrecursorMass) || (precursorMass < minPrecursorMass)) {
                    continue;
                }
                TreeMap<Float, Float> plMap = preSpectrumObj.preSpectrum(rawMzIntensityMap, precursorMass, precursorCharge, ms2Tolerance);
                if (plMap.size() <= minPeakNum) {
                    continue;
                }
                spectrumEntry = new SpectrumEntry(scanNum, precursorMz, precursorMass, precursorCharge, plMap);
            }

            numSpectrumMap.put(scanNum, spectrumEntry);
        }

        System.setOut(originalStream);
    }

    public Map<Integer, SpectrumEntry> returnNumSpectrumMap() {
        return numSpectrumMap;
    }
}
