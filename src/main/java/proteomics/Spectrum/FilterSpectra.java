package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.jmzreader.model.Spectrum;

import java.util.*;

public class FilterSpectra {

    private static final Logger logger = LoggerFactory.getLogger(FilterSpectra.class);

    private final MassScan[] scanNumArray;

    public FilterSpectra(JMzReader spectraReader, Map<String, String> parameterMap, Map<String, Float> massTable) {
        int minMs1Charge = Integer.valueOf(parameterMap.get("min_ms1_charge"));
        int maxMs1Charge = Integer.valueOf(parameterMap.get("max_ms1_charge"));
        float minPrecursorMass =  Float.valueOf(parameterMap.get("min_precursor_mass"));
        float maxPrecursorMass = Float.valueOf(parameterMap.get("max_precursor_mass"));
        int minPeakNum = Integer.valueOf(parameterMap.get("min_peak_num"));
        List<MassScan> tempMassScan = new LinkedList<>();
        Iterator<Spectrum> spectraIterator = spectraReader.getSpectrumIterator();

        while (spectraIterator.hasNext()) {
            Spectrum spectrum = spectraIterator.next();

            if (spectrum.getMsLevel() != 2) {
                continue;
            }

            if (spectrum.getPrecursorCharge() == null) {
                logger.warn("Scan {} doesn't have charge information. Skip.", spectrum.getId());
                continue;
            }
            int precursorCharge = spectrum.getPrecursorCharge();

            if ((precursorCharge < minMs1Charge) || (precursorCharge > maxMs1Charge)) {
                continue;
            }

            float precursorMz = spectrum.getPrecursorMZ().floatValue();
            float precursorMass = precursorMz * precursorCharge - precursorCharge * massTable.get("PROTON");
            if ((precursorMass > maxPrecursorMass) || (precursorMass < minPrecursorMass)) {
                continue;
            }

            if (spectrum.getPeakList().size() < minPeakNum) {
                logger.debug("Scan {} doesn't contain enough peak number ({}). Skip.", spectrum.getId(), minPeakNum);
                continue;
            }

            tempMassScan.add(new MassScan(precursorMass, spectrum.getId()));
        }

        Collections.sort(tempMassScan);
        scanNumArray = tempMassScan.toArray(new MassScan[tempMassScan.size()]);
    }

    public MassScan[] getScanNumArray() {
        return scanNumArray;
    }

    public class MassScan implements Comparable<MassScan> {

        final float mass;
        final String scanId;

        MassScan(float mass, String scanId) {
            this.mass = mass;
            this.scanId = scanId;
        }

        public int compareTo(MassScan other) {
            if (other.mass > mass) {
                return -1;
            } else if (other.mass < mass) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}