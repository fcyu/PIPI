package proteomics.Spectrum;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class IsotopeDistribution {

    private static final Logger logger = LoggerFactory.getLogger(IsotopeDistribution.class);
    private static final Pattern pattern1 = Pattern.compile("^([A-Za-z]+)\\s+([0-9]+)");
    private static final Pattern pattern2 = Pattern.compile("^([0-9\\.]+)\\s+([0-9\\.]+)");

    private static final double averagineAverageMass = 111.1254;
    private static final double averagineMonoMass = 111.05429999675;
    private static final double averagineC = 4.9384;
    private static final double averagineH = 7.7583;
    private static final double averagineN = 1.3577;
    private static final double averagineO = 1.4773;
    private static final double averagineS = 0.0417;

    private static final double CMono = 12;
    private static final double HMono = 1.0078246;
    private static final double NMono = 14.0030732;
    private static final double OMono = 15.9949141;
    private static final double SMono = 31.972070;

    private final Map<String, Peak[]> elemMap = new HashMap<>();
    private final double limit;

    public static void main(String[] args) {
        IsotopeDistribution isotopeDistribution = new IsotopeDistribution(0);

        try {
            Map<String, Integer> formMap = new HashMap<>();
            formMap.put("H", 1);
            formMap.put("C", 2);
            formMap.put("O", 3);
            formMap.put("N", 4);
            List<Peak> result = isotopeDistribution.calculate(formMap);
            int i = 0;
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    public IsotopeDistribution(double limit) {
        this.limit = limit;

        // reading ISOTOPE.DAT
        InputStream inputStream = getClass().getClassLoader().getResourceAsStream("ISOTOPE.DAT");
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!line.isEmpty()) {
                    Matcher matcher1 = pattern1.matcher(line);
                    if (matcher1.find()) {
                        String elem = matcher1.group(1).trim();
                        if (elemMap.containsKey(elem)) {
                            throw new Exception("Something is wrong in reading ISOTOPE.DAT (" + line + ").");
                        } else {
                            int peakNum = Integer.valueOf(matcher1.group(2).trim());
                            Peak[] peakArray = new Peak[peakNum];
                            for (int i = 0; i < peakNum; ++i) {
                                Matcher matcher2 = pattern2.matcher(reader.readLine().trim());
                                if (matcher2.find()) {
                                    peakArray[i] = new Peak(Double.valueOf(matcher2.group(1)), Double.valueOf(matcher2.group(2)));
                                } else {
                                    throw new Exception("Something is wrong in reading ISOTOPE.DAT (" + line + ").");
                                }
                            }
                            elemMap.put(elem, peakArray);
                        }
                    } else {
                        throw new Exception("Something is wrong in reading ISOTOPE.DAT (" + line + ").");
                    }
                }
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    Map<String, Integer> getElementMapFromMonoMass(double mass) {
        double unit = mass / averagineMonoMass;
        int C = (int) Math.round(averagineC * unit);
        int N = (int) Math.round(averagineN * unit);
        int O = (int) Math.round(averagineO * unit);
        int S = (int) Math.round(averagineS * unit);
        int H = (int) Math.round((mass - C * CMono - N * NMono - O * OMono - S * SMono) / HMono);

        Map<String, Integer> elementMap = new HashMap<>();
        if (C > 0) {
            elementMap.put("C", C);
        }
        if (N > 0) {
            elementMap.put("N", N);
        }
        if (O > 0) {
            elementMap.put("O", O);
        }
        if (S > 0) {
            elementMap.put("S", S);
        }
        if (H > 0) {
            elementMap.put("H", H);
        }
        return elementMap;
    }

    List<Peak> calculate(Map<String, Integer> formMap) throws Exception {
        List<Peak> result = new ArrayList<>();
        result.add(new Peak(0, 1));
        for (String element : formMap.keySet()) {
            List<List<Peak>> sal = new ArrayList<>();
            sal.add(Arrays.asList(elemMap.get(element)));
            int n = formMap.get(element);
            int j = 0;
            while (n > 0) {
                if (j == sal.size()) {
                    sal.add(convolute(sal.get(j - 1), sal.get(j - 1)));
                    prune(sal.get(j), limit);
                }
                if ((n & 1) != 0) {
                    result = convolute(result, sal.get(j));
                    prune(result, limit);
                }
                n >>= 1;
                ++j;
            }
        }

        prune(result, 1e-6f); // prune the final result to a reasonable precision.
        return result;
    }

    private List<Peak> convolute(List<Peak> g, List<Peak> f) {
        List<Peak> h = new ArrayList<>();
        int gN = g.size();
        int fN = f.size();
        if (gN != 0 && fN != 0) {
            for (int k = 0; k < gN + fN - 1; ++k) {
                double sumWeight = 0;
                double sumMass = 0;
                int start = Math.max(0, k - fN + 1);
                int end = Math.min(gN - 1, k);
                for (int i = start; i <= end; ++i) {
                    double weight = g.get(i).realArea * f.get(k - i).realArea;
                    double mass = g.get(i).mass + f.get(k - i).mass;
                    sumWeight += weight;
                    sumMass += weight * mass;
                }
                Peak peak;
                if (sumWeight == 0) {
                    peak = new Peak(-1, sumWeight);
                } else {
                    peak = new Peak(sumMass / sumWeight, sumWeight);
                }
                h.add(peak);
            }
        }
        return h;
    }

    private void prune(List<Peak> f, double limit) {
        List<Peak> toBeDeleteList = new LinkedList<>();
        for (int i = 0; i < f.size(); ++i) {
            if (f.get(i).realArea <= limit) {
                toBeDeleteList.add(f.get(i));
            }
        }
        f.removeAll(toBeDeleteList);
    }


    public class Peak {

        public final double mass;
        public final double realArea;

        public Peak(double mass, double realArea) {
            this.mass = mass;
            this.realArea = realArea;
        }
    }
}
