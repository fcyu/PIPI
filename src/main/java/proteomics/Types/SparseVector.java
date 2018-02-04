package proteomics.Types;

import java.util.*;

public class SparseVector {

    private Map<Integer, Double> sparseVector = new HashMap<>();

    public SparseVector(Map<Integer, Double> sparseVector) {
        for (int i : sparseVector.keySet()) {
            sparseVector.put(i, sparseVector.get(i));
        }
    }

    public SparseVector() {}

    public void add(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            if (sparseVector.containsKey(i)) {
                sparseVector.put(i, sparseVector.get(i) + v);
            } else {
                sparseVector.put(i, v);
            }
        }
    }

    public void put(int i, double v) {
        if (Math.abs(v) > 1e-6) {
            sparseVector.put(i, v);
        }
    }

    public double get(int i) {
        if (sparseVector.containsKey(i)) {
            return sparseVector.get(i);
        } else {
            return 0;
        }
    }

    public Set<Integer> idxSet() {
        return sparseVector.keySet();
    }

    public Double[] getValues() {
        return sparseVector.values().toArray(new Double[sparseVector.size()]);
    }

    public double getMaxValue() {
        List<Double> intensityList = new ArrayList<>(sparseVector.values());
        Collections.sort(intensityList, Collections.reverseOrder());
        return intensityList.get(0);
    }

    public double getMinValue() {
        List<Double> intensityList = new ArrayList<>(sparseVector.values());
        Collections.sort(intensityList);
        return intensityList.get(0);
    }

    public double norm2square() {
        double output = 0;
        for (double v : sparseVector.values()) {
            output += v * v;
        }
        return output;
    }

    public double dot(SparseVector other) {
        double output = 0;
        Map<Integer, Double> otherVector = other.sparseVector;
        Set<Integer> intersectedKeys = new HashSet<>(sparseVector.keySet());
        intersectedKeys.retainAll(otherVector.keySet());
        for (int i : intersectedKeys) {
            output += sparseVector.get(i) * otherVector.get(i);
        }
        return output;
    }

    public Map<Integer, Double> getVectorMap() {
        return sparseVector;
    }

    public Set<Integer> getNonzeroIdx() {
        return sparseVector.keySet();
    }

    public boolean isNonzero(int i) {
        return get(i) != 0;
    }
}
