package proteomics.Types;

import java.util.*;

public class SparseVector {

    private Map<Integer, Float> sparseVector = new HashMap<>();

    public SparseVector(Map<Integer, Float> sparseVector) {
        for (int i : sparseVector.keySet()) {
            sparseVector.put(i, sparseVector.get(i));
        }
    }

    public SparseVector() {}

    public void add(int i, float v) {
        if (Math.abs(v) > 1e-6) {
            if (sparseVector.containsKey(i)) {
                sparseVector.put(i, sparseVector.get(i) + v);
            } else {
                sparseVector.put(i, v);
            }
        }
    }

    public void put(int i, float v) {
        if (Math.abs(v) > 1e-6) {
            sparseVector.put(i, v);
        }
    }

    public float get(int i) {
        if (sparseVector.containsKey(i)) {
            return sparseVector.get(i);
        } else {
            return 0;
        }
    }

    public Set<Integer> idxSet() {
        return sparseVector.keySet();
    }

    public Float[] getValues() {
        return sparseVector.values().toArray(new Float[sparseVector.size()]);
    }

    public float getMaxValue() {
        List<Float> intensityList = new LinkedList<>(sparseVector.values());
        Collections.sort(intensityList, Collections.reverseOrder());
        return intensityList.get(0);
    }

    public float getMinValue() {
        List<Float> intensityList = new LinkedList<>(sparseVector.values());
        Collections.sort(intensityList);
        return intensityList.get(0);
    }

    public double norm2square() {
        float output = 0;
        for (float v : sparseVector.values()) {
            output += v * v;
        }
        return output;
    }

    public double dot(SparseVector other) {
        double output = 0;
        Map<Integer, Float> otherVector = other.sparseVector;
        for (int i : sparseVector.keySet()) {
            if (otherVector.containsKey(i)) {
                output += sparseVector.get(i) * otherVector.get(i);
            }
        }
        return output;
    }

    public Map<Integer, Float> getVectorMap() {
        return sparseVector;
    }

    public Set<Integer> getNonzeroIdx() {
        return sparseVector.keySet();
    }

    public boolean isNonzero(int i) {
        return get(i) != 0;
    }
}
