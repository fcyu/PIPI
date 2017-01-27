package proteomics.Types;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class SparseBooleanVector {

    private Set<Integer> sparseVector = new HashSet<>();

    public SparseBooleanVector(Set<Integer> sparseVector) {
        this.sparseVector = sparseVector;
    }

    public SparseBooleanVector() {}

    public void put(int idx) {
        sparseVector.add(idx);
    }

    public double norm2square() {
        return sparseVector.size();
    }

    public double dot(SparseVector other) {
        double output = 0;
        Map<Integer, Float> otherVector = other.getVectorMap();
        Set<Integer> intersectedKeys = new HashSet<>(sparseVector);
        intersectedKeys.retainAll(otherVector.keySet());
        for (int i : intersectedKeys) {
            output += otherVector.get(i);
        }
        return output;
    }

    public double dot(SparseBooleanVector other) {
        Set<Integer> intersectedKeys = new HashSet<>(sparseVector);
        intersectedKeys.retainAll(other.sparseVector);
        return intersectedKeys.size();
    }

    public SparseBooleanVector deepCopy() {
        SparseBooleanVector outputVector = new SparseBooleanVector();
        for (int idx : this.sparseVector) {
            outputVector.put(idx);
        }
        return outputVector;
    }

    public boolean isZero(int idx) {
        return !sparseVector.contains(idx);
    }

    public void delete(int idx) {
        if (sparseVector.contains(idx)) {
            sparseVector.remove(idx);
        }
    }

    public int getNonZeroNum() {
        return sparseVector.size();
    }

    public Integer[] getNonZeroIdxes() {
        Integer[] outputArray = sparseVector.toArray(new Integer[sparseVector.size()]);
        Arrays.sort(outputArray);
        return outputArray;
    }
}
