package proteomics.Types;

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
        for (int i : sparseVector) {
            if (otherVector.containsKey(i)) {
                output += otherVector.get(i);
            }
        }
        return output;
    }

    public double dot(SparseBooleanVector other) {
        double dotProduct = 0;
        for (int k : other.sparseVector) {
            if (sparseVector.contains(k)) {
                ++dotProduct;
            }
        }
        return dotProduct;
    }
}
