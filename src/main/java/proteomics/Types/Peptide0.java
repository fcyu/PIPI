package proteomics.Types;


import java.util.Set;

public class Peptide0 {

    public final String sequence;
    public final float mass;
    public final SparseBooleanVector code;
    public final boolean isTarget;
    public final Set<String> proteins;
    public final char leftFlank;
    public final char rightFlank;

    public Peptide0(String sequence, float mass, SparseBooleanVector code, boolean isTarget, Set<String> proteins, char leftFlank, char rightFlank) {
        this.sequence = sequence;
        this.mass = mass;
        this.code = code;
        this.isTarget = isTarget;
        this.proteins = proteins;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
    }
}
