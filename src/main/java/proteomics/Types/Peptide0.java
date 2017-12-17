package proteomics.Types;


import java.util.Set;

public class Peptide0 {

    public final SparseBooleanVector code;
    public final boolean isTarget;
    public final Set<String> proteins;
    public final char leftFlank;
    public final char rightFlank;

    public Peptide0(SparseBooleanVector code, boolean isTarget, Set<String> proteins, char leftFlank, char rightFlank) {
        this.code = code;
        this.isTarget = isTarget;
        this.proteins = proteins;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
    }
}
