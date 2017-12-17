package proteomics.Types;


public class Peptide0 {

    public final SparseBooleanVector code;
    public final boolean isTarget;
    public final String[] proteins;
    public final char leftFlank;
    public final char rightFlank;

    public Peptide0(SparseBooleanVector code, boolean isTarget, String[] proteins, char leftFlank, char rightFlank) {
        this.code = code;
        this.isTarget = isTarget;
        this.proteins = proteins;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
    }
}
