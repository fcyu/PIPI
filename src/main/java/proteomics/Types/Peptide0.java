package proteomics.Types;

public class Peptide0 {

    public final float peptideMass;
    public final SparseBooleanVector code;
    public final double codeNormSquare;
    public final boolean isTarget;

    public Peptide0(float peptideMass, SparseBooleanVector code, double codeNormSquare, boolean isTarget) {
        this.peptideMass = peptideMass;
        this.code = code;
        this.codeNormSquare = codeNormSquare;
        this.isTarget = isTarget;
    }
}
