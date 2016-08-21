package proteomics.Types;

public class PercolatorEntry {

    public final double percolatorScore;
    public final String qValue;
    public final String PEP;

    public PercolatorEntry(double percolatorScore, String qValue, String PEP) {
        this.percolatorScore = percolatorScore;
        this.qValue = qValue;
        this.PEP = PEP;
    }
}
