package proteomics.Types;

import java.util.Comparator;
import java.util.TreeSet;

public class PeptidePTMPattern {

    public final String ptmFreePeptide;

    private TreeSet<Peptide> peptideTreeSet = new TreeSet<>(Comparator.reverseOrder());

    public PeptidePTMPattern(String ptmFreePeptide) {
        this.ptmFreePeptide = ptmFreePeptide;
    }

    public void update(Peptide peptide) {
        if (peptideTreeSet.size() < 5) {
            peptideTreeSet.add(peptide);
        } else if (peptideTreeSet.last().compareTo(peptide) < 0) {
            peptideTreeSet.pollLast();
            peptideTreeSet.add(peptide);
        }
    }

    public TreeSet<Peptide> getPeptideTreeSet() {
        return peptideTreeSet;
    }
}
