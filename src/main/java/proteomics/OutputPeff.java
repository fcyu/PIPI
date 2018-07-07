package proteomics;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Parameter.Parameter;
import ProteomicsLibrary.*;
import ProteomicsLibrary.Types.AA;
import ProteomicsLibrary.Utilities;
import proteomics.Types.ModEntry;
import proteomics.Types.Peptide0;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Stream;


public class OutputPeff {

    private static Logger logger = LoggerFactory.getLogger(OutputPeff.class);
    private static double aScoreT = 13;

    public static void main(String[] args) {
        String parameterPath = "D:\\Dropbox\\Results\\PIPI\\TCGA_Ovarian_Cancer\\parameter.def";
        String pipiDirPath = "D:\\Dropbox\\Results\\PIPI\\TCGA_Ovarian_Cancer\\";
        String outputPath = "D:\\Dropbox\\Results\\PIPI\\TCGA_Ovarian_Cancer\\pipi.peff";
        try {
            new OutputPeff(parameterPath, pipiDirPath, outputPath);
        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }
        logger.info("Done!");
    }

    private OutputPeff(String parameterPath, String pipiDirPath, String outputPath) throws Exception {
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        BuildIndex buildIndex = new BuildIndex(parameterMap, "N14", false, false, true);
        MassTool massTool = buildIndex.returnMassTool();
        Map<Character, Double> massTable = massTool.getMassTable();
        Map<String, Peptide0> peptideProtein0Map = buildIndex.getPeptide0Map();
        DbTool dbTool = buildIndex.getDbTool();
        Map<String, String> proteinSequenceMap = dbTool.getProteinSequenceMap();
        Map<Character, Double> fixModMap = buildIndex.returnFixModMap();

        // recode all variable and fix modification so that these won't be included in the PEFF.
        Multimap<Character, Double> siteVarFixModMap = HashMultimap.create();
        for (String k : parameterMap.keySet()) {
            if (k.startsWith("mod") && !parameterMap.get(k).startsWith("0.0@")) {
                String[] tempArray = parameterMap.get(k).split("@");
                siteVarFixModMap.put(tempArray[1].trim().charAt(0), Double.valueOf(tempArray[0].trim()));
            } else if (k.startsWith("Nterm") && !parameterMap.get(k).contentEquals("0.0")) {
                String[] tempArray = parameterMap.get(k).split("#")[0].trim().split(",");
                for (String s : tempArray) {
                    siteVarFixModMap.put('n', Double.valueOf(s.trim()));
                }
            } else if (k.startsWith("Cterm") && !parameterMap.get(k).contentEquals("0.0")) {
                String[] tempArray = parameterMap.get(k).split("#")[0].trim().split(",");
                for (String s : tempArray) {
                    siteVarFixModMap.put('c', Double.valueOf(s.trim()));
                }
            }
        }
        for (char site : fixModMap.keySet()) {
            if (Math.abs(fixModMap.get(site)) > 0.01) {
                siteVarFixModMap.put(site, fixModMap.get(site));
            }
        }

        // record identified PTM-containing peptides
        Set<String> peptideSet = new HashSet<>();
        Stream<Path> pathStream = Files.find(Paths.get(pipiDirPath), 9999, (p, q) -> q.isRegularFile() && p.getFileName().toString().endsWith("pipi.csv"));
        Iterator<Path> pathIterator = pathStream.iterator();
        int processedFileNum = 0;
        while (pathIterator.hasNext()) {
            Path pipiFilePath = pathIterator.next();
            logger.info("Analyzing {}...", pipiFilePath.toString());
            String line;
            BufferedReader reader = new BufferedReader(new FileReader(pipiFilePath.toFile()));
            ++processedFileNum;
            while ((line = reader.readLine()) != null) {
                if (!line.isEmpty() && !line.startsWith("scan_num")) {
                    String[] parts = Utilities.csvSplitPattern.split(line.trim());
                    if (Double.valueOf(parts[11]) <= 0.01 && !parts[6].trim().contentEquals("-") && Double.valueOf(parts[6]) >= aScoreT) {
                        String peptide = parts[1].trim();
                        peptide = eliminateKnownMods(peptide, siteVarFixModMap);
                        if (peptide.contains("(")) {
                            peptideSet.add(peptide);
                        }
                    }
                }
            }
        }

        logger.info("Processed {} files.", processedFileNum);

        if (!peptideSet.isEmpty()) {
            // generate protein location mod table and protein location AAS table
            Multimap<Character, ModEntry> siteModMap = InferPTM.readUnimodAndGenerateAAS(-1000, 1000);
            Table<String, Integer, Set<String>> proteinLocationModTable = HashBasedTable.create();
            Table<String, Integer, Set<Character>> proteinLocationAASTable = HashBasedTable.create();
            for (String peptide : peptideSet) {
                AA[] aaArray = MassTool.seqToAAList(peptide);
                for (String protein : peptideProtein0Map.get(DbTool.getPtmFreePeptide(peptide)).proteins) {
                    if (proteinSequenceMap.containsKey(protein)) {
                        String proteinSequence = proteinSequenceMap.get(protein);
                        Set<Integer> peptideLocationSet = DbTool.findPeptideLocation(proteinSequence, peptide, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1")); // FixMe: Only consider the first enzyme if the users specify two enzymes.
                        for (int peptideLocation : peptideLocationSet) {
                            for (int i = 1; i < aaArray.length - 1; ++i) {
                                AA aa = aaArray[i];
                                if (Math.abs(aa.ptmDeltaMass) > 0) {
                                    String ptmName = getModString(aa.aa, aa.ptmDeltaMass, siteModMap);
                                    Character newAA = getAAS(aa.aa, aa.ptmDeltaMass, massTable);
                                    if (ptmName != null) {
                                        if (proteinLocationModTable.contains(protein, peptideLocation + i)) {
                                            proteinLocationModTable.get(protein, peptideLocation + i).add(ptmName);
                                        } else {
                                            Set<String> tempSet = new HashSet<>();
                                            tempSet.add(ptmName);
                                            proteinLocationModTable.put(protein, peptideLocation + i, tempSet);
                                        }
                                    } else if (newAA != null) { // we prefer mod to make sure that the amino acid substitutions do not have ambiguous.
                                        if (proteinLocationAASTable.contains(protein, peptideLocation + i)) {
                                            proteinLocationAASTable.get(protein, peptideLocation + i).add(newAA);
                                        } else {
                                            Set<Character> tempSet = new HashSet<>();
                                            tempSet.add(newAA);
                                            proteinLocationAASTable.put(protein, peptideLocation + i, tempSet);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // write PEFF
            Map<String, String> proteinAnnotationMap = dbTool.getProteinAnnotateMap();
            BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath));
            writer.write("# PEFF 1.0\n");
            writer.write("# //\n");
            writer.write(String.format(Locale.US, "# DbName=%s\n", (Paths.get(outputPath)).getFileName().toString()));
            writer.write(String.format(Locale.US, "# DbSource=pipi processed %s\n", parameterMap.get("db")));
            writer.write("# Prefix=\n");
            writer.write(String.format(Locale.US, "# NumberOfEntries=%d\n", proteinSequenceMap.size()));
            writer.write(("# SequenceType=AA\n"));
            writer.write(("# Decoy=false\n"));
            writer.write("# //\n");
            for (String protein : proteinAnnotationMap.keySet()) {
                writer.write(String.format(Locale.US, ">%s \\DbUniqueId=%s \\PName=%s \\Length=%d", protein, protein, proteinAnnotationMap.get(protein), proteinSequenceMap.get(protein).length()));
                if (proteinLocationModTable.containsRow(protein)) {
                    writer.write(" \\ModResUnimod=");
                    for (int location : proteinLocationModTable.row(protein).keySet()) {
                        for (String ptmString : proteinLocationModTable.get(protein, location)) {
                            writer.write(String.format(Locale.US, "(%d|%s)", location, ptmString));
                        }
                    }
                }
                if (proteinLocationAASTable.containsRow(protein)) {
                    writer.write(" \\VariantSimple=");
                    for (int location : proteinLocationAASTable.row(protein).keySet()) {
                        for (char aas : proteinLocationAASTable.get(protein, location)) {
                            writer.write(String.format(Locale.US, "(%d|%c)", location, aas));
                        }
                    }
                }
                writer.write("\n");
                writer.write(proteinSequenceMap.get(protein) + "\n");
            }
            writer.close();
        } else {
            logger.warn("There is no useful PTM-containint peptides.");
        }
    }

    private String eliminateKnownMods(String peptide, Multimap<Character, Double> siteVarFixModMap) { // caution: we assume that the delta mass precision is .001
        String outputPeptide = peptide;
        for (char site : siteVarFixModMap.keySet()) {
            for (double ptmDeltaMass : siteVarFixModMap.get(site)) {
                outputPeptide = outputPeptide.replaceAll(String.format(Locale.US, "%c\\(%.3f\\)", site, ptmDeltaMass), String.valueOf(site));
            }
        }
        return outputPeptide;
    }

    private static String getModString(char site, double ptmDeltaMass, Multimap<Character, ModEntry> siteModMap) {
        if (siteModMap.containsKey(site)) {
            for (ModEntry modEntry : siteModMap.get(site)) {
                if (Math.abs(modEntry.mass - ptmDeltaMass) < 0.01) {
                    return String.format(Locale.US, "%s|%s", modEntry.id, modEntry.name); // we randomly report one mod since there is no way to distinguish mods with the same mass.
                }
            }
        } else {
            logger.error("There is no site {} in the modification map.", site);
        }
        return null;
    }

    private static Character getAAS(char site, double ptmDeltaMass, Map<Character, Double> massTable) {
        double newMass = massTable.get(site) + ptmDeltaMass;
        for (char newAA : massTable.keySet()) {
            if (newAA >= 'A' && newAA <= 'Z' && Math.abs(massTable.get(newAA) - newMass) < 0.01) { // there are n, c, #, and $ in massTable.
                return newAA; // caution: we don't distinguish I and L
            }
        }
        return null;
    }
}
