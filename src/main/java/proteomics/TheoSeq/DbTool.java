package proteomics.TheoSeq;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class DbTool {

    private static final Logger logger = LoggerFactory.getLogger(DbTool.class);

    private Map<String, String> proSeqMap = new HashMap<>();
    private Map<String, String> proAnnotateMap = new HashMap<>();

    public DbTool(String dbName, String databaseType) throws IOException {
        String id = "";
        String annotate;
        StringBuilder seq = new StringBuilder(99999);

        boolean newPro = true;

        Pattern headerPattern;
        if (databaseType.contentEquals("TAIR")) {
            headerPattern = Pattern.compile("^>([^\\s]+)[\\s|]+(.+)$");
        } else if (databaseType.contentEquals("UniProt") || databaseType.contentEquals("SwissProt")) {
            headerPattern = Pattern.compile("^>[^|]+\\|(.+)\\|(.+)$");
        } else if (databaseType.contentEquals("contaminants")) {
            headerPattern = Pattern.compile("^>([^ ]+) (.+)$");
        } else {
            headerPattern = null;
            logger.error("Incorrect database type ({}) in the parameter file.", databaseType);
            System.exit(1);
        }

        BufferedReader dbReader;
        if (databaseType.contentEquals("contaminants")) {
            InputStream inputStream = getClass().getClassLoader().getResourceAsStream("contaminants.fasta");
            dbReader = new BufferedReader(new InputStreamReader(inputStream));
        } else {
            dbReader = new BufferedReader(new FileReader(dbName));
        }
        String line;
        while ((line = dbReader.readLine()) != null) {
            line = line.trim();
            Matcher headMatcher = headerPattern.matcher(line);
            if (headMatcher.matches()) {
                // This line is a header
                if (!newPro) {
                    // This isn't the first protein
                    proSeqMap.put(id, seq.toString());
                }
                id = headMatcher.group(1).trim();
                annotate = headMatcher.group(2).trim();
                proAnnotateMap.put(id, annotate);
                newPro = true;
            } else if (!line.isEmpty()) {
                // This line is a body
                if (newPro) {
                    seq = new StringBuilder(99999);
                    seq.append(line);
                    newPro = false;
                } else {
                    seq.append(line);
                }
            }
        }
        dbReader.close();
        // Last protein
        proSeqMap.put(id, seq.toString());
    }

    public Map<String, String> returnSeqMap() {
        return proSeqMap;
    }

    public Map<String, String> returnAnnotateMap() {
        return proAnnotateMap;
    }

    public Set<Integer> findPeptideLocation(String proteinId, String peptide) {
        peptide = peptide.trim().replaceAll("[^A-Z]+", "");
        Set<Integer> output = new HashSet<>();
        int idx = proSeqMap.get(proteinId).indexOf(peptide);
        while (idx >= 0) {
            output.add(idx);
            idx = proSeqMap.get(proteinId).indexOf(peptide, idx + 1);
        }
        if (!output.isEmpty()) {
            return output;
        } else {
            throw new NullPointerException(String.format(Locale.US, "Cannot find the peptide %s from the protein %s.", peptide, proteinId));
        }
    }

    public static Set<String> reduceProteinIdSet(Set<String> input) {
        if (input.size() == 1) {
            return input;
        } else {
            Map<String, Integer> tempMap = new HashMap<>();
            for (String s : input) {
                String[] tempArray = s.split("\\.");
                if (tempMap.containsKey(tempArray[0])) {
                    if (tempMap.get(tempArray[0]) > Integer.valueOf(tempArray[1])) {
                        tempMap.put(tempArray[0], Integer.valueOf(tempArray[1]));
                    }
                } else {
                    tempMap.put(tempArray[0], Integer.valueOf(tempArray[1]));
                }
            }
            Set<String> output = new HashSet<>();
            for (String s : tempMap.keySet()) {
                output.add(s + "." + tempMap.get(s));
            }
            return output;
        }
    }
}
