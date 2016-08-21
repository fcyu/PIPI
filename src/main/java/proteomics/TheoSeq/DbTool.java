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

    public DbTool(String dbName) {
        String id = "";
        String annotate;
        String seq = "";

        boolean newPro = true;

        Pattern headerPattern = Pattern.compile(">([^\\s]*)(.*)");

        try (BufferedReader dbReader = new BufferedReader(new FileReader(dbName))) {
            String line;
            while ((line = dbReader.readLine()) != null) {
                line = line.trim();
                Matcher headMatcher = headerPattern.matcher(line);
                if (headMatcher.matches()) {
                    // This line is a header
                    if (!newPro) {
                        // This isn't the first protein
                        proSeqMap.put(id, seq);
                    }
                    id = headMatcher.group(1).trim();
                    annotate = headMatcher.group(2).trim();
                    proAnnotateMap.put(id, annotate);
                    newPro = true;
                } else if (!line.isEmpty()) {
                    // This line is a body
                    if (newPro) {
                        seq = line;
                        newPro = false;
                    } else {
                        seq += line;
                    }
                }
            }
            // Last protein
            proSeqMap.put(id, seq);
        } catch (IOException | PatternSyntaxException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
    }

    public Map<String, String> returnSeqMap() {
        return proSeqMap;
    }

    public Map<String, String> returnAnnotateMap() {
        return proAnnotateMap;
    }
}
