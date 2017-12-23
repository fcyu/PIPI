package proteomics.Parameter;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.PIPI;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class Parameter {

    private static final Logger logger = LoggerFactory.getLogger(Parameter.class);
    private static final Pattern commentLinePattern = Pattern.compile("^#.*");
    private static final Pattern linePattern = Pattern.compile("([^#]+)=([^#]+)#*.*");
    private static final Pattern enzymePattern = Pattern.compile("(.+)\\s+([01])\\s+([A-Z]+)\\s+([A-Z\\-]+)");

    private Map<String, String> parameterMap = new LinkedHashMap<>();

    public Parameter(String parameterFile) {
        try (BufferedReader parameterReader = new BufferedReader(new FileReader(parameterFile))) {
            String line = parameterReader.readLine().trim();
            if (!line.contentEquals("# " + PIPI.versionStr)) {
                logger.error("The parameter file version ({}) is not compatible with current PIPI version ({}).", line.substring(2), PIPI.versionStr);
                System.exit(1);
            }
            while ((line = parameterReader.readLine()) != null) {
                line = line.trim();
                Matcher commentLineMatcher = commentLinePattern.matcher(line);
                if (!commentLineMatcher.matches()) {
                    // This is not a comment line
                    Matcher lineMatcher = linePattern.matcher(line);
                    if (lineMatcher.matches()) {
                        String parameterName = lineMatcher.group(1).trim();
                        String parameterValue = lineMatcher.group(2).trim();
                        parameterMap.put(parameterName, parameterValue);
                    } else {
                        Matcher enzymeMatcher = enzymePattern.matcher(line);
                        if (enzymeMatcher.matches()) {
                            parameterMap.put("enzyme_name", enzymeMatcher.group(1).trim());
                            parameterMap.put("cleavage_from_c_term", enzymeMatcher.group(2).trim());
                            parameterMap.put("cleavage_site", enzymeMatcher.group(3).trim());
                            parameterMap.put("protection_site", enzymeMatcher.group(4).trim());
                        }
                    }
                }
            }
        } catch (IOException | IllegalStateException ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
            System.exit(1);
        }
    }


    public Map<String, String> returnParameterMap() {
        return parameterMap;
    }
}