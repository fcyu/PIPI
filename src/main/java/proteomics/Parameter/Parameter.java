/*
 * Copyright 2016-2019 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package proteomics.Parameter;

import proteomics.PIPI;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class Parameter {

    private static final Pattern commentLinePattern = Pattern.compile("^#.*");
    private static final Pattern linePattern = Pattern.compile("([^#]+)=([^#]+)#*.*");
    private static final Pattern enzymePattern = Pattern.compile("(.+)\\s+([01])\\s+([A-Z]+)\\s+([A-Z\\-]+)");

    private Map<String, String> parameterMap = new LinkedHashMap<>();

    public Parameter(String parameterFile) throws IOException {
        BufferedReader parameterReader = new BufferedReader(new FileReader(parameterFile));
        String line = parameterReader.readLine().trim();
        if (!line.contentEquals("# " + PIPI.versionStr)) {
            throw new IOException(String.format(Locale.US, "The parameter file version (%s) is not compatible with current PIPI version (%s).", line.substring(2), PIPI.versionStr));
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
                        if (parameterMap.containsKey("enzyme_name_1")) {
                            parameterMap.put("enzyme_name_2", enzymeMatcher.group(1).trim());
                            parameterMap.put("is_from_C_term_2", enzymeMatcher.group(2).trim());
                            parameterMap.put("cleavage_site_2", enzymeMatcher.group(3).trim());
                            parameterMap.put("protection_site_2", enzymeMatcher.group(4).trim());
                        } else {
                            parameterMap.put("enzyme_name_1", enzymeMatcher.group(1).trim());
                            parameterMap.put("is_from_C_term_1", enzymeMatcher.group(2).trim());
                            parameterMap.put("cleavage_site_1", enzymeMatcher.group(3).trim());
                            parameterMap.put("protection_site_1", enzymeMatcher.group(4).trim());
                        }
                    }
                }
            }
        }
        parameterReader.close();
    }


    public Map<String, String> returnParameterMap() {
        return parameterMap;
    }
}