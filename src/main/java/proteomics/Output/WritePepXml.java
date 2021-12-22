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

package proteomics.Output;

import ProteomicsLibrary.Types.Coordinate;
import proteomics.PIPI;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.AA;
import proteomics.Types.*;

import java.io.*;
import java.sql.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Date;

public class WritePepXml {

    private final String outputPath;
    private final String baseName;
    private final String rawDataType;
    private final Map<String, String> parameterMap;

    public WritePepXml(String outputPath, String spectraName, Map<String, String> parameterMap, Map<Character, Double> massTable, Map<Integer, PercolatorEntry> percolatorResultMap, Map<String, Peptide0> peptide0Map, Map<Character, Double> fixModMap, String sqlPath) throws IOException, SQLException {
        this.outputPath = outputPath;
        int tempIdx = spectraName.lastIndexOf('.');
        baseName = spectraName.substring(0, tempIdx);
        rawDataType = spectraName.substring(tempIdx);
        this.parameterMap = parameterMap;

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath));
        writer.write(pepxmlHeader(massTable));
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass, isDecoy, score, matchedPeakNum, otherPtmPatterns, aScore FROM spectraTable");
        while (sqlResultSet.next()) {
            int scanNum = sqlResultSet.getInt("scanNum");
                String peptide = sqlResultSet.getString("peptide");
            if ((percolatorResultMap == null || percolatorResultMap.containsKey(scanNum)) && !sqlResultSet.wasNull()) {
                String ptmFreePeptide = peptide.replaceAll("[^ncA-Z]+", "");
                Peptide0 peptide0 = peptide0Map.get(ptmFreePeptide);
                TreeSet<String> proteinIdSet = new TreeSet<>();
                for (String protein : peptide0.proteins) {
                    proteinIdSet.add(protein.trim());
                }
                double expMass = sqlResultSet.getDouble("precursorMass");
                String aScore = sqlResultSet.getString("aScore");
                PercolatorEntry percolatorEntry = null;
                if (percolatorResultMap != null) {
                    percolatorEntry = percolatorResultMap.get(scanNum);
                }
                int precursorCharge = sqlResultSet.getInt("precursorCharge");
                double theoMass = sqlResultSet.getDouble("theoMass");

                writer.write(String.format(Locale.US,
                        "\t\t<spectrum_query spectrum=\"%d\" start_scan=\"%d\" end_scan=\"%d\" precursor_neutral_mass=\"%f\" assumed_charge=\"%d\" index=\"%d\">\r\n" +
                                "\t\t\t<search_result>\r\n" +
                                "\t\t\t\t<search_hit hit_rank=\"1\" peptide=\"%s\" peptide_prev_aa=\"%c\" peptide_next_aa=\"%c\" protein=\"%s\" num_tot_proteins=\"%d\" num_matched_ions=\"%d\" tot_num_ions=\"%d\" calc_neutral_pep_mass=\"%f\" massdiff=\"%f\" num_tol_term=\"2\">\r\n" +
                                "\t\t\t\t\t<search_score name=\"score\" value=\"%f\"/>\r\n" +
                                "\t\t\t\t\t<search_score name=\"A_score\" value=\"%s\"/>\r\n" +
                                "\t\t\t\t\t<search_score name=\"percolator_score\" value=\"%f\"/>\r\n" +
                                "\t\t\t\t\t<search_score name=\"percolator_error_prob\" value=\"%s\"/>\r\n" +
                                "\t\t\t\t\t<search_score name=\"q_value\" value=\"%s\"/>\r\n" +
                                "\t\t\t\t\t<search_score name=\"labelling\" value=\"%s\"/>\r\n", scanNum, scanNum, scanNum, expMass, precursorCharge, scanNum, ptmFreePeptide.replaceAll("[nc]+", ""), peptide0.leftFlank, peptide0.rightFlank, String.join(";", proteinIdSet), peptide0.proteins.length, sqlResultSet.getInt("matchedPeakNum"), (ptmFreePeptide.length() - 2) * 2 * Math.max(1, precursorCharge - 1), theoMass, PIPI.getMassDiff(expMass, theoMass, MassTool.C13_DIFF), sqlResultSet.getDouble("score"), aScore, percolatorEntry == null ? null : percolatorEntry.percolatorScore, percolatorEntry == null ? "null" : percolatorEntry.PEP, percolatorEntry == null ? "null" : percolatorEntry.qValue, sqlResultSet.getString("labelling")));

                if (!aScore.contentEquals("-")) {
                    PositionDeltaMassMap ptmMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                    AA[] aaArray = MassTool.seqToAAList(peptide);
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < aaArray.length; ++i) {
                        if (Math.abs(aaArray[i].ptmDeltaMass) > 0.5 && Math.abs(fixModMap.get(aaArray[i].aa) - aaArray[i].ptmDeltaMass) > 0.1) {
                            ptmMap.put(new Coordinate(i, i + 1), aaArray[i].ptmDeltaMass);
                            sb.append(String.format(Locale.US, "%c(%.3f)", aaArray[i].aa, aaArray[i].ptmDeltaMass));
                        } else {
                            sb.append(aaArray[i].aa);
                        }
                    }
                    if (ptmMap.containsKey(new Coordinate(0, 1)) && ptmMap.containsKey(new Coordinate(ptmMap.peptideLength - 1, ptmMap.peptideLength))) {
                        writer.write(String.format(Locale.US, "\t\t\t\t\t<modification_info modified_peptide=\"%s\" mod_nterm_mass=\"%f\" mod_cterm_mass=\"%f\">\r\n", sb.toString(), ptmMap.get(new Coordinate(0, 1)) + MassTool.PROTON, ptmMap.get(new Coordinate(ptmMap.peptideLength - 1, ptmMap.peptideLength))));
                    } else if (ptmMap.containsKey(new Coordinate(0, 1))) {
                        writer.write(String.format(Locale.US, "\t\t\t\t\t<modification_info modified_peptide=\"%s\" mod_nterm_mass=\"%f\">\r\n", sb.toString(), ptmMap.get(new Coordinate(0, 1)) + MassTool.PROTON));
                    } else if (ptmMap.containsKey(new Coordinate(ptmMap.peptideLength - 1, ptmMap.peptideLength))) {
                        writer.write(String.format(Locale.US, "\t\t\t\t\t<modification_info modified_peptide=\"%s\" mod_cterm_mass=\"%f\">\r\n", sb.toString(), ptmMap.get(new Coordinate(ptmMap.peptideLength - 1, ptmMap.peptideLength))));
                    } else {
                        writer.write(String.format(Locale.US, "\t\t\t\t\t<modification_info modified_peptide=\"%s\">\r\n", sb.toString()));
                    }
                    for (Coordinate co : ptmMap.keySet()) {
                        if (co.x != 0 && co.x != ptmMap.peptideLength - 1) {
                            writer.write(String.format(Locale.US, "\t\t\t\t\t\t<mod_aminoacid_mass position=\"%d\" mass=\"%f\"/>\r\n", co.x, massTable.get(ptmFreePeptide.charAt(co.x)) + ptmMap.get(co)));
                        }
                    }
                    writer.write("\t\t\t\t\t</modification_info>\r\n");
                }
                writer.write("\t\t\t\t</search_hit>\r\n" +
                        "\t\t\t</search_result>\r\n" +
                        " \t\t</spectrum_query>\r\n");
            }
        }
        writer.write("\t</msms_run_summary>\r\n" +
                        "</msms_pipeline_analysis>\r\n");
        writer.close();
        sqlResultSet.close();
        sqlStatement.close();
        sqlConnection.close();
    }

    private String pepxmlHeader(Map<Character, Double> massTable) {
        DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd");
        Date date = new Date();
        StringBuilder header = new StringBuilder();
        header.append(String.format(Locale.US,
                "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\r\n" +
                "<msms_pipeline_analysis date=\"%s\" xmlns=\"http://regis-web.systemsbiology.net/pepXML\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v117.xsd\" summary_xml=\"%s\">\r\n" +
                "\t<msms_run_summary base_name=\"%s\" raw_data_type=\"%s\" raw_data=\"%s\">\r\n" +
                "\t\t<sample_enzyme name=\"%s\">\r\n" +
                "\t\t\t<specificity cut=\"%s\" no_cut=\"%s\" sense=\"%c\"/>\r\n" +
                "\t\t</sample_enzyme>\r\n" +
                "\t\t<search_summary base_name=\"%s\" search_engine=\"PIPI\" search_engine_version=\"%s\" precursor_mass_type=\"monoisotopic\" fragment_mass_type=\"monoisotopic\" search_id=\"1\">\r\n" +
                "\t\t\t<search_database local_path=\"%s\" type=\"AA\"/>\r\n" +
                "\t\t\t<enzymatic_search_constraint enzyme=\"%s\" max_num_internal_cleavages=\"%s\" min_number_termini=\"2\"/>\r\n", dateFormat.format(date), outputPath, baseName, rawDataType, rawDataType, parameterMap.get("enzyme_name_1"), parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), parameterMap.get("is_from_C_term_1").contentEquals("1") ? 'C' : 'N', baseName, PIPI.versionStr, parameterMap.get("db"), parameterMap.get("enzyme_name_1"), parameterMap.get("missed_cleavage"))); // FixMe: Only has the first enzyme's information if the users specified two enzymes.
        for (String k : parameterMap.keySet()) {
            if (k.startsWith("mod") && !parameterMap.get(k).trim().startsWith("0.0")) {
                String[] parts = parameterMap.get(k).trim().split("@");
                header.append(String.format(Locale.US, "\t\t\t<aminoacid_modification aminoacid=\"%c\" massdiff=\"%s\" mass=\"%f\" variable=\"Y\"/>\r\n", parts[1].charAt(0), parts[0].trim(), massTable.get(parts[1].charAt(0)) + Double.valueOf(parts[0])));
            } else if (k.startsWith("Nterm") && Math.abs(Double.valueOf(parameterMap.get(k).trim())) > 0.5) {
                String[] parts = parameterMap.get(k).trim().split(",");
                for (String part : parts) {
                    header.append(String.format(Locale.US, "\t\t\t<terminal_modification terminus=\"N\" massdiff=\"%s\" mass=\"%f\" variable=\"Y\" protein_terminus=\"N\"/>\r\n", part.trim(), MassTool.PROTON + Double.valueOf(parts[0])));
                }
            } else if (k.startsWith("Cterm") && Math.abs(Double.valueOf(parameterMap.get(k).trim())) > 0.5) {
                String[] parts = parameterMap.get(k).trim().split(",");
                for (String part : parts) {
                    header.append(String.format(Locale.US, "\t\t\t<terminal_modification terminus=\"C\" massdiff=\"%s\" mass=\"%s\" variable=\"Y\" protein_terminus=\"N\"/>\r\n", part.trim(), part.trim()));
                }
            } else if (k.matches("[A-Z]") && Math.abs(Double.valueOf(parameterMap.get(k).trim())) > 0.5) {
                header.append(String.format(Locale.US, "\t\t\t<aminoacid_modification aminoacid=\"%c\" massdiff=\"%s\" mass=\"%f\" variable=\"N\"/>\r\n", k.charAt(0), parameterMap.get(k).trim(), massTable.get(k.charAt(0)) + Double.valueOf(parameterMap.get(k))));
            } else if (k.contentEquals("n") && Math.abs(Double.valueOf(parameterMap.get(k).trim())) > 0.5) {
                header.append(String.format(Locale.US, "\t\t\t<terminal_modification terminus=\"N\" massdiff=\"%s\" mass=\"%s\" variable=\"N\" protein_terminus=\"N\"/>\r\n", parameterMap.get(k).trim(), MassTool.PROTON + Double.valueOf(parameterMap.get(k).trim())));
            } else if (k.contentEquals("c") && Math.abs(Double.valueOf(parameterMap.get(k).trim())) > 0.5) {
                header.append(String.format(Locale.US, "\t\t\t<terminal_modification terminus=\"C\" massdiff=\"%s\" mass=\"%s\" variable=\"N\" protein_terminus=\"N\"/>\r\n", parameterMap.get(k).trim(), Double.valueOf(parameterMap.get(k).trim())));
            }
        }
        for (String k : parameterMap.keySet()) {
            if (!k.startsWith("mod") && !k.startsWith("Nterm") && !k.startsWith("Cterm")) {
                header.append(String.format(Locale.US, "\t\t\t<parameter name=\"%s\" value=\"%s\"/>\r\n", k, parameterMap.get(k).trim()));
            }
        }
        header.append("\t\t</search_summary>\r\n");
        return header.toString();
    }
}
