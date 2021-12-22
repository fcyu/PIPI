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

package proteomics.Search;


import ProteomicsLibrary.Binomial;
import ProteomicsLibrary.Score;
import proteomics.Spectrum.PreSpectra;
import proteomics.Types.Peptide;

import java.util.*;

public class CalSubscores {

    public CalSubscores(Peptide peptide, double ms2Tolerance, TreeMap<Double, Double> plMap, int precursorCharge, TreeSet<Peptide> ptmPatterns, Binomial binomial) throws Exception {
        peptide.setIonFrac(Score.calIonFraction(peptide.getIonMatrix(), precursorCharge, plMap, ms2Tolerance));
        peptide.setMatchedHighestIntensityFrac(Score.calMatchedHighestIntensityFraction(peptide.getIonMatrix(), precursorCharge, plMap, ms2Tolerance));
        peptide.setExplainedAaFrac(Score.calExplainedAAFraction(peptide.getIonMatrix(), precursorCharge, plMap, ms2Tolerance));

        // calculate A score
        if (peptide.hasVarPTM()) {
            Peptide[] tempArray = ptmPatterns.toArray(new Peptide[0]);
            peptide.setaScore(String.valueOf(Score.calAScore(plMap, PreSpectra.topN, binomial, peptide.getVarPTMs(), peptide.getIonMatrix(), tempArray.length > 1 ? tempArray[1].getVarPTMs() : null, tempArray.length > 1 ? tempArray[1].getIonMatrix() : null, ms2Tolerance, peptide.length())));
        }
    }
}
