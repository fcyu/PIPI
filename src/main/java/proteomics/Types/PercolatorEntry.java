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
