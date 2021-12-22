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


import ProteomicsLibrary.Types.SparseBooleanVector;

public class Peptide0 {

    public final SparseBooleanVector code;
    public final boolean isTarget;
    public final String[] proteins;
    public final char leftFlank;
    public final char rightFlank;

    public Peptide0(SparseBooleanVector code, boolean isTarget, String[] proteins, char leftFlank, char rightFlank) {
        this.code = code;
        this.isTarget = isTarget;
        this.proteins = proteins;
        this.leftFlank = leftFlank;
        this.rightFlank = rightFlank;
    }
}
