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


public class VarModParam {

    public final double mass;
    public final char aa;
    public final int priority; // 1 = high; 0 = low.
    public final boolean onlyProteinTerminalIfnc;

    private final int hashCode;

    public VarModParam(double mass, char aa, int priority, boolean onlyProteinTerminalIfnc) {
        this.mass = mass;
        this.aa = aa;
        this.priority = priority;
        this.onlyProteinTerminalIfnc = onlyProteinTerminalIfnc;

        String toString = Math.round(mass * 1000) + "@" + aa;
        hashCode = toString.hashCode();
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof VarModParam) {
            VarModParam temp = (VarModParam) other;
            return temp.hashCode == hashCode;
        } else {
            return false;
        }
    }
}
