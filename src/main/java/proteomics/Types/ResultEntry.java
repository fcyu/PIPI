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

public class ResultEntry implements Comparable<ResultEntry> {

    public final double score;
    public final String peptide;
    private final int hashCode;
    private final boolean isDecoy;

    public ResultEntry(double score, String peptide, boolean isDecoy) {
        this.score = score;
        this.peptide = peptide;
        String toString = peptide + "-" + score;
        hashCode = toString.hashCode();
        this.isDecoy = isDecoy;
    }

    public int compareTo(ResultEntry other) {
        return Double.compare(score, other.score);
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof ResultEntry) {
            ResultEntry temp = (ResultEntry) other;
            return ((temp.peptide.contentEquals(peptide)) && (temp.score == score));
        } else {
            return false;
        }
    }
}