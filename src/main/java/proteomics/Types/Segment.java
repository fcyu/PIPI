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


public class Segment implements Comparable<Segment> {

    private final String segmentString;
    private final int hashCode;

    public Segment(String segmentString) {
        String temp = reverseString(segmentString);
        // Compare the string and its reversed version. Kept the smaller one.
        int compareResult = segmentString.compareTo(temp);
        if (compareResult > 0) {
            this.segmentString = temp;
        } else {
            this.segmentString = segmentString;
        }

        hashCode = segmentString.hashCode();
    }

    public int length() {
        return segmentString.length();
    }

    public  boolean equals(Object other) {
        return (other instanceof Segment) && this.segmentString.contentEquals(((Segment) other).segmentString);
    }

    public String toString() {
        return segmentString;
    }

    public int hashCode() {
        return hashCode;
    }

    public int compareTo(Segment other) {
        return this.segmentString.compareTo(other.segmentString);
    }

    private String reverseString(String seq) {
        return new StringBuilder(seq).reverse().toString();
    }
}
