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


import ProteomicsLibrary.Types.Coordinate;

import java.util.Locale;
import java.util.TreeMap;

public class PositionDeltaMassMap extends TreeMap<Coordinate, Double> {

    public final int peptideLength;

    public PositionDeltaMassMap(int peptideLength) {
        super();
        this.peptideLength = peptideLength;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder(1000);
        for (Coordinate co : this.keySet()) {
            sb.append(String.format(Locale.US, "%.3f", this.get(co)));
            sb.append("@");
            if ((co.x == 0) || (co.x == 1)) {
                sb.append("([01]-");
                sb.append(co.y);
                sb.append(")");
            } else if ((co.y == peptideLength) || (co.y == peptideLength - 1)){
                sb.append("(");
                sb.append(co.x);
                sb.append("-[");
                sb.append(peptideLength - 1);
                sb.append(peptideLength);
                sb.append("])");
            } else {
                sb.append(co.toString());
            }
            sb.append(";");
        }
        return sb.toString();
    }

    public int hashCode() {
        return this.toString().hashCode();
    }

    public boolean equals(Object other) {
        if (other instanceof PositionDeltaMassMap) {
            PositionDeltaMassMap temp = (PositionDeltaMassMap) other;
            return temp.hashCode() == this.hashCode();
        } else {
            return false;
        }
    }

    public PositionDeltaMassMap clone() {
        super.clone();
        PositionDeltaMassMap other = new PositionDeltaMassMap(peptideLength);
        other.clear();
        for (Coordinate co : this.keySet()) {
            other.put(co, this.get(co));
        }
        return other;
    }
}