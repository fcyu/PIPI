# PIPI
PTM-Invariant Peptide Identification.

## How to use it?
Requirement: 
- Java 1.8 or later.
- With [Percolator](https://github.com/percolator/percolator/releases) installed.

Usage:
```
java -Xmx25g -jar PIPI.jar <parameter_file> <spectra_file>
```
- ```<parameter_file>```: Parameter file. Can be downloaded along with PIPI.
- ```<spectra_file>```: Spectra file (mzXML or mgf).

example: ```java -Xmx25g -jar PIPI.jar parameter.def data.mzXML```

## Cite
Yu, Fengchao, Ning Li, and Weichuan Yu. "PIPI: PTM-Invariant Peptide Identification Using Coding Method" (under review).
