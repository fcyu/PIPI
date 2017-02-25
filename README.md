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
Yu, F., Li, N., & Yu, W. (2016). PIPI: PTM-Invariant Peptide Identification Using Coding Method. Journal of Proteome Research, 15(12), 4423-4435.
```
@article{yu2016pipi,
  title={PIPI: PTM-Invariant Peptide Identification Using Coding Method},
  author={Yu, Fengchao and Li, Ning and Yu, Weichuan},
  journal={Journal of Proteome Research},
  volume={15},
  number={12},
  pages={4423--4435},
  year={2016},
  publisher={ACS Publications}
}
```
