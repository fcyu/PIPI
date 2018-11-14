# PIPI
PTM-Invariant Peptide Identification. An open search tool.

## Executable file:
Download `PIPI.zip` from https://github.com/fcyu/PIPI/releases/latest.

## How to use it?
Requirement: 
- Java 1.8.
- With [Percolator](https://github.com/percolator/percolator/releases) installed.

Usage:
```
java -Xmx64g -jar PIPI.jar <parameter_file> <spectra_file>
```
- ```<parameter_file>```: Parameter file. Can be downloaded along with PIPI. There are detailed explanations for the parameters in the file.
- ```<spectra_file>```: Spectra file (mzXML or mgf).

example: ```java -Xmx25g -jar PIPI.jar parameter.def data.mzXML```

## An example of the result file
| scan_num | peptide                                | charge | theo_mass | exp_mass | abs_ppm  | A_score  | protein_ID                                | score    | delta_C_n | percolator_score | posterior_error_prob | q_value  | other_PTM_patterns                                                                                                                                                                      | MGF_title | labelling | isotope_correction | MS1_pearson_correlation_coefficient |
|----------|----------------------------------------|--------|-----------|----------|----------|----------|-------------------------------------------|----------|-----------|------------------|----------------------|----------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|-----------|--------------------|-------------------------------------|
| 1436     | nESEDKPEIEDVGSDEE(129.043)EE(-0.948)Kc | 3      | 2320.008  | 2320.011 | 1.361049 | 13.02612 | NP_001017963.2;NP_005339.3;XP_011535020.1 | 5.106351 | 0.000121  | 2.69477          | 7.21E-06             | 0.002299 | nESEDKPEIEDVGSDEEE(129.043)E(-0.948)Kc-4.9165;nESEDKPEIEDVGSDEE(129.043)E(-0.948)EKc-4.8888;nESEDKPEIEDVGSDE(129.043)EEE(-0.948)Kc-4.8815;nESEDKPEIEDVGSDEE(-0.948)EE(129.043)Kc-4.6991 |           | N14       | 0                  | 0.897093                            |

### Explanations:
1. `scan_num`: Scan number (start from 1) of the spectrum.
2. `peptide`: Identified peptide sequence.
3. `charge`: Precursor charge.
4. `theo_mass`: Theoretical neutral mass of the identified peptide sequence.
5. `exp_mass`: Observed neutral mass of the spectrum.
6. `abs_ppm`: Absolute PPM between `theo_mass` and `exp_mass`.
7. `A_score`: AScore based on the paper "[A probability-based approach for high-throughput protein phosphorylation analysis and site localization](https://doi.org/10.1038/nbt1240)". In the original algorithm, it only considers phosphorylation. PIPI considers all modifications.
8. `protein_ID`: Protein accession from the fasta file.
9. `score`: Identification score. It's similar to XCorr.
10. `delta_c_n`: The relative difference between the top score and the second best score.
11. `percolator_score`: Percolator score reported by Percolator.
12. `posterior_error_prob`: Percolator error probability reported by Percolator.
13. `q_value`: q-value (used as FDR in the community) estimated by Percolator.
14. `other_PTM_patterns`: The second best to the fifth best modification patterns corresponding to the identified peptide sequence. There is an identification score following each pattern.
15. `MGF_title`: MGF title if the input spectral file is in MGF format.
16. `labelling`: N14 or N15 labelling.
17. `isotope_correction`: The number of C13 isotope correction. Only applicable to mzXML format.
18. `MS1_pearson_correlation_coefficient`: The similarity between the observed MS1 isotope pattern and the theoretical isotope pattern. The range is from 0 to 1. Only applicable to mzXML format.

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
  publisher={ACS Publications},
  doi={10.1021/acs.jproteome.6b00485}
}
```
