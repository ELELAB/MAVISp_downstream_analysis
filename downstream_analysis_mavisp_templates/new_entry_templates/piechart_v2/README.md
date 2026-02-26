## MAVISp pie chart
Last updated: 25/02/2026

### DESCRIPTION
The script `piechart_v2.py` takes as input a MAVISp CSV table and a MAVISp internal dictionary file
for ClinVar records, and generates a pie chart summarizing variant classifications.

It supports multiple ClinVar classification types via `--clinvar_class_type`:
```
  - aggregated (default)
  - germline
  - oncogenicity
  - clinical_impact
```
The script checks that the required ClinVar columns for the selected mode are present in the CSV,
then applies the provided dictionary mapping before plotting.

### REQUIREMENTS
```
  - python >= 3.7
  - pandas
  - matplotlib
  - seaborn
  - argparse
  - argcomplete
  ```

### INPUTS
  - MAVISp CSV (aggregated or “new style” entry CSV containing the selected ClinVar columns)
  - ClinVar internal dictionary file (tab-separated), matching the selected `--clinvar_class_type`

### USAGE
  python piechart_v2.py [`-h`] `-i` INPUT `-d` DICTIONARY `-o` OUT_NAME `-cct` {aggregated, germline, oncogenicity, clinical_impact}

### ARGUMENTS
  `-i`, `--input`
      Input MAVISp CSV.

  `-d`, `--dictionary`
      Input MAVISp internal dictionary file corresponding to the selected ClinVar classification type.

  `-o`, `--out_name`
      Output plot filename (e.g. .pdf or .png).

  `-cct`, `--clinvar_class_type`
      ClinVar classification type to use.
      Choices: aggregated, germline, oncogenicity, clinical_impact
      Default: aggregated (suitable for old style entries)

### NOTES
  - Provide the correct dictionary file for the selected `--clinvar_class_type`; the script validates
    this by checking the first cell of the dictionary file.

### EXAMPLES
EXAMPLES

**OLD-STYLE entry** command line:

1. Aggregated (default)
  python piechart_v2.py -i oldmavisp.csv -d clinvar_interpretation_internal_dictionary.txt -o pie.pdf


**NEW-STYLE entry** command line:

1. Germline classification
  python piechart_v2.py -i newmavisp.csv -cct germline -d clinvar_interpretation_internal_dictionary.txt -o pie.pdf

2. Oncogenicity classification
  python piechart_v2.py -i newmavisp.csv -cct oncogenicity -d oncodict.csv -o pie.pdf

3. Clinical impact classification
  python piechart_v2.py -i newmavisp.csv -cct clinical_impact -d clinicaldict.csv -o pie.pdf

### DICTIONARY FILE
  The internal MAVISp ClinVar dictionary is available in the official MAVISp GitHub repository:
  https://github.com/ELELAB/MAVISp/blob/main/mavisp/data/clinvar_interpretation_internal_dictionary.txt