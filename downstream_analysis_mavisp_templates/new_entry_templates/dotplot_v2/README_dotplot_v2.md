# MAVISp dot plot v2

*Last updated*: 25/02/26

## Description

The script `dot_plot_v2.py` takes in input the MAVISp aggregated table and returns a dot plot plus a summary of the effect of the mutations on the single modules of the MAVISp framework, dividing them into structural/functional effects.

Additionally, it writes a CSV file (`mechanistic_indicators_out.csv`) containing mutations with at least one identified module effect and a consolidated `MAVISp Effects` column. Optionally, this table can be filtered to keep only mutations predicted as damaging/pathogenic according to a selected VEP (AlphaMissense, REVEL, GEMME, EVE, or none) via `-vep`, and can be further restricted to mutations that are LoF or GoF for DeMaSk via `-lgof`.

Compared to `dot_plot.py`, `dot_plot_v2.py` extends ClinVar handling by supporting alternative ClinVar interpretation columns through `-cct/--clinvar_class_type`.

## Requirements

- Python >= 3.7
- matplotlib
- numpy
- pandas

## Inputs

- MAVISp CSV
- dictionary.csv|oncodict.csv|clinicaldict.csv (ClinVar annotation dictionary; required only when `-pltC`/`-colC` is used, and must match the selected ClinVar class type)

## Usage

`python dot_plot_v2.py [-h] -i INPUT [-v CLINVAR_DICTIONARY] [-o OUTPUT] [-m MUTATIONS [MUTATIONS ...]] [-r RESIDUES [RESIDUES ...]] [-R REVEL_THRESHOLD] [-D DEMASK_THRESHOLD] [-G GEMME_THRESHOLD] [-x X_LIM] [-f FIGSIZE FIGSIZE] [-pltR] [-pltD]
[-pltC {all,uncertain,benign,likely_benign,pathogenic,likely_pathogenic,conflicting,oncogenic,likely_oncogenic,strong,potential,unknown} [{all,uncertain,benign,likely_benign,pathogenic,likely_pathogenic,conflicting,oncogenic,likely_oncogenic,strong,potential,unknown} ...]]
[-colC] [-pltS {saturation,cosmic,cbioportal} [{saturation,cosmic,cbioportal} ...]] [-vep {none,alphamissense,revel,gemme,eve}] [-lgof] -cct {aggregated,germline,oncogenicity,clinical_impact}`

- `-i/--input` : input MAVISp aggregated CSV (required)
- `-v/--clinvar-dictionary` : path to ClinVar dictionary file. Only needed when using `-pltC` or `-colC` (default: dictionary.csv)
- `-o/--output` : plot(s) output base filename (default: dot_plot)
- `-m/--mutations` : selected mutations to be plotted to be provided comma-separated (e.g., A4G,F55K). Mutually exclusive with `-r`
- `-r/--residues` : selected residues to be plotted to be provided comma-separated (e.g., 4,55). Mutually exclusive with `-m`
- `-R/--revel_threshold` : REVEL score threshold (default: 0.5)
- `-D/--demask_threshold` : threshold to classify a mutation according to the DeMaSk delta-fitness score (default: 0.25)
- `-G/--gemme_threshold` : threshold to classify a mutation according to the GEMME score (default: -3.0)
- `-x/--x_lim` : x-axis limit (number of mutations per panel before splitting across multiple plots) (default: 50)
- `-f/--figsize` : figure size (default: 14 5). Suggested for ~50 mutations on x-axis and 7–8 labels on y-axis
- `-pltR/--plot_Revel` : include REVEL classifications in the plot (default: None)
- `-pltD/--plot_Demask` : include DeMaSk predicted consequence (LoF/GoF) for mutations satisfying the `-D` threshold (default: None)
- `-pltC/--plot_Clinvar` : select mutations to be plotted according to the ClinVar interpretation (e.g. pathogenic, uncertain); multiple categories can be included. Requires a dictionary matching `-cct` (default: None)
- `-colC/--color_Clinvar` : color x-tick labels according to the ClinVar interpretation. Requires a dictionary matching `-cct` (default: None)
- `-pltS/--plot_Source` : select mutations to be plotted according to source (saturation, cosmic, cbioportal); multiple sources can be included (default: None)
- `-vep/--vep-filter` : filter `mechanistic_indicators_out.csv` to keep only mutations predicted as damaging/pathogenic by the selected VEP (alphamissense, revel, gemme, eve, none) (default: none)
- `-lgof/--vep-filter-lgof` : when used with `-vep`, further restrict `mechanistic_indicators_out.csv` to DeMaSk LoF/GoF only (default: None)
- `-cct/--clinvar_class_type` : selects which ClinVar interpretation column is used: {aggregated,germline,oncogenicity,clinical_impact} (default is aggregated and can be applied only on old-style entries)

## ClinVar dictionaries and `-cct`

When using `-pltC` or `-colC`, the correct dictionary must be supplied via `-v` depending on `-cct`:
- aggregated or germline → dictionary.csv
- oncogenicity → somatic oncodict.csv
- clinical_impact → somatic clinicaldict.csv

## Example

Please, see the `examples_v2` directory and the `do.sh` scripts within.

## Output

The script returns 3 types of output:
- dot plots (`dot_plot` and additional numbered files if needed). Each dot is colored according to the type of effect as described in the legend. Plots are saved as `pdf` by default; `png` files are also produced when any of `-m`, `-r`, `-pltS`, or `-pltC` is used.
- a summary file `log.txt` containing useful information to summarize the effects of each module (GitBook compatible)
- a CSV file `mechanistic_indicators_out.csv` containing mutations with identified module effect and the consolidated `MAVISp Effects` column. This file honours all filters applied through `-m`, `-r`, `-pltS`, `-pltC`, `-vep`, and `-lgof`.

N.B. If `-m`/`-r`/`-pltS`/`-pltC` are used to delineate a specific subset, both the plots and `mechanistic_indicators_out.csv` will be restricted to that subset. The `-vep`/`-lgof` flags only affect `mechanistic_indicators_out.csv`; the plots always show the currently selected subset of mutations.