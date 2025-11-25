# MAVISp dot plot
*Last updated*: 25/11/25

## Description

`dot_plot.py` ingests a MAVISp aggregated CSV and, optionally, a ClinVar dictionary.  
It produces:
- publication-ready dot plots that summarise structural and functional module effects per mutation;
- a `log.txt` which contains a summary of mutational effects
- `mechanistic_indicators_out.csv`, which lists mutations with at least one MAVISp module effect together with the aggregated `MAVISp Effects` column. Optionally, this table is filtered by keeping only mutations that are damaging or pathogenic according to VEPs. This can be turned on with option `-vep` (arguments AlphaMissense, REVEL, GEMME, EVE or none) and can be further restricted to mutations that are LoF or GoF for DeMaSk via `-lgof`.

## Requirements

- Python >= 3.7
- pandas
- numpy
- matplotlib
- argcomplete (optional, for shell auto-completion)

## Inputs 

- MAVISp CSV (the full aggregated table; required)
- `dictionary.csv` (ClinVar annotation dictionary; only required when `-pltC/--plot_Clinvar` or `-colC/--color_Clinvar` is used)

## Usage

```
python dot_plot.py [-h] -i INPUT [-v CLINVAR_DICTIONARY] [-o OUTPUT]
                   [-m MUTATIONS [MUTATIONS ...]] [-r RESIDUES [RESIDUES ...]]
                   [-R REVEL_THRESHOLD] [-D DEMASK_THRESHOLD]
                   [-G GEMME_THRESHOLD] [-x X_LIM] [-f FIGSIZE FIGSIZE]
                   [-pltR] [-pltD]
                   [-pltC {all,uncertain,benign,likely_benign,pathogenic,likely_pathogenic,conflicting}
                          [{all,uncertain,benign,likely_benign,pathogenic,likely_pathogenic,conflicting} ...]]
                   [-colC]
                   [-pltS {saturation,cosmic,cbioportal} [{saturation,cosmic,cbioportal} ...]]
                   [-vep {none,alphamissense,revel,gemme,eve}]
                   [-lgof]
```

- `-i/--input`: MAVISp CSV to process (required).
- `-v/--clinvar-dictionary`: path to `dictionary.csv`. Only needed when plotting or colouring ClinVar categories (default: `dictionary.csv`).
- `-o/--output`: base name for dot-plot files (default: `dot_plot`).
- `-m/--mutations`: comma-separated mutations to display (e.g. `A4G,F55K`). Mutually exclusive with `-r`.
- `-r/--residues`: comma-separated residue positions to display (e.g. `4,55`). Mutually exclusive with `-m`.
- `-R/--revel_threshold`: REVEL pathogenic threshold (default: `0.5`).
- `-D/--demask_threshold`: DeMaSk delta-fitness threshold for LoF/GoF calls (default: `0.25`).
- `-G/--gemme_threshold`: GEMME threshold (default: `-3.0`).
- `-x/--x_lim`: number of mutations per panel before splitting across multiple figures (default: `50`).
- `-f/--figsize`: figure width and height (default: `14 5`). The default works well for ~50 mutations and 7–8 labels.
- `-pltR/--plot_Revel`: add REVEL classifications to the dot plot.
- `-pltD/--plot_Demask`: add DeMaSk predicted consequence (LoF/GoF) for mutations meeting the `-D` threshold.
- `-pltC/--plot_Clinvar`: filter to specific ClinVar categories (e.g. `pathogenic uncertain`). Requires `dictionary.csv`.
- `-colC/--color_Clinvar`: colour the x-axis labels according to ClinVar categories. Requires `dictionary.csv`.
- `-pltS/--plot_Source`: filter mutations by source (`saturation`, `cosmic`, `cbioportal`). Multiple sources can be provided; filters are additive with `-pltC`.
- `-vep/--vep-filter`: restrict `mechanistic_indicators_out.csv` to mutations predicted as pathogenic by the selected VEP. Choices are `alphamissense`, `revel`, `gemme`, `eve`, or `none` (default). Supplying `-vep` without an argument defaults to `alphamissense`.
- `-lgof/--vep-filter-lgof`: when set, only keep entries classified as DeMaSk LoF or GoF in `mechanistic_indicators_out.csv`. By default, this filtering is not performed.

## Example

See the `example` directory and the accompanying `do.sh` script for a minimal end-to-end run.

## Output

Running the script produces:
- `dot_plot.pdf` (and additional numbered PDFs if more mutations exceed `-x`). PNGs are also written when any of `-m`, `-r`, `-pltS`, or `-pltC` is used.
- `log.txt`, summarising how many variants satisfy each classifier (REVEL, GEMME, DeMaSk, EVE, AlphaMissense) and providing module-level counts.
- `mechanistic_indicators_out.csv`, containing the filtered subset of mutations with at least one module effect and the consolidated `MAVISp Effects` column. This file honours all filters applied through `-m`, `-r`, `-pltS`, `-pltC`, `-vep`, and `-lgof`.

Notes:
- `-pltS` and `-pltC` filters are additive; the resulting dataset must satisfy both when both options are provided.
- Any subset defined via `-m`, `-r`, `-pltS`, or `-pltC` applies to both the plots and the mechanistic output table.
- The `-vep`/`-lgof` flags only affect `mechanistic_indicators_out.csv`; the plots always show the currently selected subset of mutations.

## Testing

- Install the runtime requirements plus `pytest` (e.g. `pip install -r requirements.txt pytest`).
- From the `dot_plot/` directory run `pytest` to execute both unit tests (helper functions) and integration tests that drive the CLI end-to-end.
- The fixtures in `tests/data/simple_mode` and `tests/data/ensemble_mode` are trimmed copies of the example datasets together with golden `mechanistic_indicators_out.csv` files; update those goldens if behaviour intentionally changes.
- To run only the fast tests: `pytest tests/unit`. To re-run CLI tests: `pytest -m integration`.
