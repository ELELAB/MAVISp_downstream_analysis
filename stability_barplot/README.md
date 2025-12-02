# stability_barplot.py

## Affiliation

Cancer Structural Biology, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark

Cancer Systems Biology, Health and Technology Department, Section for Bioinformatics, 2800, Lyngby, Denmark

## Introduction

The stability_barplot.py script is part of the downstream analysis pipeline within the MAVISp framework. It visualizes ΔΔG values predicted by Rosetta, RaSP, and FoldX, as provided in the MAVISp output CSV files in both simple and ensemble modes.

The script automatically parses the input folder, distinguishing between simple and ensemble files, and generates two separate outputs accordingly. For each mode, it produces bar plots displaying the predicted ΔΔG values along with their associated standard deviations. Dedicated bar plots are generated for each method and for each ensemble or structural source (e.g., AlphaFold models, PDBs, molecular dynamics clusters, etc.).

Additionally, for every ensemble method or structural source, the script creates a combined bar plot that aggregates the ΔΔG predictions from all available methods (RaSP, Rosetta, FoldX) for each mutation. This combined plot can be customized by providing, via the appropriate flag, a specific set of columns to be aggregated together.

## Requirements
The pipeline uses the following Python packages:

- pandas
- seaborn
- re
- statistics
- argparse
- numpy
- os

## Script explanation

The stability_barplot.py script takes as input a folder containing the MAVISp CSV output files. The files can be either in simple or ensemble mode. If both types are present, the script automatically distinguishes between them and, for each mode, concatenates the corresponding CSVs after adding the name of the associated protein.

For each method (RaSP, FoldX, Rosetta) and for each structural source in simple mode (e.g., AlphaFold, PDB) or ensemble method (e.g., MD, MD_25, clusters), the script extracts the stability columns to be plotted along with the associated standard deviation columns. It generates bar plots for each method and structural source or ensemble method, showing the ΔΔG values with the computed standard deviations.

Additionally, for each ensemble method or structural source, the script produces a summary plot aggregating the predictions from all available methods (RaSP, FoldX, Rosetta) for every mutation for that specific structural source or ensemble method.

The data frame can be filtered using a custom file to select mutations of interest, or with flags to include only mutations predicted as destabilizing according to the FoldX/RaSP consensus or FoldX/Rosetta consensus (see the MAVISp paper for details). Optional flags can also be specified to customize bar colors, layout, and label sizes.


## Input structure

The stability_barplot.py necessitates the folder contianing the mavisp csv  files and must contai for every stability column the associated stdv columns (-i flag). The following flags can be specified in order to customize the barplot:

- **-i, --input_file** The input file containing the Rosetta, MutateX, RaSP column along with the corresponding standard deviations (output file from std.py).

- **-o, --output_folder** Path to the folder where output figures and CSVs will be saved.

- **-v, --vertical_plot** To create a vertical plot instead of the default horizontal one 

- **-m, --mutation_list** Optional CSV file with columns 'protein,Mutation' used to filter the dataset before plotting.

- **-g", "--grouped_columns** Specific columns to include in grouped methods plots. Default: None (all methods).

- **-dp, --only_destabilizing** Plot only the destabilizing mutations according to Rosetta (This will not influence the plotting of RaSP ΔΔG). In this case, the "Stability classification...(Rosetta/FoldX5)" column is required in the input file

- **-dr, --only_destabilizing** Plot only the destabilizing mutations according to RaSP (This will not influence the plotting of Rosetta ΔΔG). In this case, the "Stability classification...(RaSP/FoldX5)" column is required in the input file

- **-p, --palette** Color palette to automatically assign colors to methods (FoldX, Rosetta, RaSP).
                    The selected palette will generate N distinct colors based on the number of methods.
                    Available palettes (matplotlib):
                    - tab10
                    - tab20
                    - Set1, Set2, Set3
                    - Pastel1, Pastel2
                    - Dark2
                    - Accent
                    - Paired
                    - viridis, plasma, inferno, magma, cividis
                    - coolwarm, bwr, seismic
                    - turbo
                    Default: tab10

- **-c, --color_accordingly_AF** Color the mutation on the x-axis according to the plDDT score of the AF model used for the Rosetta and MuateX calculation. In this case  also the "AlphaFold2 model pLDDT score" or a column with the plDDt score is needed in the input file.

- **-ch, ----chunk_size** Number of mutations to show in each plot. This flag is recommended in case the number of mutations to plot is higher than 40. The script will create as many plots as the chunks of mutations formed based on the provided number. It's recommended not to exceed 40 mutations per plot and to provide a number that is divisible by the number of mutations to plot.

- **-he, --figure height** Heigh of the final figure. Default 215 mm. 

- **-w, --figure width** Width of the final figure. Default 180 mm.

- **-bw, --bars width** Bars width. Default. 0.45.

- **-xf, --x_font_labels** Fontsize of the x labels (Mutations on the x-axis): Default 7.

- **-xw, --weight_of_xlabels** Weight of the x-labels: Default regular.

- **-xt, --x_title_font** Fontsize of the x-title ("Mutation"). Default: 15.

- **-yt, --t_title_font** Fontsize of the y-title (ΔΔG [kcal/mol]). Default: 15

The mutation list must be containing the following column:

**protein** Protein name 
**mutation** mutation to select 

```
protein,Mutation
POLD1,L474P
POLE,S297F
POLE,L424V
POLE,P436R
POLE,P436S
POLE,E830G
```

## Output
The plots are organized in subfolders that indicate the MAVISp mode of the dataframe, the structural source (for simple mode) or the ensemble type (for ensemble mode), the stability method, and the grouped methods showing the ΔΔG predictions from every structural source or ensemble for all stability methods.

```
output/
├── ensemble
│   ├── md
│   │   ├── FoldX5
│   │   ├── grouped_methods
│   │   ├── RaSP
│   │   └── Rosetta Cartddg2020
│   └── md_25
│       ├── FoldX5
│       ├── grouped_methods
│       └── RaSP
└── simple
    └── alphafold
        ├── FoldX5
        ├── grouped_methods
        └── RaSP
```


## Usage

```
module load python/3.10/modulefile
python new_script.py -i input_files/ -o output --chunk_size 6 --grouped_columns "Stability (FoldX5, kcal/mol) [md_25]"  "Stability (RaSP, kcal/mol) [md_25]" "Stability (Rosetta Cartddg2020, kcal/mol) [md]" -m mutation_list.csv -he 21 -w 18 -xt 15 -xf 10 -yt 15 -dp -dr -v
```







