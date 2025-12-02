module load python
python ../stability_barplot.py -i input_files/ -o output --chunk_size 2 --grouped_columns "Stability (FoldX5, kcal/mol) [md_25]"  "Stability (RaSP, kcal/mol) [md_25]" "Stability (Rosetta Cartddg2020, kcal/mol) [md]" -m ../mutation_list.csv -he 21 -w 18 -xt 15 -xf 10 -yt 15  -p Pastel1
