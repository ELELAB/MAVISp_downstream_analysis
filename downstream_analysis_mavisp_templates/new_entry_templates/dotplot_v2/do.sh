module load python/3.10/modulefile
csv=$1
# change the csv path and the number of mutations you want to plot on the x-axis
# suggested not to use more than 50, otherwise you need to chage the figuresize;
# for more information about available options: python dot_plot.py -h

# Germline classification dotplot
python dot_plot_v2.py -i $csv -x  50 -pltD -pltR -colC -vep -lgof -cct germline -v dictionary.csv -o dot_plot_germline
mv log.txt log_germline.txt

python lolliplot.py -i mechanistic_indicators_out.csv -x 50 -s
python filter_pLDDT.py -i $csv -m mechanistic_indicators_out.csv -t 70

# Oncogenicity dotplot
python dot_plot_v2.py -i $csv -x  50 -pltD -pltR -colC -vep -lgof -cct oncogenicity -v oncodict.csv -o dot_plot_oncogenicity
mv log.txt log_oncogenicity.txt

# Clinical impact dotplot
python dot_plot_v2.py -i $csv -x  50 -pltD -pltR -colC -vep -lgof -cct clinical_impact -v clinicaldict.csv -o dot_plot_clinical_impact
mv log.txt log_clinical_impact.txt
