module load python/3.10/modulefile
csv=$1
# change the csv path and the number of mutations you want to plot on the x-axis
# suggested not to use more than 50, otherwise you need to chage the figuresize;
# for more information about available options: python dot_plot.py -h
python dot_plot_v2.py -i $csv -x  50  -pltD -pltR -pltC benign likely_benign pathogenic likely_pathogenic -vep none -cct germline -v dictionary.csv  
python lolliplot.py -i mechanistic_indicators_out.csv -x 50 -s
python filter_pLDDT.py -i $csv -m mechanistic_indicators_out.csv -t 70
