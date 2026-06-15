module load python/3.10/modulefile
csv=$1
# change the csv path and the number of mutations you want to plot on the x-axis
# suggested not to use more than 50, otherwise you need to chage the figuresize;
# for more information about available options: python dot_plot.py -h

# Germline classification results:
python dot_plot_v2.py -i $csv -x  50  -pltD -pltR -pltC benign likely_benign pathogenic likely_pathogenic -vep none -cct germline -v dictionary.csv -o dot_plot_germline
python lolliplot.py -i mechanistic_indicators_out_germline.csv -x 50 -s -o lolliplot_germline
python filter_pLDDT.py -i $csv -m mechanistic_indicators_out_germline.csv -t 70 -o mech_plddt_germline

# Oncogenicity classification results:
python dot_plot_v2.py -i $csv -x  50  -pltD -pltR -pltC benign likely_benign oncogenic likely_oncogenic -vep none -cct oncogenicity -v oncodict.csv -o dot_plot_oncogenicity 
python lolliplot.py -i mechanistic_indicators_out_oncogenicity.csv -x 50 -s -o lolliplot_oncogenicity
python filter_pLDDT.py -i $csv -m mechanistic_indicators_out_oncogenicity.csv -t 70 -o mech_plddt_oncogenicity

# Clinical impact classification results:
python dot_plot_v2.py -i $csv -x  50  -pltD -pltR -pltC benign potential strong -vep none -cct clinical_impact -v clinicaldict.csv -o dot_plot_clinical_impact
python lolliplot.py -i mechanistic_indicators_out_clinical_impact.csv -x 50 -s -o lolliplot_clinical_impact
python filter_pLDDT.py -i $csv -m mechanistic_indicators_out_clinical_impact.csv -t 70 -o mech_plddt_clinical_impact
