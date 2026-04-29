module load python/3.10/modulefile
csv=$1
# change the csv path and the number of mutations you want to plot on the x-axis
# suggested not to use more than 50, otherwise you need to chage the figuresize;
# for more information about available options: python dot_plot.py -h
# Germline VUS classification results:
python dot_plot_v2.py -i $csv -x  50  -pltD -pltR -pltC conflicting uncertain -vep -lgof -cct germline -v dictionary.csv -o dot_plot_vus_germline
mv mechanistic_indicators_out.csv mechanistic_indicators_out_vus_germline.csv
mv log.txt log_vus_germline.txt
python lolliplot.py -i mechanistic_indicators_out_vus_germline.csv -x 50 -s
mv lolliplot.pdf lolliplot_vus_germline.pdf
python filter_pLDDT.py -i $csv -m mechanistic_indicators_out_vus_germline.csv -t 70 -o mech_plddt_vus_germline

# Oncogenicity VUS classification results:
python dot_plot_v2.py -i $csv -x  50  -pltD -pltR -pltC conflicting uncertain -vep -lgof -cct oncogenicity -v oncodict.csv -o dot_plot_vus_oncogenicity
mv mechanistic_indicators_out.csv mechanistic_indicators_out_vus_oncogenicity.csv
mv log.txt log_vus_oncogenicity.txt
python lolliplot.py -i mechanistic_indicators_out_vus_oncogenicity.csv -x 50
mv lolliplot.pdf lolliplot_vus_oncogenicity.pdf
python filter_pLDDT.py -i $csv -m mechanistic_indicators_out_vus_oncogenicity.csv -t 70 -o mech_plddt_vus_oncogenicity

# Clinical impact VUS classification results:
python dot_plot_v2.py -i $csv -x  50  -pltD -pltR -pltC conflicting unknown -vep -lgof -cct clinical_impact -v clinicaldict.csv -o dot_plot_vus_clinical_impact
mv mechanistic_indicators_out.csv mechanistic_indicators_out_vus_clinical_impact.csv
mv log.txt log_vus_clinical_impact.txt
python lolliplot.py -i mechanistic_indicators_out_vus_clinical_impact.csv -x 50
mv lolliplot.pdf lolliplot_vus_clinical_impact.pdf
python filter_pLDDT.py -i $csv -m mechanistic_indicators_out_vus_clinical_impact.csv -t 70 -o mech_plddt_vus_clinical_impact