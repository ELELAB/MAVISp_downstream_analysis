module load python/3.10/modulefile
csv=$1
# change the csv path and the number of mutations you want to plot on the x-axis
# suggested not to use more than 50, otherwise you need to chage the figuresize;
# for more information about available options: python dot_plot.py -h

# Germline classification piecharts:
python piechart_v2.py -i $csv -d clinvar_interpretation_internal_dictionary.txt -o pie_germline.pdf -cct germline
python piechart_v2.py -i $csv -d clinvar_interpretation_internal_dictionary.txt -o pie_germline.png -cct germline

# Oncogenicity classification piecharts:
python piechart_v2.py -i $csv -d oncodict.csv -o pie_oncogenicity.pdf -cct oncogenicity
python piechart_v2.py -i $csv -d oncodict.csv -o pie_oncogenicity.png -cct oncogenicity

# Clinical impact classification piecharts:
python piechart_v2.py -i $csv -d clinicaldict.csv -o pie_clinical_impact.pdf -cct clinical_impact
python piechart_v2.py -i $csv -d clinicaldict.csv -o pie_clinical_impact.png -cct clinical_impact

