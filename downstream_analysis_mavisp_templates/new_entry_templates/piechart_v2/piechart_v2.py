import pandas as pd
import re
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import argcomplete
sns.set()

def get_clinvar_columns(df, clinvar_class_type: str):
    split = { "aggregated": ("ClinVar Interpretation", "ClinVar Review Status", "#ClinVar"),
              "germline": ("ClinVar Germline Interpretation", "ClinVar Germline Review Status", "#ClinVar"),
              "oncogenicity": ("ClinVar Oncogenicity Interpretation", "ClinVar Oncogenicity Review Status", "#Oncogenicity"),
              "clinical_impact": ("ClinVar Clinical Impact Interpretation", "ClinVar Clinical Impact Review Status", "#Clinical Impact")}
    pair = split[clinvar_class_type]
    if all(col in df.columns for col in pair[:2]):
        return pair
    if clinvar_class_type == "aggregated":
        raise ValueError( "ClinVar mode 'aggregated' requested (default), but aggregated columns are missing. "
                          "This looks like a 'new style' entry csv. "
                          "Use -cct germline|oncogenicity|clinical_impact.")
    raise ValueError(f"ClinVar mode '{clinvar_class_type}' requested, but the corresponding columns are missing.")

def piechart(df,out,clinvar_cols):

    # Select the "ClinVar Interpretation" column
    clinvar_interpretation = df[clinvar_cols[0]]
    # Count the occurrences of each value
    value_counts = clinvar_interpretation.value_counts()

    # Create a pie chart
    plt.figure(figsize=(12, 8))
    #colors = ['lightcoral', 'lightblue', 'lightgreen']  # Define your desired colors
    plt.pie(value_counts, startangle=140,autopct='')

    # Create a legend with percentages
    legend_labels = [f'{label} ({value_counts[label] / sum(value_counts) * 100:.1f}%)' for label in value_counts.index]
    plt.legend(legend_labels, title=clinvar_cols[0], loc="best",fontsize=18)
    plt.axis('equal')
    plt.tight_layout()
    # Show the pie chart
    plt.savefig(out,dpi=300)

if __name__ == "__main__":

    # Add arguments required to the script to a argparse.ArgumentParser instance.
    description = "Plot of the pie chart on the ClinVar intepretation " \
                  "column on mavisp aggregate csv. " 
    parser = argparse.ArgumentParser(description = description)

    i_helpstr = "Input: MAVISp aggregated csv."
    parser.add_argument("-i", "--input",
                        action = "store",
                        type = str,
                        help = i_helpstr,
                        required = True)

    clinvarclasstype_helpstr = f"ClinVar classification type to use for plotting/coloring. " \
                               f"Choices:aggregated, germline, oncogenicity, clinical_impact."
    parser.add_argument("-cct", "--clinvar_class_type",
                        default = "aggregated",
                        choices = ["aggregated", "germline", "oncogenicity", "clinical_impact"],
                        help = clinvarclasstype_helpstr)

    d_helpstr = "Input: MAVISp dictionary file corresponding to selected clinvar classification type"
    parser.add_argument("-d", "--dictionary",
                        action = "store",
                        type = str,
                        help = d_helpstr,
                        required = True)    

    o_helpstr = "Input: MAVISp aggregated csv."
    parser.add_argument("-o", "--out_name",
                        action = "store",
                        type = str,
                        help = i_helpstr,
                        required = True) 
    
    argcomplete.autocomplete(parser)
    args = parser.parse_args()   
    
    # Read csv
    df = pd.read_csv(args.input)
    try:
        clinvar_cols = get_clinvar_columns(df, args.clinvar_class_type)
    except ValueError as e:
        print(f"ERROR: {e}")
        exit(1)
    df[clinvar_cols[0]] = df[clinvar_cols[0]].astype("string").replace(r",\s+", ",", regex=True)
   
    # Check if correct dict was provided:
    first_cell = pd.read_csv(args.dictionary, sep="\t", header=None).iat[0, 0]
    if str(first_cell).strip() != clinvar_cols[2]:
        print(f"ERROR: Wrong ClinVar dictionary for --clinvar_class_type {args.clinvar_class_type}")
        exit(1)

    # Read  dictionary file and convert it to dictionary
    clinvar_dict = pd.read_csv(args.dictionary,sep="\t")
    clinvar_dict = clinvar_dict.set_index(clinvar_cols[2]).T.to_dict('records')
    clinvar_dict[0] = {k: (v.strip() if isinstance(v, str) else v) for k, v in clinvar_dict[0].items()}

    # Replace with Mavisp internal dictionary
    df[clinvar_cols[0]] = df[clinvar_cols[0]].replace(clinvar_dict[0])
    
    # Plot
    piechart(df,args.out_name,clinvar_cols)
