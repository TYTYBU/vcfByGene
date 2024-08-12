import pandas as pd
import numpy as np
from pandas.api.types import CategoricalDtype
from random import shuffle
import argparse

def get_variants(gene, dir_path="selected_genes/parsed_selected_genes_filtered_variants_5nt_flank"):
    """retrieves annotated and parsed variant file"""

    path = f'/mnt/project/{dir_path}/{gene}_parsed.csv'

    return(pd.read_csv(path))

def get_carriers(gene, dir_path="selected_genes/selected_genes_var_patient_5nt_flank", ukb_af=False):
    """retrieve variant-patient mapping file"""

    path = f'/mnt/project/{dir_path}/{gene}.ssv'

    if ukb_af:
        carriers = pd.read_csv(path,
                        sep=" ",
                        names=["Chrom", "na", "Pos", "Ref", "Alt", "ukb_af", "Carriers"],
                        usecols=lambda x: x != 'na')
    else:
        carriers = pd.read_csv(path,
                        sep=" ",
                        names=["Chrom", "na", "Pos", "Ref", "Alt", "Carriers"],
                        usecols=lambda x: x != 'na')
    carriers.loc[:, "Carriers"] = carriers["Carriers"].str.strip("|")
    carriers.loc[:, "Carriers"]= carriers["Carriers"].str.split("|", expand = False)
    carriers = carriers.explode("Carriers").rename(columns={"Carriers": "Carrier"})
    carriers = carriers[carriers['Carrier'].apply(lambda x: not x.startswith('W'))] # remove withdrawn
    carriers["Chrom"] = carriers["Chrom"].str.extract('(\d+)').astype("int64")

    return(carriers)

def get_joined(gene, ukb_af, **paths):
    """retrieves both variant and mapping files and joins them together by variant"""

    variants = get_variants(gene, paths["var_path"]) if "var_path" in paths else get_variants(gene)
    carriers = get_carriers(gene, paths["car_path"], ukb_af) if "car_path" in paths else get_carriers(gene, ukb_af=ukb_af)
    ann_carriers = carriers.join(variants.set_index(["Chrom", "Pos", "Ref", "Alt"]), on=["Chrom", "Pos", "Ref", "Alt"])

    return(ann_carriers)


def get_select_carriers(ann_carriers):
    """
    accepts a joined carrier-variant df.
    filters by consequence and creates a string to represent amino acid substitutions when present
    """

    ann_carriers = ann_carriers.loc[(ann_carriers["consequence"].isin(["Missense", "Synonymous"])) |
                                   ((ann_carriers["consequence"] == "Deleterious") & (ann_carriers["LoF_confidence"] == "HC"))]
    
    ann_carriers._is_copy = False # gets rid of SettingWithCopyWarning
    ann_carriers.loc[ann_carriers.consequence == "Deleterious", "VEST4_score"] = 1 # set vest4 to max for del variants

    if ann_carriers.loc[:, "AA_pos"].dtypes != "object":
        ann_carriers.loc[:,"AA_pos"] = ann_carriers["AA_pos"].astype('Int64').astype('str')
    ann_carriers.loc[ann_carriers.AA_ref == "Stop", "AA_ref"] = "*"
    ann_carriers.loc[ann_carriers.AA_alt == "Stop", "AA_alt"] = "*"
    ann_carriers.loc[:, "AA_name"] = ann_carriers["AA_ref"] + ann_carriers["AA_pos"] + ann_carriers["AA_alt"]
    ann_carriers.loc[ann_carriers.consequence == "Synonymous", "AA_name"] = np.NaN
    ann_carriers.loc[(ann_carriers.consequence == "Deleterious") & (ann_carriers.coding_position.isna()), "AA_name"] = ann_carriers.loc[(ann_carriers.consequence == "Deleterious") & (ann_carriers.coding_position.isna()), "Name"]

    return(ann_carriers)

def prep_for_pivot(ann_carriers):
    """accepts a joined and filtered patient-variant df and fills in missing values to allow pandas pivot"""

    ann_carriers.loc[ann_carriers["CADD"].isnull(), "CADD"] = "dummy"
    ann_carriers.loc[ann_carriers["VEST4_score"].isnull(), "VEST4_score"] = "dummy"

    if "REVEL_score" in ann_carriers.columns:
        ann_carriers.loc[ann_carriers["REVEL_score"].isnull(), "REVEL_score"] = "dummy"

    if "ukb_af" in ann_carriers.columns:
        ann_carriers.loc[ann_carriers["ukb_af"].isnull(), "ukb_af"] = "dummy"

    ann_carriers.loc[ann_carriers["coding_position"].isnull(), "coding_position"] = "dummy"
    ann_carriers.loc[ann_carriers["AA_name"].isnull(), "AA_name"] = "dummy"
    ann_carriers.loc[:, "val"] = 1

    return(ann_carriers)

def get_most_severe(ann_carriers):
    """
    accepts a filtered patient-variant df.
    selects and retains only the most severe variant per patient. returns this df.
    """

    # order variant consequence by general severity (deleterious > missense > synonymous)
    consequence_cats = CategoricalDtype(categories=["Synonymous", "Missense", "Deleterious"], ordered=True)
    ann_carriers.loc[:, "consequence"] = ann_carriers["consequence"].astype(consequence_cats)

    # choose maximum AF for each variant, where possible
    ann_carriers.loc[:, "max_AF"] = ann_carriers.loc[:, ["allele_frequency", "gnomADe_MAX_AF", "gnomADg_MAX_AF"]].max(axis=1, skipna=True)

    # order deleterious variance confidence levels (high > low)
    del_confidence_cats = CategoricalDtype(categories=["LC", "HC"], ordered=True)
    ann_carriers.loc[:, "LoF_confidence"] = ann_carriers["LoF_confidence"].astype(del_confidence_cats)

    # create a random tie breaker index
    shuffled_index = list(range(0,len(ann_carriers)))
    shuffle(shuffled_index)
    ann_carriers.loc[:, "tie_breaker"] = shuffled_index

    # group by eid, sort by consequence > CADD > confidence (for del. vars) > AF > tiebreaker
    # select top variant for each group
    ann_carriers_grouped = ann_carriers.groupby("Carrier", group_keys=False).apply(pd.DataFrame.sort_values, ["consequence", "CADD", "LoF_confidence", "max_AF", "tie_breaker"], ascending=[False, False, False, True, True])
    ann_carriers_grouped = ann_carriers_grouped.drop_duplicates("Carrier", keep="first")

    # important!! return categorical variables to normal objects so they don't clog pivot memory
    ann_carriers_grouped.loc[:, "consequence"] = ann_carriers_grouped["consequence"].astype(object)
    ann_carriers_grouped.loc[:, "LoF_confidence"] = ann_carriers_grouped["LoF_confidence"].astype(object)

    return(ann_carriers_grouped)

def reshape_consequence(ann_carriers, rename_id=False, extra_cols=None):
    """
    accepts filtered and prepped patient-variant df
    returns a reshaped df that separates consequences into unique columns
    """
    
    if ("REVEL_score" in ann_carriers.columns) and ("ukb_af" in ann_carriers.columns):
        ann_carriers = ann_carriers.pivot_table(index=["Carrier", "Name", "coding_position", "ukb_af", "CADD", "REVEL_score", "VEST4_score", "AA_name"], values="val", columns="consequence", fill_value=0)
    elif "REVEL_score" in ann_carriers.columns:
        ann_carriers = ann_carriers.pivot_table(index=["Carrier", "Name", "coding_position", "CADD", "REVEL_score", "VEST4_score", "AA_name"], values="val", columns="consequence", fill_value=0)
    elif "ukb_af" in ann_carriers.columns:
        ann_carriers = ann_carriers.pivot_table(index=["Carrier", "Name", "coding_position", "ukb_af", "CADD", "VEST4_score", "AA_name"], values="val", columns="consequence", fill_value=0)
    else:
        ann_carriers = ann_carriers.pivot_table(index=["Carrier", "Name", "coding_position", "CADD", "VEST4_score", "AA_name"], values="val", columns="consequence", fill_value=0)
    
    ann_carriers = ann_carriers.reset_index().replace('dummy',np.nan)

    new_names={"Carrier": "eid", "Name": "variant_id"}

    # this variable is set to true whenever get_most_severe was run before this function
    if rename_id:
        new_names["Name"] = "most_severe_variant"

    selected_cols=["Carrier", "Name", "Synonymous", "Missense", "Deleterious", "coding_position", "AA_name"]
    if extra_cols:
            selected_cols += extra_cols
            
    ann_carriers = ann_carriers.reindex(columns=selected_cols)
    ann_carriers = ann_carriers.rename(columns=new_names)
    ann_carriers.columns.name = None

    return(ann_carriers)

def patient_var_mappings(gene, most_severe=False, extra_cols=None, **paths):
    """
    accepts gene name and some combination of var_path and car_path directories.
    returns a pandas df with one row per patient detailing their most severe variant
    """
    ukb_af = ("ukb_af" in extra_cols)
    annotated_carriers = get_joined(gene, ukb_af=ukb_af, **paths)
    annotated_carriers = get_select_carriers(annotated_carriers)
    annotated_carriers = prep_for_pivot(annotated_carriers)

    if most_severe:
        annotated_carriers = get_most_severe(annotated_carriers)

    return(reshape_consequence(annotated_carriers, rename_id=most_severe, extra_cols=extra_cols))

def parse_args():
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument("gene", help="symbol of gene to access", type=str)
    parser.add_argument("-v", "--var_path", help="directory path in which annotated variant file is stored", type=str)
    parser.add_argument("-c", "--car_path", help="directory path in which variant-patient mapping file is stored", type=str)
    parser.add_argument("-s", "--most_severe", help="whether to only include most severe variant per patient in output",
                        action="store_true")
    parser.add_argument('-e', '--extra_cols', help="comma delimited list of additional VEP columns to include in output", type=lambda s: [str(item) for item in s.split(',')])
    args = parser.parse_args()
    return args


def main():
    inputs = parse_args()
    file_name=f'{inputs.gene}.csv'
    patient_var_mappings(**vars(inputs)).to_csv(file_name, index=False)

if __name__ == "__main__":
    main()
