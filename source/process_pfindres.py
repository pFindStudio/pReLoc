import pandas as pd
import numpy as np


def process_pFindspectra(pFindspectra_df, param):
    with open(param.AutoCustomInputPath, "w") as oc:
        oc.write("spec_name\tseq\tmod\tscore\n")
        for index, row in pFindspectra_df.iterrows():
            spec_name = row["File_Name"]
            pep = row["Sequence"].replace("I", "L")
            mod = row["Modification"]
            if isinstance(mod, float) and np.isnan(mod):
                mod = ""
            score = -np.log10(float(row["Final_Score"]))
            oc.write(spec_name + "\t" + pep + "\t" + mod + "\t" + str(score) + "\n")
    return


def write_findres(param):
    pfinddf = pd.read_csv(param.InputPath, sep="\t")
    pfinddf["Modification"] = pfinddf["Modification"].fillna("")
    pfinddf["pSite2_flag"] = False

    for index, row in pfinddf.iterrows():
        mod = row["Modification"]
        pep = row["Sequence"].replace("I", "L")
        for x in param.scoreModName:
            if x in mod:
                pfinddf.at[index, "pSite2_flag"] = True
    pfinddf.to_csv(param.InputPath, sep="\t", index=False)
    pfinddf = pfinddf[pfinddf["pSite2_flag"] == True]
    process_pFindspectra(pfinddf, param)
    return


def process_pfindres(param):
    if param.score_isSite:
        pfinddf = pd.read_csv(param.InputPath, sep="\t")
        process_pFindspectra(pfinddf, param)
    else:
        write_findres(param)
