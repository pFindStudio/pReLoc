import pandas as pd
import numpy as np


def read_summary(summarypath, topk):
    m_str_lst = []
    m_mass_lst = []
    with open(summarypath, "r") as rsp:
        slines = rsp.readlines()
        for line in slines[1:]:
            line_lst = line.split("\t")
            m_str = line_lst[1]
            m_mass = float(line_lst[2].split(" ")[0])
            m_str_lst.append(m_str)
            m_mass_lst.append(m_mass)
    if len(m_str_lst) <= topk:
        return m_str_lst, m_mass_lst
    else:
        return m_str_lst[:topk], m_mass_lst[:topk]


def write_findres(pfindspectrafp, summarypath, topk):
    m_str_lst, m_mass_lst = read_summary(summarypath, topk)
    pfinddf = pd.read_csv(pfindspectrafp, sep="\t")
    pfinddf["Modification"] = pfinddf["Modification"].fillna("")
    pfinddf["pSite2_flag"] = False

    for index, row in pfinddf.iterrows():
        mod = row["Modification"]
        pep = row["Sequence"].replace("I", "L")
        delta_mass = row[" Accurate modification mass "]
        if "PFIND_DELTA" in mod:
            for i in range(len(m_mass_lst)):
                tol_50 = m_mass_lst[i] * 50 * 1e-6
                if (
                    abs(delta_mass - m_mass_lst[i]) <= tol_50
                    and abs(delta_mass - m_mass_lst[i]) <= 0.05
                ):
                    modinfo_lst = mod.split(";")[:-1]
                    new_mod = ""
                    for modinfo in modinfo_lst:
                        if "PFIND_DELTA" not in modinfo:
                            new_mod += modinfo + ";"
                        else:
                            modindex = int(modinfo.split(",")[0]) - 1
                            if modindex != -1 and modindex != len(pep):
                                new_mod += (
                                    str(modindex + 1)
                                    + ","
                                    + m_str_lst[i]
                                    + "["
                                    + pep[modindex]
                                    + "];"
                                )
                            elif modindex == -1:
                                new_mod += (
                                    str(modindex + 1)
                                    + ","
                                    + m_str_lst[i]
                                    + "[ProteinN-term];"
                                )
                            else:
                                new_mod += (
                                    str(modindex + 1)
                                    + ","
                                    + m_str_lst[i]
                                    + "[ProteinC-term];"
                                )
                    pfinddf.at[index, "Modification"] = new_mod
                    pfinddf.at[index, "pSite2_flag"] = True
                    break

    # Keep the original search result reusable by both XGB and CL. Persisting
    # converted modifications here makes a second run discard every converted row.
    return pfinddf


def get_deltamass(param):
    pfind_df = write_findres(param.InputPath, param.InputSummary, param.pChemmodnum)
    # pfind_df = pfind_df[pfind_df["Modification"].str.contains("PFIND_DELTA", na=False)]
    pfind_df = pfind_df[pfind_df["pSite2_flag"] == True]

    # mgf_dict = read_mgf_file(res_set, wmgffp)
    with open(param.AutoCustomInputPath, "w") as oc:
        oc.write("spec_name\tseq\tmod\tscore\n")
        for index, row in pfind_df.iterrows():
            spec_name = row["File_Name"]
            pep = row["Sequence"].replace("I", "L")
            mod = row["Modification"]
            score = -np.log10(float(row["Final_Score"]))
            oc.write(spec_name + "\t" + pep + "\t" + mod + "\t" + str(score) + "\n")
    return
