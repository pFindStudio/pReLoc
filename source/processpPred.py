from resources import resource_path
import pandas as pd
import numpy as np
import subprocess
import os
import copy
from scipy.stats import spearmanr, pearsonr
from sklearn.metrics.pairwise import cosine_similarity
from itertools import islice


def keep_first_77_lines(input_file_path, output_file_path):
    with open(input_file_path, "r") as file:
        lines = list(islice(file, 77))

    with open(output_file_path, "w") as file:
        file.writelines(lines)


def getallfile_ft(true_spectra, filefolder, ms2_allfp, rt_allfp):
    df = pd.read_csv(true_spectra, sep="\t")
    df = df[["File_Name", "Sequence", "Modification"]].rename(
        columns={"File_Name": "spec_name", "Sequence": "pep", "Modification": "mod"}
    )
    filename_set = set(df["spec_name"].str.split(".").str[0])

    with open(ms2_allfp, "a") as maf:
        maf.truncate(0)
    with open(rt_allfp, "a") as raf:
        raf.truncate(0)

    i = 0
    for filename in filename_set:
        print(f"Processing file: {filename}")
        # key_fp = filefolder + filename + ".key"
        ms2_fp = filefolder + filename + ".ms2data"
        rt_fp = filefolder + filename + ".rtdata"

        # Process the .ms2data file.
        with open(ms2_fp, "r") as mf:
            with open(ms2_allfp, "a") as maf:
                while True:
                    ms2line = mf.readline()
                    if ms2line.startswith("1"):
                        if i == 0:
                            maf.write(ms2line)
                        continue
                    if not ms2line:
                        break
                    maf.write(ms2line)

        # Reopen the key file to process the .rtdata file.
        with open(rt_fp, "r") as rf:
            with open(rt_allfp, "a") as raf:
                while True:
                    rtline = rf.readline()
                    if rtline.startswith("seq"):
                        if i == 0:
                            raf.write(rtline)
                        continue
                    if not rtline:
                        break
                    raf.write(rtline)
        i += 1
    return


def finetune_pPred(param):
    # finetune MS2 RT
    ms2_model_fp = param.ftmodelfolder + "\\" + param.ftmodelname + "_ms2.pt"
    rt_model_fp = param.ftmodelfolder + "\\" + param.ftmodelname + "_rt.pt"
    keep_first_77_lines(param.pFindResfp, "pPred_data\\pFind77.spectra")

    command_GenPredData_finetune = [
        "GenPredData.exe",
        "--pfind_res",
        "pPred_data\\pFind77.spectra",
        "--qvalue_threshold",
        "0.001",
        "--precursor_tol",
        "10.0",
        "--output_path",
        "pPred_data\\finetune",
        "--ms_path",
        param.rawfolder,
        "--instrument",
        param.instrument,
        "--nce",
        str(param.nce) + "," + str(param.nce) + "," + str(param.nce),
        "--output_all",
    ]
    result = subprocess.run(
        command_GenPredData_finetune, capture_output=True, text=True
    )
    print("GenPredData.exe get finetune data:", result.stdout)
    print("GenPredData.exe get finetune data:", result.stderr)

    getallfile_ft(
        "pPred_data\\pFind77.spectra",
        "pPred_data\\finetune\\",
        "pPred_data\\pFind77.ms2data",
        "pPred_data\\pFind77.rtdata",
    )

    if not os.path.exists(rt_model_fp):
        command_pPred_finetuneRT = [
            "pPred.exe",
            "-rt",
            "--tune",
            "-m",
            "pPred_data\\mainmodel\\pPredRT_main.pt",
            "-n",
            param.ftmodelname,
            "-om",
            param.ftmodelfolder,
            "-d",
            "pPred_data\\pFind77.rtdata",
            "-i",
            param.rtparam,
            "-l",
            "pPred_data\\log.log",
        ]
        resultRT = subprocess.run(
            command_pPred_finetuneRT, capture_output=True, text=True
        )
        print("pPred.exe finetune RT model:", resultRT.stdout)
        print("pPred.exe finetune RT model:", resultRT.stderr)

    if not os.path.exists(ms2_model_fp):
        pho_flag = False
        for x in param.scoreModName:
            if "Phospho" in x:
                pho_flag = True
                break
        if not pho_flag:
            command_pPred_finetuneMS2 = [
                "pPred.exe",
                "-ms2",
                "--tune",
                "-m",
                "pPred_data\\mainmodel\\pPredMS2_main.pt",
                "-n",
                param.ftmodelname,
                "-om",
                param.ftmodelfolder,
                "-d",
                "pPred_data\\pFind77.ms2data",
                "-i",
                param.ms2param,
                "-l",
                "pPred_data\\log.log",
            ]
        else:
            command_pPred_finetuneMS2 = [
                "pPred.exe",
                "-ms2",
                "--tune",
                "-m",
                "pPred_data\\mainmodel\\pPredMS2_phos_0908.pt",
                "-n",
                param.ftmodelname,
                "-om",
                param.ftmodelfolder,
                "-d",
                "pPred_data\\pFind77.ms2data",
                "-i",
                param.ms2param,
                "-l",
                "pPred_data\\log.log",
            ]
        resultMS2 = subprocess.run(
            command_pPred_finetuneMS2, capture_output=True, text=True
        )
        print("pPred.exe finetune MS2 model:", resultMS2.stdout)
        print("pPred.exe finetune MS2 model:", resultMS2.stderr)


def get_pFindRes(param):
    df = pd.read_csv(param.pSiteResPath, sep="\t")
    df_extracted = df[["pep", "mod"]].rename(columns={"pep": "seq"})
    df_extracted.to_csv("pPred_data\\ori_pred.pepdata", sep="\t", index=False)

    df_extracted_can = df[["bm25_pep", "bm25_mod"]].rename(
        columns={"bm25_pep": "seq", "bm25_mod": "mod"}
    )
    df_extracted_can.to_csv("pPred_data\\can_pred.pepdata", sep="\t", index=False)

    title = "File_Name\tScan_No\tExp.MH+\tCharge\tQ-value\tSequence\tCalc.MH+\tMass_Shift(Exp.-Calc.)\tRaw_Score\tFinal_Score\tModification\tSpecificity\tProteins\tPositions\tLabel\tTarget/Decoy\tMiss.Clv.Sites\tAvg.Frag.Mass.Shift\tOthers\n"
    df["mod"].fillna(" ", inplace=True)
    with open("pPred_data\\true.spectra", "w") as ts:
        ts.write(title)
        for index, row in df.iterrows():
            spec_name = row["spec_name"]
            pep = row["pep"]
            mod = row["mod"]
            charge = spec_name.split(".")[3]
            scan_num = spec_name.split(".")[1]
            row_str = (
                spec_name
                + "\t"
                + scan_num
                + "\t\t"
                + charge
                + "\t0\t"
                + pep
                + "\t \t0\t \t \t"
                + mod
                + "\t \t \t \t \t \t \t \t \n"
            )
            ts.write(row_str)
    return


def predict_pPred(param):
    command_pPred_predictoriRT = [
        "pPred.exe",
        "-rt",
        "--pred",
        "-m",
        param.ftmodelfolder + "\\" + param.ftmodelname + "_rt.pt",
        "-d",
        "pPred_data\\ori_pred.pepdata",
        "-op",
        "pPred_data\\ori_pred.pepdata.rt.pred",
        "-i",
        param.rtparam,
        "-l",
        "pPred_data\\log.log",
    ]
    result_predictoriRT = subprocess.run(
        command_pPred_predictoriRT, capture_output=True, text=True
    )
    print("pPred.exe predict original pep RT:", result_predictoriRT.stdout)
    print("pPred.exe predict original pep RT:", result_predictoriRT.stderr)

    command_pPred_predictoriMS2 = [
        "pPred.exe",
        "-ms2",
        "--pred",
        "-m",
        param.ftmodelfolder + "\\" + param.ftmodelname + "_ms2.pt",
        "-d",
        "pPred_data\\ori_pred.pepdata",
        "-op",
        "pPred_data\\ori_pred.pepdata.ms2.pred",
        "-i",
        param.ms2param,
        "-l",
        "pPred_data\\log.log",
    ]
    result_predictoriMS2 = subprocess.run(
        command_pPred_predictoriMS2, capture_output=True, text=True
    )
    print("pPred.exe predict original pep MS2:", result_predictoriMS2.stdout)
    print("pPred.exe predict original pep MS2:", result_predictoriMS2.stderr)

    command_pPred_predictcanRT = [
        "pPred.exe",
        "-rt",
        "--pred",
        "-m",
        param.ftmodelfolder + "\\" + param.ftmodelname + "_rt.pt",
        "-d",
        "pPred_data\\can_pred.pepdata",
        "-op",
        "pPred_data\\can_pred.pepdata.rt.pred",
        "-i",
        param.rtparam,
        "-l",
        "pPred_data\\log.log",
    ]
    result_predictcanRT = subprocess.run(
        command_pPred_predictcanRT, capture_output=True, text=True
    )
    print("pPred.exe predict candidate pep RT:", result_predictcanRT.stdout)
    print("pPred.exe predict candidate pep RT:", result_predictcanRT.stderr)

    command_pPred_predictcanMS2 = [
        "pPred.exe",
        "-ms2",
        "--pred",
        "-m",
        param.ftmodelfolder + "\\" + param.ftmodelname + "_ms2.pt",
        "-d",
        "pPred_data\\can_pred.pepdata",
        "-op",
        "pPred_data\\can_pred.pepdata.ms2.pred",
        "-i",
        param.ms2param,
        "-l",
        "pPred_data\\log.log",
    ]
    result_predictcanMS2 = subprocess.run(
        command_pPred_predictcanMS2, capture_output=True, text=True
    )
    print("pPred.exe predict candidate pep MS2:", result_predictcanMS2.stdout)
    print("pPred.exe predict candidate pep MS2:", result_predictcanMS2.stderr)


def getallfile(true_spectra, filefolder, ms2_allfp, rt_allfp):
    ms2_allfp = ms2_allfp[:-8] + "_.ms2data"
    rt_allfp = rt_allfp[:-7] + "_.rtdata"
    df = pd.read_csv(true_spectra, sep="\t")
    df = df[["File_Name", "Sequence", "Modification"]].rename(
        columns={"File_Name": "spec_name", "Sequence": "pep", "Modification": "mod"}
    )
    filename_set = set(df["spec_name"].str.split(".").str[0])
    count = 0

    with open(ms2_allfp, "a") as maf:
        maf.truncate(0)
    with open(rt_allfp, "a") as raf:
        raf.truncate(0)

    for filename in filename_set:
        print(f"Processing file: {filename}")
        key_fp = filefolder + filename + ".key"
        ms2_fp = filefolder + filename + ".ms2data"
        rt_fp = filefolder + filename + ".rtdata"

        # Process the .ms2data file.
        with open(key_fp, "r") as kf:
            with open(ms2_fp, "r") as mf:
                with open(ms2_allfp, "a") as maf:
                    while True:
                        keyline = kf.readline()
                        ms2line = mf.readline()
                        if keyline.startswith("spec"):
                            if count == 0:
                                maf.write(keyline[:-1] + "\t" + ms2line)
                            continue
                        if not keyline or not ms2line:
                            break
                        maf.write(keyline[:-1] + "\t" + ms2line)

        # Reopen the key file to process the .rtdata file.
        with open(key_fp, "r") as kf:
            with open(rt_fp, "r") as rf:
                with open(rt_allfp, "a") as raf:
                    while True:
                        keyline = kf.readline()
                        rtline = rf.readline()
                        if keyline.startswith("spec"):
                            if count == 0:
                                raf.write(keyline[:-1] + "\t" + rtline)
                                count += 1
                            continue
                        if not keyline or not rtline:
                            break
                        raf.write(keyline[:-1] + "\t" + rtline)
                        count += 1
    return count


def write_empty_lines(filepath, num_lines):
    with open(filepath, "w") as file:
        for _ in range(num_lines):
            file.write("\n")


def read_spectra(spectrafp):
    df = pd.read_csv(spectrafp, sep="\t")
    res_dict = {}
    for index, row in df.iterrows():
        key = row["File_Name"] + row["Sequence"] + row["Modification"]
        if key in res_dict.keys():
            res_dict[key].append(index)
        else:
            res_dict[key] = [index]
    return res_dict


def writefile(ms2_allfp, rt_allfp, num_lines, res_dict):
    ms2_tmpallfp = ms2_allfp[:-8] + "_.ms2data"
    rt_tmpallfp = rt_allfp[:-7] + "_.rtdata"
    ms2_res_dict = copy.deepcopy(res_dict)
    rt_res_dict = copy.deepcopy(res_dict)

    write_empty_lines(ms2_allfp, num_lines)
    write_empty_lines(rt_allfp, num_lines)

    with open(ms2_tmpallfp, "r") as ms2_tmp:
        ms2_lines = ms2_tmp.readlines()
        ms2_lines_tmp = ms2_lines.copy()
    with open(rt_tmpallfp, "r") as rt_tmp:
        rt_lines = rt_tmp.readlines()
        rt_lines_tmp = rt_lines.copy()

    for line in ms2_lines_tmp:
        lst = line.split("\t")
        line = "\t".join(lst[1:])
        if lst[0] == "spec":
            ms2_lines[0] = line
            continue
        key = lst[0] + lst[1] + lst[2]
        if lst[2] == "":
            key = key + " "
        lines_num = ms2_res_dict[key].pop(0)
        ms2_lines[lines_num + 1] = line
    for line in rt_lines_tmp:
        lst = line.split("\t")
        line = "\t".join(lst[1:])
        if lst[0] == "spec":
            rt_lines[0] = line
            continue
        key = lst[0] + lst[1] + lst[2]
        if lst[2] == "":
            key = key + " "
        lines_num = rt_res_dict[key].pop(0)
        rt_lines[lines_num + 1] = line

    with open(ms2_allfp, "w") as ms2_all:
        ms2_all.writelines(ms2_lines)
    with open(rt_allfp, "w") as rt_all:
        rt_all.writelines(rt_lines)

    os.remove(ms2_tmpallfp)
    os.remove(rt_tmpallfp)
    return


def get_true(param):
    command_GenPredData_true = [
        "GenPredData.exe",
        "--pfind_res",
        "pPred_data\\true.spectra",
        "--qvalue_threshold",
        "1",
        "--precursor_tol",
        "100",
        "--output_path",
        "pPred_data\\true",
        "--ms_path",
        param.rawfolder,
        "--instrument",
        param.instrument,
        "--nce",
        str(param.nce) + "," + str(param.nce) + "," + str(param.nce),
        "--output_all",
    ]
    result_GenPredData_true = subprocess.run(
        command_GenPredData_true, capture_output=True, text=True
    )
    print("GenPredData.exe get true RT MS2:", result_GenPredData_true.stdout)
    print("GenPredData.exe get true RT MS2:", result_GenPredData_true.stderr)

    num_lines = getallfile(
        "pPred_data\\true.spectra",
        "pPred_data\\true\\",
        "pPred_data\\true.ms2data",
        "pPred_data\\true.rtdata",
    )
    res_dict = read_spectra("pPred_data\\true.spectra")
    writefile(
        "pPred_data\\true.ms2data", "pPred_data\\true.rtdata", num_lines, res_dict
    )


def read_ms2_pred(fp1, maxlen, start, end):
    dict = {}
    psm_count = 0
    with open(fp1, "r") as file:
        lines = file.readlines()[start:end]
        for line in lines:
            line = line.rstrip("\n")
            line_lst = line.split("\t")
            i = 0
            while i < len(line_lst) and line_lst[i] != "":
                charge = int(line_lst[i])
                num = int(line_lst[i + 1])
                x = 0
                pred_vector_2d = np.zeros((maxlen, 12), dtype=float)
                while x < num:
                    j = i + 2 + 3 * x
                    bytype = int(line_lst[j])
                    idx = int(line_lst[j + 1])
                    intensity = float(line_lst[j + 2])
                    pred_vector_2d[idx, bytype] = intensity
                    x += 1
                i = j + 3
                if psm_count in dict.keys():
                    dict[psm_count][charge] = pred_vector_2d
                else:
                    dict[psm_count] = {}
                    dict[psm_count][charge] = pred_vector_2d
            psm_count += 1
    return dict


def read_ms2_true(fp1, maxlen, start, end):
    dict = {}
    i = 0
    psm_count = 0
    with open(fp1, "r") as file:
        lines = file.readlines()[start:end]
        for line in lines:
            line = line.rstrip("\n")
            line_lst = line.split("\t")
            charge = int(line_lst[2])
            line_lst = line_lst[4:]
            if len(line_lst) % 3 != 0:
                ValueError("ERROR!")
            i = 0
            pred_vector_2d = np.zeros((maxlen, 12), dtype=float)
            while i < len(line_lst):
                bytype = int(line_lst[i])
                idx = int(line_lst[i + 1])
                intensity = float(line_lst[i + 2])
                pred_vector_2d[idx, bytype] = intensity
                i = i + 3
            if psm_count in dict.keys():
                dict[psm_count][charge] = pred_vector_2d
            else:
                dict[psm_count] = {}
                dict[psm_count][charge] = pred_vector_2d
            psm_count += 1
    return dict


def get_pepmaxlen(fp1):
    maxlen = 0
    maxc = 0
    with open(fp1, "r") as file:
        lines = file.readlines()[1:]
        for line in lines:
            line_lst = line.split("\t")
            pep = line_lst[0]
            charge = int(line_lst[2])
            peplen = len(pep)
            if peplen > maxlen:
                maxlen = peplen
            if charge > maxc:
                maxc = charge
    return maxlen, len(lines) + 1


def cosine_sim(v1, v2):
    v1_flat = v1.flatten()
    v2_flat = v2.flatten()
    dot_product = np.dot(v1_flat, v2_flat)
    norm_v1 = np.linalg.norm(v1_flat)
    norm_v2 = np.linalg.norm(v2_flat)
    if norm_v1 * norm_v2 == 0:
        return 0
    similarity = dot_product / (norm_v1 * norm_v2)
    # vector1 = np.array(v1).reshape(
    #     1, -1
    # Convert to a 2D array because cosine_similarity requires 2D input.
    # vector2 = np.array(v2).reshape(1, -1)

    # Calculate cosine similarity.
    # similarity = cosine_similarity(vector1, vector2)
    return similarity


def spearman_sim(v1, v2):
    v1_flat = v1.flatten()
    v2_flat = v2.flatten()
    correlation, _ = spearmanr(v1_flat, v2_flat)
    return correlation


def pearson_sim(v1, v2):
    v1_flat = v1.flatten()
    v2_flat = v2.flatten()
    correlation, _ = pearsonr(v1_flat, v2_flat)
    return correlation


def cal_sim(true_dict, pred_dict):
    cos_lst = []
    spear_lst = []
    pears_lst = []
    for key1, value_dict in true_dict.items():
        for key2, value in value_dict.items():
            cos_lst.append(cosine_sim(value, pred_dict[key1][key2]))
            spear_lst.append(spearman_sim(value, pred_dict[key1][key2]))
            pears_lst.append(pearson_sim(value, pred_dict[key1][key2]))
    return cos_lst, spear_lst, pears_lst


def compute_pPred_features(
    cos_lst,
    spear_lst,
    pears_lst,
    cos_canlst,
    spear_canlst,
    pears_canlst,
    truertfp,
    orirtfp,
    canrtfp,
    start,
    end,
):
    """Compute one site-feature batch in memory without a feature-table file."""
    cos_result = [a - b for a, b in zip(cos_lst, cos_canlst)]
    spear_result = [a - b for a, b in zip(spear_lst, spear_canlst)]
    pears_result = [a - b for a, b in zip(pears_lst, pears_canlst)]

    truert_df = pd.read_csv(truertfp, sep="\t")
    orirt_df = pd.read_csv(orirtfp, sep="\t")
    canrt_df = pd.read_csv(canrtfp, sep="\t")
    truert_lst = list(truert_df["rt"])[start - 1 : end - 1]
    orirt_lst = list(orirt_df["pred"])[start - 1 : end - 1]
    canrt_lst = list(canrt_df["pred"])[start - 1 : end - 1]
    orirt_lst_ = [abs(a - b) for a, b in zip(truert_lst, orirt_lst)]
    canrt_lst_ = [abs(a - b) for a, b in zip(truert_lst, canrt_lst)]
    rt_lst = [a - b for a, b in zip(orirt_lst_, canrt_lst_)]

    return pd.DataFrame(
        {
            "delta_cosine": cos_result,
            "delta_spearman": spear_result,
            "delta_pearson": pears_result,
            "delta_rt": rt_lst,
            "cosine": cos_lst,
            "spearman": spear_lst,
            "pearson": pears_lst,
            "rt": orirt_lst_,
        }
    )

    
def process_pPred(param):
    original_working_directory = os.getcwd()
    working_directory = resource_path("pPred-MS2RT")
    for subfolder in ("pPred_data", "pPred_data/finetune", "pPred_data/true"):
        os.makedirs(os.path.join(working_directory, subfolder), exist_ok=True)
    os.chdir(working_directory)
    try:
        ms2_model_fp = param.ftmodelfolder + "\\" + param.ftmodelname + "_ms2.pt"
        rt_model_fp = param.ftmodelfolder + "\\" + param.ftmodelname + "_rt.pt"
        if not (os.path.exists(ms2_model_fp) and os.path.exists(rt_model_fp)):
            finetune_pPred(param)
        get_pFindRes(param)
        predict_pPred(param)
        get_true(param)
        maxlen, row_num = get_pepmaxlen("pPred_data\\true.ms2data")
        boundaries = list(range(1, row_num, 10000)) + [row_num]
        feature_batches = []
        for start, end in zip(boundaries, boundaries[1:]):
            truems2_dict = read_ms2_true("pPred_data\\true.ms2data", maxlen + 10, start, end)
            orims2_dict = read_ms2_pred("pPred_data\\ori_pred.pepdata.ms2.pred", maxlen + 10, start, end)
            cos_lst, spear_lst, pears_lst = cal_sim(truems2_dict, orims2_dict)
            canms2_dict = read_ms2_pred("pPred_data\\can_pred.pepdata.ms2.pred", maxlen + 10, start, end)
            cos_canlst, spear_canlst, pears_canlst = cal_sim(truems2_dict, canms2_dict)
            feature_batches.append(compute_pPred_features(
                cos_lst, spear_lst, pears_lst, cos_canlst, spear_canlst, pears_canlst,
                "pPred_data\\true.rtdata", "pPred_data\\ori_pred.pepdata.rt.pred",
                "pPred_data\\can_pred.pepdata.rt.pred", start, end,
            ))
        columns = ["delta_cosine", "delta_spearman", "delta_pearson", "delta_rt"]
        features = (pd.concat(feature_batches, ignore_index=True)[columns]
                    if feature_batches else pd.DataFrame(columns=columns))
        candidates = pd.read_csv(param.pSiteResPath, sep="\t")
        if len(features) != len(candidates):
            raise ValueError("pPred feature rows do not align with candidate site rows")
        pd.concat([candidates, features], axis=1).to_csv(param.pSiteResPath, sep="\t", index=False)
    finally:
        os.chdir(original_working_directory)


def clear_files(folder_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path) or os.path.islink(file_path):
            try:
                os.unlink(file_path)  # Delete the file or symbolic link.
                print(f"Deleted file: {file_path}")
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")


def delfiles():
    folder_lst = [
        resource_path("pPred-MS2RT\\pPred_data"),
        resource_path("pPred-MS2RT\\pPred_data\\finetune"),
        resource_path("pPred-MS2RT\\pPred_data\\true"),
    ]
    for folder_path in folder_lst:
        clear_files(folder_path)
