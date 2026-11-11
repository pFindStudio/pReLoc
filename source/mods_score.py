import pandas as pd
import joblib
import time
from sklearn.preprocessing import MinMaxScaler
from processpPred import process_pPred, delfiles
import os
import numpy as np
import sys


def one_hot_encode(letter):
    # The 26 uppercase letters.
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    # Find the letter index in the alphabet.
    index = alphabet.index(letter)
    # Create the one-hot vector.
    one_hot_vector = [1 if char == letter else 0 for char in alphabet]
    return one_hot_vector


def read_IDfile(param, feature_mode="XGBOOST"):
    with open(param.pSiteResPath, "w") as modelinputfile:
        modelinputfile.write(
            "spec_name\tpep\tmod\tscore\tbm25_score\tion_ratio\tintensity_ratio\tdelta_bm25_score\tdelta_ion_ratio\tdelta_intensity_ratio\taa_position\taa_type\tpep_length\tbm25_pep\tbm25_mod\tion_pep\tion_mod\tintensity_pep\tintensity_mod\n"
        )
        with open(param.CandidatesPath, "r") as file:
            for line in file:
                line = line.strip()
                if line.startswith("psm"):
                    line_vec = line.split("\t")
                    # key = line_vec[0][-1] - "0"
                    spec_name = line_vec[1]
                    pep = line_vec[2]
                    mod = line_vec[3]
                    score = line_vec[4]
                    scoreModindex = []
                    if mod != "":
                        mod_list = mod.split(";")[:-1]
                        for modinfo in mod_list:
                            modindex = int(modinfo.split(",")[0]) - 1
                            modname = modinfo.split(",")[1]
                            if modname in param.scoreModName:
                                if modindex == -1:
                                    scoreModindex.append(0)
                                elif "C-term" in modname:
                                    scoreModindex.append(modindex - 1)
                                else:
                                    scoreModindex.append(modindex)

                elif line.startswith("s") and scoreModindex:
                    line_vec = line.split("=")
                    index = int(line_vec[0][1:])
                    if index in scoreModindex:
                        value_vec = line_vec[1].split("\t")
                        feature_vec = value_vec[:9]
                        # feature_vec = list(map(float, feature_vec))
                        if feature_mode == "CL":
                            aa_position = str(int(feature_vec[-3]))
                            aa_type = chr(int(feature_vec[-2]))
                            pep_length = int(feature_vec[-1])
                        else:
                            aa_position = str(int(feature_vec[-3]) / int(feature_vec[-1]))
                            aa_type = one_hot_encode(chr(int(feature_vec[-2])))
                            if int(feature_vec[-1]) <= 5:
                                pep_length = [1, 0, 0, 0, 0]
                            elif int(feature_vec[-1]) <= 10:
                                pep_length = [0, 1, 0, 0, 0]
                            elif int(feature_vec[-1]) <= 15:
                                pep_length = [0, 0, 1, 0, 0]
                            elif int(feature_vec[-1]) <= 20:
                                pep_length = [0, 0, 0, 1, 0]
                            else:
                                pep_length = [0, 0, 0, 0, 1]

                        if len(value_vec) == 14:
                            intensity_mod = ""
                        elif len(value_vec) == 15:
                            intensity_mod = value_vec[14]
                        else:
                            KeyError("ERROR!")

                        modelinputfile.write(
                            spec_name
                            + "\t"
                            + pep
                            + "\t"
                            + mod
                            + "\t"
                            + score
                            + "\t"
                            + "\t".join(feature_vec[:-3])
                            + "\t"
                            + aa_position
                            + "\t"
                            + str(aa_type)
                            + "\t"
                            + str(pep_length)
                            + "\t"
                            + value_vec[9]
                            + "\t"
                            + value_vec[10]
                            + "\t"
                            + value_vec[11]
                            + "\t"
                            + value_vec[12]
                            + "\t"
                            + value_vec[13]
                            + "\t"
                            + intensity_mod
                            + "\n"
                        )
    if os.path.exists(param.CandidatesPath):
        os.remove(param.CandidatesPath)
        os.remove(param.AutoCustomInputPath)
        if os.path.abspath(param.CparamFilePath) != os.path.abspath(param.paramFilePath):
            os.remove(param.CparamFilePath)
    return

def resource_path(relative_path):
    try:
        # PyInstaller stores the bundled resource directory in _MEIPASS.
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.dirname(os.path.abspath(__file__))
    
    return os.path.join(base_path, relative_path)

def modScore_XGB(param):
    a = time.time()
    read_IDfile(param)

    load_model = joblib.load(resource_path("XGBoost_final.pkl"))
    test_df = pd.read_csv(param.pSiteResPath, sep="\t")
    test_df = test_df.drop_duplicates(keep="first")

    aa_type = test_df["aa_type"].apply(lambda x: [int(i) for i in x[1:-1].split(", ")])
    pep_length = test_df["pep_length"].apply(
        lambda x: [int(i) for i in x[1:-1].split(", ")]
    )
    df_test = pd.concat(
        [
            test_df,
            pd.get_dummies(aa_type.apply(pd.Series), prefix="aa_type"),
        ],
        axis=1,
    )
    df_test = pd.concat(
        [
            df_test,
            pd.get_dummies(pep_length.apply(pd.Series), prefix="pep_length"),
        ],
        axis=1,
    )
    df_test = df_test.drop(["aa_type", "pep_length"], axis=1)
    selected_columns = list(range(4, 11)) + list(range(-31, 0))
    test_X = df_test.iloc[:, selected_columns].values

    a_ = time.time()
    test_preds = load_model.predict_proba(test_X)[:, 1]
    b_ = time.time()
    test_df["pSite2_score"] = test_preds

    finalize_mod_scores(test_df, param)
    b = time.time()
    with open(param.TimePath, "a") as tf:
        tf.write(
            "\nPredict Sites Score:\nfull time: "
            + str(b - a)
            + " s\npredict time: "
            + str(b_ - a_)
            + " s"
        )


def finalize_mod_scores(test_df, param):
    test_df["pSite2_score_modsite"] = test_df.groupby(
        ["spec_name", "pep", "mod", "score"]
    )["pSite2_score"].transform(lambda x: ",".join(x.astype(str)))

    x_namelst = []
    for x in param.scoreModName:
        if x.split("[")[0] not in x_namelst:
            x_namelst.append(x.split("[")[0])
    test_df["scoremodsite_num"] = test_df["mod"].apply(
        lambda mod_str: sum(mod_str.count(name) for name in x_namelst)
    )

    # scaler = MinMaxScaler()
    # test_df["normalized_pSite2_score"] = test_df.groupby(["spec_name", "pep", "mod"])[
    #     "pSite2_score"
    # ].transform(
    #     lambda x: (
    #         scaler.fit_transform(x.values.reshape(-1, 1)).flatten() if len(x) > 1 else x
    #     )
    # )

    test_df = (
        test_df.groupby(
            [
                "spec_name",
                "pep",
                "mod",
                "score",
                "pSite2_score_modsite",
                "bm25_score",
                "scoremodsite_num",
            ]
        )["pSite2_score"]
        .agg("min")
        .reset_index()
    )
    # test_df = (
    #     test_df.groupby(
    #         ["spec_name", "pep", "mod", "score", "pSite2_score_modsite", "bm25_score"]
    #     )["pSite2_score"]
    #     .agg("mean")
    #     .reset_index()
    # )
    # Harmonic mean.
    # test_df = (
    #     test_df.groupby(
    #         ["spec_name", "pep", "mod", "score", "pSite2_score_modsite", "bm25_score"]
    #     )["pSite2_score"]
    #     .agg(lambda x: len(x) / sum(1 / i for i in x))
    #     .reset_index()
    # )
    # grouped = test_df.groupby(["spec_name", "pep"])["pSite2_score"].transform("sum")
    # Calculate normalized values.
    # test_df["normalized_pSite2_score"] = test_df["pSite2_score"] / grouped
    if param.isFLRest:
        test_df.rename(columns={"score": "target_label"}, inplace=True)
    # Legacy debugging output for development builds.
    # filepath_o = param.pSiteResPath[:-4] + "-original.tsv"
    # test_df.to_csv(filepath_o, sep="\t", index=False)
    final_process(test_df, param)


def modScore_XGB_pPred(param):
    a = time.time()
    read_IDfile(param)
    process_pPred(param)
    load_model = joblib.load(resource_path("XGBoost_pPred.pkl"))

    test_df = pd.read_csv(param.pSiteResPath, sep="\t")
    test_df = test_df.groupby(["spec_name", "pep", "mod"]).first().reset_index()

    aa_type = test_df["aa_type"].apply(lambda x: [int(i) for i in x[1:-1].split(", ")])
    pep_length = test_df["pep_length"].apply(
        lambda x: [int(i) for i in x[1:-1].split(", ")]
    )
    df_test = pd.concat(
        [
            test_df,
            pd.get_dummies(aa_type.apply(pd.Series), prefix="aa_type"),
        ],
        axis=1,
    )
    df_test = pd.concat(
        [
            df_test,
            pd.get_dummies(pep_length.apply(pd.Series), prefix="pep_length"),
        ],
        axis=1,
    )
    df_test = df_test.drop(["aa_type", "pep_length"], axis=1)
    selected_columns = list(range(4, 11)) + list(range(19, 23)) + list(range(-31, 0))
    test_X = df_test.iloc[:, selected_columns].values

    a_ = time.time()
    test_preds = load_model.predict_proba(test_X)[:, 1]
    b_ = time.time()
    test_df["pSite2_score"] = test_preds

    test_df["pSite2_score_modsite"] = test_df.groupby(
        ["spec_name", "pep", "mod", "score"]
    )["pSite2_score"].transform(lambda x: ",".join(x.astype(str)))

    for x in param.scoreModName:
        test_df["scoremodsite_num"] = test_df["mod"].str.count(x)

    test_df = (
        test_df.groupby(
            [
                "spec_name",
                "pep",
                "mod",
                "score",
                "pSite2_score_modsite",
                "bm25_score",
                "scoremodsite_num",
            ]
        )["pSite2_score"]
        .agg("min")
        .reset_index()
    )

    if param.isFLRest:
        test_df.rename(columns={"score": "target_label"}, inplace=True)
    final_process(test_df, param)
    # delfiles()
    b = time.time()
    with open(param.TimePath, "a") as tf:
        tf.write(
            "\nPredict Sites Score:\nfull time: "
            + str(b - a)
            + " s\npredict time: "
            + str(b_ - a_)
            + " s"
        )


# def phomodScore_XGB(param):
#     a = time.time()
#     read_IDfile(param)

#     load_model = joblib.load("XGBoost_pho.pkl")

#     test_df = pd.read_csv(param.pSiteResPath, sep="\t")
#     test_df = test_df.groupby(["spec_name", "pep", "mod"]).first().reset_index()

#     aa_type = test_df["aa_type"].apply(lambda x: [int(i) for i in x[1:-1].split(", ")])
#     pep_length = test_df["pep_length"].apply(
#         lambda x: [int(i) for i in x[1:-1].split(", ")]
#     )
#     df_test = pd.concat(
#         [
#             test_df,
#             pd.get_dummies(aa_type.apply(pd.Series), prefix="aa_type"),
#         ],
#         axis=1,
#     )
#     df_test = pd.concat(
#         [
#             df_test,
#             pd.get_dummies(pep_length.apply(pd.Series), prefix="pep_length"),
#         ],
#         axis=1,
#     )
#     df_test = df_test.drop(["aa_type", "pep_length"], axis=1)
#     selected_columns = list(range(4, 11)) + list(range(-31, 0))
#     test_X = df_test.iloc[:, selected_columns].values

#     a_ = time.time()
#     test_preds = load_model.predict_proba(test_X)[:, 1]
#     b_ = time.time()
#     test_df["pSite2_score"] = test_preds

#     test_df["pSite2_score_modsite"] = test_df.groupby(
#         ["spec_name", "pep", "mod", "score"]
#     )["pSite2_score"].transform(lambda x: ",".join(x.astype(str)))

#     for x in param.scoreModName:
#         test_df["scoremodsite_num"] = test_df["mod"].str.count(x)

#     test_df = (
#         test_df.groupby(
#             [
#                 "spec_name",
#                 "pep",
#                 "mod",
#                 "score",
#                 "pSite2_score_modsite",
#                 "bm25_score",
#                 "scoremodsite_num",
#             ]
#         )["pSite2_score"]
#         .agg("min")
#         .reset_index()
#     )

#     if param.isFLRest:
#         test_df.rename(columns={"score": "target_label"}, inplace=True)
#     final_process(test_df, param)
#     # test_df.to_csv(param.pSiteResPath, sep="\t", index=False)
#     b = time.time()
#     with open(param.TimePath, "a") as tf:
#         tf.write(
#             "\nPredict Sites Score:\nfull time: "
#             + str(b - a)
#             + " s\npredict time: "
#             + str(b_ - a_)
#             + " s"
#         )


def softmax(x):
    e_x = np.exp(x - np.max(x))
    return e_x / e_x.sum()


def sitepep_bm25(test_df, stopthreshold, maxiternum, param):
    if param.isFLRest:
        test_df = (
            test_df.sort_values("target_label", ascending=False)
            .groupby(["spec_name", "pep", "mod"])
            .head(1)
        )
        test_df = test_df.sort_values(
            by=["spec_name", "target_label"], ascending=[True, False]
        )

    test_df["bm25_score_softmax"] = test_df.groupby("spec_name")[
        "bm25_score"
    ].transform(lambda x: softmax(x))
    if "sitepep_bm25" not in test_df.columns:
        test_df["product"] = test_df["bm25_score_softmax"] * test_df["pSite2_score"]
        grouped = test_df.groupby(["spec_name"]).agg({"product": "sum"}).reset_index()
        test_df = pd.merge(test_df, grouped, on=["spec_name"], suffixes=("", "_sum"))
        res = (test_df["product"]) / test_df["product_sum"]
        test_df["sitepep_bm25"] = res
        del test_df["product"]
        del test_df["product_sum"]

    test_df.fillna({"sitepep_bm25": 0}, inplace=True)
    test_df["sitepep_bm25_1"] = test_df["sitepep_bm25"]

    mean_df = test_df.loc[test_df.groupby(["spec_name"])["sitepep_bm25_1"].idxmax()]
    mean0 = mean_df["sitepep_bm25_1"].values.mean()
    mean1 = mean0 + stopthreshold + 1
    num = 1

    while abs(mean1 - mean0) > stopthreshold and num < maxiternum:
        num += 1
        if num != 2:
            mean0 = mean1
        test_df["product"] = (
            test_df["bm25_score_softmax"] * test_df["sitepep_bm25_" + str(num - 1)]
        )
        grouped = test_df.groupby(["spec_name"]).agg({"product": "sum"}).reset_index()
        test_df = pd.merge(test_df, grouped, on=["spec_name"], suffixes=("", "_sum"))
        res = (test_df["product"]) / test_df["product_sum"]
        test_df["sitepep_bm25_" + str(num)] = res
        del test_df["product"]
        del test_df["product_sum"]
        del test_df["sitepep_bm25_" + str(num - 1)]
        test_df.fillna({"sitepep_bm25_" + str(num): 0}, inplace=True)
        mean_df = test_df.loc[
            test_df.groupby(["spec_name"])["sitepep_bm25_" + str(num)].idxmax()
        ]
        mean1 = mean_df["sitepep_bm25_" + str(num)].values.mean()
    del test_df["sitepep_bm25"]
    del test_df["sitepep_bm25_" + str(num)]
    test_df["sitepep_bm25_iter"] = res
    test_df.fillna({"sitepep_bm25_iter": 0}, inplace=True)
    return test_df


def smooth_list(lst):
    max_right = lst[-1]
    for i in range(len(lst) - 2, -1, -1):
        if lst[i] > max_right:
            lst[i] = max_right
        else:
            max_right = lst[i]
    return lst


def find_first_greater(arr, tar):
    left, right = 0, len(arr) - 1
    while left <= right:
        mid = (left + right) // 2
        if arr[mid] > tar:
            right = mid - 1
        else:
            left = mid + 1
    return left - 1


def FLRthreshold(df_tail, df, FLRthresh, param):
    df["sitepep_bm25_iter_normalized"] = (
        df["sitepep_bm25_iter"] - df["sitepep_bm25_iter"].min()
    ) / (df["sitepep_bm25_iter"].max() - df["sitepep_bm25_iter"].min())

    df = df.sort_values(
        by=["sitepep_bm25_iter_normalized", "target_label"],
        ascending=[False, False],
    )
    target_label = list(df["target_label"])

    count_true = 0
    FP_estimate = 0
    FLR_estimate_list = []
    count = 0
    PSMs_num = []
    for i in range(len(target_label)):
        count += 1
        if target_label[i] != -1:
            count_true += 1
            if FLR_estimate_list:
                FLR_estimate_list.append(FLR_estimate_list[-1])
            else:
                FLR_estimate_list.append(0)
            PSMs_num.append(count)
        if target_label[i] == -1:
            FP_estimate += 1
            FLR_estimate = FP_estimate / count_true
            FLR_estimate_list.append(FLR_estimate)
            PSMs_num.append(count)

    df_tail = df_tail.sort_values(
        by=["sitepep_bm25_iter", "target_label"],
        ascending=[False, False],
    )
    df_tail = df_tail.tail(int(0.75 * len(df_tail)))
    decoy_num = len(df_tail[df_tail["target_label"] == -1])
    target_num = len(df_tail[df_tail["target_label"] != -1])
    c = (target_num) / decoy_num
    FLR_estimate_list = [x * 100 * c for x in FLR_estimate_list]
    FLR_estimate_list = smooth_list(FLR_estimate_list)

    index = find_first_greater(FLR_estimate_list, FLRthresh * 100)
    df["FLR"] = FLR_estimate_list
    flr_df = df.copy()
    filtered_df = df.head(PSMs_num[index])
    filtered_df = filtered_df[filtered_df["target_label"] != -1].copy()
    del filtered_df["sitepep_bm25_iter_normalized"]
    return filtered_df, flr_df


def restore_original_sequences(result_df, param):
    if result_df.empty:
        return result_df
    if param.ispFind:
        input_df = pd.read_csv(
            param.InputPath, sep="\t", usecols=["File_Name", "Sequence"]
        ).rename(columns={"File_Name": "spec_name", "Sequence": "seq"})
    else:
        input_df = pd.read_csv(
            param.InputPath, sep="\t", usecols=["spec_name", "seq"]
        )
    sequence_map = input_df.drop_duplicates("spec_name").set_index("spec_name")["seq"]
    result_df["seq"] = result_df["spec_name"].map(sequence_map)
    return result_df


def final_process(output_df, param):
    if output_df is None:
        result_columns = [
            "spec_name",
            "seq",
            "scoremodsite_num",
            "ori_mod",
            "ori_pepScore",
            "ori_pReLocsiteScore",
            "ori_pReLocScore",
            "mod",
            "pepScore",
            "pReLocsiteScore",
            "pReLocScore",
        ]
        main_result_columns = result_columns.copy()
        if param.isFLRest:
            main_result_columns.append("FLR")
        empty_result = pd.DataFrame(columns=main_result_columns)
        empty_result.to_csv(param.pSiteResPath, sep="\t", index=False)
        if param.isFLRest:
            pd.DataFrame(columns=result_columns).to_csv(
                param.pSiteResPath[:-4] + "-filtered.tsv", sep="\t", index=False
            )
        return
    output_df = sitepep_bm25(output_df, 0.001, 100, param)
    # output_df.to_csv(param.pSiteResPath[:-4] + "-all.tsv", sep="\t", index=False)
    if not param.isFLRest:
        output_df0 = output_df.loc[output_df["score"] != -1].copy()
        del output_df0["score"]
    else:
        output_df0 = output_df.loc[output_df["target_label"] == 1].copy()
        del output_df0["target_label"]
    del output_df0["pSite2_score"]

    output_df0.rename(
        columns={
            "mod": "ori_mod",
            "bm25_score": "ori_pepScore",
            "sitepep_bm25_iter": "ori_pReLocScore",
            "pSite2_score_modsite": "ori_pReLocsiteScore",
        },
        inplace=True,
    )

    output_df1 = output_df.loc[
        output_df.groupby(["spec_name"])["sitepep_bm25_iter"].idxmax()
    ].copy()

    if param.isFLRest:
        filtered_df, competition_df = FLRthreshold(
            output_df, output_df1, param.FLRthresh, param
        )
        target_flr_df = competition_df.loc[
            competition_df["target_label"] != -1, ["spec_name", "FLR"]
        ]
        target_df = output_df.loc[output_df["target_label"] != -1]
        output_df1 = target_df.loc[
            target_df.groupby("spec_name")["sitepep_bm25_iter"].idxmax()
        ].copy()
        output_df1 = pd.merge(
            output_df1,
            target_flr_df,
            on="spec_name",
            how="left",
            validate="one_to_one",
        )
        del output_df1["target_label"]
        del filtered_df["target_label"]
        del filtered_df["pSite2_score"]
        filtered_df.rename(
            columns={
                "bm25_score": "pepScore",
                "sitepep_bm25_iter": "pReLocScore",
                "pSite2_score_modsite": "pReLocsiteScore",
            },
            inplace=True,
        )
        mergedfiltered_df = pd.merge(
            output_df0,
            filtered_df,
            on=["spec_name", "pep", "scoremodsite_num"],
            how="inner",
        )
        mergedfiltered_df = mergedfiltered_df[
            [
                "spec_name",
                "pep",
                "scoremodsite_num",
                "ori_mod",
                "ori_pepScore",
                "ori_pReLocsiteScore",
                "ori_pReLocScore",
                "mod",
                "pepScore",
                "pReLocsiteScore",
                "pReLocScore",
                "FLR",
            ]
        ]
        mergedfiltered_df = mergedfiltered_df.rename(columns={"pep": "seq"})
        mergedfiltered_df = restore_original_sequences(mergedfiltered_df, param)
        mergedfiltered_df.sort_values(
            ["FLR", "pReLocScore"],
            ascending=[True, False],
            kind="stable",
            inplace=True,
        )
        del mergedfiltered_df["FLR"]
        filepath2 = param.pSiteResPath[:-4] + "-filtered.tsv"
        mergedfiltered_df.to_csv(filepath2, sep="\t", index=False)
    else:
        del output_df1["score"]

    del output_df1["pSite2_score"]
    # del output_df1["sitepep_bm25_iter_normalized"]
    output_df1.rename(
        columns={
            "bm25_score": "pepScore",
            "sitepep_bm25_iter": "pReLocScore",
            "pSite2_score_modsite": "pReLocsiteScore",
        },
        inplace=True,
    )
    merged_df = pd.merge(
        output_df0,
        output_df1,
        on=["spec_name", "pep", "scoremodsite_num"],
        how="inner",
    )
    result_columns = [
        "spec_name",
        "pep",
        "scoremodsite_num",
        "ori_mod",
        "ori_pepScore",
        "ori_pReLocsiteScore",
        "ori_pReLocScore",
        "mod",
        "pepScore",
        "pReLocsiteScore",
        "pReLocScore",
    ]
    if param.isFLRest:
        result_columns.append("FLR")
    merged_df = merged_df[result_columns]
    merged_df = merged_df.rename(columns={"pep": "seq"})
    merged_df = restore_original_sequences(merged_df, param)
    if param.isFLRest:
        merged_df.sort_values(
            ["FLR", "pReLocScore"],
            ascending=[True, False],
            kind="stable",
            inplace=True,
        )
    merged_df.to_csv(param.pSiteResPath, sep="\t", index=False)
