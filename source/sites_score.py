import pandas as pd
import joblib
import time
from processpPred import process_pPred, delfiles
import os
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
                elif line.startswith("s"):
                    line_vec = line.split("=")
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
    return

def resource_path(relative_path):
    try:
        # PyInstaller stores the bundled resource directory in _MEIPASS.
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.dirname(os.path.abspath(__file__))
    
    return os.path.join(base_path, relative_path)

def siteScore_XGB(param):
    a = time.time()
    read_IDfile(param)

    load_model = joblib.load(resource_path("XGBoost_final.pkl"))

    test_df = pd.read_csv(param.pSiteResPath, sep="\t")

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
    columns_to_drop = [
        "bm25_score",
        "ion_ratio",
        "intensity_ratio",
        "delta_bm25_score",
        "delta_ion_ratio",
        "delta_intensity_ratio",
        "aa_position",
        "aa_type",
        "pep_length",
        "bm25_pep",
        "bm25_mod",
        "ion_pep",
        "ion_mod",
        "intensity_pep",
        "intensity_mod",
    ]
    test_df.drop(columns=columns_to_drop, inplace=True)
    test_df = test_df.groupby(
        ["spec_name", "pep", "mod"], as_index=False, dropna=False
    ).agg({"pSite2_score": lambda x: ",".join(map(str, x))})
    test_df.rename(columns={"pep": "seq"}, inplace=True)
    test_df.to_csv(param.pSiteResPath, sep="\t", index=False)
    b = time.time()
    with open(param.TimePath, "a") as tf:
        tf.write(
            "\nPredict Sites Score:\nfull time: "
            + str(b - a)
            + " s\npredict time: "
            + str(b_ - a_)
            + " s"
        )


def siteScore_XGB_pPred(param):
    a = time.time()
    read_IDfile(param)
    process_pPred(param)

    load_model = joblib.load(resource_path("XGBoost_pPred.pkl"))

    load_model1 = joblib.load(resource_path("XGBoost.pkl"))

    test_df = pd.read_csv(param.pSiteResPath, sep="\t")

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

    # no pPred
    selected_columns = list(range(4, 11)) + list(range(-31, 0))
    test_X = df_test.iloc[:, selected_columns].values

    a_ = time.time()
    test_preds = load_model1.predict_proba(test_X)[:, 1]
    b_ = time.time()
    test_df["pSite2_score"] = test_preds

    # pPred
    selected_columns = list(range(4, 11)) + list(range(19, 23)) + list(range(-31, 0))
    test_X = df_test.iloc[:, selected_columns].values

    a_ = time.time()
    test_preds = load_model.predict_proba(test_X)[:, 1]
    b_ = time.time()
    test_df["pSite2_score_pPred"] = test_preds
    columns_to_drop = [
        # "bm25_score",
        # "ion_ratio",
        # "intensity_ratio",
        "delta_bm25_score",
        "delta_ion_ratio",
        "delta_intensity_ratio",
        "aa_position",
        "aa_type",
        "pep_length",
        "bm25_pep",
        "bm25_mod",
        "ion_pep",
        "ion_mod",
        "intensity_pep",
        "intensity_mod",
        "delta_cosine",
        "delta_spearman",
        "delta_pearson",
        "delta_rt",
    ]
    test_df.drop(columns=columns_to_drop, inplace=True)
    test_df = test_df.groupby(
        ["spec_name", "pep", "mod", "bm25_score", "ion_ratio", "intensity_ratio"],
        as_index=False,
        dropna=False,
    ).agg(
        {
            "pSite2_score": lambda x: ",".join(map(str, x)),
            "pSite2_score_pPred": lambda x: ",".join(map(str, x)),
        }
    )
    test_df.rename(columns={"pep": "seq"}, inplace=True)
    test_df.to_csv(param.pSiteResPath, sep="\t", index=False)
    delfiles()
    b = time.time()
    with open(param.TimePath, "a") as tf:
        tf.write(
            "\nPredict Sites Score:\nfull time: "
            + str(b - a)
            + " s\npredict time: "
            + str(b_ - a_)
            + " s"
        )
