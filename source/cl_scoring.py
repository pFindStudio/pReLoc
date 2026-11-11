"""Adapt CL site predictions to the existing pReLoc result pipeline."""

import time

import numpy as np
import pandas as pd

from resources import resource_path


def predict_cl_sites(features, param):
    # Keep torch, ESM and their checkpoints out of the XGBoost startup path.
    from cl import main_predict

    features = features.copy().reset_index(drop=True)
    features["can_bm25_score"] = features["bm25_score"] - features["delta_bm25_score"]
    features["can_ion_ratio"] = features["ion_ratio"] - features["delta_ion_ratio"]
    features["can_intensity_ratio"] = (
        features["intensity_ratio"] - features["delta_intensity_ratio"]
    )
    config = {
        "model_path": resource_path(getattr(param, "cl_model_path", "best_model.pth")),
        "esm_path": resource_path(
            getattr(param, "cl_esm_path", "esm2_t33_650M_UR50D.pt")
        ),
        "device": getattr(param, "cl_device", "auto"),
        "batch_size": getattr(param, "cl_batch_size", 1),
        "num_workers": 0,
    }
    predictions = main_predict(features, config=config)
    if len(predictions) != len(features) or not predictions.index.equals(features.index):
        raise ValueError("CL predictions do not align with the candidate site rows")
    probabilities = predictions["pred_prob"].to_numpy(dtype=float)
    if not np.isfinite(probabilities).all() or ((probabilities < 0) | (probabilities > 1)).any():
        raise ValueError("CL returned invalid site probabilities")
    return probabilities


def _write_timing(param, start, prediction_start, prediction_end):
    with open(param.TimePath, "a", encoding="utf-8") as output:
        output.write(
            f"\nPredict Sites Score (CL):\nfull time: {time.time() - start} s"
            f"\npredict time: {prediction_end - prediction_start} s"
        )


def modScore_CL(param):
    from mods_score import final_process, finalize_mod_scores, read_IDfile

    start = time.time()
    read_IDfile(param, feature_mode="CL")
    features = pd.read_csv(param.pSiteResPath, sep="\t").drop_duplicates().reset_index(drop=True)
    if features.empty:
        final_process(None, param)
        return
    prediction_start = time.time()
    features["pSite2_score"] = predict_cl_sites(features, param)
    prediction_end = time.time()
    # Preserve each row's original target/decoy label; never merge duplicate site keys.
    finalize_mod_scores(features, param)
    _write_timing(param, start, prediction_start, prediction_end)


def siteScore_CL(param):
    from sites_score import read_IDfile

    start = time.time()
    read_IDfile(param, feature_mode="CL")
    features = pd.read_csv(param.pSiteResPath, sep="\t")
    prediction_start = time.time()
    if features.empty:
        features["pSite2_score"] = pd.Series(dtype=float)
    else:
        features["pSite2_score"] = predict_cl_sites(features, param)
    prediction_end = time.time()
    result = features.groupby(
        ["spec_name", "pep", "mod"], as_index=False, dropna=False
    ).agg({"pSite2_score": lambda values: ",".join(map(str, values))})
    result.rename(columns={"pep": "seq"}, inplace=True)
    result.to_csv(param.pSiteResPath, sep="\t", index=False)
    _write_timing(param, start, prediction_start, prediction_end)
