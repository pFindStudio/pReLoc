import argparse
from contextlib import contextmanager
import sys
from myparam import PARAM
from process_OtherRes import init_processOtherRes
from sites_score import siteScore_XGB, siteScore_XGB_pPred
from mods_score import final_process, modScore_XGB, modScore_XGB_pPred
import subprocess
import os
import time
import logging
import traceback
import tempfile
from write_param import (
    writepChemparam,
    update_pfind_modifications,
    fill_missing_mod_masses,
)


def setup_logging():
    log_dir = os.path.dirname(os.path.realpath(__file__))  # Application directory.
    log_file = os.path.join(log_dir, "pReLoc.log")  # Log file path.

    # Write logs to both the file and console.
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(log_file, mode="a", encoding="utf-8"),
            logging.StreamHandler(),
        ],
    )


def execute_command(command, expected_output=None):
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    while True:
        output = process.stdout.readline()
        if output == "" and process.poll() is not None:
            break
        if output:
            print(output.strip())
    _, error = process.communicate()
    if process.returncode == 0:
        if expected_output and not os.path.exists(expected_output):
            raise FileNotFoundError(
                f"Candidate generation exited successfully but did not create: {expected_output}"
            )
        print("Generate Candidates successfully!")
    else:
        raise RuntimeError(
            f"Candidate generation failed with exit code {process.returncode}: {error}"
        )

def resource_path(relative_path):
    try:
        # PyInstaller stores the bundled resource directory in _MEIPASS.
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.dirname(os.path.abspath(__file__))
    
    return os.path.join(base_path, relative_path)


@contextmanager
def generated_parameter_cleanup(param):
    """Remove an internally generated parameter file without deleting user input."""
    original_path = os.path.normcase(os.path.realpath(param.paramFilePath))
    generated_path = os.path.abspath(param.CparamFilePath)
    owns_generated_file = os.path.normcase(os.path.realpath(generated_path)) != original_path
    try:
        yield
    finally:
        if owns_generated_file and os.path.isfile(generated_path):
            os.remove(generated_path)


@contextmanager
def candidate_parameter_file(param):
    # Legacy phospho engines misread comments containing '=' as modification entries.
    # Give the native process a comment-free copy without changing the user's settings.
    with open(param.CparamFilePath, encoding="utf-8-sig") as source:
        lines = [line for line in source if not line.lstrip().startswith("#")]
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", prefix="pReLoc-native-", suffix=".txt",
        dir=os.path.dirname(os.path.abspath(param.CandidatesPath)), delete=False,
    ) as temporary:
        temporary.writelines(lines)
        native_path = temporary.name
    try:
        yield native_path
    finally:
        os.remove(native_path)

def configure_model(param, model_override=None):
    model = (model_override or getattr(param, "score_model", "XGBOOST")).upper()
    if model not in {"XGBOOST", "CL"}:
        raise ValueError("model must be XGBoost or CL")
    if model == "CL" and param.ispPred:
        raise ValueError("CL does not support ispPred=true; use ispPred=false or model=XGBoost")
    param.score_model = model
    return param


def score_results(param):
    model = getattr(param, "score_model", "XGBOOST")
    if model == "CL":
        from cl_scoring import modScore_CL, siteScore_CL

        scorer = siteScore_CL if param.score_isSite else modScore_CL
    elif param.ispPred:
        scorer = siteScore_XGB_pPred if param.score_isSite else modScore_XGB_pPred
    else:
        scorer = siteScore_XGB if param.score_isSite else modScore_XGB
    scorer(param)


def main(paramfilepath=None, model_override=None):
    a = time.time()
    if paramfilepath is None:
        paramfilepath = input("Input the path of the parameter file:")
    paramfilepath = paramfilepath.replace('"', "")
    paramfilepath = os.path.abspath(paramfilepath)
    if not os.path.exists(paramfilepath):
        print(
            f"Error: The file {paramfilepath} does not exist. Please check if the path is correct."
        )
    else:
        param = configure_model(PARAM(paramfilepath), model_override)
        with generated_parameter_cleanup(param):
            if param.ispFindOpen:
                if not param.ispFind:
                    raise ValueError("ispfindOpen=true requires ispfind=true and a pFind result file")
                update_pfind_modifications(paramfilepath, param.InputPath)
            fill_missing_mod_masses(
                paramfilepath, resource_path(os.path.join("local_denovo", "modification.ini"))
            )
            param = configure_model(PARAM(paramfilepath), model_override)
            if param.ispChem:
                writepChemparam(paramfilepath, param.CparamFilePath, param.InputSummary)
                param = configure_model(PARAM(param.CparamFilePath), model_override)
            print(f"Scoring model: {param.score_model}")
            # Preprocess external identification formats.
            # Write the converted input to AutoCustomInput.tsv.
            has_relocation_candidates = init_processOtherRes(param)

            param1 = param
            if param.score_isMod and not has_relocation_candidates:
                final_process(None, param1)
            else:
                pho_flag = False
                for x in param.scoreModName:
                    if "Phospho" in x:
                        pho_flag = True
                # Generate candidate peptides and baseline features with the native engine.
                original_working_directory = os.getcwd()
                if pho_flag:
                    executablepath = r".\loc_denovo_pho.exe"
                else:
                    executablepath = r".\loc_denovo.exe"
                working_directory = resource_path("local_denovo")
                os.chdir(working_directory)
                print("Start Generating Candidates ...\n")
                try:
                    with candidate_parameter_file(param) as native_path:
                        execute_command([executablepath, native_path], param.CandidatesPath)
                finally:
                    os.chdir(original_working_directory)

                param1 = configure_model(PARAM(param.CparamFilePath), model_override)
                score_results(param1)

            b = time.time()
            with open(param1.TimePath, "a") as tf:
                tf.write("\n\nFull Process:" + str(b - a) + " s\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="pReLoc with XGBoost or CL scoring")
    parser.add_argument("param_file", nargs="?", help="parameter file; prompts when omitted")
    parser.add_argument("--model", type=str.upper, choices=["XGBOOST", "CL"],
                        help="override [score] model for this run")
    args = parser.parse_args()
    setup_logging()
    try:
        main(args.param_file, args.model)
    except Exception as e:
        logging.error("Running error\n", exc_info=True)
        print(f"ERROR:\n{e}")
        traceback.print_exc()
        sys.exit(1)
    finally:
        if args.param_file is None:
            input("Program ends, press Enter to exit...")
