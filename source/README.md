# pReLoc source distribution

This directory contains the Python source for pReLoc, its scoring checkpoints, a 100-PSM example, and the files needed to build a standalone Windows application. pReLoc calls local de novo for candidate generation and optionally calls pPred for retention-time and MS/MS prediction features. Both native components are supplied as Windows executables with their runtime assets.

## Install and run the example

Use 64-bit Windows and install the Microsoft Visual C++ x64 Redistributable if it is missing. Source validation used Windows 11 and Python 3.13.11. The direct dependency versions from that environment are in [requirements-verified.txt](requirements-verified.txt).

To download just the source distribution, run the following in PowerShell with Git and Git LFS installed:

```powershell
git lfs install
$env:GIT_LFS_SKIP_SMUDGE = '1'
git clone --sparse https://github.com/pFindStudio/pReLoc.git
cd pReLoc
git sparse-checkout set source
Remove-Item Env:GIT_LFS_SKIP_SMUDGE
git lfs pull --include="source/**"
cd source
python -m pip install -r requirements.txt
python run_demo.py --model XGBoost
```

For the basic XGBoost demo, you can omit the optional pPred runtime from the LFS download by adding `--exclude="source/pPred-MS2RT/**"` to `git lfs pull`. The demo resolves all input paths and writes `results/xgboost/pReLocRes.tsv`, `pReLocRes-filtered.tsv`, `DirectRes.tsv`, and `Time.txt`. The main result should contain 100 rows plus the header. Both ordinary XGBoost and CL can run on a CPU.

## Run CL

Install the CL dependencies and download the official ESM-2 650M checkpoint into this directory. The checkpoint is 2,604,537,549 bytes and exceeds this repository's single-file Git LFS limit. The bundled contact-regression checkpoint belongs beside it.

```powershell
python -m pip install -r requirements-cl.txt
curl.exe --fail --location --retry 3 --output esm2_t33_650M_UR50D.pt https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt
if ($LASTEXITCODE -ne 0) { throw 'ESM2 download failed' }
$expected = 'ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c'
if ((Get-FileHash -Algorithm SHA256 .\esm2_t33_650M_UR50D.pt).Hash -ne $expected) { throw 'ESM2 checksum mismatch' }
python run_demo.py --model CL
```

The CL demo writes to `results/cl/`. A previous CPU reference run took approximately six minutes for the 100-PSM example. CL uses `best_model.pth` for scoring and ESM-2 for sequence embeddings.

## Run your own data

Copy [param.txt](param.txt), set absolute paths to your input and output directories, and select exactly one of `isSite=true` and `isMod=true`. Create the output directory before running the application. For custom identifications, use a TSV with the same structure as [examples/data_input.tsv](examples/data_input.tsv), and match `spec_name` to the MGF `TITLE`. For pFind input, set `ispfind=true` and provide the pFind `.spectra` file and its MGF data. For pChem, use `isFLR=0` and supply the corresponding pChem summary file.

```powershell
python pReLoc.py "C:\analysis\param.txt"
python pReLoc.py "C:\analysis\param.txt" --model CL
```

### Optional pPred features

Use `model=XGBoost`, `ispPred=true`, and `isFLR=1` (or `isFLR=0` to disable FLR filtering). Set `pFindResfp` to the pFind result file, `rawfolder` to the matching parsed mass-spectrometry inputs (`ms1`/`ms2` or `pf1`/`pf2` and their indexes), and `modelfolder` and `modelname` to the fine-tuned RT/MS2 checkpoint location and basename. The three base checkpoints in `pPred-MS2RT/pPred_data/mainmodel/` are used for fine-tuning when the matching fine-tuned checkpoints are absent. Keep `pPred-MS2RT/` together and writable. Its native runtime includes CUDA 11 and cuDNN 8 libraries.

The repository's [pPred example](../examples/mod_xgbpPred) supplies matched data and `demofinetune` checkpoints. For that example, use `instrument=2`, `nce=30`, `rtparam=10_130`, and `ms2param=0_30_30_30_2_4`. The main scoring input and MGF can use the two files in `source/examples/`. If you used the sparse checkout above, add the example with `git sparse-checkout add examples/mod_xgbpPred` and retrieve its data with `git lfs pull --include="examples/mod_xgbpPred/**"`.

## Build the standalone Windows application

Download all `source/` LFS assets and the ESM-2 checkpoint before building. From this directory:

```powershell
python -m pip install -r requirements-build.txt
.\build_exe.bat
```

Keep `dist/pReLoc/pReLoc.exe` and its `_internal/` directory together. The build uses PyInstaller's directory format to keep the large weights outside the executable.

## Models and native components

| Component | Files |
| --- | --- |
| pReLoc XGBoost scoring | `XGBoost_final.pkl`, `XGBoost_pPred.pkl`, `XGBoost.pkl` |
| pReLoc CL scoring | `best_model.pth`, the official ESM-2 checkpoint above, and `esm2_t33_650M_UR50D-contact-regression.pt` |
| local de novo candidate generation | `local_denovo/loc_denovo.exe`, `loc_denovo_pho.exe`, and the adjacent `.ini` files |
| pPred prediction and feature extraction | `pPred-MS2RT/pPred.exe`, `GenPredData.exe`, their runtime libraries and configuration files, and three base checkpoints in `pPred_data/mainmodel/` |

`processpPred.py` is pReLoc's adapter. It invokes the external executables and prepares the prediction features for pReLoc scoring.
