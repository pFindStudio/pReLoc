# pReLoc

A framework for site-level confidence assessment in proteomics, supporting amino-acid site scoring and modification-site localization with XGBoost or contrastive learning (CL).

This Windows distribution combines both scoring backends in **`pReLoc.exe`**.

This repository provides **compiled standalone software** together with a **small real demonstration dataset** (`examples/data_input.tsv` and `examples/dataMGF.mgf`, 100 peptide-spectrum matches, or **PSMs**). Each PSM links an identified peptide to a spectrum.

## Start here

1. **Download the application** using the [installation guide](#installation-guide), including its Git LFS files. Keep `pReLoc.exe` and `_internal/` together.
2. **Choose a workflow** from the [template table](#choose-a-workflow). For a first run, use `examples/pReLocparam_mod_xgb.txt` with the small example dataset.
3. **Edit a copy of the parameter file:** update the input paths and `OutputFolderPath` for your computer, and create the output folder. The [demo setup](#demo-modification-localization) can prepare these paths automatically.
4. **If using CL**, first [download the ESM2 model](#download-the-esm2-model-for-cl) into `_internal/`. Ordinary XGBoost runs do not need this download.
5. **Double-click `pReLoc.exe`**, enter the full path to your edited parameter file, and press Enter. See [Run pReLoc](#run-preloc) for both double-click and command-line instructions.
6. **Open the configured output folder.** Start with `pReLocRes.tsv`; see [Read the results](#read-the-results) for the difference between scored, FLR-filtered, and directly assigned results.

## System requirements and tested versions

### Operating system and hardware

| Item | Requirement / verified configuration |
| --- | --- |
| Operating system | 64-bit Windows. This distribution was tested on **Windows 11 Pro 25H2, build 26200.9711** on 11 November 2026. |
| Reference desktop | **Intel Core i7-12700**, 12 cores / 20 logical processors, approximately **64 GB RAM**. This is the machine used for the demo measurements below, not a minimum specification. |
| CPU/GPU | The ordinary XGBoost and CL demos run on CPU. **No GPU or other non-standard hardware is required for these demos.** The supplied PyTorch build is CPU-only; setting `clDevice=cuda` does not add CUDA support to it. |
| Memory | The minimum RAM requirement has not been established. CL loads a 650-million-parameter ESM2 model and needs substantially more memory than XGBoost. Use the tested desktop configuration as the reproducibility reference; reduce `clBatchSize` if memory is insufficient. |
| Storage | Approximately **11.5 GiB** of repository files, plus the separately downloaded **2.43 GiB** ESM2 model for CL (about **13.9 GiB / 14.95 GB** total). Allow about **35 GB free space** for a Git/LFS checkout, its cache, the model, and demo outputs; larger analyses need additional space. |
| Optional pPred | The separate native pPred runtime includes CUDA 11/cuDNN 8 libraries. |

### Software dependencies and version numbers

The executable includes its Windows runtime: **no separate Python environment, `pip install`, compiler, or CUDA toolkit installation is needed for the ordinary XGBoost/CL demos**. XGBoost's models are included. **CL additionally requires the ESM2 model download described below.** Keep the entire `_internal/` directory. Git and Git LFS are needed for the repository download procedure below, rather than for executing an already downloaded copy.

Versions were read from the packaged runtime and its distribution metadata, not from an unrelated Python environment on the test machine:

| Bundled runtime / analysis dependency | Version |
| --- | --- |
| CPython | 3.13.11 (64-bit) |
| NumPy | 2.4.3 |
| pandas | 3.0.1 |
| SciPy | 1.17.1 |
| scikit-learn | 1.8.0 |
| XGBoost | 3.2.0 |
| PyTorch / torchvision | 2.11.0+cpu / 0.26.0+cpu |
| fair-esm | 2.0.0 |
| Transformers / tokenizers | 4.48.1 / 0.21.4 |
| huggingface-hub / safetensors | 0.36.2 / 0.7.0 |
| Matplotlib / Pillow | 3.10.8 / 12.3.0 |
| NetworkX / SymPy / mpmath | 3.6.1 / 1.14.0 / 1.3.0 |
| joblib / threadpoolctl | 1.5.3 / 3.6.0 |
| python-dateutil | 2.9.0.post0 |

Additional bundled distributions whose metadata is retained include:

| Dependency | Version | Dependency | Version |
| --- | --- | --- | --- |
| attrs | 26.1.0 | cryptography | 46.0.3 |
| filelock | 3.25.2 | Jinja2 | 3.1.6 |
| MarkupSafe | 3.0.3 | packaging | 25.0 |
| prompt-toolkit | 3.0.52 | pydantic | 2.12.4 |
| PyYAML | 6.0.3 | regex | 2026.4.4 |
| requests | 2.32.5 | rich | 14.2.0 |
| sentencepiece | 0.2.1 | tqdm | 4.67.1 |
| certifi | 2025.11.12 | urllib3 | 2.6.1 |
| charset-normalizer | 3.4.4 | lxml | 6.1.3 |

Native dependencies are also supplied with the package, including the Microsoft C/C++ runtime (root `msvcp140.dll` / `vcruntime140.dll`: **14.44.35208.0**), native candidate-generation executables, pFindSDK, and the separate pPred runtime. The pPred CUDA runtime DLL reports file version `6,14,11,11080`, and its cuDNN DLL reports `6,14,11,6050` (Windows file-version identifiers). Some native components do not expose a release version; use the exact bundled files rather than substituting system libraries. Their complete file/version inventory is recorded in [BUNDLED_DEPENDENCIES.tsv](BUNDLED_DEPENDENCIES.tsv), including explicit `not reported` entries where version metadata is absent. This inventory also records the packaged Python distributions; it is not a list of packages to install manually.

The demo validation covers this exact package, not every possible combination of the dependency versions above. The application does not declare a semantic release version. The tested `pReLoc.exe` SHA-256 is:

```text
708dc35c675a7b19138db955e0d16bafce7e35844d21621cf1f97cd5b8d4c176
```

## Installation guide

1. Install Git for Windows and [Git LFS](https://git-lfs.com/) if needed. The repository preparation used **Git 2.53.0.windows.2** and **Git LFS 3.7.1**; these are recorded tool versions, not minimum required versions.
2. In PowerShell, download the repository files:

   ```powershell
   git lfs install
   git clone https://github.com/pFindStudio/pReLoc.git
   cd pReLoc
   git lfs pull
   git lfs fsck
   ```

3. Keep `pReLoc.exe`, `_internal/`, and `examples/` together. Use a writable installation directory because the application and optional pPred tools create working/log files. No installer or compilation step follows the download.
4. Verify the executable starts:

   ```powershell
   .\pReLoc.exe --help
   ```

Large repository assets are stored with Git LFS. An executable or model that contains `version https://git-lfs.github.com/spec/v1` is only a pointer: complete `git lfs pull` before running. Do not assume a source ZIP contains the actual LFS objects. The ESM2 model exceeds this repository's per-file LFS limit and is intentionally excluded from Git; **`git lfs pull` does not download it**. Follow the additional step below before running CL.

### Download the ESM2 model for CL

Download the **ESM-2, 33-layer, 650-million-parameter UR50D** checkpoint from the [official Meta/FAIR ESM model list](https://github.com/facebookresearch/esm#available-models). Use this exact checkpoint; other ESM2 sizes and Hugging Face `.bin`/`.safetensors` exports are not drop-in replacements for this package.

| Item | Required value |
| --- | --- |
| File to download | **`esm2_t33_650M_UR50D.pt`** |
| Official direct download | [Download esm2_t33_650M_UR50D.pt](https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt) |
| Destination | **`<pReLoc installation folder>\_internal\esm2_t33_650M_UR50D.pt`** |
| Example destination | `C:\tools\pReLoc\_internal\esm2_t33_650M_UR50D.pt` |
| Exact file size | **2,604,537,549 bytes** (approximately 2.43 GiB) |
| SHA-256 of the checkpoint used in the demo tests | `ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c` |

Open PowerShell in the folder containing `pReLoc.exe` and download the file directly into `_internal/`:

```powershell
curl.exe --fail --location --retry 3 --output ".\_internal\esm2_t33_650M_UR50D.pt" "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt"
if ($LASTEXITCODE -ne 0) { throw 'ESM2 download failed; do not run CL yet.' }
$model = Get-Item -LiteralPath '.\_internal\esm2_t33_650M_UR50D.pt'
if ($model.Length -ne 2604537549) { throw 'ESM2 file size is incorrect.' }
$expected = 'ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c'
$actual = (Get-FileHash -LiteralPath $model.FullName -Algorithm SHA256).Hash
if ($actual -ne $expected) { throw 'ESM2 checksum mismatch; download the official checkpoint again.' }
```

If the direct download is interrupted or blocked by your network, download the same linked file in a browser or on another machine, copy it to the destination above, and run the size/checksum checks. Do not save an HTML error page or rename a different checkpoint to this filename.

Keep these three files together:

```text
pReLoc/
  pReLoc.exe
  _internal/
    best_model.pth                               (included)
    esm2_t33_650M_UR50D-contact-regression.pt      (included)
    esm2_t33_650M_UR50D.pt                        (download separately)
```

The small **`esm2_t33_650M_UR50D-contact-regression.pt`** companion is already included in Git LFS. The [ESM local-file loader](https://github.com/facebookresearch/esm/blob/main/esm/pretrained.py) expects it next to the main checkpoint. If it is missing, first run `git lfs pull`; its [official recovery download](https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt) must also be saved in `_internal/` under that exact filename. Keep the supplied pReLoc-specific `best_model.pth`; it is a different model and must not be replaced by ESM2 weights.

The templates already specify the following. A relative `clEsmPath` resolves inside `_internal/`, so **do not add another `_internal\` prefix** to this setting:

```ini
[score]
clModelPath=best_model.pth
clEsmPath=esm2_t33_650M_UR50D.pt
```

After this one-time download, CL loads the checkpoint locally and can run offline. The ordinary XGBoost workflow does not require this ESM2 checkpoint.

### Typical installation time

Allow approximately **5–10 minutes for local file preparation after the download** on a desktop comparable to the reference machine. Copying the complete 7,745-file package to a separate local directory took **5 minutes 8 seconds** in the validation setup. This is a measured local copy, not a measured fresh network installation or an installer benchmark.

Network transfer time is additional. For 14.95 GB including the separate ESM2 download, the ideal transfer time is about **20 minutes at a sustained 100 Mbit/s**, or **100 minutes at 20 Mbit/s**, before protocol overhead, GitHub/LFS throttling, antivirus scanning, and disk work. The ESM2 file alone takes about **3.5 minutes at 100 Mbit/s** under ideal conditions. Actual download speeds can differ substantially. Complete the ESM2 download before timing the CL demo; download time is not part of its reported runtime.

## Package contents

| Path | Purpose |
| --- | --- |
| `pReLoc.exe` | Combined XGBoost / CL application. |
| `_internal/` | Python runtime, included scoring models, native candidate generators and pPred dependencies; also the destination for the separately downloaded ESM2 weights. |
| `examples/data_input.tsv` | The 100-PSM custom identification example. |
| `examples/dataMGF.mgf` | Matching spectra for the basic examples. |
| `examples/pReLocparam_*.txt` | Editable workflow parameter files. |
| `examples/mod_*`, `examples/site_*` | Example output folders; the pPred folder also contains its input data and fine-tuned models. |
| `examples/pChem/pChem_preparedata/` | pChem summary, pFind identification table and spectra. |
| `USAGE.txt` | Short command reference. |

Use 64-bit Windows and keep the complete application directory, including `_internal/`, together. A separate Python installation is not required. This packaged CL runtime uses CPU PyTorch and reads its ESM2 weights locally. Allow several minutes for a complete CL example on CPU, plus disk space for analysis results.

## Run pReLoc

### Option 1: double-click the executable

You can run pReLoc without opening PowerShell:

1. Prepare and save your parameter file. Check its input paths, selected model, and output folder. Create the output folder if it does not exist.
2. In Windows File Explorer, **double-click `pReLoc.exe`**. A console window opens.
3. When asked for the parameter-file path, paste its **full path**, for example `C:\tools\pReLoc\examples\my_mod_xgb.txt`, and press **Enter**. The file must be your saved configuration; the example path is not created automatically.
4. Wait for the analysis to finish. Check any messages in the console, then open the folder specified by `OutputFolderPath` to view the results.
5. Press **Enter** at the final prompt to close the window.

The backend is selected by `[score] model=XGBoost` or `model=CL` in the parameter file. CL requires the separately downloaded ESM2 checkpoint before starting.

### Option 2: use the command line

Open PowerShell in the folder containing `pReLoc.exe`:

```powershell
.\pReLoc.exe "C:\analysis\my_parameters.txt"
```

You can override the model for one run or display help:

```powershell
.\pReLoc.exe "C:\analysis\my_parameters.txt" --model CL
.\pReLoc.exe "C:\analysis\my_parameters.txt" --model XGBoost
.\pReLoc.exe --help
```

`--model` takes priority over the parameter-file setting. Older parameter files without a `model` entry use XGBoost. Command-line runs with a parameter-file argument exit automatically; a failed run returns a nonzero exit code. Starting the executable without that argument uses the same interactive prompts as double-clicking it.

## Demo: modification localization

The small real dataset is `examples/data_input.tsv` (100 identifications) with matching `examples/dataMGF.mgf` spectra. The following setup prepares both demos and redirects their outputs away from the supplied reference results. Open PowerShell in this application's root directory and run:

```powershell
$root = (Get-Location).Path
$originalRoot = 'C:\Users\Jiaxiang\Desktop\pReLoc_all'
foreach ($backend in @('xgb', 'cl')) {
    $out = Join-Path $root "results\demo_$backend"
    New-Item -ItemType Directory -Force $out | Out-Null
    $text = Get-Content ".\examples\pReLocparam_mod_$backend.txt" -Raw
    $text = $text.Replace($originalRoot, $root)
    $text = [regex]::Replace($text, '(?m)^OutputFolderPath=[^\r\n]*', "OutputFolderPath=$out")
    $text = $text.Replace('clDevice=auto', 'clDevice=cpu')
    $config = Join-Path $root "examples\demo_mod_$backend.txt"
    [IO.File]::WriteAllText($config, $text, [Text.UTF8Encoding]::new($false))
}
```

Run the XGBoost demo:

```powershell
.\pReLoc.exe ".\examples\demo_mod_xgb.txt"
```

For CL, first complete **Download the ESM2 model for CL** above, then run:

```powershell
.\pReLoc.exe ".\examples\demo_mod_cl.txt"
```

The original templates contain absolute paths from the packaging machine. The preparation commands above replace that prefix with your checkout path and create new output directories. If you edit parameters manually, update every enabled input/model path and `OutputFolderPath` yourself. Re-running a configuration can replace its prior results.

To run either prepared demo by double-clicking, start `pReLoc.exe` and enter the full path to `examples/demo_mod_xgb.txt` or `examples/demo_mod_cl.txt` in your installation folder.

### Expected demo output and run time

Results are written to `results/demo_xgb/` and `results/demo_cl/`. Expect `pReLocRes.tsv`, `pReLocRes-filtered.tsv`, `DirectRes.tsv`, and `Time.txt`. Both example configurations perform modification localization with `isFLR=2`, `FLR=0.01`, and `ispPred=false`. The CL measurement uses `clDevice=cpu` and `clBatchSize=8`.

The table below reports a single validation run on the reference desktop described in **System requirements**, on 11 November 2026. Wall time includes executable startup; `Full Process` is the program's own timing field and excludes some startup overhead. These are practical expectations for this small dataset, not guaranteed performance for all computers or production datasets.

| Backend | End-to-end wall time | `Full Process` | Main result rows | Rows after 1% FLR filtering | `DirectRes.tsv` rows |
| --- | --- | --- | --- | --- | --- |
| XGBoost, CPU | 16.0 s | 4.5 s | 100 | 17 | 0 |
| CL, CPU, batch size 8 | 137.7 s (about 2 min 18 s) | 135.4 s | 100 | 13 | 0 |

For planning on a similar desktop, allow roughly **15–30 seconds for XGBoost** and **2–5 minutes for CL** for this demo. These ranges are guidance based on the measured run and the packaged example timings, rather than a multi-machine benchmark.

Successful command-line runs return exit code `0`. `Time.txt` must include `Full Process`, and the main result should contain the localization/FLR columns described below. `DirectRes.tsv` contains only its header in this demo: none of the input PSMs has just one allowed placement of the target modifications. Row counts exclude the header; the filtered result can differ between the XGBoost and CL models. Floating-point results may vary across machines.

The demonstration does not validate the optional pPred or pChem workflows. Their templates and usage are included below, but no runtime claim is made for those workflows.

## Instructions for use on your own data

1. Select the appropriate workflow template from the table below and make a working copy. Keep your input data and original parameter file as the reproducibility record.
2. Export corresponding spectra as MGF and prepare either a pFind identification table or the custom TSV format described in **Prepare your inputs**. Match the spectrum identifiers exactly.
3. Set `MGFPath`, `isMGFfolder`, `Input`, and `ispfind`; enable pChem/pPred only for a matching workflow. Replace every path used by that workflow with a path on your own machine.
4. Set the backend (`model=XGBoost` or `model=CL`), one scoring mode, modification definitions, target modifications, and an appropriate FLR setting. Create a separate writable `OutputFolderPath` for each analysis.
5. Run the configuration and check the exit code:

   ```powershell
   .\pReLoc.exe "C:\analysis\my_parameters.txt"
   $LASTEXITCODE
   ```

6. Check `Time.txt`, inspect the main and filtered result tables, and retain the parameter file with those results. Report the application revision, backend, modifications, and filtering settings with the analysis. The **Read the results** section defines the output fields.

## Choose a workflow

Use exactly one of `isSite=true` and `isMod=true`. Update the selected template's paths and create the configured output directory before running.

| Workflow | Template in `examples/` | Important settings |
| --- | --- | --- |
| Modification localization, XGBoost | `pReLocparam_mod_xgb.txt` | `model=XGBoost`, `isMod=true`, `ispPred=false`, `isFLR=2` |
| Modification localization, CL | `pReLocparam_mod_cl.txt` | `model=CL`, `isMod=true`, `ispPred=false`, `isFLR=2` |
| Site scoring, XGBoost | `pReLocparam_site_xgb.txt` | `model=XGBoost`, `isSite=true`, `isFLR=0` |
| Site scoring, CL | `pReLocparam_site_cl.txt` | `model=CL`, `isSite=true`, `ispPred=false`, `isFLR=0` |
| pChem localization, XGBoost | `pReLocparam_pChem.txt` | `model=XGBoost`, `ispChem=true`, `ispfind=true`, `isFLR=0` |
| pChem localization, CL | `pReLocparam_pChem_cl.txt` | `model=CL`, `ispChem=true`, `ispfind=true`, `isFLR=0` |
| Localization with pPred features | `pReLocparam_mod_xgbpPred.txt` | `model=XGBoost`, `ispPred=true`, `isFLR=1` |

CL does not support `ispPred=true`. Its settings are `clDevice`, `clBatchSize`, `clModelPath` and `clEsmPath`. The supplied examples use `clDevice=auto` and `clBatchSize=8`. Relative scoring-model paths resolve inside the application resources, so keep the default checkpoint names unless using compatible replacement weights.

## Prepare your inputs

### Spectra

Set `MGFPath` to an MGF file with `isMGFfolder=false`, or to a folder of MGF files with `isMGFfolder=true`. Spectrum identifiers must match the MGF `TITLE` values. The additional raw-data inputs required for pPred are configured separately.

### Custom identification table

Use `ispfind=false` and a tab-separated table following `examples/data_input.tsv`:

| Column | Meaning |
| --- | --- |
| `spec_name` | Spectrum identifier matching an MGF title. |
| `seq` | Peptide amino-acid sequence. |
| `mod` | Modification annotations, for example `3,Oxidation[M];6,Phospho[S];`. |
| `score` | Original identification score. |

The demonstration table includes additional metadata columns that can be retained. Residue positions are 1-based; position 0 can represent an N-terminal modification. Use modification names matching those in the configuration.

### pFind identification table

Set `ispfind=true` and point `Input` to a pFind result file, such as `pFind-Filtered.spectra`. Set `MGFPath` to its corresponding spectra. With `ispfindOpen=true`, modification sections are populated from the pFind identifications. Unknown modification masses may be requested and written back to the parameter file.

## Configure modifications and FLR

List residue modifications under `[mods]`, terminal modifications under `[NTerm_mods]` and `[CTerm_mods]`, and modifications to localize under `[scoreModName]`. For example:

```ini
[mods]
mod1=Carbamidomethyl[C]
mod2=Oxidation[M]
mod3=Phospho[S]
mod4=Phospho[T]
mod5=Phospho[Y]

[NTerm_mods]
Nmod1=Acetyl[ProteinN-term]

[CTerm_mods]

[scoreModName]
scoremod1=Phospho[S]
scoremod2=Phospho[T]
scoremod3=Phospho[Y]
```

The `[FLR]` section applies to modification localization:

| Setting | Meaning |
| --- | --- |
| `isFLR=0` | Disable FLR filtering. |
| `isFLR=1` | Generate decoys by relocating amino acids. Use this strategy with the supplied pPred workflow. |
| `isFLR=2` | Generate decoys by relocating modification mass. Used by the ordinary XGBoost / CL localization examples. |
| `FLR=0.01` | Filter at 1% FLR. The main result's `FLR` column is expressed in percent. |

Use `isFLR=0` for pChem and site scoring. With pPred, `isFLR=2` generates modifications such as `Phospho[A]` that the native predictor does not support; use `isFLR=1` or disable FLR.

## pChem workflow

Start with the appropriate pChem template. Supply the MGF folder, pFind table, pChem summary and number of top summary modifications to consider:

```ini
[path_info]
isMGFfolder=true
MGFPath=C:\analysis\pChem_preparedata
ispfind=true
Input=C:\analysis\pChem_preparedata\pFind-Filtered.spectra
ispChem=true
pChemmodnum=1
InputSummary=C:\analysis\pChem_preparedata\pChem.summary
```

Keep `ispPred=false` and `isFLR=0`, select XGBoost or CL, and configure an existing output folder. These lines are an excerpt; retain all other required template sections. Modification conversion is performed in memory and preserves the input `.spectra` file, so the same input can be reused for both backends.

## Optional pPred features

Enable pPred only with XGBoost. The main scoring input can be the custom table, while `pFindResfp` separately supplies pFind identifications for prediction.

| Setting in `[pPred]` | Purpose |
| --- | --- |
| `ispPred=true` | Enable prediction-derived features. |
| `pFindResfp` | Matching pFind identification file. |
| `rawfolder` | Parsed raw-data folder with matching `ms1`/`ms2` or `pf1`/`pf2` and index files. |
| `modelfolder` | Directory for fine-tuned RT/MS2 models. |
| `modelname` | Shared basename of `<modelname>_rt.pt` and `<modelname>_ms2.pt`. |
| `instrument` | Instrument code; retain the template's documented mapping. |
| `nce` | Normalized collision energy. |
| `rtparam` | Retention-time range as `minRT_maxRT`. |
| `ms2param` | `0_nce_nce_nce_minCharge_maxCharge`. |

The supplied example's raw-data files and fine-tuned models are in `examples/mod_xgbpPred/`. Use data-specific models for your own analyses. The pPred tools and native libraries are bundled under `_internal/pPred-MS2RT/`; their working files are kept there rather than in the main result folder. The package includes RT/MS2 base models under its `pPred_data/mainmodel/` directory for fine-tuning.

Prediction-derived features are calculated and combined in memory. `Candidates_pPred.tsv` is not created in the analysis output directory.

## Read the results

For modification localization, pReLoc separates records that need model scoring from records whose modification positions are already fixed by the configured rules. The output files have different meanings:

| File | Contents |
| --- | --- |
| `pReLocRes.tsv` | Model-scored results. In localization mode, these are PSMs with alternative allowed placements of the target modifications. |
| `pReLocRes-filtered.tsv` | Scored localization results that pass the configured FLR threshold, when filtering is enabled. |
| `DirectRes.tsv` | PSMs with only one allowed placement of the target modifications under the configured rules. These are written directly, without model scoring or FLR estimation. |
| `Time.txt` | Stage timings; a completed run includes `Full Process`. |

### What does `DirectRes.tsv` mean?

The program checks the allowed modification sites and the number of target modifications for each PSM, taking account of positions already occupied by other modifications. If these rules leave **only one possible assignment**, there are no alternative positions for the localization model to compare, so that PSM is written to `DirectRes.tsv`.

For example, suppose a peptide has **one phosphorylation**, the configuration allows phosphorylation on **S, T, or Y**, and the peptide has **only one available S/T/Y residue**. That residue is the only allowed location, so the PSM is reported directly. One possible assignment can also involve multiple modified residues.

`DirectRes.tsv` contains `spec_name`, `seq`, `ori_mod`, `mod`, and `scoremodsite_num`. It has **no model-derived confidence score or estimated FLR**. Keep these directly assigned records separate when reporting results that passed an FLR threshold; an absent FLR does not mean an estimated FLR of zero. If every input PSM has only one allowed assignment, the main result table contains only its header.

### Columns in the scored results

Site-scoring results contain `spec_name`, `seq`, `mod` and `pSite2_score`; multiple site scores are comma-separated. Localization results include:

| Columns | Meaning |
| --- | --- |
| `spec_name`, `seq` | Spectrum identifier and peptide sequence. |
| `scoremodsite_num` | Number of localized target modifications. |
| `ori_mod`, `mod` | Original and selected modification assignments. |
| `ori_pepScore`, `pepScore` | Original and selected peptide scores. |
| `ori_pReLocsiteScore`, `pReLocsiteScore` | Original and selected site scores. |
| `ori_pReLocScore`, `pReLocScore` | Original and selected localization scores. |
| `FLR` | Estimated FLR in percent, when enabled. |

When a decoy wins the competition, the retained best target has no assigned FLR and does not enter the filtered result. XGBoost and CL values are not calibrated to the same score scale.

## Troubleshooting

- **An input file cannot be found:** update all paths for the enabled workflow, using absolute paths. Create the output directory and check write access.
- **The program waits at startup or exit:** provide the parameter-file path at the prompt, or pass it as a command-line argument for unattended runs.
- **A CL checkpoint cannot be found:** download `esm2_t33_650M_UR50D.pt` into `_internal/` using the instructions above. Keep the included contact-regression companion and `best_model.pth`, check file sizes/hashes, and retain the default checkpoint names. `git lfs pull` alone does not fetch the excluded ESM2 checkpoint.
- **pPred reports unknown modifications:** use `isFLR=1` or `isFLR=0`, and confirm that the configured modifications are supported by its predictor.
- **pPred reports a DLL or model error:** retain its full bundled dependency directory and check the raw-data and model paths.
- **Only intermediate candidate features are present:** inspect the error and `Time.txt`; candidate generation alone does not establish successful scoring.
- **A `paramtmp.txt` left by an older version is present:** the updated package removes known old generated copies and cleans new copies at run completion. Keep your editable input parameters outside these generated-file locations.

For a reproducible issue report, include the workflow, parameter file with private paths redacted, console error and a small representative input.
