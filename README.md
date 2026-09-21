# pReLoc

A framework for site-level confidence assessment in proteomics, with workflows for amino-acid site scoring and modification-site localization.

This repository distributes Windows executables, example inputs and results, and the optional pPred prediction package. Some configuration comments, messages, and output columns retain the earlier name **pSite2**.

## Download

Install Git and [Git LFS](https://git-lfs.com/), then run the following in PowerShell:

```powershell
git lfs install
git clone https://github.com/pFindStudio/pReLoc.git
cd pReLoc
git lfs pull
```

The executables, model weights, and large spectrum files are stored with Git LFS. A file containing `version https://git-lfs.github.com/spec/v1` is a pointer, not the actual executable or data; run `git lfs pull` to retrieve its contents. Access to the repository is required if it is private.

The full distribution contains approximately **9.5 GiB** of files. Allow additional disk space for the Git LFS cache, executable extraction, and analysis outputs.

## Package contents

| Path | Purpose |
| --- | --- |
| `pReLoc_XGB.exe` | XGBoost scoring executable; supports the optional pPred workflow. |
| `pReLoc_CL.exe` | CL scoring executable; use with `ispPred=false`. |
| `examples/data_input.tsv` | Small custom identification table used by the basic examples. |
| `examples/dataMGF.mgf` | Matching spectra for the basic examples. |
| `examples/pReLocparam_*.txt` | Editable parameter templates. |
| `examples/mod_*`, `examples/site_*` | Example outputs and, for the pPred example, supporting data/models. |
| `examples/pChem/pChem_preparedata/` | pChem summary, pFind results, and MGF input. |
| `pPred-20241111/` | Optional prediction executables, DLLs, model weights, and working data. |

Use Windows for the supplied `.exe` programs. No source-build or Linux/macOS execution workflow is provided in this distribution. Keep the pPred directory and its bundled dependencies together when using pPred.

## Quick start: modification localization with XGBoost

This example uses the small bundled dataset, without pPred.

1. Open PowerShell in the repository root. Copy the template and create an output directory:

   ```powershell
   Copy-Item .\examples\pReLocparam_mod_xgb.txt .\examples\my_mod_xgb.txt
   New-Item -ItemType Directory -Force .\results\mod_xgb
   notepad .\examples\my_mod_xgb.txt
   ```

2. Update the paths in the copied file. For example, if the repository is at `C:\tools\pReLoc`, use the following values in their existing sections:

   ```ini
   [path_info]
   isMGFfolder=false
   MGFPath=C:\tools\pReLoc\examples\dataMGF.mgf
   ispfind=false
   Input=C:\tools\pReLoc\examples\data_input.tsv
   ispChem=false

   [score]
   isSite=false
   isMod=true

   [pPred]
   ispPred=false

   [FLR]
   isFLR=2
   FLR=0.01

   [output_info]
   OutputFolderPath=C:\tools\pReLoc\results\mod_xgb
   ```

   This is an excerpt, not a complete parameter file. Retain the other template settings, including the modification sections. Replace `C:\tools\pReLoc` with your actual installation path. Use absolute paths, and create the output directory before running.

3. Start the executable:

   ```powershell
   .\pReLoc_XGB.exe
   ```

4. When the program displays `Input the path of the parameter file:`, enter the full path to your edited file, for example:

   ```text
   C:\tools\pReLoc\examples\my_mod_xgb.txt
   ```

   The supplied executables read the parameter path interactively. Passing it as a command-line argument does not replace this prompt. At completion, press Enter when prompted to exit.

5. Inspect `pSiteRes.tsv`, `pSiteRes-filtered.tsv`, and `Time.txt` in your output directory. Check the console for errors as well as checking that the result files were produced.

The XGBoost modification workflow above was checked with the bundled 100-row input and `isFLR=2`, `FLR=0.01`: it produced 100 result rows and 11 filtered rows. These counts describe this example, not a general expected result for other datasets. Floating-point scores can vary across environments.

## Choose a workflow

Copy the corresponding template, update its paths, create a separate output directory, and supply the edited configuration at the executable's prompt.

| Workflow | Executable | Template in `examples/` | Key settings |
| --- | --- | --- | --- |
| Modification localization, XGBoost | `pReLoc_XGB.exe` | `pReLocparam_mod_xgb.txt` | `isMod=true`, `isSite=false`, `ispPred=false` |
| Modification localization, CL | `pReLoc_CL.exe` | `pReLocparam_mod_cl.txt` | `isMod=true`, `isSite=false`, `ispPred=false` |
| Site scoring, XGBoost | `pReLoc_XGB.exe` | `pReLocparam_site_xgb.txt` | `isSite=true`, `isMod=false` |
| Site scoring, CL | `pReLoc_CL.exe` | `pReLocparam_site_cl.txt` | `isSite=true`, `isMod=false` |
| pChem localization, XGBoost | `pReLoc_XGB.exe` | `pReLocparam_pChem.txt` | `ispChem=true`, `ispfind=true`, `isFLR=0` |
| pChem localization, CL | `pReLoc_CL.exe` | `pReLocparam_pChem_cl.txt` | `ispChem=true`, `ispfind=true`, `isFLR=0` |
| Localization with pPred features | `pReLoc_XGB.exe` | `pReLocparam_mod_xgbpPred.txt` | `ispPred=true`; additional pFind/raw-data inputs required |

Enable one scoring mode at a time. Existing example output folders are reference results; use a new output folder for your own runs.

## Prepare your inputs

### Spectra

Set `MGFPath` to an MGF file with `isMGFfolder=false`, or to a directory containing MGF files with `isMGFfolder=true`. The scoring input uses MGF spectra; the additional pPred raw-data inputs are configured separately.

Spectrum identifiers in the identification table must match the MGF `TITLE` values. For example, identifiers in the bundled data look like:

```text
20150615_HepG2_Phos_control_70ug_1.10493.10493.2.3.dta
```

### Custom identification table

Use `ispfind=false` and set `Input` to a tab-delimited file. Follow [examples/data_input.tsv](examples/data_input.tsv). The fields read by the custom-input workflow are:

| Column | Meaning |
| --- | --- |
| `spec_name` | Spectrum identifier matching an MGF title. |
| `seq` | Peptide amino-acid sequence. |
| `mod` | Modification annotations such as `3,Oxidation[M];6,Phospho[S];`. |
| `score` | Original identification score. |

The bundled table also includes `raw_name` and `phospho_count`; retain its structure when adapting the example. Residue positions in the examples are 1-based. A protein N-terminal acetylation is represented as `0,Acetyl[ProteinN-term];`. Match modification names exactly to the configured names.

For example, peptide `YGMGTSVER` has the modification string `3,Oxidation[M];6,Phospho[S];` in the example input.

### pFind results

Set `ispfind=true` and point `Input` to the pFind result file, such as `pFind-Filtered.spectra`. Provide its corresponding MGF spectra through `MGFPath`.

## Configure modifications and FLR filtering

List residue modifications under `[mods]`, terminal modifications under `[NTerm_mods]` and `[CTerm_mods]`, and modifications to localize under `[scoreModName]`. For ordinary modification localization, the target modifications must also appear in `[mods]`.

For the phosphorylation example:

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

Some shipped templates repeat the `mod3` key for S/T/Y. The example above uses distinct keys and was used in the quick-start check.

The `[FLR]` section controls false localization rate filtering for modification localization (`isMod=true`):

| Parameter | Meaning |
| --- | --- |
| `isFLR=0` | Disable FLR filtering. |
| `isFLR=1` | Use the amino-acid relocation decoy strategy described in the template. |
| `isFLR=2` | Use the modification-mass relocation decoy strategy recommended by the template. |
| `FLR=0.01` | Set the filtering threshold to 0.01 (1%). |

The supplied pChem workflow does not support `isFLR=1` or `2`; keep `isFLR=0`. The site-scoring templates also use `isFLR=0`.

## pChem workflow

Start from a pChem template and update these paths:

```ini
[path_info]
isMGFfolder=true
MGFPath=C:\tools\pReLoc\examples\pChem\pChem_preparedata
ispfind=true
Input=C:\tools\pReLoc\examples\pChem\pChem_preparedata\pFind-Filtered.spectra
ispChem=true
pChemmodnum=1
InputSummary=C:\tools\pReLoc\examples\pChem\pChem_preparedata\pChem.summary
```

Keep `ispPred=false` and `isFLR=0`, set your output directory, and run the matching executable. The example uses `pChemmodnum=1`; its summary's first modification is `PFIND_DELTA_334`. Use the summary and identification files from the same pChem analysis when working with your own data.

## Optional pPred features

pPred features are supported by the XGBoost executable. The template states that this workflow requires pFind results and does not support pChem. The main scoring table can still be the custom table used by the shipped example: `pFindResfp` separately supplies the pFind results needed by pPred.

Use `pReLocparam_mod_xgbpPred.txt` and update all of the following:

| Setting in `[pPred]` | Purpose |
| --- | --- |
| `ispPred=true` | Enable prediction-derived features. |
| `pFindResfp` | pFind result file used by pPred. |
| `rawfolder` | Matching parsed MS data directory (`ms1`/`ms2` or `pf1`/`pf2`, with accompanying index files where present). |
| `modelfolder` | Directory for the fine-tuned model files. |
| `modelname` | Model basename; the supplied example uses `demofinetune`. |
| `instrument` | `0`: QE; `1`: Astral; `2`: QE HF; `3`: Lumos; `4`: Exploris; `5`: Thermo-others; `6`: TIMS. |
| `nce` | Normalized collision energy; example: `30`. |
| `rtparam` | Retention-time range in minutes, `minRT_maxRT`; example: `10_130`. |
| `ms2param` | Format documented by the template: `0_nce_nce_nce_minCharge_maxCharge`; example: `0_30_30_30_2_4`. |

The example data and fine-tuned models are under `examples/mod_xgbpPred/`. Run from the repository root so the adjacent `pPred-20241111` directory is available. Ensure the model and output directories exist and are writable. Use model names/directories appropriate to your dataset; the bundled fine-tuned models belong to the supplied example.

The pPred package includes CUDA 11/cuDNN 8 libraries. A supported GPU/driver matrix and a CPU-only pPred procedure are not specified in the repository. The basic XGBoost quick-start check did not exercise pPred. Also, the checked-in `examples/mod_xgbpPred/pSiteRes.tsv` contains candidate-feature columns, and its `Time.txt` records candidate generation only; treat these as intermediate example artifacts, not evidence of a completed pPred run.

## Read the results

| Output | Contents |
| --- | --- |
| `pSiteRes.tsv` | Main scoring/localization result; schema depends on the workflow and completion stage. |
| `pSiteRes-filtered.tsv` | FLR-filtered localization results when filtering is enabled and completes. |
| `Time.txt` | Timing records for processing stages. A completed basic run includes `Full Process`. |

Completed **site-scoring** examples contain `spec_name`, `seq`, `mod`, and `pSite2_score`. The last field is a comma-separated list of site scores in peptide sequence order.

Completed **modification-localization** examples contain the following columns:

| Column(s) | Description |
| --- | --- |
| `spec_name`, `seq` | Spectrum identifier and peptide sequence. |
| `scoremodsite_num` | Number of target modification sites. |
| `ori_mod`, `mod` | Original and selected modification assignments. |
| `ori_pepScore`, `pepScore` | Peptide scores for the original and selected assignments. |
| `ori_pSite2siteScore`, `pSite2siteScore` | Site scores for the original and selected assignments. |
| `ori_pSite2Score`, `pSite2Score` | Localization scores for the original and selected assignments. |

Preserve these original column names when parsing results. An individual score is not the same field as the dataset-level `FLR` filtering threshold.

## Troubleshooting

- **The parameter file or data cannot be found:** replace the original author's absolute paths in the template with your own. Check `Input`, `MGFPath`, `OutputFolderPath`, and all paths for enabled optional workflows.
- **`AutoCustomInput.tsv` cannot be created:** create the configured output directory first and confirm that it is writable.
- **The executable is only a few hundred bytes or contains LFS pointer text:** run `git lfs pull` from the clone and resolve any download/authentication error before launching it.
- **The process waits at startup or exit:** enter the parameter-file path at the first prompt and press Enter at the final prompt. A `--help` argument does not bypass interactive input in this build.
- **Results are missing or still contain candidate-feature columns:** inspect the console error and `Time.txt`; candidate generation alone does not mean scoring completed. The program may catch an error and still reach its final exit prompt, so an exit code alone is insufficient.
- **pPred fails to load dependencies or run predictions:** check the complete `pPred-20241111` package, matching pFind/raw-data inputs, model paths, and the reported CUDA/DLL error. Use `ispPred=false` when following the basic non-pPred example.

For reproducible reports, include the executable used, your parameter file with private paths redacted if necessary, the console error, and a small representative input when opening an [issue](https://github.com/pFindStudio/pReLoc/issues).
