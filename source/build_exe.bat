@echo off
setlocal
pushd "%~dp0" || exit /b 1
if not defined PYTHON set "PYTHON=python"
set "OUTPUT_DIR=%~dp0dist"
if not "%~1"=="" set "OUTPUT_DIR=%~f1"
echo Building pReLoc with XGBoost and CL support...
"%PYTHON%" -c "import PyInstaller, numpy, pandas, scipy, sklearn, joblib, xgboost, torch, esm, tqdm; from pathlib import Path; files=['XGBoost_final.pkl','XGBoost_pPred.pkl','XGBoost.pkl','best_model.pth','esm2_t33_650M_UR50D.pt','esm2_t33_650M_UR50D-contact-regression.pt','local_denovo/loc_denovo.exe','local_denovo/loc_denovo_pho.exe','pPred-MS2RT/pPred.exe','pPred-MS2RT/GenPredData.exe']; missing=[f for f in files if not Path(f).is_file()]; assert not missing, 'Missing build assets: '+repr(missing); print('Build dependencies and XGB/CL/native assets: OK')"
if errorlevel 1 goto :failed

REM ESM exceeds 2 GB. Keep model weights and native tools in an onedir bundle.
"%PYTHON%" -m PyInstaller --noconfirm --clean --distpath "%OUTPUT_DIR%" --workpath "%~dp0build" "%~dp0pReLoc.spec"
if errorlevel 1 goto :failed
"%OUTPUT_DIR%\pReLoc\pReLoc.exe" --help
if errorlevel 1 goto :failed

echo Build complete: %OUTPUT_DIR%\pReLoc\pReLoc.exe
echo Keep the entire pReLoc folder, including _internal, together.
echo Select model=XGBoost or model=CL in the parameter file [score] section.
popd
exit /b 0

:failed
echo Build FAILED. See the error above.
popd
exit /b 1
