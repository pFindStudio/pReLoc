# -*- mode: python ; coding: utf-8 -*-
from pathlib import Path
from PyInstaller.utils.hooks import collect_all

root = Path(SPECPATH)
datas = []
binaries = []
hiddenimports = ["cl", "cl_scoring", "torch", "esm", "xgboost.sklearn", "sklearn"]
for package in ("xgboost", "esm"):
    package_data, package_binaries, package_imports = collect_all(package)
    datas += package_data
    binaries += package_binaries
    hiddenimports += package_imports
for filename in (
    "XGBoost_final.pkl", "XGBoost_pPred.pkl", "XGBoost.pkl", "best_model.pth",
    "esm2_t33_650M_UR50D.pt", "esm2_t33_650M_UR50D-contact-regression.pt",
):
    if not (root / filename).is_file():
        raise FileNotFoundError(f"Missing model asset: {root / filename}")
    datas.append((str(root / filename), "."))
for folder in ("local_denovo", "pPred-MS2RT"):
    for path in (root / folder).rglob("*"):
        relative = path.relative_to(root)
        if (path.is_file()
                and path.suffix.lower() not in {".pdb", ".exp", ".lib", ".zip"}
                and not any(part in {".vs", "__pycache__"} for part in relative.parts)
                and not (folder == "local_denovo" and "source" in relative.parts)
                and ("pPred_data" not in relative.parts or "mainmodel" in relative.parts)):
            datas.append((str(path), str(relative.parent)))
a = Analysis(
    [str(root / "pReLoc.py")], pathex=[str(root)], binaries=binaries, datas=datas,
    hiddenimports=hiddenimports, hookspath=[], hooksconfig={}, runtime_hooks=[],
    excludes=[], noarchive=False, optimize=0,
)
pyz = PYZ(a.pure)
# Keep the 2.6 GB ESM checkpoint outside the executable archive.
exe = EXE(
    pyz, a.scripts, [], exclude_binaries=True, name="pReLoc", debug=False,
    bootloader_ignore_signals=False, strip=False, upx=False, console=True,
)
coll = COLLECT(exe, a.binaries, a.datas, strip=False, upx=False, name="pReLoc")
