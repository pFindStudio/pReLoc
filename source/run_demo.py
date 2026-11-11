"""Run the supplied demonstration with absolute paths resolved on this computer."""
import argparse
from pathlib import Path

from pReLoc import main


def run():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=str.upper, choices=["XGBOOST", "CL"], default="XGBOOST")
    parser.add_argument("--output", type=Path, help="output directory (default: results/<model>)")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    output = (args.output or root / "results" / args.model.lower()).resolve()
    output.mkdir(parents=True, exist_ok=True)
    replacements = {
        "MGFPath": str(root / "examples/dataMGF.mgf"),
        "Input": str(root / "examples/data_input.tsv"),
        "OutputFolderPath": str(output), "model": args.model,
    }
    lines = []
    for line in (root / "param.txt").read_text(encoding="utf-8-sig").splitlines():
        key = line.partition("=")[0]
        lines.append(f"{key}={replacements[key]}" if key in replacements else line)
    parameter_file = output / "param.txt"
    parameter_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    main(str(parameter_file))
    print(f"Results: {output}")


if __name__ == "__main__":
    run()
