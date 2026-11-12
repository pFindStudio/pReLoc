import pandas as pd
from sklearn.model_selection import train_test_split


def split_by_peptide(data, random_state=42):
    """Group by pep sequence, keeping all spectra, sites, and modifications in one split."""
    if "pep" not in data.columns:
        raise ValueError("The required peptide sequence column 'pep' is missing")
    if data["pep"].isna().any() or data["pep"].str.strip().eq("").any():
        raise ValueError("The 'pep' column contains missing or empty peptide sequences")

    peptides = data["pep"].unique()

    # Split unique peptides into initial training and test sets at approximately 10:1.
    train_peptides, test_peptides = train_test_split(
        peptides, test_size=1 / 11, random_state=random_state
    )
    # Reserve approximately 1/10 of the initial training peptides for validation.
    train_peptides, val_peptides = train_test_split(
        train_peptides, test_size=0.1, random_state=random_state
    )

    train_data = data.loc[data["pep"].isin(train_peptides)].copy()
    val_data = data.loc[data["pep"].isin(val_peptides)].copy()
    test_data = data.loc[data["pep"].isin(test_peptides)].copy()
    return train_data, val_data, test_data


def balance_train_data(train_data, random_state=42):
    """Downsample training rows so label=1 and label=0 have equal counts."""
    if "label" not in train_data.columns:
        raise ValueError("The required label column 'label' is missing")
    if not train_data["label"].isin([0, 1]).all():
        raise ValueError("The 'label' column must contain only 0 and 1, with no missing values")

    positive_samples = train_data.loc[train_data["label"].eq(1)]
    negative_samples = train_data.loc[train_data["label"].eq(0)]
    sample_count = min(len(positive_samples), len(negative_samples))
    if sample_count == 0:
        raise ValueError("The training set must contain both positive and negative samples for balancing")

    # Keep all minority samples and downsample the majority without replacement.
    # Discarded samples are not moved to other splits.
    balanced_data = pd.concat(
        [
            positive_samples.sample(n=sample_count, random_state=random_state),
            negative_samples.sample(n=sample_count, random_state=random_state),
        ]
    )
    return balanced_data.sample(frac=1, random_state=random_state).copy()


def main():
    data = pd.read_csv("pLocate_Phosites.tsv", sep="\t")
    train_data, val_data, test_data = split_by_peptide(data)
    # Balance only the training split after grouping peptides.
    # Preserve the original validation and test distributions.
    train_data = balance_train_data(train_data)

    train_data.to_csv("train_pLocate_Phosites.tsv", sep="\t", index=False)
    val_data.to_csv("val_pLocate_Phosites.tsv", sep="\t", index=False)
    test_data.to_csv("test_pLocate_Phosites.tsv", sep="\t", index=False)

    for name, subset in (("Train", train_data), ("Validation", val_data), ("Test", test_data)):
        print(
            f"{name}: peptides={subset['pep'].nunique()}, "
            f"spectra={subset['spec_name'].nunique()}, samples={len(subset)}, "
            f"positive={subset['label'].eq(1).sum()}, negative={subset['label'].eq(0).sum()}"
        )


if __name__ == "__main__":
    main()
