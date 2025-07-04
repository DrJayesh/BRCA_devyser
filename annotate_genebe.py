# pip install genebe pandas openpyxl
import os
import sys
import pandas as pd
import genebe as gnb
from datetime import date

# 1) Either set credentials once…
# gnb.set_credentials(username="tgarg2@jhu.edu", api_key="ak-xo8hVORNe1bHu3oqei8PKtxiI")

# …or skip set_credentials and pass them on each call by setting use_netrc=False below.

def _build_variant_df(df: pd.DataFrame) -> pd.DataFrame:
    """Return a DataFrame with columns chr,pos,ref,alt from a parsed Excel file."""
    return pd.DataFrame({
        "chr": df["CHROM"].astype(str).str.replace("chr", "", case=False, regex=False),
        # Treat ``pos`` as a string to avoid dtype mismatches when merging with
        # stored annotations that may have been saved with a different type.
        "pos": df["POS"].astype(str),
        "ref": df["REF"].astype(str),
        "alt": df["ALT"].astype(str),
    })

# ROOT_ANNOTATION_FILE stays in script directory
ROOT_ANNOTATION_FILE = os.path.join(os.path.dirname(__file__), "annotated_variants.xlsx")


def _load_local_annotations() -> pd.DataFrame:
    """Return the local annotation DataFrame if present, else an empty DataFrame."""
    if os.path.exists(ROOT_ANNOTATION_FILE):
        df = pd.read_excel(ROOT_ANNOTATION_FILE, engine="openpyxl")
        # Ensure merge columns use consistent dtypes
        for col in ["chr", "pos", "ref", "alt"]:
            if col in df.columns:
                df[col] = df[col].astype(str)
        # Drop any ``chr`` prefix so the column matches the format produced by
        # ``_build_variant_df``.
        if "chr" in df.columns:
            df["chr"] = df["chr"].str.replace("chr", "", case=False, regex=False)
        return df
    return pd.DataFrame()


def _save_local_annotations(df: pd.DataFrame) -> None:
    """Write ``df`` to the persistent annotation Excel file."""
    df.to_excel(ROOT_ANNOTATION_FILE, index=False, engine="openpyxl")


def annotate_folder(input_folder: str) -> None:
    """Annotate all Excel files in ``input_folder`` using a single set of unique
    variants and write the annotated files sorted by AF into a ``gnb_anno``
    directory in the parent of ``input_folder``."""

    excel_files = [f for f in os.listdir(input_folder) if f.lower().endswith((".xlsx", ".xls"))]
    if not excel_files:
        print("No Excel files found for annotation.")
        return

    dfs = {}
    all_var_dfs = []

    # Collect variant rows from each file
    for fn in excel_files:
        path = os.path.join(input_folder, fn)
        df = pd.read_excel(path, sheet_name=0, engine="openpyxl")
        var_df = _build_variant_df(df)
        dfs[fn] = (df, var_df)
        all_var_dfs.append(var_df)

    unique_vars = pd.concat(all_var_dfs, ignore_index=True).drop_duplicates().reset_index(drop=True)
    local_ann = _load_local_annotations()

    merge_cols = ["chr", "pos", "ref", "alt"]
    annotated_unique = pd.DataFrame()
    to_annotate = unique_vars.copy()

    if not local_ann.empty:
        # Merge to bring in existing annotations
        annotated_unique = unique_vars.merge(local_ann, on=merge_cols, how="left")
        # Identify only those rows where all annotation columns are missing
        ann_cols = local_ann.columns.difference(merge_cols)
        to_annotate = annotated_unique[ann_cols].isna().all(axis=1)
        to_annotate = annotated_unique.loc[to_annotate, merge_cols]
        # Keep the previously annotated rows
        annotated_unique = annotated_unique.dropna(subset=ann_cols, how="all")

    if not to_annotate.empty:
        new_ann = gnb.annotate(
            to_annotate,
            genome="hg19",
            use_ensembl=False,
            use_refseq=True,
            flatten_consequences=True,
            output_format="dataframe",
            use_netrc=False,
            username="tgarg2@jhu.edu",
            api_key="ak-xo8hVORNe1bHu3oqei8PKtxiI",
            batch_size=500,
        )
        # Normalize dtypes of merge columns returned by GeneBe
        for col in merge_cols:
            new_ann[col] = new_ann[col].astype(str)
        new_ann["chr"] = new_ann["chr"].str.replace("chr", "", case=False, regex=False)
        new_ann["FIRST_ANNOTATED"] = date.today().isoformat()

        # Concatenate new annotations with existing, drop duplicates
        local_ann = pd.concat([local_ann, new_ann], ignore_index=True)
        local_ann.drop_duplicates(subset=merge_cols, keep="first", inplace=True)
        _save_local_annotations(local_ann)

        # Update annotated_unique by concatenation
        annotated_unique = pd.concat([annotated_unique, new_ann], ignore_index=True)
        annotated_unique.drop_duplicates(subset=merge_cols, inplace=True)
    elif annotated_unique.empty:
        # No local_ann and no new to annotate: just keep unique_vars
        annotated_unique = unique_vars.copy()

    # Prepare output folder
    output_folder = os.path.join(os.path.dirname(os.path.abspath(input_folder)), "gnb_anno")
    os.makedirs(output_folder, exist_ok=True)

    for fn in excel_files:
        df, var_df = dfs[fn]
        # Merge back annotations
        annotated_df = var_df.merge(annotated_unique, on=merge_cols, how="left")
        merged = pd.concat(
            [df.reset_index(drop=True), annotated_df.drop(columns=merge_cols).reset_index(drop=True)],
            axis=1,
        )
        # Sort by AF (leave ambiguous column as-is)
        merged.sort_values(by="AF", ascending=False, inplace=True)
        out_path = os.path.join(output_folder, fn)
        merged.to_excel(out_path, index=False, engine="openpyxl")
        print(f"Annotated {fn} → {out_path}")


def main():
    folder = os.getcwd()
    if len(sys.argv) > 1:
        folder = sys.argv[1]
    annotate_folder(folder)


if __name__ == '__main__':
    main()
