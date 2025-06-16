# pip install genebe pandas openpyxl
import os
import sys
import pandas as pd
import genebe as gnb

# 1) Either set credentials once…
# gnb.set_credentials(username="tgarg2@jhu.edu", api_key="ak-xo8hVORNe1bHu3oqei8PKtxiI")

# …or skip set_credentials and pass them on each call by setting use_netrc=False below.

def _build_variant_df(df: pd.DataFrame) -> pd.DataFrame:
    """Return a DataFrame with columns chr,pos,ref,alt from a parsed Excel file."""
    return pd.DataFrame({
        "chr": df["CHROM"].astype(str).str.replace("chr", "", case=False),
        "pos": df["POS"].astype(int),
        "ref": df["REF"].astype(str),
        "alt": df["ALT"].astype(str),
    })


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

    # Annotate the unique variants via GeneBe
    annotated_unique = gnb.annotate(
        unique_vars,
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

    output_folder = os.path.join(os.path.dirname(os.path.abspath(input_folder)), "gnb_anno")
    os.makedirs(output_folder, exist_ok=True)

    for fn in excel_files:
        df, var_df = dfs[fn]
        annotated_df = var_df.merge(annotated_unique, on=["chr", "pos", "ref", "alt"], how="left")
        merged = pd.concat(
            [df.reset_index(drop=True), annotated_df.drop(columns=["chr", "pos", "ref", "alt"]).reset_index(drop=True)],
            axis=1,
        )
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
