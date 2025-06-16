# pip install genebe pandas openpyxl
import os
import pandas as pd
import genebe as gnb

# 1) Either set credentials once…
# gnb.set_credentials(username="tgarg2@jhu.edu", api_key="ak-xo8hVORNe1bHu3oqei8PKtxiI")

# …or skip set_credentials and pass them on each call by setting use_netrc=False below.

def annotate_excel_file(input_path, output_path):
    # Read your sheet
    df = pd.read_excel(input_path, sheet_name=0, engine='openpyxl')
    
    # Build the minimal variant DataFrame
    df_vars = pd.DataFrame({
        'chr': df['CHROM'].astype(str).str.replace('chr', '', case=False),
        'pos': df['POS'].astype(int),
        'ref': df['REF'].astype(str),
        'alt': df['ALT'].astype(str),
    })
    
    # Annotate, requesting a DataFrame back
    annotated_df = gnb.annotate(
        df_vars,
        genome='hg19',                   # Match your build
        use_ensembl=False,
        use_refseq=True,
        flatten_consequences=True,
        output_format='dataframe',       # <-- ensure a DataFrame is returned
        use_netrc=False,                 # <-- disable netrc so username/api_key are honored
        username="tgarg2@jhu.edu",       # optional if you didn’t call set_credentials
        api_key="ak-xo8hVORNe1bHu3oqei8PKtxiI",
        batch_size=500
    )
    
    # Side-by-side concat (same row order) instead of merge
    merged = pd.concat([df.reset_index(drop=True), annotated_df.reset_index(drop=True)], axis=1)
    
    # Write it out
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    merged.to_excel(output_path, index=False, engine='openpyxl')

def main():
    input_folder = '.'
    output_folder = 'gnb_anno'
    for fn in os.listdir(input_folder):
        if fn.lower().endswith(('.xlsx', '.xls')):
            in_f  = os.path.join(input_folder, fn)
            out_f = os.path.join(output_folder, fn)
            annotate_excel_file(in_f, out_f)
            print(f"Annotated {fn} → {out_f}")

if __name__ == '__main__':
    main()
