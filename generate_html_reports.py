#!/usr/bin/env python3
"""
Script: generate_html_reports.py
--------------------------------
Prompts the user for an Excel (.xlsx) file and a Jinja2-based HTML template,
then generates one HTML report per row using the 'Sample_ID' column for filenames.
For each variant (from 'HGVSG'), it fetches the GeneBe summary and makes it available
in the template as {{ Summary }}.

Requirements:
  pip install pandas jinja2 openpyxl requests beautifulsoup4

Usage:
  python generate_html_reports.py
"""

import os
import pandas as pd
import requests
from bs4 import BeautifulSoup
from jinja2 import Template

BASE_URL = "https://genebe.net/variant/hg19/"

def fetch_variant_summary(variant_id: str) -> str:
    """
    Given a variant identifier like 'chr17-41246537-TT-T',
    fetch the GeneBe page and return the detailed summary text.
    """
    url = BASE_URL + variant_id
    resp = requests.get(url)
    resp.raise_for_status()  # abort on HTTP errors

    soup = BeautifulSoup(resp.text, "html.parser")
    heading = soup.find("h2", string="Variant summary")
    if heading is None:
        raise RuntimeError(f"Couldn't find the summary heading at {url}")

    paras = heading.find_next_siblings("p", limit=2)
    if len(paras) < 2:
        raise RuntimeError(f"Couldn't locate summary paragraphs at {url}")

    # Return the detailed description paragraph
    return paras[1].get_text(strip=True)


def main():
    # Prompt user for Excel file path
    excel_path = input("Enter path to the Excel (.xlsx) file: ").strip()
    if not os.path.isfile(excel_path):
        print(f"Error: File '{excel_path}' not found.")
        return

    # Load data
    df = pd.read_excel(excel_path, engine='openpyxl')

    # Prompt user for template file path (default: template.html)
    template_path = input("Enter path to the HTML template [template.html]: ").strip() or 'template.html'
    if not os.path.isfile(template_path):
        print(f"Error: Template '{template_path}' not found.")
        return

    # Read and compile template
    with open(template_path, 'r', encoding='utf-8') as f:
        template = Template(f.read())

    # Create output directory
    out_dir = 'reports'
    os.makedirs(out_dir, exist_ok=True)

    # Iterate rows and generate reports
    for idx, row in df.iterrows():
        # Build context from row
        context = row.to_dict()

        # Fetch and add the GeneBe summary
        hgvs_id = row.get('HGVSG')
        if pd.isna(hgvs_id):
            summary = "(no HGVSG provided)"
        else:
            try:
                summary = fetch_variant_summary(str(hgvs_id))
            except Exception as e:
                summary = f"(error fetching summary: {e})"
        context['Summary'] = summary

        # Render HTML
        html_content = template.render(**context)

        # Build filename from 'Sample_ID'
        if 'Sample_ID' not in row:
            print("Error: 'Sample_ID' column not found in the data.")
            return
        safe_id = str(row['Sample_ID']).replace(' ', '_')
        filename = f"{safe_id}.html"
        out_path = os.path.join(out_dir, filename)

        # Write the HTML report
        with open(out_path, 'w', encoding='utf-8') as out_file:
            out_file.write(html_content)

        print(f"Generated report: {out_path}")

if __name__ == '__main__':
    main()
