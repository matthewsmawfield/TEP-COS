#!/usr/bin/env python3
"""
SDSS SkyServer DR18 Schema Explorer

Dumps full table/column metadata from SDSS SkyServer for cosmological analysis.
Outputs:
- JSON schema dumps (data/sdss_schema/tables.json, columns.json)
- CSV reports (data/sdss_schema/*.csv)
- Markdown summary of high-value TEP probe columns
"""

import requests
import json
import csv
import os
import time
from datetime import datetime

SDSS_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data", "sdss_schema")
os.makedirs(OUTPUT_DIR, exist_ok=True)

TEP_KEYWORDS = [
    "time", "redshift", "velocity", "dispersion", "magnitude", "flux",
    "age", "metallicity", "distance", "luminosity", "mass", "shear",
    "peculiar", "environment", "density", "alpha", "iron", "mgb",
    "d4000", "hbeta", "hdelta", "lick", "sfr", "ssfr", "bpt"
]


def run_sdss_query(sql, max_retries=3, delay=1.0):
    """Execute SQL query against SDSS SkyServer using GET."""
    for attempt in range(max_retries):
        try:
            response = requests.get(
                SDSS_URL,
                params={"cmd": sql, "format": "json"},
                timeout=300
            )
            if response.status_code == 200:
                data = response.json()
                if data and len(data) > 0 and "Rows" in data[0]:
                    return data[0]["Rows"]
                return []
            else:
                print(f"  HTTP {response.status_code}")
        except requests.exceptions.Timeout:
            print(f"  Timeout (attempt {attempt + 1})")
        except Exception as e:
            print(f"  Error: {e}")
        if attempt < max_retries - 1:
            time.sleep(delay)
    return []


def discover_tables_via_dbobjects():
    """Try SDSS's DBObjects metadata table."""
    sql = """
    SELECT name, type, description
    FROM DBObjects
    WHERE type IN ('U', 'V')
    ORDER BY name
    """
    return run_sdss_query(sql)


def discover_tables_via_sysobjects():
    """Fallback: use sys.objects."""
    sql = """
    SELECT name, type_desc AS type
    FROM sys.objects
    WHERE type IN ('U', 'V')
    ORDER BY name
    """
    return run_sdss_query(sql)


def get_columns_via_dbcolumns(table_name):
    """Try SDSS's DBColumns metadata table."""
    sql = f"""
    SELECT field AS column_name, type AS data_type, description, unit
    FROM DBColumns
    WHERE tableName = '{table_name}'
    ORDER BY id
    """
    return run_sdss_query(sql)


def get_columns_via_information_schema(table_name):
    """Fallback: use INFORMATION_SCHEMA."""
    sql = f"""
    SELECT
        column_name,
        data_type,
        character_maximum_length AS char_len,
        numeric_precision,
        numeric_scale,
        is_nullable
    FROM INFORMATION_SCHEMA.COLUMNS
    WHERE table_name = '{table_name}'
    ORDER BY ordinal_position
    """
    return run_sdss_query(sql)


def identify_tep_columns(columns, table_name):
    """Flag columns potentially useful for TEP analysis."""
    tep_cols = []
    for col in columns:
        col_name = col.get("column_name", "").lower()
        col_desc = str(col.get("description", "")).lower()
        for kw in TEP_KEYWORDS:
            if kw in col_name or kw in col_desc:
                col["table_name"] = table_name
                col["tep_relevance"] = kw
                tep_cols.append(col)
                break
    return tep_cols


def main():
    print("=" * 70)
    print("SDSS SKYSERVER DR18 SCHEMA EXPLORER")
    print("=" * 70)
    print(f"Started: {datetime.now().isoformat()}")
    print(f"Output dir: {OUTPUT_DIR}")

    # Step 1: Discover tables
    print("\n[1/4] Discovering tables...")
    tables = discover_tables_via_dbobjects()
    method = "DBObjects"
    if not tables:
        print("  DBObjects not available, trying sys.objects...")
        tables = discover_tables_via_sysobjects()
        method = "sys.objects"
    
    if not tables:
        print("  ERROR: Could not retrieve table list from SDSS.")
        return
    
    print(f"  Found {len(tables)} tables/views via {method}")
    
    # Save tables
    tables_path = os.path.join(OUTPUT_DIR, "tables.json")
    with open(tables_path, "w") as f:
        json.dump(tables, f, indent=2)
    print(f"  Saved: {tables_path}")

    # Step 2: Get columns for each table
    print("\n[2/4] Discovering columns (this may take a few minutes)...")
    all_columns = {}
    tep_columns = []
    use_dbcolumns = True
    
    # Test which column method works
    test_cols = get_columns_via_dbcolumns("SpecObj")
    if not test_cols:
        print("  DBColumns not available, using INFORMATION_SCHEMA...")
        use_dbcolumns = False
    
    for i, table in enumerate(tables):
        table_name = table["name"]
        if i % 50 == 0:
            print(f"  Processing table {i+1}/{len(tables)}: {table_name}...")
        
        if use_dbcolumns:
            cols = get_columns_via_dbcolumns(table_name)
        else:
            cols = get_columns_via_information_schema(table_name)
        
        if cols:
            all_columns[table_name] = cols
            tep_cols = identify_tep_columns(cols, table_name)
            tep_columns.extend(tep_cols)
        
        # Rate limit
        time.sleep(0.1)
    
    # Save all columns
    columns_path = os.path.join(OUTPUT_DIR, "columns.json")
    with open(columns_path, "w") as f:
        json.dump(all_columns, f, indent=2)
    print(f"  Saved: {columns_path}")
    
    # Save TEP-relevant columns
    tep_path = os.path.join(OUTPUT_DIR, "tep_columns.json")
    with open(tep_path, "w") as f:
        json.dump(tep_columns, f, indent=2)
    print(f"  Found {len(tep_columns)} TEP-relevant columns")
    print(f"  Saved: {tep_path}")

    # Step 3: Generate CSV summary
    print("\n[3/4] Generating CSV reports...")
    
    # Tables CSV
    tables_csv = os.path.join(OUTPUT_DIR, "tables.csv")
    with open(tables_csv, "w", newline="") as f:
        if tables:
            writer = csv.DictWriter(f, fieldnames=tables[0].keys())
            writer.writeheader()
            writer.writerows(tables)
    print(f"  Saved: {tables_csv}")
    
    # TEP columns CSV
    tep_csv = os.path.join(OUTPUT_DIR, "tep_columns.csv")
    with open(tep_csv, "w", newline="") as f:
        fieldnames = ["table_name", "column_name", "data_type", "description", "tep_relevance"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(tep_columns)
    print(f"  Saved: {tep_csv}")

    # Step 4: Generate markdown report
    print("\n[4/4] Generating markdown report...")
    report_path = os.path.join(OUTPUT_DIR, "tep_columns_report.md")
    with open(report_path, "w") as f:
        f.write("# SDSS SkyServer Columns Relevant for TEP Cosmology Tests\n\n")
        f.write(f"Generated: {datetime.now().isoformat()}\n\n")
        f.write(f"Total tables: {len(tables)}\n")
        f.write(f"Total TEP-relevant columns: {len(tep_columns)}\n\n")
        
        # Group by table
        by_table = {}
        for col in tep_columns:
            tbl = col.get("table_name", "unknown")
            if tbl not in by_table:
                by_table[tbl] = []
            by_table[tbl].append(col)
        
        f.write("## Tables with TEP-Relevant Columns\n\n")
        for tbl in sorted(by_table.keys()):
            cols = by_table[tbl]
            f.write(f"### `{tbl}` ({len(cols)} columns)\n\n")
            f.write("| Column | Type | TEP Relevance | Description |\n")
            f.write("|--------|------|---------------|-------------|\n")
            for col in cols:
                desc = str(col.get("description", ""))[:60].replace("|", "/")
                f.write(f"| `{col.get('column_name', '')}` | {col.get('data_type', '')} | {col.get('tep_relevance', '')} | {desc} |\n")
            f.write("\n")
    
    print(f"  Saved: {report_path}")

    # Summary
    print("\n" + "=" * 70)
    print("SCHEMA EXPLORATION COMPLETE")
    print("=" * 70)
    print(f"Tables discovered: {len(tables)}")
    print(f"Tables with columns: {len(all_columns)}")
    print(f"TEP-relevant columns: {len(tep_columns)}")
    print(f"\nOutput files:")
    print(f"  - {tables_path}")
    print(f"  - {columns_path}")
    print(f"  - {tep_path}")
    print(f"  - {tables_csv}")
    print(f"  - {tep_csv}")
    print(f"  - {report_path}")


if __name__ == "__main__":
    main()
