#!/usr/bin/env python3
"""PDF processing for TEP-COS.

Compresses a PDF (Ghostscript) and embeds metadata (ExifTool; Ghostscript pdfmark fallback).

Usage:
    python3 scripts/utils/process_pdf.py <input_pdf> [--quality ebook|printer|prepress|default]

Notes:
- Requires `gs` (Ghostscript) for compression.
- Uses `exiftool` if available for robust metadata; otherwise falls back to Ghostscript pdfmark.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional


def _load_json(path: Path) -> Optional[Dict]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _read_text(path: Path) -> Optional[str]:
    try:
        return path.read_text(encoding="utf-8")
    except Exception:
        return None


def _load_results_data(project_root: Path) -> Dict[str, object]:
    """Load key analysis results from output JSON files.
    
    Returns a dictionary with key findings for PDF metadata embedding.
    """
    results_dir = project_root / "results" / "outputs"
    data: Dict[str, object] = {}
    
    # Load hierarchical density scaling results
    density_results = _load_json(results_dir / "step_5_33_hierarchical_density_results.json")
    if density_results:
        data["density_slope"] = density_results.get("model_b_mixed_slope")
        data["density_slope_error"] = density_results.get("model_b_mixed_error")
        data["density_p_value"] = density_results.get("model_b_mixed_p")
        data["newtonian_predicted"] = density_results.get("newtonian_predicted_slope")
        data["rejection_sigma"] = density_results.get("rejection_sigma")
        data["suppression_factor"] = density_results.get("suppression_factor")
    
    # Load binary control results
    binary_results = _load_json(results_dir / "step_5_36_integrated_binary_control.json")
    if binary_results:
        gc = binary_results.get("gc_binary_analysis", {})
        field = binary_results.get("field_binary_control", {})
        data["gc_binary_n"] = gc.get("n_binary")
        data["gc_isolated_n"] = gc.get("n_isolated")
        data["binary_inversion_dex"] = gc.get("diff_dex")
        data["binary_inversion_p"] = gc.get("mw_p")
        data["field_binary_n"] = field.get("n_binary")
        data["field_isolated_n"] = field.get("n_isolated")
    
    # Load sample size results
    sample_results = _load_json(results_dir / "step_5_27_hybrid_maximum_analysis.json")
    if sample_results:
        sizes = sample_results.get("sample_sizes", {})
        data["total_gc_msps"] = sizes.get("gc_total")
        data["total_field_msps"] = sizes.get("field_total")
        data["total_msps"] = sizes.get("total")
    
    # Load CMC Gold Standard results
    cmc_results = _load_json(results_dir / "step_5_50_cmc_gold_standard.json")
    if cmc_results:
        data["cmc_n_pulsars"] = cmc_results.get("n_cmc_pulsars")
        data["cmc_verdict"] = cmc_results.get("verdict", {}).get("overall_verdict")
        tests = cmc_results.get("tests", {})
        raw_excess = tests.get("raw_excess", {})
        data["cmc_raw_excess_sigma"] = raw_excess.get("sigma")
        density_scaling = tests.get("density_scaling", {})
        data["cmc_density_sigma"] = density_scaling.get("sigma_difference")
    
    return data


def _parse_citation_cff(cff_text: str) -> Dict[str, object]:
    """Minimal, dependency-free parser for the subset of CFF we use.

    This is intentionally conservative: it only extracts a few top-level scalar fields plus
    the keywords list and the first author name.
    """

    data: Dict[str, object] = {}

    def set_scalar(key: str, value: str) -> None:
        value = value.strip().strip('"').strip("'")
        if value:
            data[key] = value

    lines = cff_text.splitlines()
    i = 0
    in_abstract = False
    abstract_lines: List[str] = []
    in_keywords = False
    keywords: List[str] = []

    first_author_family = None
    first_author_given = None

    while i < len(lines):
        line = lines[i]

        if in_abstract:
            if not line.startswith("  ") and line.strip() and ":" in line:
                in_abstract = False
                data["abstract"] = "\n".join(abstract_lines).strip()
                abstract_lines = []
                continue
            abstract_lines.append(line[2:] if line.startswith("  ") else line)
            i += 1
            continue

        if in_keywords:
            stripped = line.strip()
            if stripped.startswith("-"):
                kw = stripped[1:].strip().strip('"').strip("'")
                if kw:
                    keywords.append(kw)
                i += 1
                continue
            in_keywords = False
            continue

        stripped = line.strip()

        if stripped.startswith("abstract:") and stripped.endswith(">"):
            in_abstract = True
            i += 1
            continue

        if stripped.startswith("keywords:"):
            in_keywords = True
            i += 1
            continue

        if stripped.startswith("title:"):
            set_scalar("title", stripped.split(":", 1)[1])
            i += 1
            continue

        if stripped.startswith("doi:"):
            set_scalar("doi", stripped.split(":", 1)[1])
            i += 1
            continue

        if stripped.startswith("date-released:"):
            set_scalar("date_released", stripped.split(":", 1)[1])
            i += 1
            continue

        if stripped.startswith("version:"):
            set_scalar("version", stripped.split(":", 1)[1])
            i += 1
            continue

        if stripped.startswith("url:"):
            set_scalar("url", stripped.split(":", 1)[1])
            i += 1
            continue

        if stripped.startswith("repository-code:"):
            set_scalar("repository_code", stripped.split(":", 1)[1])
            i += 1
            continue

        if stripped.startswith("license:"):
            set_scalar("license", stripped.split(":", 1)[1])
            i += 1
            continue

        if stripped.startswith("authors:"):
            # Scan forward for first author given/family.
            j = i + 1
            while j < len(lines):
                s = lines[j].strip()
                # Authors are typically YAML list entries (e.g., "- family-names:")
                s_key = s.lstrip("- ").strip()
                if s.startswith("preferred-citation:") or (s and not lines[j].startswith(" ") and ":" in s):
                    break
                if s_key.startswith("family-names:") and first_author_family is None:
                    first_author_family = s_key.split(":", 1)[1].strip().strip('"').strip("'")
                if s_key.startswith("given-names:") and first_author_given is None:
                    first_author_given = s_key.split(":", 1)[1].strip().strip('"').strip("'")
                if first_author_family and first_author_given:
                    break
                j += 1
            i += 1
            continue

        i += 1

    if in_abstract and abstract_lines:
        data["abstract"] = "\n".join(abstract_lines).strip()

    if keywords:
        data["keywords"] = keywords

    if first_author_family or first_author_given:
        data["author"] = " ".join([p for p in [first_author_given, first_author_family] if p])

    return data


def _load_default_metadata(project_root: Path) -> Dict[str, str]:
    version_json = project_root / "VERSION.json"
    # Check site/ subdirectory first (TEP-COS structure), then project root
    citation_cff = project_root / "site" / "CITATION.cff"
    if not citation_cff.exists():
        citation_cff = project_root / "CITATION.cff"

    version = {}
    try:
        version = json.loads(version_json.read_text(encoding="utf-8"))
    except Exception:
        version = {}

    cff_text = _read_text(citation_cff) or ""
    cff = _parse_citation_cff(cff_text) if cff_text else {}

    title = str(cff.get("title") or "The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars")
    author = str(cff.get("author") or "")
    doi = str(cff.get("doi") or "")
    date_released = str(cff.get("date_released") or "")
    codename = str(version.get("codename") or "")
    vnum = str(version.get("version") or cff.get("version") or "")

    abstract = str(cff.get("abstract") or "")

    keywords_list = cff.get("keywords") if isinstance(cff.get("keywords"), list) else []
    keywords = "; ".join([str(k) for k in keywords_list])
    if codename and vnum:
        keywords = f"{keywords}; {codename} v{vnum}" if keywords else f"{codename} v{vnum}"

    repo = str(cff.get("repository_code") or cff.get("url") or "")

    # ExifTool expects PDF date format like YYYY:MM:DD HH:MM:SS
    creation_date = ""
    if date_released:
        try:
            yyyy, mm, dd = date_released.split("-")
            creation_date = f"{yyyy}:{mm}:{dd} 00:00:00"
        except Exception:
            creation_date = ""

    subject_parts = []
    if abstract:
        subject_parts.append(" ".join(abstract.split()))
    if doi:
        subject_parts.append(f"DOI: {doi}")
    if repo:
        subject_parts.append(f"Code: {repo}")
    subject = " ".join(subject_parts).strip()

    license_str = str(cff.get("license") or "CC-BY-4.0")

    producer = "TEP-COS Research Project"
    if codename and vnum:
        producer = f"TEP-COS Research Project ({codename} v{vnum})"

    metadata: Dict[str, str] = {
        "Title": title,
        "Author": author,
        "Creator": author,
        "Producer": producer,
        "Subject": subject,
        "Keywords": keywords,
        "Copyright": f"Creative Commons Attribution 4.0 International License (CC BY 4.0)" if "CC-BY" in license_str.upper() else license_str,
    }

    if creation_date:
        metadata["CreationDate"] = creation_date
        metadata["ModifyDate"] = creation_date

    # Optional-but-useful fields (ExifTool recognizes XMP*: and some PDF keys too, but keep conservative)
    if doi:
        metadata["Identifier"] = doi
    
    # Load and embed key results data
    results_data = _load_results_data(project_root)
    if results_data:
        # Add sample sizes
        if results_data.get("total_msps"):
            metadata["XMP-TEPCOS:TotalMSPs"] = str(results_data["total_msps"])
        if results_data.get("total_gc_msps"):
            metadata["XMP-TEPCOS:GCMSPs"] = str(results_data["total_gc_msps"])
        if results_data.get("total_field_msps"):
            metadata["XMP-TEPCOS:FieldMSPs"] = str(results_data["total_field_msps"])
        
        # Add density scaling results
        if results_data.get("density_slope"):
            metadata["XMP-TEPCOS:DensitySlope"] = f"{results_data['density_slope']:.3f}"
        if results_data.get("rejection_sigma"):
            metadata["XMP-TEPCOS:RejectionSigma"] = f"{results_data['rejection_sigma']:.2f}"
        if results_data.get("suppression_factor"):
            metadata["XMP-TEPCOS:SuppressionFactor"] = f"{results_data['suppression_factor']:.3f}"
        
        # Add binary inversion results
        if results_data.get("binary_inversion_dex"):
            metadata["XMP-TEPCOS:BinaryInversionDex"] = f"{results_data['binary_inversion_dex']:.3f}"
        if results_data.get("binary_inversion_p"):
            metadata["XMP-TEPCOS:BinaryInversionP"] = f"{results_data['binary_inversion_p']:.4f}"
        
        # Add CMC Gold Standard results
        if results_data.get("cmc_n_pulsars"):
            metadata["XMP-TEPCOS:CMCNPulsars"] = str(results_data["cmc_n_pulsars"])
        if results_data.get("cmc_verdict"):
            metadata["XMP-TEPCOS:CMCVerdict"] = str(results_data["cmc_verdict"])
        if results_data.get("cmc_raw_excess_sigma"):
            metadata["XMP-TEPCOS:CMCRawExcessSigma"] = f"{results_data['cmc_raw_excess_sigma']:.1f}"
        if results_data.get("cmc_density_sigma"):
            metadata["XMP-TEPCOS:CMCDensitySigma"] = f"{results_data['cmc_density_sigma']:.2f}"
    
    return metadata


def compress_pdf(input_path: str, output_path: str, quality: str = "ebook") -> Dict[str, float]:
    quality_settings = {
        "screen": "/screen",
        "ebook": "/ebook",
        "printer": "/printer",
        "prepress": "/prepress",
        "default": "/default",
    }

    if quality not in quality_settings:
        raise ValueError(f"Quality must be one of: {', '.join(quality_settings.keys())}")

    original_size = os.path.getsize(input_path)

    cmd = [
        "gs",
        "-sDEVICE=pdfwrite",
        "-dCompatibilityLevel=1.4",
        f"-dPDFSETTINGS={quality_settings[quality]}",
        "-dNOPAUSE",
        "-dQUIET",
        "-dBATCH",
        f"-sOutputFile={output_path}",
        input_path,
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
        compressed_size = os.path.getsize(output_path)
        reduction = ((original_size - compressed_size) / original_size) * 100 if original_size else 0.0
        return {
            "original_mb": original_size / (1024 * 1024),
            "compressed_mb": compressed_size / (1024 * 1024),
            "reduction_pct": reduction,
        }
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Ghostscript compression failed: {e.stderr.decode(errors='replace')}")


def embed_metadata(pdf_path: str, metadata: Dict[str, str]) -> None:
    cmd = ["exiftool"]
    for key, value in metadata.items():
        if value is None:
            continue
        v = str(value).strip()
        if not v:
            continue
        cmd.extend([f"-{key}={v}"])

    cmd.extend(["-overwrite_original", pdf_path])

    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("  ⚠ exiftool failed or not found, falling back to Ghostscript pdfmark...")
        embed_metadata_gs(pdf_path, metadata)


def embed_metadata_gs(pdf_path: str, metadata: Dict[str, str]) -> None:
    def escape_ps(s: str) -> str:
        return s.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".ps", delete=False, encoding="utf-8") as f:
        meta_str = ""
        for key, value in metadata.items():
            v = str(value).strip()
            if not v:
                continue
            meta_str += f"/{key} ({escape_ps(v)}) "
        f.write(f"[ {meta_str} /DOCINFO pdfmark")
        pdfmark_path = f.name

    output_path = f"{pdf_path}.tmp"
    cmd = [
        "gs",
        "-sDEVICE=pdfwrite",
        "-dNOPAUSE",
        "-dBATCH",
        "-dQUIET",
        f"-sOutputFile={output_path}",
        pdf_path,
        pdfmark_path,
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
        os.replace(output_path, pdf_path)
    finally:
        try:
            if os.path.exists(pdfmark_path):
                os.unlink(pdfmark_path)
        except Exception:
            pass


def verify_metadata(pdf_path: str, fields: List[str]) -> Optional[str]:
    cmd = ["exiftool"] + [f"-{f}" for f in fields] + [pdf_path]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return result.stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Compress and embed metadata into a TEP-COS PDF")
    parser.add_argument("input_pdf", help="Path to input PDF")
    parser.add_argument(
        "--quality",
        choices=["screen", "ebook", "printer", "prepress", "default"],
        default="ebook",
        help="Ghostscript compression quality (default: ebook)",
    )
    parser.add_argument("--doi", default=None, help="Override DOI metadata")
    parser.add_argument("--title", default=None, help="Override Title metadata")
    parser.add_argument("--author", default=None, help="Override Author metadata")
    parser.add_argument("--url", default=None, help="Override URL (added into Subject)")

    args = parser.parse_args()

    input_path = Path(args.input_pdf).resolve()
    if not input_path.exists():
        print(f"Error: file not found: {input_path}")
        return 1

    project_root = Path(__file__).resolve().parents[2]
    metadata = _load_default_metadata(project_root)
    
    # Load results data for display
    results_data = _load_results_data(project_root)

    if args.title:
        metadata["Title"] = args.title
    if args.author:
        metadata["Author"] = args.author
        metadata["Creator"] = args.author
    if args.doi:
        metadata["Identifier"] = args.doi
        if metadata.get("Subject"):
            metadata["Subject"] = f"{metadata['Subject']} DOI: {args.doi}".strip()
        else:
            metadata["Subject"] = f"DOI: {args.doi}".strip()

    if args.url:
        if metadata.get("Subject"):
            metadata["Subject"] = f"{metadata['Subject']} URL: {args.url}".strip()
        else:
            metadata["Subject"] = f"URL: {args.url}".strip()

    print(f"Processing TEP-COS PDF: {input_path}")
    print(f"Quality: {args.quality}")
    print()
    
    # Display key results data if available
    if results_data:
        print("Key Results Data:")
        if results_data.get("total_msps"):
            print(f"  Sample: {results_data.get('total_gc_msps', 'N/A')} GC + {results_data.get('total_field_msps', 'N/A')} Field = {results_data['total_msps']} MSPs")
        if results_data.get("rejection_sigma"):
            print(f"  Density Scaling: {results_data['rejection_sigma']:.2f}σ rejection of Newtonian prediction")
        if results_data.get("binary_inversion_dex"):
            print(f"  Binary Inversion: {results_data['binary_inversion_dex']:.3f} dex (p={results_data.get('binary_inversion_p', 'N/A'):.4f})")
        if results_data.get("cmc_n_pulsars"):
            print(f"  CMC Analysis: {results_data['cmc_n_pulsars']:,} synthetic pulsars")
        if results_data.get("cmc_raw_excess_sigma"):
            print(f"  CMC Discrepancy: {results_data['cmc_raw_excess_sigma']:.1f}σ (raw excess)")
        if results_data.get("cmc_verdict"):
            print(f"  CMC Verdict: {results_data['cmc_verdict']}")
        print()

    print("Step 1: Compressing PDF...")
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        stats = compress_pdf(str(input_path), tmp_path, args.quality)
        os.replace(tmp_path, str(input_path))
        print(f"  Original:   {stats['original_mb']:.2f} MB")
        print(f"  Compressed: {stats['compressed_mb']:.2f} MB")
        print(f"  Reduction:  {stats['reduction_pct']:.1f}%")
        print()
    except Exception as e:
        try:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)
        except Exception:
            pass
        print(f"Error during compression: {e}")
        return 1

    print("Step 2: Embedding metadata...")
    embed_metadata(str(input_path), metadata)
    print("  Metadata embedded")
    print()

    print("Step 3: Verifying metadata (if exiftool available)...")
    verification = verify_metadata(str(input_path), ["Title", "Author", "Subject", "Keywords", "Creator", "Producer", "Copyright"])
    if verification:
        print("  ✓ Metadata verified")
        print()
        print(verification)
    else:
        print("  ⚠ Verification skipped (exiftool not available).")

    print()
    print(f"✓ Done: {input_path}")
    print(f"  Final size: {os.path.getsize(input_path) / (1024 * 1024):.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
