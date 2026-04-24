#!/usr/bin/env python3
"""
Build 10-TEP-COS-v0.5-Caracas.md from site/components HTML files.
Converts HTML to Markdown and concatenates in order.
"""

from pathlib import Path
import re

def html_to_md(html_content):
    """Simple HTML to Markdown conversion."""
    md = html_content
    
    # Remove HTML comments
    md = re.sub(r'<!--.*?-->', '', md, flags=re.DOTALL)
    
    # Convert headers
    md = re.sub(r'<h1>(.*?)</h1>', r'# \1', md)
    md = re.sub(r'<h2>(.*?)</h2>', r'## \1', md)
    md = re.sub(r'<h3>(.*?)</h3>', r'### \1', md)
    md = re.sub(r'<h4>(.*?)</h4>', r'#### \1', md)
    
    # Convert bold and italic
    md = re.sub(r'<strong>(.*?)</strong>', r'**\1**', md)
    md = re.sub(r'<em>(.*?)</em>', r'*\1*', md)
    md = re.sub(r'<b>(.*?)</b>', r'**\1**', md)
    md = re.sub(r'<i>(.*?)</i>', r'*\1*', md)
    
    # Convert code
    md = re.sub(r'<code>(.*?)</code>', r'`\1`', md)
    
    # Convert paragraphs
    md = re.sub(r'<p>(.*?)</p>', r'\1\n\n', md, flags=re.DOTALL)
    
    # Convert lists
    md = re.sub(r'<ul>(.*?)</ul>', r'\1', md, flags=re.DOTALL)
    md = re.sub(r'<ol>(.*?)</ol>', r'\1', md, flags=re.DOTALL)
    md = re.sub(r'<li>(.*?)</li>', r'- \1\n', md, flags=re.DOTALL)
    
    # Convert tables (simplified)
    md = re.sub(r'<table[^>]*>(.*?)</table>', r'\1', md, flags=re.DOTALL)
    md = re.sub(r'<thead>(.*?)</thead>', r'\1', md, flags=re.DOTALL)
    md = re.sub(r'<tbody>(.*?)</tbody>', r'\1', md, flags=re.DOTALL)
    md = re.sub(r'<tr>(.*?)</tr>', r'\1\n', md, flags=re.DOTALL)
    md = re.sub(r'<th>(.*?)</th>', r'| \1 ', md)
    md = re.sub(r'<td>(.*?)</td>', r'| \1 ', md)
    
    # Convert divs with classes to markdown sections
    md = re.sub(r'<div[^>]*>(.*?)</div>', r'\1', md, flags=re.DOTALL)
    
    # Clean up empty lines
    md = re.sub(r'\n{3,}', '\n\n', md)
    
    # Clean up HTML tags we didn't handle
    md = re.sub(r'<[^>]+>', '', md)
    
    return md.strip()

def build_manuscript():
    """Build manuscript from components."""
    repo_root = Path(__file__).resolve().parents[1]
    components_dir = repo_root / "site" / "components"
    output_file = repo_root / "10-TEP-COS-v0.5-Caracas.md"
    
    # Component files in order
    component_files = [
        "0_abstract.html",
        "1_introduction.html",
        "2_theory.html",
        "3_pulsars.html",
        "4_discussion.html",
        "5_conclusions.html",
        "6_references.html",
        "7_appendix.html",
        "8_reproducibility.html",
    ]
    
    manuscript_parts = []
    
    for component_file in component_files:
        component_path = components_dir / component_file
        if component_path.exists():
            with open(component_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            md_content = html_to_md(html_content)
            manuscript_parts.append(md_content)
            manuscript_parts.append('\n\n---\n\n')  # Section separator
    
    # Add Jakarta-style header
    from datetime import datetime
    header = """# The Temporal Equivalence Principle: Suppressed Density Scaling in Globular Cluster Pulsars
**Matthew Lukin Smawfield**
Version: v0.5 (Caracas)
First published: 9 January 2026 · Last updated: {date}
DOI: 10.5281/zenodo.18165798

---

""".format(date=datetime.now().strftime('%Y-%m-%d'))
    
    # Write manuscript with header
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(header)
        f.write('\n'.join(manuscript_parts))
    
    print(f"Manuscript built: {output_file}")
    return output_file

if __name__ == "__main__":
    build_manuscript()
