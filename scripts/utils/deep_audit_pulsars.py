import os
import ast
import re
import json
import glob

def find_python_files(directory):
    return glob.glob(os.path.join(directory, "**/*.py"), recursive=True)

def audit_file(filepath):
    issues = []
    warnings = []
    
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 1. Check for manual fabrications (e.g. lists of floats that look like data, but not np.array or from file)
    # Simple heuristic: looking for arrays with many numbers that might be fake data.
    if re.search(r'\[\s*\d+\.\d+\s*,\s*\d+\.\d+\s*,\s*\d+\.\d+', content):
        warnings.append("Possible hardcoded data array found (check for fabrication).")
        
    # 2. Check for data cherry-picking
    if '.drop(' in content or '.dropna(' in content or '.query(' in content or '!= ' in content:
        # Check if it drops specific clusters
        if re.search(r"\!=\s*['\"](?:Terzan|47 Tuc|NGC|M)", content):
            issues.append("Data cherry-picking detected: excluding specific clusters.")
            
    # 3. Check for statistical rigor (hardcoded p-values, t-stats)
    if re.search(r'p(?:_value)?\s*=\s*0\.\d+', content):
        if not re.search(r'p(?:_value)?\s*=\s*(?:scipy|stats|ttest|mannwhitney)', content):
            warnings.append("Possible hardcoded p-value found.")
            
    if re.search(r't(?:_stat|statistic)?\s*=\s*-?\d+\.\d+', content):
        warnings.append("Possible hardcoded t-statistic found.")
        
    # 4. Consistency of magic numbers
    # rho_intra = 0.3 is the standard from manuscript, let's see if it's altered
    rho_matches = re.findall(r'rho(?:_intra)?\s*=\s*(\d+\.\d+)', content)
    for match in rho_matches:
        if match != '0.3' and match != '0.30':
            issues.append(f"Inconsistent rho_intra value: {match} (expected 0.3)")
            
    # 5. Fixed random seeds (Reproducibility)
    if 'random' in content or 'sample' in content or 'bootstrap' in content.lower():
        if 'np.random.seed' not in content and 'random.seed' not in content:
            warnings.append("Random process might missing a fixed seed.")
            
    return issues, warnings

def run_audit():
    script_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'scripts', 'steps')
    files = find_python_files(script_dir)
    
    report = {
        "files_scanned": len(files),
        "issues_found": [],
        "warnings_found": []
    }
    
    for file in files:
        issues, warnings = audit_file(file)
        if issues:
            report['issues_found'].append({"file": os.path.basename(file), "issues": issues})
        if warnings:
            report['warnings_found'].append({"file": os.path.basename(file), "warnings": warnings})
            
    print(json.dumps(report, indent=2))
    
    output_path = os.path.join(os.path.dirname(script_dir), 'results', 'outputs', 'deep_audit_codebase_report.json')
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)

if __name__ == '__main__':
    run_audit()
