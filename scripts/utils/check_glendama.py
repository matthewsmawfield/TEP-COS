#!/usr/bin/env python3
import requests
import re
import os
from pathlib import Path

def download_glendama_q0957():
    """
    Attempt to scrape Q0957+561 data links from GLENDAMA database.
    Since the main page is a form, we need to inspect how it submits.
    Usually GLENDAMA links are static files if we can find the directory.
    
    Target URL: https://grupos.unican.es/glendama/database/
    """
    print("Checking GLENDAMA structure...")
    
    # Try common directory patterns for data
    base_urls = [
        "https://grupos.unican.es/glendama/database/objects/QSO_B0957+561/",
        "https://grupos.unican.es/glendama/database/data/",
        "https://grupos.unican.es/glendama/database/lightcurves/"
    ]
    
    for url in base_urls:
        try:
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                print(f"Found accessible directory: {url}")
                # Look for .dat, .txt, .rdb files
                files = re.findall(r'href="([^"]+\.(?:dat|txt|rdb|lc))"', r.text)
                if files:
                    print(f"Found files: {files}")
                    return True
        except Exception as e:
            print(f"Failed {url}: {e}")
            
    print("Could not find direct data directory. Interactive form submission required.")
    return False

if __name__ == "__main__":
    download_glendama_q0957()
