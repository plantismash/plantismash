#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011,2012 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Script to download and install Pfam and ClusterBlast databases (latest from Zenodo)."""

import argparse
import ctypes
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.error
import urllib.parse
import urllib.request
from os import path

# ------------------ PFAM (static) ------------------
PFAM_URL = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz"
PFAM_CHECKSUM_SHA256 = "b29bc2c54db8090531df0361a781b8d7397f60ebedc0c36a16e7d45e999cc329"

# ------------------ ClusterBlast (latest via Zenodo concept DOI) ------------------
CLUSTERBLAST_CONCEPT_RECID = "16927684"
CLUSTERBLAST_FILENAME = "clusterblast.tar.gz"

if sys.platform in ("win32", "darwin"):
    os.environ["EXEC"] = os.getcwd() + "\\exec"
    os.environ["PATH"] = os.pathsep + os.environ["EXEC"] + os.pathsep + os.environ["PATH"]


# ------------------ Helpers ------------------
def execute(commands, input=None):
    stdin_redir = subprocess.PIPE if input is not None else None
    proc = subprocess.Popen(commands, stdin=stdin_redir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate(input=input)
    return out, err, proc.returncode


def get_free_space(folder):
    if platform.system() == "Windows":
        free_bytes = ctypes.c_ulonglong(0)
        ctypes.windll.kernel32.GetDiskFreeSpaceExW(
            ctypes.c_wchar_p(folder), None, None, ctypes.pointer(free_bytes)
        )
        return free_bytes.value
    else:
        s = os.statvfs(folder)
        return s.f_bfree * s.f_frsize


def download_file(url, filename):
    print(f"Downloading large file {path.basename(filename)}. Please be patient...")
    try:
        req = urllib.request.urlopen(url)
    except urllib.error.URLError:
        print("ERROR: File not found on server.\nPlease check your internet connection.")
        sys.exit(1)

    CHUNK = 128 * 1024
    with open(filename, "wb") as fp:
        while True:
            chunk = req.read(CHUNK)
            if not chunk:
                break
            fp.write(chunk)
    print(f"Downloading {path.basename(filename)} finished successfully.")
    return filename


def checksum_file(filename, algo="sha256", chunksize=2 ** 20):
    algo = algo.lower()
    if algo == "md5":
        h = hashlib.md5()
    elif algo in ("sha256", "sha-256"):
        h = hashlib.sha256()
    else:
        print(f"ERROR: Unsupported checksum algorithm: {algo}")
        sys.exit(1)
    with open(filename, "rb") as fh:
        for chunk in iter(lambda: fh.read(chunksize), b""):
            h.update(chunk)
    return h.hexdigest()


def unzip_file_gzip(filename):
    newfilename, _ = path.splitext(filename)
    with gzip.open(filename, "rb") as zf, open(newfilename, "wb") as out:
        CHUNK = 128 * 1024
        while True:
            chunk = zf.read(CHUNK)
            if not chunk:
                break
            out.write(chunk)
    print(f"Extraction of {path.basename(filename)} finished successfully.")
    return newfilename


def compile_hmm(filename):
    out, err, code = execute(["hmmpress", "-f", filename])
    if code != 0:
        sys.stdout.write(out.decode(errors="ignore"))
        sys.stderr.write(err.decode(errors="ignore"))
        print("ERROR: hmmpress failed.")
        sys.exit(code)


def delete_file(filename):
    try:
        os.remove(filename)
    except OSError:
        pass


# ------------------ Zenodo ------------------
def zenodo_latest_file_info(concept_recid, target_filename=None):
    api_url = f"https://zenodo.org/api/records/{concept_recid}/versions/latest"
    with urllib.request.urlopen(api_url) as resp:
        data = json.loads(resp.read().decode("utf-8"))

    files = data.get("files", [])
    if not files:
        print("ERROR: No files found in latest Zenodo record.")
        sys.exit(1)

    chosen = None
    if target_filename:
        for f in files:
            if f.get("key") == target_filename:
                chosen = f
                break
    if chosen is None:
        for f in files:
            if f.get("key", "").endswith(".tar.gz"):
                chosen = f
                break
    if chosen is None:
        print("ERROR: Could not find a .tar.gz artifact in latest Zenodo record.")
        sys.exit(1)

    ch = chosen.get("checksum", "")
    algo, hexdigest = ("md5", ch)
    if ":" in ch:
        algo, hexdigest = ch.split(":", 1)

    links = chosen.get("links", {}) or {}
    download_url = links.get("self") or links.get("download")
    if not download_url:
        rec_id = data.get("id")
        key = chosen.get("key")
        download_url = f"https://zenodo.org/api/records/{rec_id}/files/{urllib.parse.quote(key)}?download=1"

    rec_id = data.get("id")
    rec_ver = (data.get("metadata", {}) or {}).get("version", "unknown")
    size_bytes = int(chosen.get("size") or 0)

    return download_url, algo.lower(), hexdigest, rec_id, rec_ver, size_bytes


def download_with_verify(url, filename, algo, expected_hex):
    needs_download = True
    if path.exists(filename):
        have = checksum_file(filename, algo=algo)
        if have == expected_hex:
            needs_download = False
        else:
            print(f"{algo.upper()} mismatch: expected {expected_hex}, got {have}")
    if needs_download:
        download_file(url, filename)
        have = checksum_file(filename, algo=algo)
        if have != expected_hex:
            print(f"Error downloading {filename}, {algo.upper()} mismatch. Expected {expected_hex}, got {have}.")
            sys.exit(1)


# ------------------ Domain-specific steps ------------------
def ensure_dirs():
    base = path.join(os.getcwd(), "antismash", "generic_modules")
    dirs = {
        "fullhmmer": path.join(base, "fullhmmer"),
        "clusterblast": path.join(base, "clusterblast"),
        "smcogs": path.join(base, "smcogs"),
    }
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    return dirs


def download_pfam(fullhmmer_dir):
    gz_path = path.join(fullhmmer_dir, PFAM_URL.rpartition("/")[2])
    if not path.exists(gz_path):
        download_file(PFAM_URL, gz_path)

    have = checksum_file(gz_path, algo="sha256")
    if have != PFAM_CHECKSUM_SHA256:
        print(f"Error downloading {gz_path}, SHA256 mismatch. Expected {PFAM_CHECKSUM_SHA256}, got {have}.")
        sys.exit(1)

    hmm_path = unzip_file_gzip(gz_path)
    compile_hmm(hmm_path)
    delete_file(gz_path)


def download_clusterblast(cluster_dir, overwrite="ask"):
    url, algo, hexdigest, rec_id, rec_ver, _ = zenodo_latest_file_info(
        CLUSTERBLAST_CONCEPT_RECID, target_filename=CLUSTERBLAST_FILENAME
    )
    print(f"Latest ClusterBlast: Zenodo record {rec_id} (version: {rec_ver})")

    tarball = path.join(cluster_dir, CLUSTERBLAST_FILENAME)
    download_with_verify(url, tarball, algo=algo, expected_hex=hexdigest)

    # confirm overwrite
    if overwrite == "ask" and any(path.isfile(path.join(cluster_dir, f)) for f in os.listdir(cluster_dir)):
        reply = input(f"ClusterBlast files already exist in {cluster_dir}. Overwrite? [y/N]: ").strip().lower()
        if reply not in ("y", "yes"):
            print("Keeping existing ClusterBlast database. Skipping extraction.")
            return
    elif overwrite == "no":
        print("Keeping existing ClusterBlast database. Skipping extraction.")
        return

    with tarfile.open(tarball, "r:gz") as tar:
        tar.extractall(path=cluster_dir)

    print(f"Extraction of {path.basename(tarball)} finished successfully.")
    delete_file(tarball)


# ------------------ Main ------------------
def main():
    parser = argparse.ArgumentParser(description="Download and install Pfam and ClusterBlast databases.")
    parser.add_argument(
        "--overwrite",
        choices=["ask", "yes", "no"],
        default="ask",
        help="Whether to overwrite existing ClusterBlast database files (default: ask).",
    )
    args = parser.parse_args()

    dirs = ensure_dirs()
    download_pfam(dirs["fullhmmer"])
    download_clusterblast(dirs["clusterblast"], overwrite=args.overwrite)
    compile_hmm(path.join(dirs["smcogs"], "smcogs.hmm"))


if __name__ == "__main__":
    main()