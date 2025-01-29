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

"""Script to download and install Pfam and ClusterBlast databases"""

import urllib.request, urllib.error, urllib.parse
import tarfile
import gzip
import backports.lzma as lzma
import hashlib
import os
import sys
import subprocess
import platform
import ctypes
from os import path

PFAM_URL = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam27.0/Pfam-A.hmm.gz"
PFAM_CHECKSUM = "b29bc2c54db8090531df0361a781b8d7397f60ebedc0c36a16e7d45e999cc329"
# CLUSTERBLAST_URL = "https://bitbucket.org/antismash/antismash/downloads/clusterblast_dmnd07.tar.xz"
# CLUSTERBLAST_CHECKSUM = "0b4911ee3e30bc2fbe00b85fdd100ec7389c835c22527586d7b93ab661a0a57b"

if sys.platform == ('win32') or sys.platform == ('darwin'):
    os.environ['EXEC'] = os.getcwd() + "\exec"
    os.environ['PATH'] = os.pathsep + os.environ['EXEC'] + os.pathsep + os.environ['PATH']


# pylint: disable=redefined-builtin
def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    proc = subprocess.Popen(commands, stdin=stdin_redir,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    out, err = proc.communicate(input=input)
    retcode = proc.returncode
    return out, err, retcode
# pylint: enable=redefined-builtin


def get_remote_filesize(url):
    try:
        usock = urllib.request.urlopen(url)
        dbfilesize = usock.info().get('Content-Length')
        if dbfilesize is None:
            dbfilesize = 0
    except urllib.error.URLError:
        dbfilesize = 0
    dbfilesize = float(int(dbfilesize))  # db file size in bytes
    return dbfilesize


def get_free_space(folder):
    """ Return folder/drive free space (in bytes)"""
    if platform.system() == 'Windows':
        free_bytes = ctypes.c_ulonglong(0)
        ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(folder), None, None, ctypes.pointer(free_bytes))
        return free_bytes.value
    else:
        return os.statvfs(folder).f_bfree * os.statvfs(folder).f_frsize


def check_diskspace(file_url):
    """ Check if sufficient disk space is available"""
    dbfilesize = int(get_remote_filesize(file_url))
    free_space = int(get_free_space("."))
    if free_space < dbfilesize:
        print('ERROR: Insufficient disk space available (required: %d, free: %d).' % (dbfilesize, free_space))
        sys.exit()


def download_file(url, filename):
    """ Download a file"""
    print("Downloading large file %s. Please be patient..." % (path.basename(filename)))
    try:
        req = urllib.request.urlopen(url)
    except urllib.error.URLError:
        print('ERROR: File not found on server.\nPlease check your internet connection.')
        sys.exit()
    CHUNK = 128 * 1024
    with open(filename, 'wb') as fp:
        while True:
            try:
                chunk = req.read(CHUNK)
                if not chunk:
                    break
                fp.write(chunk)
            except IOError:
                print('ERROR: Download interrupted.')
                sys.exit()
    #Report download success
    print("Downloading %s finished successfully." % (path.basename(filename)))
    return filename


def checksum(filename, chunksize=2 ** 20):
    """Get the SHA256 checksum of a file"""
    sha = hashlib.sha256()
    with open(filename, 'rb') as fh:
        for chunk in iter(lambda: fh.read(chunksize), b''):
            sha.update(chunk)

    return sha.hexdigest()


def unzip_file(filename, decompressor, error_type):
    """Decompress a compressed file"""
    newfilename, _ = path.splitext(filename)
    try:
        zipfile = decompressor.open(filename, 'rb')
        chunksize = 128 * 1024
        with open(newfilename, 'wb') as fp:
            while True:
                try:
                    chunk = zipfile.read(chunksize)
                    if not chunk:
                        break
                    fp.write(chunk)
                except IOError:
                    print('ERROR: Unzipping interrupted.')
                    sys.exit()
    except error_type:
        print("ERROR: Error extracting %s. Please try to extract it manually." % (path.basename(filename)))
        return
    print("Extraction of %s finished successfully." % (path.basename(filename)))
    return newfilename


def untar_file(filename):
    """ Extract a TAR/GZ file"""
    try:
        tar = tarfile.open(filename)
        tar.extractall(path=filename.rpartition(os.sep)[0])
        tar.close()
    except tarfile.ReadError:
        print("ERROR: Error extracting %s. Please try to extract it manually." % (filename.rpartition(os.sep)[2]))
        return
    print("Extraction of %s finished successfully." % (filename.rpartition(os.sep)[2]))


def compile_pfam(filename):
    """ Compile a HMMer database with hmmpress"""
    command = ['hmmpress', '-f', filename]
    execute(command)


def delete_file(filename):
    try:
        os.remove(filename)
    except OSError:
        pass


def download_if_not_present(url, filename, sha256sum):
    """Download a file if it's not present or checksum doesn't match"""
    if path.exists(filename):
        print("Creating checksum of %s" % path.basename(filename))
        csum = checksum(filename)
        if csum == sha256sum:
            return
        else:
            print("checksum mismatch: expected %s, got %s" % (sha256sum, csum))

    download_file(url, filename)

    print("Creating checksum of %s" % path.basename(filename))
    csum = checksum(filename)
    if csum != sha256sum:
        print("Error downloading %s, sha256sum mismatch. Expected %s, got %s." % \
              (filename, sha256sum, csum))
        sys.exit(1)


def download_clusterblast():
    """Download and extract the ClusterBlast database"""
    
    CLUSTERBLAST_URL = "https://plantismash.bioinformatics.nl/clusterblast/clusterblast.tar.gz"
    CLUSTERBLAST_CHECKSUM = "1e032f1f3b297e8e924a328ef26912d324cb01deae25332933f2ba50ff7af2fd"
    
    # Check for sufficient disk space before download
    check_diskspace(CLUSTERBLAST_URL)
    
    # Define the filename for the downloaded file
    filename = path.join(os.getcwd(), "antismash", "generic_modules", "clusterblast_dmnd07.tar.gz")
    
    # Download the file if not present or if the checksum does not match
    download_if_not_present(CLUSTERBLAST_URL, filename, CLUSTERBLAST_CHECKSUM)
    
    # Extract the .gz file to its contents
    extracted_filename = unzip_file(filename, gzip, gzip.zlib.error)
    
    # Extract the contents of the tar file into the generic_modules/clusterblast directory
    clusterblast_dir = path.join(os.getcwd(), "antismash", "generic_modules", "clusterblast")
    try:
        with tarfile.open(extracted_filename, 'r') as tar:
            for member in tar.getmembers():
                # Remove any leading directories from the member's path to avoid subdirectory creation
                member.name = path.relpath(member.name, start=member.name.split('/')[0])
                tar.extract(member, path=clusterblast_dir)  # Extract directly into clusterblast directory
    except tarfile.ReadError:
        print(("ERROR: Error extracting %s. Please try to extract it manually." % extracted_filename))
        return
    
    print(("Extraction of %s finished successfully." % extracted_filename))
    
    # Clean up by deleting the original compressed files
    delete_file(filename)  # delete the .tar.gz file
    delete_file(extracted_filename)  # delete the .tar file




def main():
    #Download and compile PFAM database
    check_diskspace(PFAM_URL)
    filename = path.join(os.getcwd(), "antismash", "generic_modules", "fullhmmer",
                         PFAM_URL.rpartition('/')[2])
    download_if_not_present(PFAM_URL, filename, PFAM_CHECKSUM)
    filename = unzip_file(filename, gzip, gzip.zlib.error)
    compile_pfam(filename)
    delete_file(filename + ".gz")

    # download the clusterblast database 
    download_clusterblast()

    #Download and compile ClusterBlast database
    # check_diskspace(CLUSTERBLAST_URL)
    # filename = path.join(os.getcwd(), "antismash", "generic_modules", "clusterblast",
    #                      CLUSTERBLAST_URL.rpartition('/')[2])
    # download_if_not_present(CLUSTERBLAST_URL, filename, CLUSTERBLAST_CHECKSUM)
    # filename = unzip_file(filename, lzma, lzma.LZMAError)
    # untar_file(filename)
    # delete_file(filename)
    # delete_file(filename + ".xz")

    # hmmpress the NRPS/PKS specific databases
    # compile_pfam(path.join("antismash", "specific_modules", "nrpspks", "abmotifs.hmm"))
    # compile_pfam(path.join("antismash", "specific_modules", "nrpspks", "dockingdomains.hmm"))
    # compile_pfam(path.join("antismash", "specific_modules", "nrpspks", "ksdomains.hmm"))
    # compile_pfam(path.join("antismash", "specific_modules", "nrpspks", "nrpspksdomains.hmm"))

    # hmmpress the smcog specific database
    compile_pfam(path.join("antismash", "generic_modules", "smcogs", "smcogs.hmm"))

if __name__ == "__main__":
    main()
