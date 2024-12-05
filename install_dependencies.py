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
# Copyright (C) 2024 Elena Del Pup 
# Wageningen University & Research, NL
# Bioinformatics Group, Department of Plant Sciences 
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""Install dependencies for plantiSMASH"""

import subprocess
import os

def install_packages(packages):
    for package in packages:
        try:
            # Construct the mamba install command
            command = ['mamba', 'install', '-y', package]
            print("Installing {}...".format(package))
            # Execute the command
            subprocess.check_call(command)
            print("{} installed successfully.".format(package))
        except subprocess.CalledProcessError as e:
            print("Failed to install {}: {}".format(package, e))

if __name__ == "__main__":

    # List of packages to be installed
    packages_to_install = [
        'glimmer',
        'glimmerhmm',
        'hmmer2',
        'hmmer=3.1b2',
        'fasttree=2.1.8',
        'diamond=0.8.9=boost1.60_1',
        'muscle=3.8.31',
        'prodigal',
        'blast=2.2.31',
        'cd-hit=4.6.6',
        'xz=5.2.2',
        'libxml2=2.9.3',
        "pplacer",
        'biopython == 1.76',
        'backports.lzma',
        "ncbi-datasets-cli",
        "unzip"
    ]

    # Install the packages
    install_packages(packages_to_install)

    command = ['pip', 'install', "-r", "requirements.txt"]
    print("pip installing requirements in requirements.txt...")
    subprocess.check_call(command)

    # have put graphlan in plantismash folder ahead in case of permission issue
    # command = ["git", "clone", "https://github.com/biobakery/graphlan.git"]
    # print("Cloning graphlan repository...")
    # subprocess.check_call(command)
    command = ["pip", "install", "graphlan/."]
    print("Installing graphlan...")
    subprocess.check_call(command)

    print("python2 run_antismash.py -h")
    subprocess.check_call(["python2", "run_antismash.py", "-h"])
    print("if any errors, please use conda list to check any missing dependencies and install those manually")
    # unable use nrpspks in specific_modules from antismash, but that is not for plantismash, so it was deleted
    print("!! genefinding module cannot work, so --genefinding none is default now!!")