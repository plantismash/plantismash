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

"""Configuration handling for antiSMASH

"""

import sys
from os import path
import ConfigParser
from argparse import Namespace

_config = None
_basedir = path.dirname(path.abspath(__file__))
_default_name = 'default.cfg'
_sys_name = sys.platform + '.cfg'
_user_file_name = path.expanduser('~/.antismash.cfg')
_instance_file_name = 'instance.cfg'

def load_config(namespace):
    """Load config from a default and system-specific config file and
    add it to a namespace object, but don't overwrite existing settings
    """
    default_file = path.join(_basedir, _default_name)
    sys_file = path.join(_basedir, _sys_name)
    instance_file = path.join(_basedir, _instance_file_name)

    # load generic configuration settins
    config = ConfigParser.ConfigParser()
    with open(default_file, 'r') as fp:
        config.readfp(fp)

    # load system-specific config file if available
    # also load .antismash.cfg from the user's home dir
    # and last, overriding all the other settings, instance.cfg
    config.read([sys_file, _user_file_name, instance_file])

    for s in config.sections():
        if s not in namespace:
            namespace.__dict__[s] = Namespace()
        for key, value in config.items(s):
            key = key.replace('-', '_')
            if key not in namespace.__dict__[s]:
                namespace.__dict__[s].__dict__[key] = value

    # settings from the [DEFAULT] section go to the global namespace
    for key, value in config.items('DEFAULT'):
        key = key.replace('-', '_')
        if key not in namespace:
            namespace.__dict__[key] = value

def set_config(namespace):
    """Set a namespace object to be the global configuration"""
    global _config
    _config = namespace

def get_config():
    """Get the global configuration"""
    return _config
