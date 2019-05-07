# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""
Plant cyclopeptides repetitive sequences detection module
"""
from antismash import utils
from .specific_analysis import specific_analysis
from .html_output import generate_details_div, generate_sidepanel, will_handle

name = "plant_cyclopeptides"

# The tuple is the name of the binary and whether it is an optional requirement
_required_binaries = [
    ]

short_description = name.capitalize()

def check_prereqs():
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)
    return failure_messages

__all__ = [ check_prereqs, specific_analysis, generate_sidepanel, generate_details_div, will_handle ]