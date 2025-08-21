class Signature(object):
    """Secondary metabolite signature"""
    def __init__(self, name, _type, description, cutoff, path):
        self.name = name
        self.type = _type
        self.description = description
        self.cutoff = cutoff
        self.path = path

from antismash.generic_modules import (
        fullhmmer,
#        fullhmmer_dblookup,
        genefinding,
        hmm_detection,
        clusterblast,
        subclusterblast,
        knownclusterblast,
        smcogs,
        coexpress,
        # gff_parser,
        tfbs_finder
    )

def check_prereqs(options):
    failure_msgs = []

    if options.full_hmmer or options.inclusive:
        failure_msgs.extend(fullhmmer.check_prereqs(options))

    failure_msgs.extend(genefinding.check_prereqs(options))
    failure_msgs.extend(hmm_detection.check_prereqs())

    if options.smcogs:
        failure_msgs.extend(smcogs.check_prereqs(options))

    if options.clusterblast:
        failure_msgs.extend(clusterblast.check_prereqs(options))

    if options.subclusterblast:
        failure_msgs.extend(subclusterblast.check_prereqs(options))

    if options.knownclusterblast:
        failure_msgs.extend(knownclusterblast.check_prereqs(options))

    if options.coexpress:
        failure_msgs.extend(coexpress.check_prereqs(options))

    if options.tfbs:
        failure_msgs.extend(tfbs_finder.check_prereqs(options))

    return failure_msgs
