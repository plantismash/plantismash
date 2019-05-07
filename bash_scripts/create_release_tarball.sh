#!/bin/bash
# Create an antiSMASH release tarball

# turn on paranoia
set -o errexit
set -o nounset

if [ $# -lt 2 ]; then
    echo "Usage $0 <version_number> <output directory>"
    exit 2
fi

# more readable variable names
readonly VERSION=$1
readonly OUTPUT_DIR=$(readlink -f $2)
readonly RELEASE_TAG="antismash-$(echo $VERSION | tr . -)"
readonly PROG_VERSION="antismash-${VERSION}"
readonly CWD=$(pwd)

# get a temporary directory
TMPDIR=$(mktemp -d)

function cleanup {
    rm -rf "$TMPDIR"
}
trap cleanup EXIT

# build the tarball
echo "Creating ${PROG_VERSION}.tar.gz in ${OUTPUT_DIR} for ${RELEASE_TAG}"
cd $TMPDIR
git clone git@bitbucket.org:antismash/antismash.git antismash-release
cd antismash-release
git archive --format=tar.gz --prefix="${PROG_VERSION}/" $RELEASE_TAG > "${OUTPUT_DIR}/${PROG_VERSION}.tar.gz"
cd $CWD
