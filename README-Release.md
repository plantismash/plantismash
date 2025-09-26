PlantiSMASH release process
=========================

This file documents the release process of plantiSMASH.


Preparation work
----------------

Make sure all relevant pull requests have been merged, and check if any
showstopper bugs are still open. Showstopper bugs mainly include regressions
from previous versions.


What version number will the new release get?
---------------------------------------------

PlantiSMASH is using [semantic versioning](http://semver.org/). Unfortunately,
semantic versioning focuses on libaries and other things that present APIs,
making it an imperfect match for antiSMASH.

Reflecting the basic ideas behind semantic versioning, we should consider our
command line the "API" of antiSMASH for now, as antiSMASH is certainly called
from many in-house scripts. As a result, we should be increasing the MINOR
version if we add new features that have additional command line options, and
increase the MAJOR version when we remove command line options.

We should also increase the MINOR version if we add new secondary metabolite
classes to be detected. This rule wasn't followed previously, but should be
followed for future releases.


Tag the actual release
----------------------

* Make sure the `CONTRIBUTORS` file is up-to-date. Update and commit if
  necessary.
* Update the version number in `antismash/__init__.py` and commit. You can check what is the current specified version by running `python run_antismash.py -h`. This will automatically update the version specified in the 
* Grab a list of detailed changes using `git shortlog <last-release-tag>.. >
  changes.txt`
* Edit `changes.txt` to add some release notes on top.
* Tag the release using `git tag -s -F changes.txt
  plantismash-<MAJOR>-<MINOR>-<PATCH>`, filling in the correct values for MAJOR,
  MINOR and PATCH, of course.
* The releases are automatically tracked in [Zenodo](https://zenodo.org)
* Make sure to tag the right version of the changelog for the detection rules in the script `antismash/output_modules/html/generator.py` function `add_overview_entry`. E.g., "https://plantismash.github.io/documentation/changelog/2.0/#supported-cluster-types-version-2"



Post-push work
--------------

Now that the git repository contains the new release, a couple of further steps
are required.

### Create release tarball

This is simply done by calling the appropriate script:
`./bash_scripts/create_release_tarball.sh <VERSION> <OUTPUT_DIR>`

and verify it looks sane `tar tf <OUTPUT_DIR>/antismash-<VERSION>.tar.gz`

### Upload the tarball to dl.secondarymetabolites.org for testing

`scp <OUTPUT_DIR>/antismash-<VERSION>.tar.gz
root@dl.secondarymetabolites.org:/vol/downloads/`

### Update and test `install_deb.sh`

From the `antismash/installer.git` repository on bitbucket, update the
`install_deb.sh` file with the current version number. Then, use `vagrant up`
for the respective distro to check the installer script.

Connect your browser to the [output page](http://localhost:8080/) to check if
the output looks sane.

Commit and push your changes to the `install_deb.sh` file if everything is fine.

### Upload tarball and installer script to Bitbucket

This is best done in a web browser via the Bitbucket web UI.

### Update the docker images

From the `antismash/docker.git` repository on bitbucket, update the `Dockerfile`
for the `standalone` image. Bump the version number of the file and adjust the
`ANTISMASH_VERSION` variable.

Then, rebuild the docker image using `docker build -t antismash/standalone .`,
tag the version using `docker tag antismash/standalone
antismash/standalone:<VERSION>`. Test the docker image. If everything works,
push it to the Docker Hub by using `docker push antismash/standalone && docker
push antismash/standalone:<VERSION>`.

Then, repeat the process for the other container builds.

### Send the release announcement to the mailing list

You can reuse the contents of changes.txt for this. Make sure there's a link to
the [antiSMASH download
page](http://antismash.secondarymetabolites.org/dowload/) in the email.

### Add a notification message on the antiSMASH website

Use the `smashctl notice` tool for this.
