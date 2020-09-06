## Overview

This document describes the setup process to allow you to run jobs on remote machines directly from your local machine.

- on SHARCnet (we are using `graham.sharcnet.ca`)
- on "yushan" (`jdserv.mcmaster.ca`)
- on `earnserv[12]@mcmaster.ca`

**Fixme**: earnserv1 is dead

## set up your remote account

### SHARCnet (`graham.mcmaster.ca`)

You need an active Compute Canada account, which then allows you to request a SHARCnet account, which needs to be verified by a PI. My (BMB) CC identifier is `cze-501`, in case that's useful. Once you get your account, I believe you have to sit through a training webinar before you can be [certified](https://www.sharcnet.ca/my/systems/user_certification) to use more than 8 CPUs at a time or run jobs > 24 hours.  (After that the cap goes up to 256 CPUs.)

### yushan (`jdserv.mcmaster.ca`)

Ask Jonathan Dushoff for access.

### earnserv (`earnserv[12].mcmaster.ca`)

Ask David Earn for access. (If you have access, your username and password are shared with the math & stats (`ms.mcmaster.ca`) server.)

## If necessary, set up an RSA key pair: see [here](https://www.ssh.com/ssh/copy-id)

Use your own judgment/degree of paranoia about whether to password-protect your RSA key and/or set up a local storage daemon (so you can have the security of a password-protected key without the hassle of entering your password multiple times per session)

## Copy SSH public key (from local machine to remote machine) [one time]

You need this for passwordless communication from your local machine to the remote machine

```
ssh ssh-copy-id -i ~/.ssh/rsa_id.pub USERNAME@REMOTE
```

## Set up VPN (install once on local machine; run every session on local machine)

- McMaster servers (yushan, earnserv) aren't accessible beyond the campus network unless you are running VPN software. If you want to query/run jobs directly from your off-campus machine you'll need to have a VPN connection
- For Linux machines you can download a tarball from [here](https://softdepot.mcmaster.ca/anyconnect-linux64-4.6.04056-predeploy-k9.tar.gz) (you need to enter your McMaster username and password); unpack the tarball somewhere; and run `vpn_install.sh` as root (!)
- using the Linux VPN, you will be prompted twice about dodgy-looking certificates ("Connect Anyway" rather than "Keep Me Safe") in addition to the username/password prompt

## Store github credentials (on remote machine)

In order to be able to pull onto the remote machine (**don't do this if you're paranoid/careful**)

In your home directory:

```
cat >~/.netrc <<EOF
machine github.com login YOURGITHUBNAME password YOURGITHUBPASSWORD
EOF
chmod +600 ~/.netrc  ## make file write-only
```

- there should be a way to store credentials, but I haven't figured it out
- might also be a way to interactively enter credentials?
- two alternatives: 

* create a new SSH key pair that you use to communicate from the remote machine to Github (this is the best/most principled approach); [add the new SSH key to your Github account](https://help.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account)
* copy your **private** key onto the remote machine. I definitely wouldn't do this onto SHARCnet, but it might be OK on yushan?  This also allows communication within yushan

## Optional (yushan only, one time): set up passwordless communication from head node to servers

* if you want to be able to communicate from yushan's head to the workers, you need to be able to passwordlessly communicate. You can do this by setting up a new RSA key pair (since the workers and the head share a file system, the public key will automatically be accessible to the workers).
* However, this is only necessary if you want to be able to run commands like `uptime` (see the `yup` alias) on the workers. This, in turn, is only necessary if we have a mix of scheduled and directly run jobs on yushan. Otherwise, `qsub`/`qstat`/`qstat -g c` (and their remote aliases `ymake`/`yq`/`yqc`) should be sufficient for communicating with the workers.

## Clone the repo (on remote machine; one time)

```
git clone https://github.com/mac-theobio/PHAC_covid.git
```

(I normally use `git@` rather than `https://`, but (?) can't use `git@` unless we also set up an SSH *private* key for talking to GitHub ...)

## create personal library (on remote machine; one time)

```
## get major version
mver=`R --version | grep "R version" | awk '{print $3}' | cut -c1-3`
mkdir -p ~/R/x86_64-pc-linux-gnu-library/$mver
```

* This is hardcoded, which is too bad, but R finds it automatically (don't need to set libPaths). The alternative is to put it wherever you like and make sure to run `.libPaths("your_pkg_dir")` in every R session ...
* This step is unnecessary if you've already established a personal library from previous use, or if the remote admin is willing to install/update all the necessary packages in the system library
* FIXME: there may be a more direct/less hacky way to extract the major version and/or the expected location of the personal library from R ... ??

## load modules (on graham, every session)

This is built into the `snmake` alias, so you don't actually need it unless you're working interactively (e.g. for debugging, or for the step below).

```
module load nixpkgs/16.09  gcc/8.3.0; module load r/4.0.0
```

## Installing packages (on remote machine; once)

On graham, this has to be done on the head node (which should be OK, it's not *too* intensive) *or* we have to download all the relevant tarballs and install locally. The latter is a pain in the butt because there isn't/I don't know of a good way to handle dependencies/ordering properly.

- run `setup.R` (e.g. via `R CMD BATCH setup.R`, or we can put it in the Makefile): this (1) installs the `remotes` package; (2) installs `bbmle` and `McMasterPandemic` from Github; (3) installs the rest of the packages listed in the `pkgs` file

(Let's hope we don't need Stan! (Would need to install dependencies; download tarball; set up a batch job with enough memory & time to install it))

## setting up paths (on yushan, one time)

Non-shell logins work a little differently from shell (interactive) logins, and on yushan they don't work quite as expected. Create a `.bashrc` file in your home directory if you don't have one already, and add the following lines to it:

```
export PATH=$PATH:/usr/local/sge/2011.11p1/bin/linux-x64
export SGE_ROOT=/usr/local/sge/2011.11p1
```

to test this, try `ssh USERNAME@yushan qstat -g c` from your local machine (with VPN running)

(Presumably a configuration issue RHPCS could fix centrally ... ?)

## make cache on remote machine

(To make `wrapR` make rules work) If Dropbox isn't set up on the remote machine you need:

```
ssh REMOTEHOST mkdir $PHACDIR/cache
```

## Local configuration (on local machine)

(From here on I'm assuming you're using `bash` as your shell.)

Create the `batch_setup` file that defines your username and the location of your `PHAC_covid` directory on the remote machines (yushan/SHARCnet/earnserv), as well as the maximum number of cores to use for local runs. Mine looks like this:

```
## default user name
export USER=`whoami`
export WD=PHAC_covid
## username and location for yushan (change defaults if necessary)
export YUSHAN_USER=bbolker  ## override default
export YUSHAN_WD=$WD
export SN_USER=$USER
export SN_WD=$WD
export EARNSERV_USER=$USER
export EARNSERV_WD=$WD
export earnserv_MAXCORES=10
## max number of cores to use locally
export MAXCORES=5
```

## Running jobs (on local machine)

- locally, run `. batch_aliases.mk` in the shell (with the repo as your working directory). You should then have access to
    - `yq`, `snq`: query the batch queue on yushan or graham (`yq` lists all jobs; `snq`, just your own)
	- `eload x`: print load average on `earnservx` where `x`=1 or 2
    - `ypull`, `snpull`, `epull`: pull the `PHAC_covid` repo remotely
    - `ymake`, `snmake`, `emake x`: run a make rule remotely (`x` as above)
- within your R script:
    - `source("batchtools.R")` (or, if using `wrapR`, make your target depend on `batchtools.Rout`)
    - run `batch_setup()` to set up the proper batch environment (works so far for yushan, graham, earnserv1 and 2, or local parallel computation). **note**: for graham, the current version assumes a "wall clock time" of 6 hours (per job) by default, you can change this by specifying the wall clock time **in seconds**: `batch_setup(walltime=4*3600)` for a 4-hour wall clock time. On earnserv or locally, you can also specify a `workers=` argument to specify the maximum number of cores to use. (The default is set by your `batch_setup` script.)
    - run `library("furrr"); library("future.batchtools")`
    - use one of the `future_map*` functions to set up your jobs (`future_map()` to return a list, `future_map_dfr()` to combine returned data frames into one big data frame)

## testing/troubleshooting

- `yup`
- `yq`
- `ypull`
- `ymake batch_test.Rout`

## Files

- this file (`batch_notes.md`)
- `batch_setup`: user-level configuration (remote usernames, working directories, resource limitations, etc.)
- `batch_aliases`: utilities for remote communication and spawning `make` jobs remotely
- `batch_test.R`: test code (e.g. `ymake batch_test.Rout`)
- `batchtools.R`: defines `batch_setup()`

## To do/issues

- clean up / DRY
- make rules and/or aliases for rsync'ing caches
- it's OK to run _short_ jobs on the head (login) node of `graham` (e.g. installing packages), and anything that needs network access needs to be run on the head node, but we have to be **very** careful not to accidentally run anything big directly on the head node
- can we share an R package directory across users on sharcnet?

---

## further notes/comments

The current system is set up to use `furrr`/`future.batchtools` as an interface to whatever batch system is running on a given server (i.e. SLURM scheduler [SHARCnet], SGE scheduler [yushan], nothing [earnserv]). The following components need to fit together:

- remote execution via `ssh`: special incantations are needed to make sure that jobs can be started remotely and run in the background (i.e., won't hang when the process that created them quits)
- `make`/`makestuff`/wrapR
- scheduler (e.g. `qsub`): some schedulers need to be passed a *script* (don't see a way to pass a simple command line), with command line arguments or a template file to determine job characteristics

There a variety of ways one can interact with the scheduler to submit a bunch of jobs:

- via `for` loop in a shell script
- *job arrays*, e.g. for [SGE](https://arc.leeds.ac.uk/using-the-systems/why-have-a-scheduler/advanced-sge-task-arrays/) or [SLURM](https://slurm.schedmd.com/job_array.html) (the arguments to individual jobs take numeric values). Each scheduler passes the information about which job in an array is being run to the downstream script a little differently ...
- via R, using `furrr` or `foreach` (via the `doFuture` package)
