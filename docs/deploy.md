# Running fastq2EZbakR

This page provides instructions for running fastq2EZbakR. The Quickstart section includes a code block that succinctly lays out all of the required steps. The Setup section goes into more detail about each step touched upon in the Quickstart.

## Quickstart

All of the steps necessary to deploy fastq2EZbakR are discussed in great detail below. Here, I will present a super succinct description of what needs to be done, with all necessary code included:

``` bash
### 
# PREREQUISITES: INSTALL MAMBA or CONDA AND GIT (only need to do once per system)
###

# CREATE ENVIRONMENT (only need to do once per system)
mamba create -c conda-forge -c bioconda --name deploy_snakemake snakemake snakedeploy

# CREATE AND NAVIGATE TO WORKING DIRECTORY (only need to do once per dataset)
mkdir path/to/working/directory
cd path/to/working/directory

# DEPLOY PIPELINE TO YOUR WORKING DIRECTORY (only need to do once per dataset)
conda activate deploy_snakemake
snakedeploy deploy-workflow https://github.com/isaacvock/fastq2EZbakR.git . --branch main

###
# EDIT CONFIG FILE (need to do once for each new dataset)
###

# RUN PIPELINE

# See [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html) for details on all of the configurable parameters
snakemake --cores all --use-conda --rerun-triggers mtime --keep-going
```

## Detailed instructions

There are 4 steps required to get up and running with THE_Aligner

1. [Install conda (or mamba) on your system](#conda). This is the package manager that THE_Aligner uses to make setting up the necessary dependencies a breeze.
1. [Deploy workflow](#deploy) with [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html)
1. [Edit the config file](#config) (located in config/ directory of deployed/cloned repo) to your liking
1. [Run it!](#run)

The remaining documentation on this page will describe each of these steps in greater detail and point you to additional documentation that might be useful.

### Install conda (or mamba)<a name="conda"></a>
[Conda](https://docs.conda.io/projects/conda/en/latest/index.html) is a package/environment management system. [Mamba](https://mamba.readthedocs.io/en/latest/) is a newer, faster, C++ reimplementation of conda. While often associated with Python package management, lots of software, including all of the fastq2EZbakR pipeline dependencies, can be installed with these package managers. They have pretty much the same syntax and can do the same things, so I highly suggest using Mamba in place of Conda whenever possible. 

**NOTE**: As of version 22.11 of conda, you can now use the mamba solver to get the speed advantages of Mamba. Thus, another option is to follow the steps outlined [here](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community) to use this solver in your conda installation. If you do this though, you will need to add `--conda-frontend conda` to your calls to `snakemake`, as Snakemake throws an error if it can't find `mamba`.

One way to install Mamba is to first install Conda following the instructions at [this link](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). Then you can call:

``` bash
conda install -n base -c conda-forge mamba
```
to install Mamba.

A second strategy would be to install Mambaforge, which is similar to something called Miniconda but uses Mamba instead of Conda. I will reproduce the instructions to install Mambaforge below, as this is probably the easiest way to get started with the necessary installation of Mamba. These instructions come from the [Snakemake Getting Started tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html), so go to that link if you'd like to see the full original details:

* For Linux users with a 64-bit system, run these two lines of code from the terminal:

``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh -o Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```
* For Mac users with x86_64 architecture: 
``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh -o Mambaforge-MacOSX-x86_64.sh
bash Mambaforge-MacOSX-x86_64.sh
```
* And for Mac users with ARM/M1 architecture:
``` bash
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh -o Mambaforge-MacOSX-arm64.sh
bash Mambaforge-MacOSX-arm64.sh
```

When asked this question:
``` bash
Do you wish the installer to preprend the install location to PATH ...? [yes|no]
```
answer with `yes`. Prepending to PATH means that after closing your current terminal and opening a new one, you can call the `mamba` (or `conda`) command to install software packages and create isolated environments. We'll be using this in the next step.

### Deploy workflow<a name="deploy"></a>

Snakemake pipelines can be deployed using the tool [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html). This is often more convenient than cloning the full repository locally. To get started with Snakedeploy, you first need to create a simple conda environment with Snakemake and Snakedeploy:


``` bash
mamba create -c conda-forge -c bioconda --name deploy_snakemake snakemake snakedeploy
```

Next, create a directory that you want to run THE_Aligner in (I'll refer to it as `workdir`) and move into it:
``` bash
mkdir workdir
cd workdir
```

Now, activate the `deploy_snakemake` environment and deploy the workflow as follows:

``` bash
conda activate deploy_snakemake
snakedeploy deploy-workflow https://github.com/isaacvock/fastq2EZbakR.git . --branch main
```

`snakedeploy deploy-workflow https://github.com/isaacvock/fastq2EZbakR.git` copies the content of the `config` directory in fastq2EZbakR's Github repo into the directoy specified (`.`, which means current directory, i.e., `workdir` in this example). It also creates a directory called `workflow` that contains a singular Snakefile that instructs Snakemake to use the workflow hosted on the main branch (that is what `--branch main` determines) of fastq2EZbakR's Github repo. `--branch main` can be replaced with any other existing branch.

### Edit the config file<a name="config"></a>

In the `config/` directory you will find a file named `config.yaml`. If you open it in a text editor, you will see several parameters which you can alter to your heart's content. See the [Configuration](configuration.md) section for a detailed description of all parameters present in fastq2EZbakR's config.yaml file.

### Run it!<a name="run"></a>

Once steps 1-3 are complete, fastq2EZbakR can be run from the directory you deployed the workflow to as follows:

``` bash
snakemake --cores all --use-conda --rerun-triggers mtime
```

The `--rerun-triggers mtime` addition is a suggestion that will prevent the pipline from rerunning certain steps who's output already exists and who's input has not been modified since the last run. See [this post](https://github.com/snakemake/snakemake/issues/1694) for a discussion as to why this is necessary as of Snakemake version 7.8.0.

Some additional parameters that you might find useful include:

* `--show-failed-logs`: When you include this, the log files for any rules that fail will be print to the screen. This can make it easier to figure out which steps went wrong and to quickly check the error messages.
* `--keep-going`: Including this will make Snakemake continue running independent rules even if one rule fails. Snakemake doesn't do this by default because the idea is if something went wrong in a rule, there could be upstream problems lurking that will plague all downstream rules. The `--keep-going` option can be useful in pipelines like fastq2EZbakR though, where there are two independent rules at the end of the pipeline (makecB and maketdf) that don't need each other to complete successfully for the other to complete, and where one of the rules (maketdf) can fail due to a number of reasons (not enough available RAM, IGVtools bugs, etc.) that will not impact the other rule (makecB).

There are **A LOT** of adjustable parameters that you can play with when running a Snakemake pipeline. I would point you to the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html) 
for the details on everything you can change when running the pipeline.

