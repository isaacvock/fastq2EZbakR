# Running Snakemake pipelines on HPC systems with Slurm

This page provides instructions for running fastq2EZbakR on HPC clusters that use a Slurm scheduler (such as Yale's clusters). Some of the information here is also presented in the [general deployment documentation](deploy.md), to make this a completely self-contained tutorial for Simon lab members and other Slurm users.

## Important note

Version 8 of Snakemake completely changed how Snakemake submits jobs to clusters. The main walkthrough below describes the root `yale_profile/` configuration, which works with Snakemake versions < 8.0.0. If you are using Snakemake >= 8.0.0, `yale_profile` now also includes a `yale_profile/v8plus/` subdirectory compatible with v8+ (adapted from [this repository](https://github.com/jdblischak/smk-simple-slurm)); to use it, install `snakemake-executor-plugin-cluster-generic` alongside `snakemake` in your `deploy_snakemake` environment, add `yale_profile/v8plus` to your `PATH`, and pass `--profile yale_profile/v8plus/`. 

## Advice: set the SNAKEMAKE_CONDA_PREFIX environment variable

TLDR: Follow these steps to create an environment variable called SNAKEMAKE_CONDA_PREFIX that will specify where to save conda environments that Snakemake creates. Instructions below are specific to Yale HPC, but this is a good idea for anyone running fastq2EZbakR in a folder that is periodically deleted:

1. Create a folder in your project directory called `snakemake_envs`. Technically, you can call it whatever you want, it just has to be whatever you specify in step 5. This can be done by running `mkdir ~/project/snakemake_envs` from anywhere on the cluster.
1. Go to your home directory (where you start when you log in. If you are in some other directory, you can run `cd ~` to get there).
1. If you view all files with `ls -a`, you will see a hidden file called `.bashrc`. This is a script that gets executed everytime you log onto the cluster. Here, you can specify environment variables to be created everytime you log on.
1. Prepare to edit .bashrc (e.g., run `vi .bashrc`, and then click `i` to enter editing mode).
1. Add the following line of code at the bottom of your .bashrc file: `export SNAKEMAKE_CONDA_PREFIX=~/project/snakemake_envs`. This is just the path to the folder in your project directory that you created. 
1. You can exit vi by running `:wq`. To create this new environment variable without logging out of the cluster, run `source ~/.bashrc` from anywhere on the cluster.

Why do this? One of the most common annoyances I would run into when consistently using a Snakemake pipeline on the cluster was the need to recreate all of the pipeline environments. For context, when you run a Snakemake pipeline as described above for the first time, conda environments are created which will contain installations of all dependencies necessary for the various steps of the pipeline. By default, these conda environments will be saved in a hidden `.snakemake` directory in the directory that you launced the pipeline from. The problem is that since you should probably be running these kind of pipelines in your scratch60 folder (for Yale McCleary users), some of the installation files will get deleted in 60 days. This causes the pipeline to fail in weird ways and eventually forces you to recreate all environments by deleting the whole `.snakemake` directory.

To avoid this, you can specify a different directory for Snakemake to save conda environments it creates. The easiest way to do this is to create an environment variable called SNAKEMAKE_CONDA_PREFIX and specify this directory. To do this, you can edit your .bashrc file, found in the directory you start in when you log into the cluster. This is what the instructions above explain how to do.

## Quickstart

All of the steps necessary to deploy a pipeline are discussed in great detail below. Here, I will present a super succinct description of what needs to be done, with all necessary code included:

``` bash
### 
# PREREQUISITES: INSTALL MAMBA AND GIT (only need to do once per system; already installed on Yale HPC)
###

# CREATE ENVIRONMENT (only need to do once per system)
### system = compute system, e.g., a particular cluster
### Syntax for requesting an interactive node and loading conda 
### installation are Yale specific.
### NOTE: newest version of Yale McCleary's miniconda module
### replaces mamba with conda using the libmamba solver. Thus,
### conda rather than mamba is called below:
salloc
module load miniconda
conda create -c conda-forge -c bioconda --name deploy_snakemake 'snakemake<8.0.0' snakedeploy

# CREATE AND NAVIGATE TO WORKING DIRECTORY (only need to do once per dataset)
  # This should be in palmer_scratch or project directory if on Yale McCleary; not the home directory!!
mkdir path/to/working/directory
cd path/to/working/directory

# DEPLOY PIPELINE TO YOUR WORKING DIRECTORY (only need to do once per dataset)
conda activate deploy_snakemake
snakedeploy deploy-workflow https://github.com/isaacvock/fastq2EZbakR.git . --branch main


###
# EDIT CONFIG FILE (need to do once for each new dataset)
###

# COPY PROFILE TO OPTIMIZE DEPLOYMENT ON HPC (only need to do once)
git clone https://github.com/isaacvock/yale_profile.git
cp yale_profile/run_slurm.sh ./

###
# EDIT RUN_SLURM.SH AND PROFILE AS NECESSARY (double check before each pipeline run) 
###

# RUN PIPELINE
module purge
sbatch run_slurm.sh
```


## Detailed instructions

As YCRC has already installed miniconda on all clusters, there are 3 steps required to get up and running with fastq2EZbakR on Yale HPC:

1. [Deploy workflow](#deploy_s) with [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html)
1. [Edit the config file](#config_s) (located in config/ directory of deployed/cloned repo) to your liking
1. [Run it!](#run_s)

### Deploy workflow<a name="deploy_s"></a>

Snakemake pipelines can be deployed using the tool [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html). This is often more convenient than cloning the full repository locally. To get started with Snakedeploy, you first need to create a simple conda environment with Snakemake and Snakedeploy:


``` bash
salloc
module load miniconda
mamba create -c conda-forge -c bioconda --name deploy_snakemake 'snakemake<8.0.0' snakedeploy
```

**NOTE**: Due to major breaking changes in version 8 of Snakemake, the exact strategy for deploying fastq2EZbakR on a slurm scheduling system will only work with an earlier version of Snakemake. This is why `'snakemake<8.0.0'` is specified above.

Next, create a directory that you want to run the pipeline in (I'll refer to it as `workdir`) and move into it. This is best kept in your scratch60 directory to avoid using too much project directory disk space:
``` bash
mkdir workdir
cd workdir
```

Now, activate the `deploy_snakemake` environment and deploy the workflow as follows:

``` bash
conda activate deploy_snakemake
snakedeploy deploy-workflow <path to pipeline .git> . --branch main
```

`snakedeploy deploy-workflow https://github.com/isaacvock/<pipelien.git` copies the content of the `config` directory in the pipeline's Github repo into the directoy specified (`.`, which means current directory, i.e., `workdir` in this example). It also creates a directory called `workflow` that contains a singular Snakefile that instructs Snakemake to use the workflow hosted on the main branch (that is what `--branch main` determines) of the pipeline's Github repo. `--branch main` can be replaced with any other existing branch.

### Edit the config file<a name="config_s"></a>

In the `config/` directory you will find a file named `config.yaml`. If you open it in a text editor, you will see several parameters which you can alter to your heart's content. See the [Configuration](configuration.md) section for a detailed description of all parameters present in fastq2EZbakR's config.yaml file.

### Run it!<a name="run_s"></a>

While at this point you can run fastq2EZbakR as described in the general deployment documentation, this strategy will not make use of the immense amount of computational resources available at your fingertips as a Yale HPC user. 

To make the most of what McCleary and other clusters have to offer, it is best if you specify an in-depth [profile for Snakemake](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). If you don't want to figure out how to do this from scratch, have no fear, I have already created and made available a profile I use for running Snakemake pipelines on McCleary. To use it, clone the repository containing the profile into the directory in which you will run the pipeline:

``` bash
git clone https://github.com/isaacvock/yale_profile.git
```

You can then copy a shell script from this directory that can be used for launching the pipeline from a non-interactive node (so you don't have to sit there and watch the pipeline run):

``` bash
cp yale_profile/run_slurm.sh ./
```

fastq2EZbakR can then be run from your working directory with (first purging any loaded modules to reset the environment. Not strictly necessary but will eliminate an inconsequential warning message at a later step):

``` bash
module purge
sbatch run_slurm.sh
```

## Additional details

### Description of the run script<a name="runscript"></a>

If you check out `run_slurm.sh`, you will see the following simple shell script:

``` bash
#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=PROseq
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@yale.edu

module purge

module load miniconda

conda activate deploy_snakemake

conda config --set channel_priority strict

snakemake --profile yale_profile/ --rerun-triggers mtime
```

The first line specifies the interpreter to use. This is a bash script, so the bash interpreter is specified.

The set of `#SBATCH` lines are resource requests. You'll recognize many of them from the discussion of the profile below. For details about all available job request options, check out [YCRC's documentation](https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/). Of unique note are: 

1. `#SBATCH --job-name=PROseq`: This profile was originally developed for use with the PROseq_etal pipeline I developed. Change the job name to whatever best describes your use case.
1. `#SBATCH --mail-type=ALL` tells the cluster to send you emails when the job starts and finishes.
1. `#SBATCH --mail-user=your.email@yale.edu` is the email address that will receive the above mentioned emails. Change this to be your email, or delete these two lines entirely if you don't want to receive any such emails.

The next three lines sets up the environment. First, all loaded modules are purged to start from a fresh slate and avoid dependency conflicts. Next, miniconda is loaded so that we can activate conda environments. Finally, the `deploy_snakemake` conda environment that we made earlier is activated.

`conda config --set channel_priority strict` ensures that installed dependencies are the exact ones requested.

Finally, Snakemake is called, telling it to use the profile present in the `yale_profile` directory. `--rerun-triggers mtime` makes it so that if you need to restart the pipeline, the criterion for whether a step needs to be rerun is that its input hasn't been modified since you last run the pipeline. This avoids some oddities that can arise when using custom scripts in the pipeline; see [this Issue](https://github.com/snakemake/snakemake/issues/1694) for more details.

### Description of the profile<a name="profile"></a>

Other than `run_slurm.sh`, what is in the `yale_profile` directory? If you peek inside, you will see three additional files: `config.yaml`, `status-sacct-robust.sh`, and a README. `config.yaml` contains all the information that Snakemake will use to request jobs on the cluster with slurm. If you check out its contents, the first section specifies the code that will be run to request a job:

``` yaml
cluster:
  mkdir -p logs/{rule} &&
  sbatch
    --partition=day
    --cpus-per-task={threads}
    --mem={resources.mem_mb}
    --job-name=smk-{rule}-{wildcards}
    --output=logs/{rule}/{rule}-{wildcards}-%j.out
    --time={resources.time}
    --parsable
```

The first line tells the cluster to create a directory that will capture any messages emitted by the cluster (e.g., information about what went wrong if a rule fails). Next is the command that is run to request a job to run each rule in. Rules are what each step of a Snakemake pipeline are called. The lines in this example are:

* `--partition=day`: name of the partition to submit your job to. On McCleary, this will usually be day. Edit this to request submission to a different partition.
* `--cpus-per-task={threads}`: number of cpus to request. `{threads}` is a Snakemake-specified parameter that will detect the number of cores provided to the pipeline. DON'T EDIT THIS!
* `--mem={resources.mem_mb}`: amount of RAM to request. `{resources.mem_mb}` will be defined later in the file. DON'T EDIT THIS!
* `--job-name=smk-{rule}-{wildcards}`: name to be given to requested job. `{rule}` is the name of the Snakemake rule the job is being requested for and `{wildcards}` is another Snakemake-specific parameter that will distinguish different instances of the same rule. For example, `{wildcards}` might be the sample ID for the fastq files being aligned, since each sample is run in parallel whenever possible. DON'T EDIT THIS!
* `--output=logs/{rule}/{rule}-{wildcards}-%j.out`: path to where messages from the cluster will be sent. DON'T EDIT THIS!
* `--time={resources.time}`: amount of runtime to request. `{resources.time}` will be defined later in the file. DON'T EDIT THIS!
* `--parsable`: allows the `status-sacct-robust.sh` to do what it does to track the status of the submitted job.

The next couple lines set the default resources to request:

``` yaml
default-resources:
  - mem_mb=70000
  - threads=20
  - time="4:00:00"
```

This will set default values for `{resources.mem_mb}`, `{threads}`, and `{resources.time}`. So unless the pipeline specifies otherwise, 70 GB of RAM, 20 cpus, and 4 hours of runtime will be requested. Some rules in fastq2EZbakR set the maximum amount of these resources you can request. For example, not all tools can use more than one cpu, so the rules implementing those tools will set the max value for `{threads}` to 1. For rules with the max number of threads allowed set, the number of cpus requested will be the mininmum of the max threads allowed and the number of cpus set in `default-resources`.

Finally, a number of additional settings are specified. These are optional parameters that are effectively passed to the call to `snakemake` in `run_slurm.sh`:

``` yaml
restart-times: 0
max-jobs-per-second: 0.055
latency-wait: 60
jobs: 200
keep-going: True
rerun-incomplete: True
printshellcmds: True
use-conda: True
cluster-status: status-sacct-robust.sh
```

All available options are documented [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html). The settings in this example file are:

* `restart-times`: Number of times to try to rerun a job that failed. Weird things can happen on clusters that may cause a job to fail, and rerunning will often resolve these problems. A far more common source of pipeline failures though is user input error, which a job restart can never solve. Therefore, I like to keep this at 0, meaning once a job fails, the pipeline will eventually. A fantastic perk of Snakemake is that once you have diagnosed the problem, you can restart the pipeline by running it exactly as described above. Snakemake won't rerun steps for which the output has been successfully generated, and will thus pick up right where you left off! 
* `max-jobs-per-second`: Max number of jobs that can be submitted by one user per second. I got this from Yale HPC documentation, so don't edit this!
* `latency-wait`: Wait given seconds if an output file of a job is not present after the job finished. This helps prevent pipeline crashes that can sometimes arise due to latency
* `jobs`: Max number of jobs to allow running in parallel
* `keep-going`: If True, than if one job fails, don't immediately kill the pipeline. Instead, let all jobs run that don't require the failed job's intended output.
* `rerun-incomplete`: If True, then if some of the output of a rule is recognized as "incomplete", rerun the rule.
* `printshellcmds`: If True, print out the shell commands that will be executed.
* `use-conda`: If True, use conda to automatically install dependencies when first running the pipeline.
* `cluster-status`: specifies path to script to keep an eye on the status of every job requested. 

I won't discuss the details of `status-sacct-robust.sh`, as you will likely never need to mess with it. More information about it can be found from where I got it and the framework for this entire profile, which is [this awesome repository](https://github.com/jdblischak/smk-simple-slurm) on simple Snakemake profiles for slurm.