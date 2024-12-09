## Release Notes

### Version 0.3.0

New types of output can now be generated, specifically designed to facilitate analyses of larger datasets.

* In the config file, you will now see a parameter called "final_output", with three sub-parameters, "cB", "cUP", and "arrow". See [output.md] for more details as to how these differ and what combination can be set to `True`.
* A bug was also fixed where additional, optional fastp parameters were not getting passed to fastp.

### Version 0.2.0

Adds a handful of quality of life improvements meant to reduce the complexity of pipeline configuration. Namely:

* Whether or not the user has provided BAM or FASTQ files is auto-detected by looking for .bam at the end of all strings specified under `samples` in the config file.
* Old or overly experimental functionality has been removed
* TEC assignment flag has been moved from a unique parameter sub-heading ("strategies") to a new sub-heading of the "features" parameter

### Version 0.1.0

First official release of fastq2EZbakR. Also represents version used for EZbakR suite preprint analyses.