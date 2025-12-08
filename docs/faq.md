## fastq2EZbakR FAQs


### Why are there a high number of T-to-C mutations in -s<sup>4</sup>U samples?

If you have high background mutation rates (e.g., lots of mutations in samples not exposed to label), the most likely culprit is improper trimming of reads. That is, technical sequences (adapters, UMIs, etc.) or high error rates known to commonly plague the ends of reads are to blame. See my advice for read trimming settings [here](tips.md) for details. In short, make sure you have explicitly specified technical sequences that need to be removed, and it is also worth adding a bit of hard clipping to remove common end-of-read artifacts.


### Why are there a low number of T-to-C mutations in the cB?

If there are very few T-to-C mutations in the final cB.csv file (e.g., if sample-wide mutation rates in +s<sup>4</sup>U samples are < 0.003), then you may have used the incorrect value for the `strandedness` parameter in the config. One way to tell if this is the case is by looking at one of the +s<sup>4</sup>U sample counts.csv files in `results/counts/` and checking for an abundance of A-to-G mutations. If this is the case, flip the value of `strandedness` to the opposite of whatever you used.

Related to the first point, a good sanity check after running the pipeline is going into R and checking the raw mutation rates as such:

```r
library(data.table)

# To unzip and read cB, also need to have R.utils package installed
cB <- fread("path/to/cB.csv.gz")

# Assess sample-wide T-to-C mutation rate in each sample
cB[,.(mutrate = sum(TC*n)/sum(nT*n), by = sample]
  # Want to see that +s4U samples has higher mutation rate than -s4U samples
```

Similarly, checking a counts.csv file for an abundance of A-to-G mutations can be done as follows:

```r
library(data.table)

counts <- fread("path/to/+s4U/counts.csv.gz")

## Check if A-to-G mutation rate is higher than T-to-C mutation rate:

# A-to-G mutation rate
sum(counts$AG)/sum(counts$nA)

# T-to-C mutation rate
sum(counts$TC)/sum(counts$nT)
```

### Lots of "__no_feature" in cB feature columns

If a large number of reads in the cB.csv.gz file have a feature (i.e., the annotated feature the read was assigned to) value of "__no_feature", this likely has a similar cause to the low T-to-C mutation problem discussed above. That is, the strandedness of the library may be improperly defined in the config. Thus, check that the `strandedness` parameter in the config file is properly set, and confirm the strandedness of your library.

