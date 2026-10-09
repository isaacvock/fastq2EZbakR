"""Synthetic fixtures that execute the production clustering/support rules."""

config.setdefault("annotation", "unused.gtf")
config.setdefault("genome", "unused.fa")
config.setdefault("strandedness", "forward")
SAMP_NAMES = ["fixture"]
NUM_SAMPS = 1


def fetch_informative_read(wildcards):
    return "unused.bam"


rule all:
    input:
        "results/threeputr_tests/clusters.ok",


include: "../workflow/rules/threeputr.smk"


# Feed synthetic pooled coverage into the real downstream rules.
ruleorder: write_cluster_fixture > merge_3pend_bg


rule write_cluster_fixture:
    output:
        "results/merge_3pend_bg/merged_3pend_{strand}.bg",
    run:
        fixtures = {
            "weighted": [(100, 105, 3), (110, 112, 4)],
            "split": [(100, 102, 3), (102, 105, 3), (110, 112, 4)],
            "gaps": [(100, 101, 1), (201, 202, 1), (303, 304, 1)],
        }
        with open(output[0], "w") as handle:
            for start, end, depth in fixtures[wildcards.strand]:
                handle.write(f"chr1\t{start}\t{end}\t{depth}\n")


rule check_cluster_fixtures:
    input:
        expand(
            "results/summarise_PAS_clusters/summarise_clusters_{strand}.bg",
            strand=["weighted", "split", "gaps"],
        ),
    output:
        "results/threeputr_tests/clusters.ok",
    params:
        distance=config.get("bedtools_cluster_distance", 100),
    run:
        def read_rows(path):
            with open(path) as handle:
                return [line.rstrip().split("\t") for line in handle]

        weighted, split, gaps = [read_rows(path) for path in input]
        assert weighted == [["chr1", "1", "100", "112", "23"]], weighted
        assert split == weighted, split
        expected = (
            [["chr1", "1", "100", "202", "2"], ["chr1", "2", "303", "304", "1"]]
            if params.distance == 100
            else [
                ["chr1", "1", "100", "101", "1"],
                ["chr1", "2", "201", "202", "1"],
                ["chr1", "3", "303", "304", "1"],
            ]
        )
        assert params.distance in (20, 100), params.distance
        assert gaps == expected, gaps
        with open(output[0], "w") as handle:
            handle.write(
                "Weighted support, interval splitting and clustering gaps pass.\n"
            )
