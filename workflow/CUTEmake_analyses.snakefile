#######################################
## CUTEmake: Snakefile for analyses  ##
#######################################

import os
import sys
from typing import List
import pathlib
import re
import math
import numpy
import pandas as pd
from itertools import chain


# ========================================================================================
#                                    HELPER FUNCTIONS
# ========================================================================================

def str2bool(value) -> bool:
    """
    Reads a configuration value as a boolean without going through eval().
    YAML already converts True/true/yes to a python bool, this only covers the cases in
    which the value reaches us as a string.
    """
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ["true", "t", "yes", "y", "1"]


def constraint_to(values: List[str]) -> str:
    """
    From a list, return a regular expression allowing each value and not other.
    ex: ["a", "b", "v"] -> (a|b|v)
    Values are escaped: sample names containing dots or brackets would otherwise be read as
    regular expressions and match more than themselves.
    """
    if isinstance(values, str):
        raise ValueError("constraint_to(): Expected a list, got str instead")
    if len(values) == 0:
        raise ValueError("constraint_to(): Expected a non-empty list")

    return "({})".format("|".join([re.escape(str(i)) for i in values]))


def split_cores(n_jobs, reserve = 0) -> int:
    """
    Splits the available cores over the number of jobs a rule will spawn, keeping at least
    one core per job.
    """
    return max(math.floor((workflow.cores - reserve) / max(int(n_jobs), 1)), 1)


def shifted_bam(sample):
    return os.path.join("01_BAM_filtered/shifted", ''.join([str(sample), "_mapq", MAPQ, "_woMT_shifted_sorted.bam"]))


def shifted_bai(sample):
    return ''.join([shifted_bam(sample), ".bai"])


def get_control_files(wildcards):
    """
    Control bam and its index, or nothing at all when the sample carries its own name in the
    control column of the sample table.
    """
    control = control_of[wildcards.TARGET]
    if control == wildcards.TARGET:
        return []
    return [shifted_bam(control), shifted_bai(control)]


def get_control_option(wildcards, flag):
    control = control_of[wildcards.TARGET]
    if control == wildcards.TARGET:
        return ""
    return ' '.join([flag, shifted_bam(control)])


def get_macs_ratio_table(wildcards):
    """
    The MACS3 --ratio table is required only for the targets that have a control and only
    when greenlist normalization feeds the peak calling.
    """
    if greenlist_for_peak_calling and (control_of[wildcards.TARGET] != wildcards.TARGET):
        return [MACS_RATIO_TABLE]
    return []


# ========================================================================================
#                                  GENERAL VARIABLES
# ========================================================================================

### working directory
home_dir = os.path.join(config["workflow_configuration"]["output_directory"], "")
shell('mkdir -p {home_dir}')
workdir: home_dir


### Directories
ANNOTATIONDIR = "00_annotations/"
MACSDIR = "04_MACS3_peaks/"
GOPEAKSDIR = "05_GoPeaks_peaks/"
SUMMARYDIR = "06_Overall_quality_and_info/"
GREENLISTDIR = os.path.join(SUMMARYDIR, "greenlist_normalization/")


### Global parameters
MAPQ = str(config["bam_features"]["MAPQ_threshold"])
genome_fasta = str(config["genomic_annotations"]["genome_fasta"])
BLACKLIST = str(config["genomic_annotations"]["blacklist"])
BIN_SIZE = str(config["normalization"]["bigWig_binSize"])

# samtools faidx writes <file>.fai next to the file it indexes, keeping the compression
# extension when the fasta is bgzip-compressed
genome_fai = ''.join([genome_fasta, ".fai"])

# harmonised annotation files, see rule harmonize_blacklist
BLACKLIST_HARMONIZED = os.path.join(ANNOTATIONDIR, "blacklist_harmonized.bed")
GREENLIST_HARMONIZED = os.path.join(ANNOTATIONDIR, "greenlist_harmonized.bed")
MACS_RATIO_TABLE = os.path.join(GREENLISTDIR, "MACS3_ratio_per_target.tsv")
GREENLIST_SIZE_FACTORS = os.path.join(GREENLISTDIR, "all_groups_greenlist_sizeFactors.tsv")


### Read filtering flags
if str2bool(config["bam_features"]["keep_only_proper_pairs"]):
    # -f 2 drops the orphan mates left behind by the MAPQ filter: they are still flagged as
    # paired but their mate is gone, and every paired-end aware tool downstream discards
    # them anyway
    proper_pair_flag = "-f 2"
else:
    proper_pair_flag = ""


# ========================================================================================
#                                    SAMPLE METADATA
# ========================================================================================

sample_metadata = pd.read_csv(str(config["workflow_configuration"]["sample_config_table"]), sep='\t+', engine='python')

# the fourth column is optional: a table without it behaves the way the workflow did before
# normalization groups existed, with every sample estimated against every other one
if sample_metadata.shape[1] >= 4:
    sample_metadata = sample_metadata.iloc[:, 0:4].set_axis(['target_id', 'control_id', 'broad', 'normalization_group'], axis=1, copy=False)
    sample_metadata['normalization_group'] = sample_metadata['normalization_group'].astype(str).str.strip()
    sample_metadata.loc[sample_metadata['normalization_group'].isin(['', 'nan', 'NA', 'None', '-']), 'normalization_group'] = "all_samples"
else:
    sample_metadata = sample_metadata.iloc[:, 0:3].set_axis(['target_id', 'control_id', 'broad'], axis=1, copy=False)
    sample_metadata['normalization_group'] = "all_samples"

sample_metadata['target_id'] = sample_metadata['target_id'].astype(str).str.strip()
sample_metadata['control_id'] = sample_metadata['control_id'].astype(str).str.strip()

if sample_metadata['target_id'].duplicated().any():
    duplicated = ', '.join(sorted(set(sample_metadata.loc[sample_metadata['target_id'].duplicated(), 'target_id'])))
    sys.exit(''.join(["\033[1;31m\n!!! Duplicated target_id in the sample table: ", duplicated, " !!!\n\n\033[0m"]))

if sample_metadata['control_id'].isin(['', 'nan', 'NA', 'None', '-']).any():
    sys.exit("\033[1;31m\n!!! Empty control_id in the sample table !!!\n"
             "    A sample without a control must carry its own name in the control_id column.\n\n\033[0m")

TARGETNAMES = list(numpy.unique(list(sample_metadata.target_id)))
INPUTNAMES = list(numpy.unique(list(sample_metadata.control_id)))
SAMPLENAMES = list(numpy.unique(TARGETNAMES + INPUTNAMES))

# peak calling mode and control assignment are known from the sample table, so they are
# resolved once here instead of being grepped out of the table inside every shell block
peak_mode = {t: ("broad" if str2bool(b) else "narrow") for t, b in zip(sample_metadata.target_id, sample_metadata.broad)}
control_of = dict(zip(sample_metadata.target_id, sample_metadata.control_id))

peak_extension = {t: ("broadPeak" if m == "broad" else "narrowPeak") for t, m in peak_mode.items()}
macs_mode_option = {t: ("--broad" if m == "broad" else "--call-summits") for t, m in peak_mode.items()}
gopeaks_mode_option = {t: ("--broad" if m == "broad" else "") for t, m in peak_mode.items()}
gopeaks_mdist = {t: ("3000" if m == "broad" else "1000") for t, m in peak_mode.items()}


### normalization groups -------------------------------------------------------------------
# Greenlist size factors are estimated by comparing every sample to the others, so the samples
# put together have to share a comparable background. Antibodies that cover very different
# fractions of the genome do not, and requiring every sample of a large cohort to carry signal
# in the same region shrinks the usable part of the greenlist very quickly. Grouping by
# antibody addresses both, at the price of factors that are comparable within a group and not
# across groups.
group_of = dict(zip(sample_metadata.target_id, sample_metadata.normalization_group))

invalid_groups = sorted({str(g) for g in group_of.values() if not re.fullmatch(r"[A-Za-z0-9._+-]+", str(g))})
if len(invalid_groups) > 0:
    sys.exit(''.join(["\033[1;31m\n!!! Invalid normalization_group name(s): ", ', '.join(invalid_groups), " !!!\n",
                      "    Group names are used as directory names: only letters, digits, '.', '_', '+' and '-'.\n\n\033[0m"]))

# a control belongs to the group of the targets it serves. MACS3 --ratio divides the size
# factor of a target by the one of its control, and two factors coming from different groups
# are not on a common scale, so the assignment has to be unambiguous.
for target, control in zip(sample_metadata.target_id, sample_metadata.control_id):
    if control == target:
        continue
    if control in group_of:
        if group_of[control] != group_of[target]:
            sys.exit(''.join(["\033[1;31m\n!!! '", str(control), "' is the control of '", str(target),
                              "' but the two sit in different normalization groups ('", str(group_of[control]),
                              "' and '", str(group_of[target]), "') !!!\n",
                              "    A target and its control have to be normalized together. Either put both\n",
                              "    targets in the same group, or give each group its own control.\n\n\033[0m"]))
    else:
        group_of[control] = group_of[target]

GROUPS = sorted(set(group_of.values()))
samples_in_group = {g: sorted([s for s in SAMPLENAMES if group_of[s] == g]) for g in GROUPS}


### generation of global wildcard_constraints
wildcard_constraints:
    SAMPLES = constraint_to(SAMPLENAMES),
    TARGET = constraint_to(TARGETNAMES),
    GROUP = constraint_to(GROUPS)


# ========================================================================================
#                               GREENLIST CONFIGURATION
# ========================================================================================
# Greenlist normalization, de Mello et al. 2024 (Brief Bioinform 25:bbad538).
# The bed files are not shipped with CUTEmake, they live at
# https://github.com/fndemello/CUT-RUN_greenlist

greenlist_config = config["normalization"].get("greenlist", {})

use_greenlist = str2bool(greenlist_config.get("use_greenlist", False))
GREENLIST_BED = str(greenlist_config.get("greenlist_bed", "")).strip()
greenlist_counting_method = str(greenlist_config.get("counting_method", "multiBamSummary")).strip()
greenlist_size_factor_method = str(greenlist_config.get("size_factor_method", "median_of_ratios")).strip()
greenlist_min_shared_regions = int(greenlist_config.get("min_shared_regions", 200))
greenlist_min_samples = int(greenlist_config.get("min_samples", 3))
greenlist_rescale = str(greenlist_config.get("rescale_factors", "geometric_mean")).strip()
greenlist_for_peak_calling = use_greenlist and str2bool(greenlist_config.get("use_for_peak_calling", True))

if use_greenlist:
    if len(GREENLIST_BED) == 0 or not os.path.exists(GREENLIST_BED):
        sys.exit("\033[1;31m\n!!! *greenlist_bed* is empty or does not exist !!!\n"
                 "    Download the bed matching your assay and genome build from\n"
                 "    https://github.com/fndemello/CUT-RUN_greenlist and set the path in the config.\n\n\033[0m")

    if greenlist_counting_method not in ["multiBamSummary", "featureCounts"]:
        sys.exit("\033[1;31m\n!!! *counting_method* must be 'multiBamSummary' or 'featureCounts' !!!\n\n\033[0m")

    if greenlist_size_factor_method not in ["median_of_ratios", "poscounts", "total_counts"]:
        sys.exit("\033[1;31m\n!!! *size_factor_method* must be 'median_of_ratios', 'poscounts' or 'total_counts' !!!\n\n\033[0m")

# size factors are estimated within a group, so the minimum sample number applies per group
# and not to the cohort as a whole. Groups below the minimum are skipped individually, the
# others still run.
greenlist_groups_run = [g for g in GROUPS if len(samples_in_group[g]) >= greenlist_min_samples] if use_greenlist else []
greenlist_groups_skipped = [g for g in GROUPS if g not in greenlist_groups_run] if use_greenlist else []

run_greenlist = use_greenlist and (len(greenlist_groups_run) > 0)
greenlist_for_peak_calling = greenlist_for_peak_calling and run_greenlist


# ========================================================================================
#                                   OUTPUT COLLECTIONS
# ========================================================================================

### Correlations and heatmaps outputs
correlation_outputs = []
if (len(SAMPLENAMES) > 1):
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA1.2_on_BigWigs_wholeGenome.pdf"))
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA2.3_on_BigWigs_wholeGenome.pdf"))
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf"))
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf"))
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf"))
    correlation_outputs.append(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf"))


### Greenlist outputs
greenlist_outputs = []
if run_greenlist:
    greenlist_outputs.append(os.path.join(GREENLISTDIR, "greenlist_blacklist_overlap.txt"))
    greenlist_outputs.append(GREENLIST_SIZE_FACTORS)

    for greenlist_group in greenlist_groups_run:
        greenlist_outputs.append(os.path.join(GREENLISTDIR, greenlist_group, "greenlist_sizeFactors.tsv"))
        greenlist_outputs.append(os.path.join(GREENLISTDIR, greenlist_group, "greenlist_sizeFactors_barplot.pdf"))
        greenlist_outputs.append(os.path.join(GREENLISTDIR, greenlist_group, "greenlist_background_consistency_boxplot.pdf"))
        greenlist_outputs += [os.path.join("03_bigWig_bamCoverage/greenlist_normalized/",
                                           ''.join([greenlist_sample, "_mapq", MAPQ, "_greenlist.normalized_bs", BIN_SIZE, ".bw"]))
                              for greenlist_sample in samples_in_group[greenlist_group]]

if use_greenlist:
    for greenlist_group in greenlist_groups_skipped:
        greenlist_outputs.append(os.path.join(GREENLISTDIR, greenlist_group, "not_generated_due_to_low_sample_number.txt"))


### Chromosome remove pattern
if (len(str(config["bam_features"]["remove_other_chromosomes_pattern"])) > 0):
    chr_remove_pattern = '^chrM|^M|' + str(config["bam_features"]["remove_other_chromosomes_pattern"])
else:
    chr_remove_pattern = '^chrM|^M'


# ========================================================================================
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ========================================================================================

# ========================================================================================
# Function to run all functions
rule AAA_initialization:
    input:
        BAM_shifted_sorted = expand(os.path.join("01_BAM_filtered/shifted", ''.join(["{sample}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])), sample = SAMPLENAMES),
        multiQC_BAM_html = os.path.join(SUMMARYDIR, "multiQC_bams/multiQC_report_BAMs.html"),
        report_pdf = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf"),
        report_ggplot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots_ggplot.version.pdf"),
        raw_bigWig = expand(os.path.join("03_bigWig_bamCoverage/raw_coverage/", ''.join(["{sample}_mapq", MAPQ, "_raw.coverage_bs", BIN_SIZE, ".bw"])), sample = SAMPLENAMES),
        normalized_bigWig = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", MAPQ, "_RPGC.normalized_bs", BIN_SIZE, ".bw"])), sample = SAMPLENAMES),
        macs_peaks = expand(os.path.join(MACSDIR, "{target}.filtered.BAMPE_peaks.blacklist.filtered.bed"), target = TARGETNAMES),
        gopeaks_peaks = expand(os.path.join(GOPEAKSDIR, "{target}.filtered.BAMPE_peaks.blacklist.filtered.bed"), target = TARGETNAMES),
        correlation_outputs = correlation_outputs,
        greenlist_outputs = greenlist_outputs,
        lorenz_plot_combined = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples_combined.pdf"),
        lorenz_plot_ggplot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf"),
        aggregated_qc_macs = os.path.join(SUMMARYDIR, "peaks_stats/all_samples_FRiP_report_MACS.peaks.tsv"),
        aggregated_qc_gopeaks = os.path.join(SUMMARYDIR, "peaks_stats/all_samples_FRiP_report_GoPeaks.peaks.tsv")
    shell:
        """
        printf '\033[1;36mPipeline ended!\\n\033[0m'
        """
# ========================================================================================


# ----------------------------------------------------------------------------------------
# Genome index
if not os.path.exists(genome_fai):
    rule generate_genome_index:
        input:
            genome = ancient(genome_fasta)
        output:
            genome_fai = ancient(genome_fai)
        threads: 1
        benchmark:
            "benchmarks/generate_genome_index/generate_genome_index---benchmark.txt"
        shell:
            """
            printf '\033[1;36mIndexing the genome fasta...\\n\033[0m'
            $CONDA_PREFIX/bin/samtools faidx {input.genome}
            printf '\033[1;36mGenome index done.\\n\033[0m'
            """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Chromosome naming harmonisation of the annotation files
#
# The blacklists distributed by ENCODE and the greenlists from de Mello et al. use UCSC
# names, while a genome downloaded from Ensembl does not. Subtracting a 'chr1' blacklist
# from a '1' peak file removes nothing at all and no tool complains, so both annotation
# files are rewritten once into whatever naming the genome index carries and every rule
# downstream reads the harmonised copies.
rule harmonize_blacklist:
    input:
        genome_fai = ancient(genome_fai),
        blacklist = ancient(BLACKLIST)
    output:
        blacklist = BLACKLIST_HARMONIZED
    params:
        annotation_dir = ANNOTATIONDIR
    threads: 1
    benchmark:
        "benchmarks/harmonize_blacklist/harmonize_blacklist---benchmark.txt"
    shell:
        """
        printf '\033[1;36mHarmonising the blacklist to the genome chromosome naming...\\n\033[0m'

        mkdir -p {params.annotation_dir}

        if cut -f 1 {input.genome_fai} | grep -q -m 1 '^chr'; then
            awk 'BEGIN{{OFS="\\t"}} !/^#/ && NF>=3 {{if ($1 !~ /^chr/) {{if ($1=="MT") {{$1="chrM"}} else {{$1="chr"$1}}}} ; print}}' {input.blacklist} > {output.blacklist}.tmp
        else
            awk 'BEGIN{{OFS="\\t"}} !/^#/ && NF>=3 {{sub(/^chr/,"",$1); if ($1=="M") {{$1="MT"}} ; print}}' {input.blacklist} > {output.blacklist}.tmp
        fi

        # contigs absent from the genome index are dropped, they would only make bedtools noisy
        awk 'BEGIN{{OFS="\\t"}} NR==FNR {{keep[$1]; next}} ($1 in keep)' {input.genome_fai} {output.blacklist}.tmp | sort -k1,1 -k2,2n > {output.blacklist}
        rm -f {output.blacklist}.tmp

        printf '\033[1;36m    %s regions kept in the harmonised blacklist\\n\033[0m' $(wc -l < {output.blacklist})
        """


if run_greenlist:
    rule harmonize_greenlist:
        input:
            genome_fai = ancient(genome_fai),
            greenlist = ancient(GREENLIST_BED),
            blacklist = BLACKLIST_HARMONIZED
        output:
            greenlist = GREENLIST_HARMONIZED,
            overlap_report = os.path.join(GREENLISTDIR, "greenlist_blacklist_overlap.txt")
        params:
            annotation_dir = ANNOTATIONDIR,
            greenlist_dir = GREENLISTDIR,
            source_bed = GREENLIST_BED
        threads: 1
        benchmark:
            "benchmarks/harmonize_greenlist/harmonize_greenlist---benchmark.txt"
        shell:
            """
            printf '\033[1;36mHarmonising the greenlist to the genome chromosome naming...\\n\033[0m'

            mkdir -p {params.annotation_dir}
            mkdir -p {params.greenlist_dir}

            if cut -f 1 {input.genome_fai} | grep -q -m 1 '^chr'; then
                awk 'BEGIN{{OFS="\\t"}} !/^#/ && NF>=3 {{if ($1 !~ /^chr/) {{if ($1=="MT") {{$1="chrM"}} else {{$1="chr"$1}}}} ; print $1,$2,$3}}' {input.greenlist} > {output.greenlist}.tmp
            else
                awk 'BEGIN{{OFS="\\t"}} !/^#/ && NF>=3 {{sub(/^chr/,"",$1); if ($1=="M") {{$1="MT"}} ; print $1,$2,$3}}' {input.greenlist} > {output.greenlist}.tmp
            fi

            awk 'BEGIN{{OFS="\\t"}} NR==FNR {{keep[$1]; next}} ($1 in keep)' {input.genome_fai} {output.greenlist}.tmp | sort -k1,1 -k2,2n > {output.greenlist}
            rm -f {output.greenlist}.tmp

            # Greenlist regions are regions of reproducible background, so a fraction of them
            # can sit inside a blacklist. They are counted here rather than removed: dropping
            # them would silently shrink the count matrix, and the blacklist is deliberately
            # not applied to the counting step.
            N_SOURCE=$(grep -v -c '^#' {params.source_bed} || true)
            N_KEPT=$(wc -l < {output.greenlist})
            N_OVERLAP=$($CONDA_PREFIX/bin/bedtools intersect -nonamecheck -u -a {output.greenlist} -b {input.blacklist} | wc -l)

            printf 'greenlist_source\\t%s\\n' "{params.source_bed}" > {output.overlap_report}
            printf 'regions_in_source\\t%s\\n' "$N_SOURCE" >> {output.overlap_report}
            printf 'regions_on_genome_contigs\\t%s\\n' "$N_KEPT" >> {output.overlap_report}
            printf 'regions_overlapping_blacklist\\t%s\\n' "$N_OVERLAP" >> {output.overlap_report}

            if [ "$N_KEPT" -eq 0 ]; then
                printf '\033[1;31mNo greenlist region matches the genome contigs, check the genome build!\\n\033[0m'
                exit 1
            fi

            printf '\033[1;36m    %s greenlist regions kept, %s of them overlap the blacklist\\n\033[0m' "$N_KEPT" "$N_OVERLAP"
            """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# BAM filtering
if not str2bool(config["bam_features"]["skip_bam_filtering"]):
    rule MAPQ_filter:
        input:
            source_bam = os.path.join(config["workflow_configuration"]["runs_directory"], ''.join(["{SAMPLES}", config["bam_features"]["bam_suffix"]]))
        output:
            bam_mapq_only = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, ".bam"]))),
            bam_mapq_only_sorted_toFix = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_toFix.bam"]))),
            bam_mapq_only_sorted_index_toFix = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_toFix.bai"]))),
            bam_mapq_only_sorted = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted.bam"]))),
            bam_mapq_only_sorted_index = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted.bai"]))),
            keep_chromosomes = temp(os.path.join("01_BAM_filtered", ''.join(["{SAMPLES}_mapq", MAPQ, "_keptChromosomes.txt"]))),
            idxstats_file = os.path.join(SUMMARYDIR, "reads_per_chromosome/{SAMPLES}_idxstats_read_per_chromosome.txt"),
            bam_mapq_only_sorted_woMT = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam"])),
            bam_mapq_only_sorted_woMT_index = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam.bai"])),
            flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_flagstat.txt"]))
        params:
            sample = "{SAMPLES}",
            MAPQ_threshold = MAPQ,
            proper_pair_flag = proper_pair_flag,
            chr_remove_pattern = chr_remove_pattern
        threads:
            split_cores(len(SAMPLENAMES))
        log:
            fixmate_log = "01_BAM_filtered/FixMateInformation_logs/{SAMPLES}_FixMateInformation.log"
        benchmark:
            "benchmarks/MAPQ_MT_filter/MAPQ_MT_filter---{SAMPLES}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: filtering MAPQ...\\n\033[0m'

            mkdir -p 01_BAM_filtered/FixMateInformation_logs
            mkdir -p 01_BAM_filtered/flagstat

            $CONDA_PREFIX/bin/samtools view -@ {threads} -h -q {params.MAPQ_threshold} {params.proper_pair_flag} {input.source_bam} -o {output.bam_mapq_only}

            $CONDA_PREFIX/bin/samtools sort -@ {threads} {output.bam_mapq_only} -o {output.bam_mapq_only_sorted_toFix}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted_toFix} {output.bam_mapq_only_sorted_index_toFix}

            $CONDA_PREFIX/bin/gatk FixMateInformation \
            --INPUT {output.bam_mapq_only_sorted_toFix} \
            --OUTPUT {output.bam_mapq_only_sorted} \
            --ASSUME_SORTED false \
            --ADD_MATE_CIGAR true \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY LENIENT &> {log.fixmate_log}

            $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} > {output.idxstats_file}

            printf '\033[1;36m{params.sample}: Removing MT from BAM...\\n\033[0m'

            # the list of contigs to keep is written to a file and expanded in a single
            # samtools call: piping it through xargs splits the command line on fragmented
            # assemblies and concatenates several BAMs, each with its own header, into one
            # truncated file
            $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} | cut -f 1 | grep -v -E '{params.chr_remove_pattern}' | grep -v '^\\*$' > {output.keep_chromosomes}

            $CONDA_PREFIX/bin/samtools view -@ {threads} -b {output.bam_mapq_only_sorted} $(tr '\\n' ' ' < {output.keep_chromosomes}) > {output.bam_mapq_only_sorted_woMT}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted_woMT} {output.bam_mapq_only_sorted_woMT_index}

            $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mapq_only_sorted_woMT} > {output.flagstat_filtered}
            """
else:
    rule bam_link__skip_filtering:
        input:
            source_bam = os.path.join(config["workflow_configuration"]["runs_directory"], ''.join(["{SAMPLES}", config["bam_features"]["bam_suffix"]]))
        output:
            bam_mdup = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam"])),
            bai_mdup = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam.bai"])),
            flagstat_filtered = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_flagstat.txt"]))
        params:
            sample = "{SAMPLES}"
        threads:
            split_cores(len(SAMPLENAMES))
        benchmark:
            "benchmarks/bam_link__skip_filtering/bam_link__skip_filtering---{SAMPLES}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample} (skip filtering): linking bam, indexing and computing flagstat...\\n\033[0m'

            mkdir -p 01_BAM_filtered/flagstat

            BAM_REAL=$(realpath {input.source_bam})
            ln -s -f $BAM_REAL {output.bam_mdup}
            $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mdup} {output.bai_mdup}

            $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mdup} > {output.flagstat_filtered}
            """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Generate ATAC-like shifted bams
rule bam_shifting:
    input:
        BAM = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam"])),
        BAM_index = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam.bai"]))
    output:
        BAM_shifted_toSort = temp(os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_shifted.ToSort.bam"]))),
        BAM_shifted_sorted = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        BAM_shifted_sorted_index = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        BAM_shifted_sorted_flagstat = os.path.join("01_BAM_filtered/flagstat/", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted_flagstat.txt"]))
    params:
        sample = "{SAMPLES}",
        minFragmentLength = str(config["read_shifting"]["minFragmentLength"]),
        maxFragmentLength = str(config["read_shifting"]["maxFragmentLength"])
    threads:
        workflow.cores
    log:
        out = "01_BAM_filtered/shifted/logs/{SAMPLES}_alignmentSieve.log"
    benchmark:
        "benchmarks/bam_shifting/bam_shifting---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Shifting reads in BAM...\\n\033[0m'

        mkdir -p 01_BAM_filtered/shifted/logs

        $CONDA_PREFIX/bin/alignmentSieve \
        -p {threads} \
        --ATACshift \
        --bam {input.BAM} \
        --outFile {output.BAM_shifted_toSort} \
        --minFragmentLength {params.minFragmentLength} \
        --maxFragmentLength {params.maxFragmentLength} > {log.out} 2>&1

        printf '\033[1;36m{params.sample}: Sorting shifted BAM...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools sort -@ {threads} {output.BAM_shifted_toSort} -o {output.BAM_shifted_sorted}
        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.BAM_shifted_sorted} {output.BAM_shifted_sorted_index}

        printf '\033[1;36m{params.sample}: Getting flagstat from shifted BAM...\\n\033[0m'
        $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.BAM_shifted_sorted} > {output.BAM_shifted_sorted_flagstat}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# FastQC on BAMs
rule fastQC_BAMs:
    input:
        BAM = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam"])),
        BAM_index = os.path.join("01_BAM_filtered/unshifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT.bam.bai"]))
    output:
        html = os.path.join("02_BAM_fastQC/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_fastqc.html"])),
        zip = os.path.join("02_BAM_fastQC/", ''.join(["{SAMPLES}_mapq", MAPQ, "_sorted_woMT_fastqc.zip"]))
    params:
        fastQC_BAMs_outdir = "02_BAM_fastQC/",
        sample = "{SAMPLES}"
    threads:
        split_cores(len(SAMPLENAMES))
    benchmark:
        "benchmarks/fastQC_BAMs/fastQC_BAMs---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: Performing fastQC on filtered bam...\\n\033[0m'

        mkdir -p {params.fastQC_BAMs_outdir}

        $CONDA_PREFIX/bin/fastqc -t {threads} --outdir {params.fastQC_BAMs_outdir} {input.BAM}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform multiQC for BAMs
rule multiQC_BAMs:
    input:
        BAM_fastqc_zip = expand(os.path.join("02_BAM_fastQC/", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_fastqc.zip"])), sample = SAMPLENAMES),
        macs_peaks = expand(os.path.join(MACSDIR, "{target}.filtered.BAMPE_peaks.xls"), target = TARGETNAMES)
    output:
        multiqcReportBAM = os.path.join(SUMMARYDIR, "multiQC_bams/multiQC_report_BAMs.html")
    params:
        fastQC_BAM_reports_dir = "02_BAM_fastQC/",
        BAM_flagstat_dir = "01_BAM_filtered/flagstat/",
        macs_dir = MACSDIR,
        multiQC_BAM_outdir = os.path.join(SUMMARYDIR, "multiQC_bams/")
    log:
        out = os.path.join(SUMMARYDIR, "multiQC_bams/multiQC_report_BAMs.out"),
        err = os.path.join(SUMMARYDIR, "multiQC_bams/multiQC_report_BAMs.err")
    threads: 1
    benchmark:
        "benchmarks/multiQC_BAMs/multiQC_BAMs---benchmark.txt"
    shell:
        """
        printf '\033[1;36mGenerating multiQC report from BAM quality tests...\\n\033[0m'

        mkdir -p {params.multiQC_BAM_outdir}

        # GoPeaks has no multiqc module, only the MACS3 xls files are parsed here
        $CONDA_PREFIX/bin/multiqc -f \
        --outdir {params.multiQC_BAM_outdir} \
        -n multiQC_report_BAMs.html \
        --dirs \
        {params.fastQC_BAM_reports_dir} \
        {params.BAM_flagstat_dir} \
        {params.macs_dir} > {log.out} 2> {log.err}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Fragment size distribution
rule fragment_size_distribution:
    input:
        BAM = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        BAM_index = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        blacklist = BLACKLIST_HARMONIZED
    output:
        plot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/{SAMPLES}_fragment_size_distribution.pdf"),
        fragmentSize_RawFragmentLengths = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize/{SAMPLES}_fragmentSize_RawFragmentLengths.txt"),
        fragmentSize_metrics = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize/{SAMPLES}_fragmentSize_metrics.txt")
    params:
        build_summary_directory = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log"),
        build_summary_directory_table = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize"),
        sample = "{SAMPLES}",
        plotFormat = "pdf",
        binSize = str(config["quality_controls"]["fragmentSize_window_length"]),
        maxFragmentLength = config["bam_features"]["maxFragmentLength"]
    log:
        out = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/{SAMPLES}_fragmentSize_log.out"),
        err = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/{SAMPLES}_fragmentSize_log.err")
    threads:
        split_cores(len(SAMPLENAMES))
    benchmark:
        "benchmarks/fragment_size_distribution/fragment_size_distribution---{SAMPLES}_benchmark.txt"
    priority: -1
    shell:
        """
        printf '\033[1;36m{params.sample}: Plotting the fragment size distribution...\\n\033[0m'

        mkdir -p {params.build_summary_directory}
        mkdir -p {params.build_summary_directory_table}

        $CONDA_PREFIX/bin/bamPEFragmentSize \
        -p {threads} \
        -b {input.BAM} \
        --plotFileFormat {params.plotFormat} \
        --plotTitle {params.sample} \
        --samplesLabel {params.sample} \
        --binSize {params.binSize} \
        --maxFragmentLength {params.maxFragmentLength} \
        --blackListFileName {input.blacklist} \
        --outRawFragmentLengths {output.fragmentSize_RawFragmentLengths} \
        --table {output.fragmentSize_metrics} \
        --histogram {output.plot} > {log.out} 2> {log.err}
        """


rule fragment_size_distribution_report:
    input:
        plots = expand(os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/{sample}_fragment_size_distribution.pdf"), sample = SAMPLENAMES),
        tables = expand(os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/table_and_fragmentSize/{sample}_fragmentSize_RawFragmentLengths.txt"), sample = SAMPLENAMES)
    output:
        table_list = temp(os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/fragmentSize_table_list.txt")),
        replot_script = temp(os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/replot_script.R")),
        report_pdf = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf"),
        report_ggplot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots_ggplot.version.pdf")
    params:
        maxFragmentLength = config["bam_features"]["maxFragmentLength"]
    threads: 1
    log:
        ggplot = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/ggplot_replotting.log"),
        pdfcombine = os.path.join(SUMMARYDIR, "fragmentSizeDistribution_plots/log/pdfcombine.log")
    benchmark:
        "benchmarks/fragment_size_distribution_report/fragment_size_distribution_report---benchmark.txt"
    shell:
        """
        printf '\033[1;36mMerging fragmentSizeDistribution reports in a unique PDF...\\n\033[0m'

        # the files come from the rule input rather than from a wildcard glob, otherwise
        # leftovers of a previous run made with a different sample set end up in the report
        $CONDA_PREFIX/bin/pdfcombine {input.plots} -o {output.report_pdf} -sf &> {log.pdfcombine}

        printf '%s\\n' {input.tables} > {output.table_list}

        printf '\033[1;36mReplotting fragmentSizeDistribution reports in R (ggplot version)...\\n\033[0m'

        cat > {output.replot_script} << 'ENDOFRSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
table_list <- args[1]
out_pdf <- args[2]
max_fragment_length <- as.numeric(args[3])

files <- readLines(table_list)
tb <- do.call(rbind, lapply(files, function(x) read.delim(x, h = TRUE, skip = 1)))
n.samples <- length(unique(tb[,3]))

plot <- ggplot2::ggplot(data = tb, ggplot2::aes(x = Size, y = Occurrences, color = Sample)) +
  ggplot2::geom_smooth(method = 'loess', formula = y ~ x, span = 0.05, show.legend = FALSE, se = FALSE, color = 'navyblue', linewidth = 0.5) +
  ggplot2::xlim(c(1, max_fragment_length)) +
  ggplot2::theme_classic() +
  ggplot2::facet_wrap(~Sample, scale = 'free', ncol = floor(sqrt(n.samples))) +
  ggplot2::theme(axis.ticks = ggplot2::element_line(color = 'black'),
                 axis.text = ggplot2::element_text(color = 'black'),
                 strip.background = ggplot2::element_blank())

pdf(file = out_pdf,
    width = floor(sqrt(n.samples)) * 2.7,
    height = ceiling(n.samples / floor(sqrt(n.samples))) * 1.5)
print(plot)
invisible(dev.off())
ENDOFRSCRIPT

        $CONDA_PREFIX/bin/Rscript --vanilla {output.replot_script} {output.table_list} {output.report_ggplot} {params.maxFragmentLength} &> {log.ggplot}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Extract chromosome sizes for the conversion to bigWig
rule extract_chromosome_sizes:
    input:
        genome_fai = ancient(genome_fai)
    output:
        chrSizes = "03_bigWig_bamCoverage/.genome_chromosome_sizes.txt"
    threads: 1
    benchmark:
        "benchmarks/extract_chromosome_sizes/extract_chromosome_sizes---benchmark.txt"
    shell:
        """
        printf '\033[1;36mExtracting chromosome sizes from genome index (.fai)...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage
        cut -f 1,2 {input.genome_fai} > {output.chrSizes}
        """
# ----------------------------------------------------------------------------------------


##########################################
###            PEAK CALLING            ###
##########################################

# ----------------------------------------------------------------------------------------
# PeakCalling on shifted bams (MACS3 callpeak)
rule peakCalling_MACS3:
    input:
        target_bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        control_files = get_control_files,
        ratio_table = get_macs_ratio_table,
        blacklist = BLACKLIST_HARMONIZED
    output:
        peaksPE = os.path.join(MACSDIR, "{TARGET}.filtered.BAMPE_peaks.xls"),
        peaks_filtered = os.path.join(MACSDIR, "{TARGET}.filtered.BAMPE_peaks.blacklist.filtered.bed"),
        peaks_filtered_chr = os.path.join(MACSDIR, "{TARGET}.filtered.BAMPE_peaks.blacklist.filtered_chr.bed")
    params:
        target = "{TARGET}",
        name = "{TARGET}.filtered.BAMPE",
        macs_version = "macs3",
        outdir = MACSDIR,
        control_option = lambda wildcards: get_control_option(wildcards, "-c"),
        mode_option = lambda wildcards: macs_mode_option[wildcards.TARGET],
        peak_extension = lambda wildcards: peak_extension[wildcards.TARGET],
        ratio_table = lambda wildcards: MACS_RATIO_TABLE if (greenlist_for_peak_calling and control_of[wildcards.TARGET] != wildcards.TARGET) else "",
        genomeSize = config["genomic_annotations"]["effective_genomeSize"],
        macs_qValue_cutoff = config["peak_calling"]["macs_qValue_cutoff"]
    threads:
        split_cores(len(TARGETNAMES))
    log:
        out = os.path.join(MACSDIR, "logs/{TARGET}_MACS3_log.out"),
        err = os.path.join(MACSDIR, "logs/{TARGET}_MACS3_log.err")
    benchmark:
        "benchmarks/peakCalling_MACS3/peakCalling_MACS3---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.target}: calling peaks ({params.macs_version})...\\n\033[0m'

        mkdir -p {params.outdir}logs

        # the control assignment and the narrow/broad mode come from the sample table, which
        # is parsed once when the workflow is built. Reading them back here with grep would
        # also match sample names contained in other sample names.
        RATIO_OPT=""
        if [ -n "{params.ratio_table}" ]; then
            RATIO=$(awk -F'\\t' -v s="{params.target}" '$1==s {{print $4}}' {params.ratio_table})
            if [ -n "$RATIO" ]; then
                RATIO_OPT="--ratio $RATIO"
                printf '\033[1;36m    greenlist scaling ratio: %s\\n\033[0m' "$RATIO"
            fi
        fi

        $CONDA_PREFIX/bin/{params.macs_version} callpeak \
        -t {input.target_bam} \
        {params.control_option} \
        -f BAMPE \
        -g {params.genomeSize} \
        -q {params.macs_qValue_cutoff} \
        --keep-dup all \
        --outdir {params.outdir} \
        --name {params.name} \
        {params.mode_option} ${{RATIO_OPT}} > {log.out} 2> {log.err}

        # the blacklist is subtracted while the peaks still carry the genome chromosome
        # naming, the renamed copy is produced afterwards for genome browsers
        $CONDA_PREFIX/bin/bedtools subtract -nonamecheck \
        -a {params.outdir}{params.name}_peaks.{params.peak_extension} \
        -b {input.blacklist} > {output.peaks_filtered}

        awk 'BEGIN{{OFS="\\t"}} {{if ($1 !~ /^chr/) {{if ($1=="MT") {{$1="chrM"}} else {{$1="chr"$1}}}} ; print}}' {output.peaks_filtered} > {output.peaks_filtered_chr}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# PeakCalling on shifted bams (GoPeaks)
rule peakCalling_GoPeaks:
    input:
        target_bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        control_files = get_control_files,
        chrSizes = "03_bigWig_bamCoverage/.genome_chromosome_sizes.txt",
        blacklist = BLACKLIST_HARMONIZED
    output:
        peaksPE = os.path.join(GOPEAKSDIR, "{TARGET}.filtered.BAMPE_peaks.bed"),
        peaks_filtered = os.path.join(GOPEAKSDIR, "{TARGET}.filtered.BAMPE_peaks.blacklist.filtered.bed"),
        peaks_filtered_chr = os.path.join(GOPEAKSDIR, "{TARGET}.filtered.BAMPE_peaks.blacklist.filtered_chr.bed")
    params:
        target = "{TARGET}",
        prefix = os.path.join(GOPEAKSDIR, "{TARGET}.filtered.BAMPE"),
        outdir = GOPEAKSDIR,
        control_option = lambda wildcards: get_control_option(wildcards, "--control"),
        mode_option = lambda wildcards: gopeaks_mode_option[wildcards.TARGET],
        mdist = lambda wildcards: gopeaks_mdist[wildcards.TARGET],
        gopeaks_pValue_cutoff = config["peak_calling"]["GoPeaks_pValue_cutoff"]
    threads:
        split_cores(len(TARGETNAMES))
    log:
        out = os.path.join(GOPEAKSDIR, "logs/{TARGET}_GoPeaks_log.out"),
        err = os.path.join(GOPEAKSDIR, "logs/{TARGET}_GoPeaks_log.err")
    benchmark:
        "benchmarks/peakCalling_GoPeaks/peakCalling_GoPeaks---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.target}: calling peaks (GoPeaks)...\\n\033[0m'

        mkdir -p {params.outdir}logs

        # GoPeaks has no equivalent of the MACS3 --ratio option, so the greenlist size
        # factors do not reach this rule
        $CONDA_PREFIX/bin/gopeaks \
        --bam {input.target_bam} \
        {params.control_option} \
        --chromsize {input.chrSizes} \
        --mdist {params.mdist} \
        --pval {params.gopeaks_pValue_cutoff} \
        --prefix {params.prefix} {params.mode_option} > {log.out} 2> {log.err}

        $CONDA_PREFIX/bin/bedtools subtract -nonamecheck \
        -a {output.peaksPE} \
        -b {input.blacklist} > {output.peaks_filtered}

        awk 'BEGIN{{OFS="\\t"}} {{if ($1 !~ /^chr/) {{if ($1=="MT") {{$1="chrM"}} else {{$1="chr"$1}}}} ; print}}' {output.peaks_filtered} > {output.peaks_filtered_chr}
        """
# ----------------------------------------------------------------------------------------


##########################################
###      Normalization (bigWigs)       ###
##########################################

rule raw_bigWig:
    input:
        bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        blacklist = BLACKLIST_HARMONIZED
    output:
        raw_bigWig = os.path.join("03_bigWig_bamCoverage/raw_coverage/", ''.join(["{SAMPLES}_mapq", MAPQ, "_raw.coverage_bs", BIN_SIZE, ".bw"]))
    params:
        sample = "{SAMPLES}",
        genomeSize = config["genomic_annotations"]["effective_genomeSize"],
        ignore_for_normalization = config["genomic_annotations"]["ignore_for_normalization"],
        bw_binSize = BIN_SIZE
    threads:
        split_cores(len(SAMPLENAMES))
    log:
        out = "03_bigWig_bamCoverage/raw_coverage/logs/{SAMPLES}_raw_bamCoverage_log.out",
        err = "03_bigWig_bamCoverage/raw_coverage/logs/{SAMPLES}_raw_bamCoverage_log.err"
    benchmark:
        "benchmarks/raw_bigWig/raw_bigWig---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating raw coverage bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/raw_coverage/logs

        $CONDA_PREFIX/bin/bamCoverage \
        -b {input.bam} \
        -o {output.raw_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing None \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {input.blacklist} \
        --extendReads \
        -p {threads} > {log.out} 2> {log.err}
        """


rule RPGC_normalized_bigWig:
    input:
        bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        blacklist = BLACKLIST_HARMONIZED
    output:
        RPGC_bigWig = os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{SAMPLES}_mapq", MAPQ, "_RPGC.normalized_bs", BIN_SIZE, ".bw"]))
    params:
        sample = "{SAMPLES}",
        genomeSize = config["genomic_annotations"]["effective_genomeSize"],
        ignore_for_normalization = config["genomic_annotations"]["ignore_for_normalization"],
        bw_binSize = BIN_SIZE
    threads:
        split_cores(len(SAMPLENAMES))
    log:
        out = "03_bigWig_bamCoverage/RPGC_normalized/logs/{SAMPLES}_RPGC.normalized_log.out",
        err = "03_bigWig_bamCoverage/RPGC_normalized/logs/{SAMPLES}_RPGC.normalized_log.err"
    benchmark:
        "benchmarks/RPGC_normalized_bigWig/RPGC_normalized_bigWig---{SAMPLES}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: generating RPGC-normalized bigWig...\\n\033[0m'

        mkdir -p 03_bigWig_bamCoverage/RPGC_normalized/logs

        $CONDA_PREFIX/bin/bamCoverage \
        -b {input.bam} \
        -o {output.RPGC_bigWig} \
        --binSize {params.bw_binSize} \
        --normalizeUsing RPGC \
        --effectiveGenomeSize {params.genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --blackListFileName {input.blacklist} \
        --extendReads \
        -p {threads} > {log.out} 2> {log.err}
        """
# ------------------------------------------------------------------------------


##########################################
###      Greenlist normalization       ###
##########################################
# de Mello et al. 2024, Brief Bioinform 25:bbad538
#
# Greenlist regions carry reproducible background signal. Their read pile-up is quantified
# for every sample, a median-of-ratios estimator turns the count matrix into one size factor
# per sample, and the reciprocal of that factor scales the coverage tracks.
#
# Size factors are estimated separately for every normalization group, the fourth column of
# the sample table. Two reasons for that: the estimator compares each sample to the others,
# which only means something when they share a comparable background, and the classic
# estimator keeps only the regions covered in every sample, so a large mixed cohort shrinks
# the usable part of the greenlist very quickly. The price is that factors are comparable
# within a group and not across groups.
#
# The greenlist tracks are produced in addition to the RPGC ones, both sets stay available.

if use_greenlist and len(greenlist_groups_skipped) > 0:
    rule greenlist_low_sample_number:
        output:
            notice = os.path.join(GREENLISTDIR, "{GROUP}/not_generated_due_to_low_sample_number.txt")
        params:
            group_dir = lambda wildcards: os.path.join(GREENLISTDIR, wildcards.GROUP),
            group = "{GROUP}",
            n_samples = lambda wildcards: len(samples_in_group[wildcards.GROUP]),
            min_samples = greenlist_min_samples,
            sample_list = lambda wildcards: ', '.join(samples_in_group[wildcards.GROUP])
        threads: 1
        shell:
            """
            printf '\033[1;33mGreenlist normalization skipped for group {params.group}: %s sample(s) available, %s required.\\n\033[0m' "{params.n_samples}" "{params.min_samples}"

            mkdir -p {params.group_dir}

            printf 'Greenlist normalization was not performed for the group "%s".\\n\\n' "{params.group}" > {output.notice}
            printf 'Samples in the group: %s\\n' "{params.sample_list}" >> {output.notice}
            printf 'Samples available: %s\\n' "{params.n_samples}" >> {output.notice}
            printf 'Samples required (normalization/greenlist/min_samples): %s\\n\\n' "{params.min_samples}" >> {output.notice}
            printf 'Greenlist size factors are estimated within a normalization group by\\n' >> {output.notice}
            printf 'comparing each sample to the geometric mean of the others in the same\\n' >> {output.notice}
            printf 'group. Below the minimum the estimate carries no information and would\\n' >> {output.notice}
            printf 'return factors close to 1 that look meaningful but are not.\\n\\n' >> {output.notice}
            printf 'Either merge this group with another one in the fourth column of the\\n' >> {output.notice}
            printf 'sample table, or lower min_samples if you know what the factors will be\\n' >> {output.notice}
            printf 'worth.\\n\\n' >> {output.notice}
            printf 'The RPGC-normalized bigWig files in 03_bigWig_bamCoverage/RPGC_normalized/\\n' >> {output.notice}
            printf 'are unaffected and can be used as usual.\\n' >> {output.notice}
            """


if run_greenlist:

    if greenlist_counting_method == "multiBamSummary":
        rule greenlist_counts:
            input:
                bams = lambda wildcards: [shifted_bam(s) for s in samples_in_group[wildcards.GROUP]],
                bais = lambda wildcards: [shifted_bai(s) for s in samples_in_group[wildcards.GROUP]],
                greenlist = GREENLIST_HARMONIZED
            output:
                npz = os.path.join(GREENLISTDIR, "{GROUP}/greenlist_counts_matrix.npz"),
                raw = temp(os.path.join(GREENLISTDIR, "{GROUP}/greenlist_raw_counts.deeptools.txt")),
                counts = os.path.join(GREENLISTDIR, "{GROUP}/greenlist_raw_counts.tsv")
            params:
                group = "{GROUP}",
                log_dir = os.path.join(GREENLISTDIR, "{GROUP}/logs"),
                labels = lambda wildcards: ' '.join(samples_in_group[wildcards.GROUP])
            threads:
                max((workflow.cores - 1), 1)
            log:
                out = os.path.join(GREENLISTDIR, "{GROUP}/logs/greenlist_multiBamSummary.out"),
                err = os.path.join(GREENLISTDIR, "{GROUP}/logs/greenlist_multiBamSummary.err")
            benchmark:
                "benchmarks/greenlist_counts/greenlist_counts---multiBamSummary---{GROUP}_benchmark.txt"
            shell:
                """
                printf '\033[1;36m{params.group}: quantifying the greenlist regions (multiBamSummary)...\\n\033[0m'

                mkdir -p {params.log_dir}

                # no --blackListFileName here on purpose: greenlist regions are regions of
                # reproducible background and some of them sit inside a blacklist, removing
                # them would silently drop rows from the count matrix
                $CONDA_PREFIX/bin/multiBamSummary BED-file \
                --BED {input.greenlist} \
                --bamfiles {input.bams} \
                --labels {params.labels} \
                --extendReads \
                --centerReads \
                -p {threads} \
                -o {output.npz} \
                --outRawCounts {output.raw} > {log.out} 2> {log.err}

                # the column order of --outRawCounts follows the order of --bamfiles, so the
                # header is rewritten from the sample list rather than parsed back
                printf 'region\\t%s\\n' "$(echo '{params.labels}' | tr ' ' '\\t')" > {output.counts}
                tail -n +2 {output.raw} | awk -F'\\t' 'BEGIN{{OFS="\\t"}} {{printf "%s:%s-%s", $1, $2, $3; for (i=4; i<=NF; i++) printf "%s%d", OFS, $i; printf "\\n"}}' >> {output.counts}
                """

    else:
        rule greenlist_counts:
            input:
                bams = lambda wildcards: [shifted_bam(s) for s in samples_in_group[wildcards.GROUP]],
                bais = lambda wildcards: [shifted_bai(s) for s in samples_in_group[wildcards.GROUP]],
                greenlist = GREENLIST_HARMONIZED
            output:
                saf = temp(os.path.join(GREENLISTDIR, "{GROUP}/greenlist_regions.saf")),
                raw = temp(os.path.join(GREENLISTDIR, "{GROUP}/greenlist_raw_counts.featureCounts.txt")),
                counts = os.path.join(GREENLISTDIR, "{GROUP}/greenlist_raw_counts.tsv")
            params:
                group = "{GROUP}",
                log_dir = os.path.join(GREENLISTDIR, "{GROUP}/logs"),
                labels = lambda wildcards: ' '.join(samples_in_group[wildcards.GROUP])
            threads:
                max((workflow.cores - 1), 1)
            log:
                out = os.path.join(GREENLISTDIR, "{GROUP}/logs/greenlist_featureCounts.out"),
                err = os.path.join(GREENLISTDIR, "{GROUP}/logs/greenlist_featureCounts.err")
            benchmark:
                "benchmarks/greenlist_counts/greenlist_counts---featureCounts---{GROUP}_benchmark.txt"
            shell:
                """
                printf '\033[1;36m{params.group}: quantifying the greenlist regions (featureCounts)...\\n\033[0m'

                mkdir -p {params.log_dir}

                # SAF is 1-based and end-inclusive, BED is 0-based and end-exclusive
                printf 'GeneID\\tChr\\tStart\\tEnd\\tStrand\\n' > {output.saf}
                awk 'BEGIN{{OFS="\\t"}} {{print $1":"$2"-"$3, $1, $2+1, $3, "."}}' {input.greenlist} >> {output.saf}

                # --countReadPairs counts fragments rather than the two mates separately,
                # -B keeps only pairs with both ends aligned and -C drops chimeric pairs
                $CONDA_PREFIX/bin/featureCounts \
                -F SAF \
                -a {output.saf} \
                -p --countReadPairs -B -C \
                -T {threads} \
                -o {output.raw} \
                {input.bams} > {log.out} 2> {log.err}

                printf 'region\\t%s\\n' "$(echo '{params.labels}' | tr ' ' '\\t')" > {output.counts}
                grep -v '^#' {output.raw} | awk -F'\\t' 'BEGIN{{OFS="\\t"}} NR>1 {{printf "%s", $1; for (i=7; i<=NF; i++) printf "%s%d", OFS, $i; printf "\\n"}}' >> {output.counts}
                """


    rule greenlist_sizeFactors:
        input:
            counts = os.path.join(GREENLISTDIR, "{GROUP}/greenlist_raw_counts.tsv"),
            flagstat = lambda wildcards: [os.path.join("01_BAM_filtered/flagstat/", ''.join([s, "_mapq", MAPQ, "_woMT_shifted_sorted_flagstat.txt"])) for s in samples_in_group[wildcards.GROUP]]
        output:
            size_factors = os.path.join(GREENLISTDIR, "{GROUP}/greenlist_sizeFactors.tsv"),
            barplot = os.path.join(GREENLISTDIR, "{GROUP}/greenlist_sizeFactors_barplot.pdf"),
            boxplot = os.path.join(GREENLISTDIR, "{GROUP}/greenlist_background_consistency_boxplot.pdf"),
            fragments = temp(os.path.join(GREENLISTDIR, "{GROUP}/greenlist_total_fragments.tsv")),
            script = temp(os.path.join(GREENLISTDIR, "{GROUP}/greenlist_sizeFactors.R"))
        params:
            group = "{GROUP}",
            log_dir = os.path.join(GREENLISTDIR, "{GROUP}/logs"),
            method = greenlist_size_factor_method,
            min_shared_regions = greenlist_min_shared_regions,
            rescale = greenlist_rescale
        threads: 1
        log:
            out = os.path.join(GREENLISTDIR, "{GROUP}/logs/greenlist_sizeFactors.log")
        benchmark:
            "benchmarks/greenlist_sizeFactors/greenlist_sizeFactors---{GROUP}_benchmark.txt"
        run:
            # The R code below is written to a file as a plain string instead of being
            # embedded in a shell block: every brace of an R script would otherwise have to
            # be doubled for the snakemake formatter, which is a reliable source of silent
            # typos.
            os.makedirs(params.log_dir, exist_ok = True)

            group_samples = samples_in_group[wildcards.GROUP]

            ### total fragments per sample, read from the flagstat of the shifted bam
            with open(output.fragments, "w") as fragments_file:
                fragments_file.write("sample\tfragments\n")
                for sample, flagstat_path in zip(group_samples, input.flagstat):
                    n_fragments = 0
                    n_total = 0
                    with open(flagstat_path) as flagstat:
                        for line in flagstat:
                            if re.search(r"\bread1\b", line):
                                n_fragments = int(line.split(" ")[0])
                            elif " in total " in line:
                                n_total = int(line.split(" ")[0])
                    if n_fragments == 0:
                        n_fragments = max(int(n_total / 2), 1)
                    fragments_file.write(''.join([str(sample), "\t", str(n_fragments), "\n"]))

            ### R script
            r_script = """
### CUTEmake, greenlist size factors -------------------------------------------------------
### de Mello et al. 2024, Brief Bioinform 25:bbad538

options(scipen = 999)

args <- commandArgs(trailingOnly = TRUE)

counts_file    <- args[1]
fragments_file <- args[2]
method         <- args[3]
min_regions    <- as.integer(args[4])
rescale_mode   <- args[5]
out_table      <- args[6]
out_barplot    <- args[7]
out_boxplot    <- args[8]

counts <- as.matrix(read.delim(counts_file, header = TRUE, row.names = 1, check.names = FALSE))
mode(counts) <- "numeric"

### the classic median-of-ratios estimator drops any region that is empty in even one
### sample, so a single shallow control can empty most of the matrix on its own
shared <- counts[apply(counts, 1, function(x) all(x > 0)), , drop = FALSE]
non_zero_per_sample <- colSums(counts > 0)

if (method == "median_of_ratios" && nrow(shared) < min_regions) {
  stop(paste0("only ", nrow(shared), " of ", nrow(counts),
              " greenlist regions carry signal in every sample, the minimum is ", min_regions, ".",
              "\n  Regions with signal, per sample: ",
              paste0(names(non_zero_per_sample), "=", non_zero_per_sample, collapse = ", "),
              "\n  When one sample is far below the others the cause is sequencing depth, and",
              "\n  size_factor_method 'poscounts' handles it. When every sample is low, check the",
              "\n  genome build and the assay of the greenlist bed before changing anything."))
}

if (method == "median_of_ratios") {
  suppressPackageStartupMessages(library(DESeq2))
  size_factor <- DESeq2::estimateSizeFactorsForMatrix(shared, type = "ratio")
} else if (method == "poscounts") {
  ### geometric means taken over the positive counts only, so a region empty in one sample
  ### still contributes for the others instead of being discarded
  suppressPackageStartupMessages(library(DESeq2))
  size_factor <- DESeq2::estimateSizeFactorsForMatrix(counts, type = "poscounts")
} else {
  total <- colSums(counts)
  size_factor <- total / exp(mean(log(total)))
}

if (any(!is.finite(size_factor)) || any(size_factor <= 0)) {
  stop("non finite or non positive size factors, check the greenlist count table")
}

### DESeq2 divides the counts by the size factor while deeptools multiplies the coverage by
### --scaleFactor, so the value handed over to bamCoverage is the reciprocal
scale_factor <- 1 / size_factor

if (rescale_mode == "geometric_mean") {
  scale_factor <- scale_factor / exp(mean(log(scale_factor)))
} else if (rescale_mode == "mean") {
  scale_factor <- scale_factor / mean(scale_factor)
}

### how tightly each sample follows the cohort background: a wide spread means the greenlist
### assumption does not hold for that sample and its factor should not be trusted.
### With a sparse matrix there are too few fully covered regions to say anything, so the
### diagnostic falls back to the regions covered in at least half of the samples and treats
### the remaining zeros as missing.
if (nrow(shared) >= 20) {
  diagnostic_counts <- shared
} else {
  diagnostic_counts <- counts[rowSums(counts > 0) >= ceiling(ncol(counts) / 2), , drop = FALSE]
  diagnostic_counts[diagnostic_counts == 0] <- NA
}

log_ratio <- log(diagnostic_counts) - rowMeans(log(diagnostic_counts), na.rm = TRUE)
log_ratio_sd <- apply(log_ratio, 2, sd, na.rm = TRUE)

fragments <- read.delim(fragments_file, header = TRUE, check.names = FALSE)
rownames(fragments) <- as.character(fragments$sample)

result <- data.frame(sample = colnames(counts),
                     size_factor = as.numeric(size_factor),
                     scale_factor = as.numeric(scale_factor),
                     greenlist_signal = as.numeric(colSums(counts)),
                     total_fragments = as.numeric(fragments[colnames(counts), "fragments"]),
                     log_ratio_sd = as.numeric(log_ratio_sd[colnames(counts)]),
                     regions_with_signal = as.numeric(non_zero_per_sample[colnames(counts)]),
                     total_regions = nrow(counts),
                     shared_regions = nrow(shared),
                     stringsAsFactors = FALSE)

result$fraction_in_greenlist <- result$greenlist_signal / result$total_fragments

result$size_factor <- round(result$size_factor, 6)
result$scale_factor <- round(result$scale_factor, 6)
result$log_ratio_sd <- round(result$log_ratio_sd, 4)
result$fraction_in_greenlist <- round(result$fraction_in_greenlist, 6)

write.table(result, file = out_table, quote = FALSE, sep = "\\t", row.names = FALSE)

### plots -----------------------------------------------------------------------------------
factor_plot <- ggplot2::ggplot(result, ggplot2::aes(x = stats::reorder(sample, scale_factor), y = scale_factor)) +
  ggplot2::geom_col(fill = "#4c9a6b") +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  ggplot2::coord_flip() +
  ggplot2::xlab(NULL) +
  ggplot2::ylab("applied scale factor (1 / size factor)") +
  ggplot2::ggtitle("Greenlist normalization factors") +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text = ggplot2::element_text(color = "black"),
                 axis.ticks = ggplot2::element_line(color = "black"))

pdf(file = out_barplot, width = 7, height = max(2.5, 0.25 * nrow(result) + 1.5))
print(factor_plot)
invisible(dev.off())

ratio_long <- data.frame(sample = rep(colnames(log_ratio), each = nrow(log_ratio)),
                         log2_ratio = as.vector(log_ratio) / log(2),
                         stringsAsFactors = FALSE)

ratio_long <- ratio_long[is.finite(ratio_long$log2_ratio), ]

consistency_plot <- ggplot2::ggplot(ratio_long, ggplot2::aes(x = sample, y = log2_ratio)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  ggplot2::geom_boxplot(fill = "grey90", outlier.size = 0.2, outlier.alpha = 0.3) +
  ggplot2::xlab(NULL) +
  ggplot2::ylab("log2 (region signal / region geometric mean)") +
  ggplot2::ggtitle("Consistency of the greenlist background") +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text = ggplot2::element_text(color = "black"),
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                 axis.ticks = ggplot2::element_line(color = "black"))

pdf(file = out_boxplot, width = max(5, 0.4 * ncol(counts) + 2), height = 5)
print(consistency_plot)
invisible(dev.off())
"""

            with open(output.script, "w") as script_file:
                script_file.write(r_script)

            shell("printf '\033[1;36m{params.group}: computing greenlist size factors ({params.method})...\\n\033[0m'")

            shell("$CONDA_PREFIX/bin/Rscript --vanilla {output.script} "
                  "{input.counts} {output.fragments} {params.method} {params.min_shared_regions} "
                  "{params.rescale} {output.size_factors} {output.barplot} {output.boxplot} "
                  "&> {log.out}")


    rule greenlist_aggregate_sizeFactors:
        input:
            size_factors = [os.path.join(GREENLISTDIR, group, "greenlist_sizeFactors.tsv") for group in greenlist_groups_run]
        output:
            all_size_factors = GREENLIST_SIZE_FACTORS,
            macs_ratio = MACS_RATIO_TABLE
        params:
            groups = greenlist_groups_run
        threads: 1
        benchmark:
            "benchmarks/greenlist_aggregate_sizeFactors/greenlist_aggregate_sizeFactors---benchmark.txt"
        run:
            shell("printf '\033[1;36mCollecting the greenlist size factors of all groups...\\n\033[0m'")

            size_factor = {}
            header_written = False

            with open(output.all_size_factors, "w") as merged_file:
                for group, table_path in zip(params.groups, input.size_factors):
                    with open(table_path) as group_table:
                        header = group_table.readline().rstrip("\n")
                        if not header_written:
                            merged_file.write(''.join(["normalization_group\t", header, "\n"]))
                            header_written = True
                        for line in group_table:
                            fields = line.rstrip("\n").split("\t")
                            size_factor[fields[0]] = float(fields[1])
                            merged_file.write(''.join([str(group), "\t", line.rstrip("\n"), "\n"]))

            ### MACS3 scaling ratios, target over control. A target and its control always
            ### belong to the same group, so the two factors are on a common scale. Targets
            ### whose group was skipped simply get no row and no --ratio option.
            with open(output.macs_ratio, "w") as ratio_file:
                ratio_file.write("target_id\tcontrol_id\tsize_factor_target\tratio\tnormalization_group\n")
                for target in TARGETNAMES:
                    control = control_of[target]
                    if control == target:
                        continue
                    if (target not in size_factor) or (control not in size_factor):
                        continue
                    ratio = size_factor[target] / size_factor[control]
                    ratio_file.write(''.join([target, "\t", control, "\t",
                                              "{:.6f}".format(size_factor[target]), "\t",
                                              "{:.6f}".format(ratio), "\t",
                                              str(group_of[target]), "\n"]))


    rule greenlist_normalized_bigWig:
        input:
            bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
            bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
            size_factors = lambda wildcards: os.path.join(GREENLISTDIR, group_of[wildcards.SAMPLES], "greenlist_sizeFactors.tsv"),
            blacklist = BLACKLIST_HARMONIZED
        output:
            greenlist_bigWig = os.path.join("03_bigWig_bamCoverage/greenlist_normalized/", ''.join(["{SAMPLES}_mapq", MAPQ, "_greenlist.normalized_bs", BIN_SIZE, ".bw"]))
        params:
            sample = "{SAMPLES}",
            group = lambda wildcards: group_of[wildcards.SAMPLES],
            genomeSize = config["genomic_annotations"]["effective_genomeSize"],
            ignore_for_normalization = config["genomic_annotations"]["ignore_for_normalization"],
            bw_binSize = BIN_SIZE
        threads:
            split_cores(len(SAMPLENAMES))
        log:
            out = "03_bigWig_bamCoverage/greenlist_normalized/logs/{SAMPLES}_greenlist.normalized_log.out",
            err = "03_bigWig_bamCoverage/greenlist_normalized/logs/{SAMPLES}_greenlist.normalized_log.err"
        benchmark:
            "benchmarks/greenlist_normalized_bigWig/greenlist_normalized_bigWig---{SAMPLES}_benchmark.txt"
        shell:
            """
            printf '\033[1;36m{params.sample}: generating greenlist-normalized bigWig (group {params.group})...\\n\033[0m'

            mkdir -p 03_bigWig_bamCoverage/greenlist_normalized/logs

            SCALE=$(awk -F'\\t' -v s="{params.sample}" '$1==s {{print $3}}' {input.size_factors})

            if [ -z "$SCALE" ]; then
                printf '\033[1;31m{params.sample}: no greenlist scale factor found in %s\\n\033[0m' "{input.size_factors}"
                exit 1
            fi

            printf '\033[1;36m    scale factor: %s\\n\033[0m' "$SCALE"

            # --normalizeUsing None: the depth correction is already carried by the greenlist
            # scale factor, deeptools would otherwise multiply the two together
            $CONDA_PREFIX/bin/bamCoverage \
            -b {input.bam} \
            -o {output.greenlist_bigWig} \
            --binSize {params.bw_binSize} \
            --normalizeUsing None \
            --scaleFactor $SCALE \
            --effectiveGenomeSize {params.genomeSize} \
            --ignoreForNormalization {params.ignore_for_normalization} \
            --blackListFileName {input.blacklist} \
            --extendReads \
            -p {threads} > {log.out} 2> {log.err}
            """
# ------------------------------------------------------------------------------


##########################################
###          Quality controls          ###
##########################################

# ----------------------------------------------------------------------------------------
# Generation of matrix scores
rule multiBigwigSummary:
    input:
        norm_bw = expand(os.path.join("03_bigWig_bamCoverage/RPGC_normalized/", ''.join(["{sample}_mapq", MAPQ, "_RPGC.normalized_bs", BIN_SIZE, ".bw"])), sample = SAMPLENAMES),
        blacklist = BLACKLIST_HARMONIZED
    output:
        matrix = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/multiBigWigSummary_matrix_allSamples.npz"),
        rawCounts = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/multiBigWigSummary_matrix_allSamples.txt")
    params:
        make_directory = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/"),
        labels = ' '.join(SAMPLENAMES),
        window = config["quality_controls"]["multiBigwigSummary_binning_window_size"]
    threads:
        max((workflow.cores - 1), 1)
    benchmark:
        "benchmarks/multiBigwigSummary/multiBigwigSummary---benchmark.txt"
    shell:
        """
        printf '\033[1;36mComparing the whole signal among samples...\\n\033[0m'

        mkdir -p {params.make_directory}

        $CONDA_PREFIX/bin/multiBigwigSummary bins \
        -p {threads} \
        -b {input.norm_bw} \
        --labels {params.labels} \
        --binSize {params.window} \
        --blackListFileName {input.blacklist} \
        -o {output.matrix} \
        --outRawCounts {output.rawCounts}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Generation of samples PCA and heatmap
rule PCA_and_samples_correlation:
    input:
        matrix = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/multiBigWigSummary_matrix_allSamples.npz")
    output:
        PCA_OneTwo = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA1.2_on_BigWigs_wholeGenome.pdf"),
        PCA_TwoThree = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA2.3_on_BigWigs_wholeGenome.pdf"),
        PCA_tab = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/PCA_on_BigWigs_wholeGenome_values.txt"),
        hetamap_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf"),
        hetamap_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf"),
        matrix_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_matrix_on_BigWigs_wholeGenome_spearmanMethod.txt"),
        matrix_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_matrix_on_BigWigs_wholeGenome_pearsonMethod.txt"),
        scatterplot_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf"),
        scatterplot_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf")
    params:
        make_directory = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/"),
        heatmap_color = config["quality_controls"]["correlation_heatmap_color"]
    threads: 1
    benchmark:
        "benchmarks/PCA_and_samples_correlation/PCA_and_samples_correlation---benchmark.txt"
    shell:
        """
        printf '\033[1;36mPlotting the correlation and variability of the whole signal among samples...\\n\033[0m'

        mkdir -p {params.make_directory}

        printf '\033[1;36m    - plotting PCAs (whole genome)...\\n\033[0m'
        $CONDA_PREFIX/bin/plotPCA -in {input.matrix} -o {output.PCA_OneTwo} -T 'PCA on BigWigs (whole genome)' --plotFileFormat 'pdf' --outFileNameData {output.PCA_tab}
        $CONDA_PREFIX/bin/plotPCA -in {input.matrix} -o {output.PCA_TwoThree} -T 'PCA on BigWigs (whole genome)' --plotFileFormat 'pdf' --PCs 2 3

        printf '\033[1;36m    - plotting Spearman correlation heatmap (whole genome)...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod spearman \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Spearman correlation of BigWigs" \
        --whatToPlot heatmap \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        --outFileCorMatrix {output.matrix_spearman} \
        -o {output.hetamap_spearman}

        printf '\033[1;36m    - plotting Pearson correlation heatmap (whole genome)...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod pearson \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson correlation of BigWigs" \
        --whatToPlot heatmap \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        --outFileCorMatrix {output.matrix_pearson} \
        -o {output.hetamap_pearson}

        printf '\033[1;36m    - plotting Spearman correlation scatterplot (whole genome)...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod spearman \
        --log1p \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Spearman correlation of BigWigs - ln values" \
        --whatToPlot scatterplot \
        --plotFileFormat 'pdf' \
        -o {output.scatterplot_spearman}

        printf '\033[1;36m    - plotting Pearson correlation scatterplot (whole genome)...\\n\033[0m'
        $CONDA_PREFIX/bin/plotCorrelation -in {input.matrix} \
        --corMethod pearson \
        --log1p \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson correlation of BigWigs - ln values" \
        --whatToPlot scatterplot \
        --plotFileFormat 'pdf' \
        -o {output.scatterplot_pearson}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Generation of Lorenz curves
rule Lorenz_curve:
    input:
        BAM_sorted = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        BAM_sorted_index = os.path.join("01_BAM_filtered/shifted", ''.join(["{SAMPLES}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        blacklist = BLACKLIST_HARMONIZED
    output:
        lorenz_plot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_plots/{SAMPLES}_Lorenz_curve_deeptools.plotFingreprint.pdf"),
        lorenz_metrics = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_metrics/{SAMPLES}_Lorenz_quality.metrics_deeptools.plotFingreprint.txt"),
        lorenz_counts = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_counts/{SAMPLES}_Lorenz_raw.counts_deeptools.plotFingreprint.txt")
    params:
        labels = "{SAMPLES}",
        binSize = config["quality_controls"]["plotFingerprint"]["binSize"],
        sampledRegions = config["quality_controls"]["plotFingerprint"]["sampledRegions"],
        extra_params = config["quality_controls"]["plotFingerprint"]["extra_parameters"]
    threads:
        max(math.floor((workflow.cores - 1) / 2), 1)
    log:
        out = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/log/{SAMPLES}_deeptools_plotFingreprint.log")
    benchmark:
        "benchmarks/Lorenz_curve/Lorenz_curve---{SAMPLES}_benchmark.txt"
    priority: -10
    shell:
        """
        printf '\033[1;36m{params.labels}: plotting Lorenz curves-Fingerprint...\\n\033[0m'

        $CONDA_PREFIX/bin/plotFingerprint \
        --bamfiles {input.BAM_sorted} \
        --plotFile {output.lorenz_plot} \
        --labels {params.labels} \
        --blackListFileName {input.blacklist} \
        --binSize {params.binSize} \
        --numberOfSamples {params.sampledRegions} \
        --outQualityMetrics {output.lorenz_metrics} \
        --outRawCounts {output.lorenz_counts} \
        -p {threads} {params.extra_params} &> {log.out}
        """


rule Lorenz_curve_merge_plots:
    input:
        lorenz_plots = expand(os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_plots/{sample}_Lorenz_curve_deeptools.plotFingreprint.pdf"), sample = SAMPLENAMES),
        lorenz_counts = expand(os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_counts/{sample}_Lorenz_raw.counts_deeptools.plotFingreprint.txt"), sample = SAMPLENAMES)
    output:
        counts_list = temp(os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_counts_list.txt")),
        replot_script = temp(os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/lorenz_replot_script.R")),
        lorenz_plot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples_combined.pdf"),
        lorenz_plot_ggplot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf")
    threads: 1
    log:
        ggplot = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/log/lorenz_replotting.log"),
        pdfcombine = os.path.join(SUMMARYDIR, "LorenzCurve_plotFingreprint/log/lorenz_pdfcombine.log")
    benchmark:
        "benchmarks/Lorenz_curve/Lorenz_curve_merge_plots---benchmark.txt"
    priority: -5
    shell:
        """
        printf '\033[1;36mCombining Lorenz curves-Fingerprint for all samples...\\n\033[0m'

        $CONDA_PREFIX/bin/pdfcombine {input.lorenz_plots} -o {output.lorenz_plot} -sf &> {log.pdfcombine}

        printf '%s\\n' {input.lorenz_counts} > {output.counts_list}

        printf '\033[1;36mMaking the combined Lorenz curves-Fingerprint plot...\\n\033[0m'

        cat > {output.replot_script} << 'ENDOFRSCRIPT'
suppressPackageStartupMessages(require(dplyr))

args <- commandArgs(trailingOnly = TRUE)
counts_list <- args[1]
out_pdf <- args[2]

tables <- readLines(counts_list)

read_one <- function(x) read.delim(x, skip = 2, h = FALSE) %>%
  dplyr::rename(counts = V1) %>%
  dplyr::mutate(sample = gsub('_Lorenz_raw[.]counts_deeptools[.]plotFingreprint[.]txt', '', basename(x))) %>%
  dplyr::arrange(counts) %>%
  dplyr::mutate(cumulative_sum = cumsum(counts), rank = (1:nrow(.)) / nrow(.)) %>%
  dplyr::mutate(cumulative_sum = cumulative_sum / max(cumulative_sum))

combined_table <- do.call(rbind, lapply(tables, read_one))

plot <- ggplot2::ggplot(data = combined_table, ggplot2::aes(x = rank, y = cumulative_sum, color = sample)) +
  ggplot2::geom_line() +
  ggplot2::ggtitle('Fingerprints (Lorenz curves) all samples') +
  ggplot2::xlim(c(0, 1)) +
  ggplot2::xlab('Normalized rank') +
  ggplot2::ylab('Fraction with reference to the bin with highest coverage') +
  ggplot2::theme_classic() +
  ggplot2::theme(axis.text = ggplot2::element_text(color = 'black'),
                 axis.ticks = ggplot2::element_line(color = 'black'))

pdf(out_pdf, width = 8, height = 6.5)
print(plot)
invisible(dev.off())
ENDOFRSCRIPT

        $CONDA_PREFIX/bin/Rscript --vanilla {output.replot_script} {output.counts_list} {output.lorenz_plot_ggplot} &> {log.ggplot}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
## Compute FRiP for MACS peaks
rule MACS_peak_QC:
    input:
        target_bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        peaks = os.path.join(MACSDIR, "{TARGET}.filtered.BAMPE_peaks.blacklist.filtered.bed")
    output:
        qc = temp(os.path.join(SUMMARYDIR, "peaks_stats/{TARGET}.filtered.BAMPE_peaks.qc_MACS.txt"))
    params:
        target = "{TARGET}",
        calling_mode = lambda wildcards: peak_mode[wildcards.TARGET],
        peaks_stats_dir = os.path.join(SUMMARYDIR, "peaks_stats"),
        genomeSize = config["genomic_annotations"]["effective_genomeSize"]
    threads:
        split_cores(len(TARGETNAMES))
    benchmark:
        "benchmarks/MACS_peak_QC/MACS_peak_QC---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.target}: computing peak stats (MACS)...\\n\033[0m'

        mkdir -p {params.peaks_stats_dir}

        # the statistics describe the blacklist-filtered peak set, which is the one the user
        # is going to work with
        peak_count=$(wc -l < {input.peaks})

        # -f 64 keeps read1 only, so the counts are fragments rather than mates
        mapped_reads=$($CONDA_PREFIX/bin/samtools view -@ {threads} -c -f 64 -F 260 {input.target_bam})
        reads_in_peaks=$($CONDA_PREFIX/bin/samtools view -@ {threads} -c -f 64 -F 260 -L {input.peaks} {input.target_bam})

        peak_len=$(awk 'BEGIN{{total=0}} {{total+=$3-$2}} END{{print total+0}}' {input.peaks})

        frip=$(awk -v a="$reads_in_peaks" -v b="$mapped_reads" 'BEGIN{{if (b>0) printf "%.3f", a/b; else printf "NA"}}')
        genomecov=$(awk -v a="$peak_len" -v b="{params.genomeSize}" 'BEGIN{{if (b>0) printf "%.5f", a/b; else printf "NA"}}')

        printf '{params.target}\\t{params.calling_mode}\\t%s\\t%s\\t%s\\n' "$peak_count" "$frip" "$genomecov" > {output.qc}
        """


rule aggregate_FRiP_MACS:
    input:
        qc = expand(os.path.join(SUMMARYDIR, "peaks_stats/{target}.filtered.BAMPE_peaks.qc_MACS.txt"), target = TARGETNAMES)
    output:
        aggregated_qc = os.path.join(SUMMARYDIR, "peaks_stats/all_samples_FRiP_report_MACS.peaks.tsv")
    threads: 1
    benchmark:
        "benchmarks/aggregate_FRiP_MACS/aggregate_FRiP_MACS---allTargets_benchmark.txt"
    shell:
        """
        printf 'sample\\tcalling_mode\\tn_peaks\\tFRiP\\tfraction_genome_coverage\\n' > {output.aggregated_qc}
        cat {input.qc} >> {output.aggregated_qc}
        """
# ------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
## Compute FRiP for GoPeaks peaks
rule GoPeaks_peak_QC:
    input:
        target_bam = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam"])),
        target_bai = os.path.join("01_BAM_filtered/shifted", ''.join(["{TARGET}_mapq", MAPQ, "_woMT_shifted_sorted.bam.bai"])),
        peaks = os.path.join(GOPEAKSDIR, "{TARGET}.filtered.BAMPE_peaks.blacklist.filtered.bed")
    output:
        qc = temp(os.path.join(SUMMARYDIR, "peaks_stats/{TARGET}.filtered.BAMPE_peaks.qc_GoPeaks.txt"))
    params:
        target = "{TARGET}",
        calling_mode = lambda wildcards: peak_mode[wildcards.TARGET],
        peaks_stats_dir = os.path.join(SUMMARYDIR, "peaks_stats"),
        genomeSize = config["genomic_annotations"]["effective_genomeSize"]
    threads:
        split_cores(len(TARGETNAMES))
    benchmark:
        "benchmarks/GoPeaks_peak_QC/GoPeaks_peak_QC---{TARGET}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.target}: computing peak stats (GoPeaks)...\\n\033[0m'

        mkdir -p {params.peaks_stats_dir}

        peak_count=$(wc -l < {input.peaks})

        mapped_reads=$($CONDA_PREFIX/bin/samtools view -@ {threads} -c -f 64 -F 260 {input.target_bam})
        reads_in_peaks=$($CONDA_PREFIX/bin/samtools view -@ {threads} -c -f 64 -F 260 -L {input.peaks} {input.target_bam})

        peak_len=$(awk 'BEGIN{{total=0}} {{total+=$3-$2}} END{{print total+0}}' {input.peaks})

        frip=$(awk -v a="$reads_in_peaks" -v b="$mapped_reads" 'BEGIN{{if (b>0) printf "%.3f", a/b; else printf "NA"}}')
        genomecov=$(awk -v a="$peak_len" -v b="{params.genomeSize}" 'BEGIN{{if (b>0) printf "%.5f", a/b; else printf "NA"}}')

        printf '{params.target}\\t{params.calling_mode}\\t%s\\t%s\\t%s\\n' "$peak_count" "$frip" "$genomecov" > {output.qc}
        """


rule aggregate_FRiP_GoPeaks:
    input:
        qc = expand(os.path.join(SUMMARYDIR, "peaks_stats/{target}.filtered.BAMPE_peaks.qc_GoPeaks.txt"), target = TARGETNAMES)
    output:
        aggregated_qc = os.path.join(SUMMARYDIR, "peaks_stats/all_samples_FRiP_report_GoPeaks.peaks.tsv")
    threads: 1
    benchmark:
        "benchmarks/aggregate_FRiP_GoPeaks/aggregate_FRiP_GoPeaks---allTargets_benchmark.txt"
    shell:
        """
        printf 'sample\\tcalling_mode\\tn_peaks\\tFRiP\\tfraction_genome_coverage\\n' > {output.aggregated_qc}
        cat {input.qc} >> {output.aggregated_qc}
        """
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#                                 END PIPELINE                                 #
# ------------------------------------------------------------------------------
