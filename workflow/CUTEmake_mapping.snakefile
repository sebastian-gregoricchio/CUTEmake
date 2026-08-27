#########################################
## CUTEmake: Snakefile for DNA mapping ##
#########################################

import os
import sys
from typing import List
import pathlib
import re
import numpy
import pandas as pd
import math
from itertools import combinations


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


# ========================================================================================
#                                  GENERAL VARIABLES
# ========================================================================================

genome_fasta = str(config["genome_fasta"])
MAPQ = str(config["MAPQ_threshold"])
fastq_suffix = str(config["fastq_suffix"])
read_suffix = list(config["read_suffix"])


### working directory
home_dir = os.path.join(config["output_directory"], "")
shell('mkdir -p {home_dir}')
workdir: home_dir


### get the unique samples names and other variables
if not (os.path.exists(config["fastq_directory"])):
    sys.exit("\033[1;31m\n!!! *fastq_directory* does not exist !!!\n\n\033[0m")

# only files carrying the expected suffix are considered: md5 checksums, READMEs and index
# reads sitting in the same folder would otherwise be picked up as samples
FILENAMES = sorted([f for f in next(os.walk(config["fastq_directory"]))[2] if f.endswith(fastq_suffix)])

if len(FILENAMES) == 0:
    sys.exit(''.join(["\033[1;31m\n!!! No file ending in '", fastq_suffix, "' found in *fastq_directory* !!!\n\n\033[0m"]))

RUNNAMES = [re.sub(rf"{re.escape(fastq_suffix)}$", "", i) for i in FILENAMES]

# the alternation must be grouped: '_R1|_R2.*$' is read as '(_R1)|(_R2.*$)' and strips the
# trailing part of R2 file names only, which splits each pair into two different samples
SAMPLENAMES = list(numpy.sort(numpy.unique(
    [re.sub(rf"({re.escape(read_suffix[0])}|{re.escape(read_suffix[1])}).*$", "", i) for i in RUNNAMES]
)))

# consistency check: the trimmed files are named <sample><read_suffix>_trimmed.fastq.gz, so
# the run names must be reconstructible from the sample names. When they are not, the most
# common cause is a lane/chunk field left out of *fastq_suffix* (e.g. '_001.fastq.gz')
EXPECTED_RUNNAMES = sorted([''.join([s, r]) for s in SAMPLENAMES for r in read_suffix])
if sorted(RUNNAMES) != EXPECTED_RUNNAMES:
    sys.exit(''.join(["\033[1;31m\n!!! *fastq_suffix* and *read_suffix* do not reconstruct the fastq file names !!!\n",
                      "    found:    ", ', '.join(sorted(RUNNAMES)[0:4]), " ...\n",
                      "    expected: ", ', '.join(EXPECTED_RUNNAMES[0:4]), " ...\n",
                      "    Check that *fastq_suffix* includes everything that follows the read tag.\n\n\033[0m"]))

RUNNAMES = EXPECTED_RUNNAMES


### Chromosome remove pattern
if (len(str(config["remove_other_chromosomes_pattern"])) > 0):
    chr_remove_pattern = '^chrM|^M|' + str(config["remove_other_chromosomes_pattern"])
else:
    chr_remove_pattern = '^chrM|^M'


### Read filtering flags
if str2bool(config["keep_only_proper_pairs"]):
    # -f 2 drops the orphan mates left behind by the MAPQ filter: they are still flagged as
    # paired but their mate is gone, and every paired-end aware tool downstream discards
    # them anyway. Removing them here keeps the flagstat counts honest.
    proper_pair_flag = "-f 2"
else:
    proper_pair_flag = ""


### Optional analysis outputs
if str2bool(config["run_fastq_qc"]):
    multiqc_fastq = "03_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.html"
else:
    multiqc_fastq = []


### Generation of global wildcard_constraints
wildcard_constraints:
    SAMPLE = constraint_to(SAMPLENAMES),
    RUNS = constraint_to(RUNNAMES)


# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================
# Function to run all functions
rule AAA_initialization:
    input:
        multiqc_bam_report = "03_quality_controls/multiQC_bam_filtered/multiQC_bam_filtered.html",
        multiqc_fastq = multiqc_fastq
    shell:
        """
        printf '\033[1;36mPipeline ended!\\n\033[0m'
        """

# ========================================================================================
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ========================================================================================

### Generate bowtie index if required
if not os.path.exists("".join([config["bowtie2_idx_prefix"], ".1.bt2"])):
    # ----------------------------------------------------------------------------------------
    # Genome index generation
    rule extra_generate_genome_index:
        input:
            genome = ancient(config["genome_fasta"])
        output:
            genome_idx = "".join([config["bowtie2_idx_prefix"], ".1.bt2"])
        params:
            bowtie2_idx_prefix = config["bowtie2_idx_prefix"],
            idx_folder = os.path.abspath(os.path.join(config["bowtie2_idx_prefix"], os.pardir))
        threads:
            workflow.cores
        log:
            out = "".join([config["bowtie2_idx_prefix"], "_bowtie2.build.log"])
        shell:
            """
            printf '\033[1;36mGenerating the genome index...\\n\033[0m'

            mkdir -p {params.idx_folder}/

            $CONDA_PREFIX/bin/bowtie2-build --threads {threads} -f {input.genome} {params.bowtie2_idx_prefix} > {log.out} 2>&1
            printf '\033[1;36mGenome index done.\\n\033[0m'
            """
# ----------------------------------------------------------------------------------------



# cutadapt -------------------------------------------------------------------------------
rule cutadapt_PE:
    input:
        R1 = os.path.join(config["fastq_directory"], "".join(["{SAMPLE}", read_suffix[0], fastq_suffix])),
        R2 = os.path.join(config["fastq_directory"], "".join(["{SAMPLE}", read_suffix[1], fastq_suffix]))
    output:
        R1_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", read_suffix[0], "_trimmed.fastq.gz"])),
        R2_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", read_suffix[1], "_trimmed.fastq.gz"]))
    params:
        sample = "{SAMPLE}",
        opts = str(config["cutadapt_trimm_options"]),
        fw_adapter_sequence = str(config["fw_adapter_sequence"]),
        rv_adapter_sequence = str(config["rv_adapter_sequence"])
    log:
        out = "01_trimmed_fastq/logs/cutadapt.{SAMPLE}.out",
        err = "01_trimmed_fastq/logs/cutadapt.{SAMPLE}.err"
    threads:
        max((workflow.cores - 1), 1)
    benchmark:
        "benchmarks/cutadapt_PE/cutadapt_PE---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: reads trimming...\\n\033[0m'
        mkdir -p 01_trimmed_fastq/logs/

        $CONDA_PREFIX/bin/cutadapt \
        -j {threads} -e 0.1 -q 16 -O 3 --trim-n --minimum-length 25 \
        -a {params.fw_adapter_sequence} -A {params.rv_adapter_sequence} {params.opts} \
        -o {output.R1_trimm} -p {output.R2_trimm} {input.R1} {input.R2} > {log.out} 2> {log.err}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Reads alignment
rule bowtie2_mapping:
    input:
        R1_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", read_suffix[0], "_trimmed.fastq.gz"])),
        R2_trimm = os.path.join("01_trimmed_fastq", "".join(["{SAMPLE}", read_suffix[1], "_trimmed.fastq.gz"])),
        genome_idx = "".join([config["bowtie2_idx_prefix"], ".1.bt2"])
    output:
        SAM = temp("02_BAM/{SAMPLE}.sam")
    params:
        build_SAM = "02_BAM/",
        build_log_dir = "02_BAM/bowtie2_aln_summary/",
        bowtie2_idx_prefix = config["bowtie2_idx_prefix"],
        minFragmentLength = str(config["bowtie2_minFragmentLength"]),
        maxFragmentLength = str(config["bowtie2_maxFragmentLength"]),
        extra_params = str(config["bowtie2_extra_parameters"]),
        sample = "{SAMPLE}"
    threads:
        workflow.cores
    log:
        metrics = os.path.join("02_BAM/bowtie2_aln_summary/{SAMPLE}_bowtie2_METRICS.out"),
        summary = os.path.join("02_BAM/bowtie2_aln_summary/{SAMPLE}_bowtie2_summary.err")
    benchmark:
        "benchmarks/bowtie2_mapping/bowtie2_mapping---{SAMPLE}_benchmark.txt"
    shell:
        """
        mkdir -p {params.build_SAM}
        mkdir -p {params.build_log_dir}

        printf '\033[1;36m{params.sample}: alignment of the trimmed reads (bowtie2)...\\n\033[0m'

        # --end-to-end and --local are mutually exclusive, bowtie2 silently keeps the last
        # one it is given. The CUT&Tag protocol calls for end-to-end alignment, so --local
        # must not appear here.
        $CONDA_PREFIX/bin/bowtie2 \
        -x {params.bowtie2_idx_prefix} \
        -1 {input.R1_trimm} \
        -2 {input.R2_trimm} \
        -S {output.SAM} \
        --end-to-end \
        --very-sensitive \
        --no-mixed \
        --no-discordant \
        --phred33 \
        -I {params.minFragmentLength} \
        -X {params.maxFragmentLength} \
        {params.extra_params} \
        --met-file {log.metrics} \
        -p {threads} 2> {log.summary}
        """
# --------------------------------------------------------------------------------------------------


# samtools mapq filter -----------------------------------------------------------------------------
rule MAPQ_MT_filter:
    input:
        source_sam = "02_BAM/{SAMPLE}.sam"
    output:
        bam_mateFixed = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_mateFixed.bam"]))),
        bam_mapq_only = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, ".bam"]))),
        bam_mapq_only_sorted = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted.bam"]))),
        bam_mapq_only_sorted_index = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted.bam.bai"]))),
        keep_chromosomes = temp(os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_keptChromosomes.txt"]))),
        idxstats_file = "02_BAM/reads_per_chromosome/{SAMPLE}_idxstats_read_per_chromosome.txt",
        bam_mapq_only_sorted_woMT = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT.bam"])),
        bam_mapq_only_sorted_woMT_index = os.path.join("02_BAM", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT.bam.bai"])),
        flagstat_filtered = os.path.join("02_BAM/flagstat/", ''.join(["{SAMPLE}_mapq", MAPQ, "_sorted_woMT_flagstat.txt"]))
    params:
        sample = "{SAMPLE}",
        MAPQ_threshold = MAPQ,
        proper_pair_flag = proper_pair_flag,
        chr_remove_pattern = chr_remove_pattern
    threads:
        split_cores(len(SAMPLENAMES))
    benchmark:
        "benchmarks/MAPQ_MT_filter/MAPQ_MT_filter---{SAMPLE}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.sample}: filtering MAPQ and removing mitochondrial chromosome...\\n\033[0m'

        mkdir -p 02_BAM/reads_per_chromosome
        mkdir -p 02_BAM/flagstat

        $CONDA_PREFIX/bin/samtools fixmate -@ {threads} -m -O bam {input.source_sam} {output.bam_mateFixed}

        $CONDA_PREFIX/bin/samtools view -@ {threads} -Sb -h -q {params.MAPQ_threshold} {params.proper_pair_flag} {output.bam_mateFixed} -o {output.bam_mapq_only}

        $CONDA_PREFIX/bin/samtools sort -@ {threads} {output.bam_mapq_only} -o {output.bam_mapq_only_sorted}
        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted} {output.bam_mapq_only_sorted_index}

        $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} > {output.idxstats_file}

        printf '\033[1;36m{params.sample}: Removing MT from BAM...\\n\033[0m'

        # the list of contigs to keep is written to a file and expanded in a single samtools
        # call: piping it through xargs splits the command line on fragmented assemblies and
        # concatenates several BAMs, each with its own header, into one truncated file
        $CONDA_PREFIX/bin/samtools idxstats {output.bam_mapq_only_sorted} | cut -f 1 | grep -v -E '{params.chr_remove_pattern}' | grep -v '^\\*$' > {output.keep_chromosomes}

        $CONDA_PREFIX/bin/samtools view -@ {threads} -b {output.bam_mapq_only_sorted} $(tr '\\n' ' ' < {output.keep_chromosomes}) > {output.bam_mapq_only_sorted_woMT}
        $CONDA_PREFIX/bin/samtools index -@ {threads} -b {output.bam_mapq_only_sorted_woMT} {output.bam_mapq_only_sorted_woMT_index}

        $CONDA_PREFIX/bin/samtools flagstat -@ {threads} {output.bam_mapq_only_sorted_woMT} > {output.flagstat_filtered}
        """
# --------------------------------------------------------------------------------------------------


# fastQC on trimmed fastq ---------------------------------------------------------------------------
rule fastQC_trimmed_fastq:
    input:
        fastq_trimm = os.path.join("01_trimmed_fastq", "{RUNS}_trimmed.fastq.gz")
    output:
        fastqc_html = os.path.join("03_quality_controls/trimmed_fastq_fastqc", "{RUNS}_trimmed_fastqc.html"),
        fastqc_zip = os.path.join("03_quality_controls/trimmed_fastq_fastqc", "{RUNS}_trimmed_fastqc.zip")
    params:
        outdir = "03_quality_controls/trimmed_fastq_fastqc",
        run = "{RUNS}"
    threads:
        split_cores(len(RUNNAMES))
    benchmark:
        "benchmarks/fastQC_trimmed_fastq/fastQC_trimmed_fastq---{RUNS}_benchmark.txt"
    shell:
        """
        printf '\033[1;36m{params.run}: performing fastQC on trimmed fastq...\\n\033[0m'

        mkdir -p {params.outdir}

        # the file to process is taken from the rule input rather than from a wildcard glob:
        # a glob would also pick up leftovers of previous runs made with a different sample set
        $CONDA_PREFIX/bin/fastqc -t {threads} --outdir {params.outdir} {input.fastq_trimm}
        """


rule multiQC_trimmed_fastq:
    input:
        fastqc_zip = expand(os.path.join("03_quality_controls/trimmed_fastq_fastqc", "{runs}_trimmed_fastqc.zip"), runs = RUNNAMES)
    output:
        multiqc_fastqc_report = "03_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.html"
    params:
        fastqc_zip_dir = "03_quality_controls/trimmed_fastq_fastqc/",
        out_directory = "03_quality_controls/trimmed_fastq_multiQC/",
        cutadapt_logs = "01_trimmed_fastq/logs/",
        multiqc_fastqc_report_name = "multiQC_report_trimmed_fastq.html"
    log:
        out = "03_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.out",
        err = "03_quality_controls/trimmed_fastq_multiQC/multiQC_report_trimmed_fastq.err"
    threads: 1
    benchmark:
        "benchmarks/multiQC_trimmed_fastq/multiQC_trimmed_fastq---benchmark.txt"
    shell:
        """
        printf '\033[1;36mGenerating multiQC report for trimmed fastq...\\n\033[0m'

        mkdir -p {params.out_directory}

        $CONDA_PREFIX/bin/multiqc -f \
        -o {params.out_directory} \
        -n {params.multiqc_fastqc_report_name} \
        --dirs \
        {params.fastqc_zip_dir} \
        {params.cutadapt_logs} > {log.out} 2> {log.err}
        """
# --------------------------------------------------------------------------------------------------


# multiQC aligned bams -------------------------------------------------------------------------------
rule multiQC_bam:
    input:
        flagstat_filtered = expand(os.path.join("02_BAM/flagstat/", ''.join(["{sample}_mapq", MAPQ, "_sorted_woMT_flagstat.txt"])), sample = SAMPLENAMES)
    output:
        multiqc_bam_report = "03_quality_controls/multiQC_bam_filtered/multiQC_bam_filtered.html"
    params:
        out_directory = "03_quality_controls/multiQC_bam_filtered/",
        flagstat_dir = "02_BAM/flagstat/",
        bowtie2_summary_dir = "02_BAM/bowtie2_aln_summary/",
        multiqc_bam_report_name = "multiQC_bam_filtered.html"
    log:
        out = "03_quality_controls/multiQC_bam_filtered/multiQC_bam_filtered.out",
        err = "03_quality_controls/multiQC_bam_filtered/multiQC_bam_filtered.err"
    threads: 1
    benchmark:
        "benchmarks/multiQC_bam/multiQC_bam---benchmark.txt"
    shell:
        """
        printf '\033[1;36mPerforming multiQC on filtered bam flagstat and alignment summaries...\\n\033[0m'

        mkdir -p {params.out_directory}

        # the bowtie2 summaries are parsed as well, otherwise the alignment rates never reach
        # the report
        $CONDA_PREFIX/bin/multiqc -f \
        -o {params.out_directory} \
        -n {params.multiqc_bam_report_name} \
        --dirs \
        {params.flagstat_dir} \
        {params.bowtie2_summary_dir} > {log.out} 2> {log.err}
        """
# --------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#                                 END pipeline
# ------------------------------------------------------------------------------
