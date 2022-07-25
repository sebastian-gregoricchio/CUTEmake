from typing import List
import pathlib
import re
import os
#os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy

# Function to handle the values for the wilcards
def constraint_to(values: List[str]) -> str:
    """
    From a list, return a regular expression allowing each
    value and not other.
    ex: ["a", "b", "v"] -> (a|b|v)
    """
    if isinstance(values, str):
        raise ValueError("constraint_to(): Expected a list, got str instead")

    return "({})".format("|".join(values))

# working diirectory
workdir: config["output_directory"]


# get the unique samples names and other variables
FILENAMES = next(os.walk(config["runs_directory"]))[2]
RUNNAMES = [re.sub(rf"{config['fastq_extension']}$", "", i) for i in FILENAMES]
SAMPLENAMES = numpy.unique([re.sub(rf"{config['runs_suffix'][0]}|{config['runs_suffix'][1]}.*$", "", i) for i in RUNNAMES])




rm_duplicates = str(config["remove_duplicates"]).lower()
if (eval(str(rm_duplicates.capitalize())) == True):
    DUP="dedup"
else:
    DUP="mdup"


# define folder names
BINS=config["bigWig_binSize"]
PEAKSDIR = "05_Peaks_SEACR/"
SUMMARYDIR = "06_Overall_quality_and_info/"
GATKDIR = "07_Variant_calling/"
PEAKCALLER = "SEACR"


# generation of global wildcard_constraints
wildcard_constraints:
    RUNS=constraint_to(RUNNAMES),
    SAMPLES=constraint_to(SAMPLENAMES)


# GATK outputs
if (eval(str(config["call_variants"]))):
    bsqr_table = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_bsqr.table"), sample=SAMPLENAMES))
    dedup_BAM_bsqr = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_bsqr.bam"), sample=SAMPLENAMES))
    dedup_BAM_bsqr_index = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_bsqr.bam.bai"), sample=SAMPLENAMES))
    vcf = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_gatk.vcf.gz"), sample=SAMPLENAMES))
else:
    dedup_BAM_bsqr = []
    dedup_BAM_bsqr_index = []
    vcf = []


if (eval(str(config["call_SNPs"]))):
    snp = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_gatk-snp.vcf"), sample=SAMPLENAMES))
else:
    snp = []


if (eval(str(config["call_indels"]))):
    indels = ancient(expand(os.path.join(GATKDIR, "{sample}/{sample}_gatk.vcf.gz"), sample=SAMPLENAMES))
else:
    indels = []


correlation_outputs = []
if (len(SAMPLENAMES) > 1):
    correlation_outputs.append(ancient(os.path.join(SUMMARYDIR, "Sample_comparisons/PCA_on_BigWigs_wholeGenome.pdf"))) # PCA
    correlation_outputs.append(ancient(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf"))) # heatamap_spearman
    correlation_outputs.append(ancient(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf"))) # hetamap_pearson
    correlation_outputs.append(ancient(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf"))) # scatterplot_spearman
    correlation_outputs.append(ancient(os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf"))) # scatterplot_pearson


# ========================================================================================
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ========================================================================================



# ========================================================================================
# Function to run all funtions
rule AAA_initialization:
    input:
        # Rule A
        fastQC_raw_zip = ancient(expand(os.path.join("01_fastQC_raw", "{run}_fastqc.zip"), run=RUNNAMES)),

        # Rule B
        multiQC_raw_html = ancient("01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw.html"),

        # Rule C
        #SAM = ancient(expand(os.path.join("01b_SAM_tempFolder/", "{sample}.sam"), sample=SAMPLENAMES)),

        # Rule D
        filtBAM_sorted_woMT = ancient(expand(os.path.join("02_BAM/", "{sample}_mapQ{MAPQ}_sorted_woMT.bam"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
        filtBAM_sorted_woMT_index = ancient(expand(os.path.join("02_BAM/", "{sample}_mapQ{MAPQ}_sorted_woMT.bam.bai"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
        flagstat_unfiltered_BAM = ancient(expand(os.path.join("02_BAM/flagstat/", "{sample}_flagstat_UNfiltered_bam.txt"), sample=SAMPLENAMES)),
        flagstat_on_filtered_woMT_BAM = ancient(expand(os.path.join("02_BAM/flagstat/", "{sample}_flagstat_filtered_bam_woMT.txt"), sample=SAMPLENAMES)),

        # Rule E
        # dedup_BAM = ancient(expand(os.path.join("03_BAM_{DUP}/unshifted_bams/", "{sample}_mapQ{MAPQ}_sorted_woMT_{DUP}.bam"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
        # dedup_BAM_index = ancient(expand(os.path.join("03_BAM_{DUP}/unshifted_bams/", "{sample}_mapQ{MAPQ}_sorted_woMT_{DUP}.bai"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]))),
        # dedup_BAM_metrics = ancient(expand(os.path.join("03_BAM_{DUP}/metrics", "{sample}_metrics_woMT_{DUP}_bam.txt"), sample=SAMPLENAMES)),
        # dedup_BAM_flagstat = ancient(expand(os.path.join("03_BAM_{DUP}/flagstat/", "{sample}_flagstat_filtered_bam_woMT_{DUP}.txt"), sample=SAMPLENAMES)),

        # Rule F
        dedup_BAM_shifted_sorted = ancient(expand(os.path.join("03_BAM_{dup}/", "{sample}_mapQ{MAPQ}_woMT_{dup}_shifted_sorted.bam"), sample=SAMPLENAMES, dup=DUP, MAPQ=str(config["mapQ_cutoff"]))),
        dedup_BAM_shifted_sorted_index = ancient(expand(os.path.join("03_BAM_{dup}/", "{sample}_mapQ{MAPQ}_woMT_{dup}_shifted_sorted.bam.bai"), sample=SAMPLENAMES, dup=DUP, MAPQ=str(config["mapQ_cutoff"]))),
        dedup_BAM_flagstat_shifted_sorted = ancient(expand(os.path.join("03_BAM_{dup}/flagstat/", "{sample}_flagstat_woMT_{dup}_shifted_sorted.txt"), dup=DUP, sample=SAMPLENAMES)),

        # Rule G
        fastQC_zip_BAM = ancient(expand(os.path.join("03_BAM_{dup}/fastQC/", "{sample}_mapQ{MAPQ}_sorted_woMT_{dup}_fastqc.zip"), sample=SAMPLENAMES, dup=DUP, MAPQ=str(config["mapQ_cutoff"]))),

        # Rule H
        multiQC_BAM_html = ancient(expand("03_BAM_{dup}/fastQC/multiQC_{dup}_bams/multiQC_report_BAMs_{dup}.html", dup=DUP)),

        # Rule I
        fragmentSizePlot = expand(os.path.join("03_BAM_{dup}/fragmentSizeDistribution_plots/", "{sample}_fragment_size_distribution.pdf"), sample=SAMPLENAMES, dup=DUP),
        report_pdf = expand("03_BAM_{dup}/fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf", dup=DUP),

        # Rule J
        #scalingFactors_txt_result = "04_Normalization/scalingFactor/scalingFactor_results.txt",

        # Rule_normalization
        norm_bw = ancient(expand(os.path.join("04_Normalization/normalized_bigWigs/", "{sample}_mapQ{MAPQ}_woMT_{dup}_shifted_normalized_bs{binSize}.bw"), sample=SAMPLENAMES, dup=DUP, MAPQ=str(config["mapQ_cutoff"]), binSize=str(BINS))),

        # Rule peakCalling
        seacr_peaks = expand(os.path.join(PEAKSDIR, "{sample}_mapQ{MAPQ}_woMT_{dup}_shifted_seacr_top{AUC}.stringent.bed"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]), dup=DUP, AUC=str(config["AUC_fration"])),

        # Rule counts_summary
        #temp_file_counts = ancient(expand(os.path.join(SUMMARYDIR, "{sample}_counts_summary.temp"), sample=SAMPLENAMES)),
        summary_file = os.path.join(SUMMARYDIR, "Counts/counts_summary.txt"),
        summary_file_temp = os.path.join(SUMMARYDIR, "Counts/summary_file.temp"),

        # Rules L PCA, plotCorrelation, LorenzCurve
        correlation_outputs = correlation_outputs,
        lorenz_plot = os.path.join(SUMMARYDIR, "Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf"),

        # Rule Heatmap zScores peaks
        rawScores_hetamap_SEACR = ancient(os.path.join(SUMMARYDIR, ''.join(["Sample_comparisons/Peak_comparison/Heatmaps/Heatmap_on_log1p.rawScores_for_", PEAKCALLER, ".peaks_union_population.pdf"]))),

        # Rules gatk variant calling
        dedup_BAM_bsqr = dedup_BAM_bsqr,
        dedup_BAM_bsqr_index = dedup_BAM_bsqr_index,
        vcf = vcf,
        snp = snp,
        indels = indels

    #params:
    #    summary_file = str(os.path.join(SUMMARYDIR, "Counts/counts_summary.txt"))

    shell:
        """
        mkdir -p 01b_SAM_tempFolder
        rm -R 01b_SAM_tempFolder

        printf '\033[1;36mPipeline ended!\\n\033[0m'
        """
# ========================================================================================


# ----------------------------------------------------------------------------------------
# Perform the FastQC on raw fastq.gz
rule A_fastQC_raw:
    input:
        fastq_gz = ancient(os.path.join(config["runs_directory"], "".join(["{RUNS}", config['fastq_extension']])))
    output:
        html = os.path.join("01_fastQC_raw","{RUNS}_fastqc.html"),
        zip =  os.path.join("01_fastQC_raw","{RUNS}_fastqc.zip")
    params:
        build_fastqcDir = os.path.dirname("01_fastQC_raw/multiQC_raw/"),
        fastQC_raw_outdir = os.path.join(config["output_directory"], "01_fastQC_raw"),
        run = "{RUNS}",
        CPUs = config["fastQC_threads"]
    threads:
        config["fastQC_threads"]
    shell:
        """
        printf '\033[1;36m{params.run}: Performing fastQC on raw fastq...\\n\033[0m'

        mkdir -p {params.build_fastqcDir}
        fastqc -t {params.CPUs} --outdir {params.fastQC_raw_outdir} {input.fastq_gz}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform multiQC for raw fastq reports
rule B_multiQC_raw:
    input:
        fastqc_zip = ancient(expand(os.path.join("01_fastQC_raw", "{run}_fastqc.zip"), run=RUNNAMES))
    output:
        multiqcReportRaw = "01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw.html"
    params:
        fastqc_Raw_reports = os.path.join("01_fastQC_raw", "*.zip"),
        multiQC_raw_outdir = os.path.join(config["output_directory"], "01_fastQC_raw/multiQC_raw/")
    shell:
        """
        printf '\033[1;36mGenerating multiQC report for fatsq quality test...\\n\033[0m'
        multiqc -f --outdir {params.multiQC_raw_outdir} -n multiQC_report_fastqRaw.html {params.fastqc_Raw_reports}
        """
# ----------------------------------------------------------------------------------------

if not os.path.exists("".join([config["bowtie2_idx_prefix"], ".1.bt2"])):
    # ----------------------------------------------------------------------------------------
    # Reads alignement
    rule Cextra_generate_genome_index:
        input:
            multiqcReportRaw = "01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw.html"
        output:
            genome_fai = "".join([config["bowtie2_idx_prefix"], ".1.bt2"])
        params:
            genome = config["genome_fasta"],
            bowtie2_idx_prefix = config["bowtie2_idx_prefix"],
            idx_folder = os.path.abspath(os.path.join(config["bowtie2_idx_prefix"], os.pardir))
        threads:
            config["bowtie2_threads"]
        shell:
            """
            printf '\033[1;36mGenerating the genome index...\\n\033[0m'

            mkdir -p {params.idx_folder}/

            bowtie2-build -f {params.genome} {params.bowtie2_idx_prefix}
            printf '\033[1;36mGenome index done.\\n\033[0m'
            """
# ----------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------
# Reads alignement
rule C_bowtie2_align:
    input:
        multiqcReportRaw = ancient("01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw.html"),
        R1 = ancient(os.path.join(config["runs_directory"], "".join(["{SAMPLES}", config['runs_suffix'][0], config['fastq_extension']]))),
        R2 = ancient(os.path.join(config["runs_directory"], "".join(["{SAMPLES}", config['runs_suffix'][1], config['fastq_extension']]))),
        genome_fai = "".join([config["bowtie2_idx_prefix"], ".1.bt2"])
    output:
        SAM = os.path.join("01b_SAM_tempFolder/", "{SAMPLES}.sam")
    params:
        build_SAM = "01b_SAM_tempFolder/",
        build_log_dir = "02_BAM/bowtie2_aln_summary/",
        bowtie2_idx_prefix = config["bowtie2_idx_prefix"],
        sample = "{SAMPLES}",
        CPUs = config["bowtie2_threads"]
    threads:
        config["bowtie2_threads"]
    log:
        out = os.path.join("02_BAM/bowtie2_aln_summary/{SAMPLES}_bowtie2_METRICS.out"),
        err = os.path.join("02_BAM/bowtie2_aln_summary/{SAMPLES}_bowtie2_summary.err")
    shell:
        """
        mkdir -p {params.build_SAM}
        mkdir -p {params.build_log_dir}

        printf '\033[1;36m{params.sample}: alignment of the reads...\\n\033[0m'

        bowtie2 \
        -x {params.bowtie2_idx_prefix} \
        -1 {input.R1} \
        -2 {input.R2} \
        -S {output.SAM} \
        --end-to-end \
        --local --very-sensitive \
        --no-mixed \
        --no-discordant \
        --phred33 \
        -I 10 \
        -X 700 \
        --met-file {log.out} \
        -p {params.CPUs} 2> {log.err}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# SAM filtering for mapping quality and BAM generation | BAM MT-reads removal
rule D_sam_to_bam:
    input:
        SAM = ancient(os.path.join("01b_SAM_tempFolder/", "{SAMPLES}.sam"))
    output:
        #filtBAM_toFix = temp(os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_toFix.bam"]))),
        filtBAM = temp(os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), ".bam"]))),
        filtBAM_sorted = temp(os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted.bam"]))),
        filtBAM_sorted_index = temp(os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted.bam.bai"]))),
        flagstat_on_unfiltered_BAM = os.path.join("02_BAM/flagstat/", "{SAMPLES}_flagstat_UNfiltered_bam.txt"),
        filtBAM_sorted_woMT = os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT.bam"])),
        filtBAM_sorted_woMT_index = os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT.bam.bai"])),
        flagstat_on_filtered_woMT_BAM = os.path.join("02_BAM/flagstat/", "{SAMPLES}_flagstat_filtered_bam_woMT.txt")
    params:
        build_BAM = os.path.dirname("02_BAM/flagstat/"),
        mapq_cutoff = str(config["mapQ_cutoff"]),
        sample = "{SAMPLES}",
        CPUs = config["SAMtools_threads"]
    threads:
        config["SAMtools_threads"]
    shell:
        """
        mkdir -p {params.build_BAM}

        printf '\033[1;36m{params.sample}: filtering SAM and generate BAM...\\n\033[0m'
        samtools view -@ {params.CPUs} -b -q {params.mapq_cutoff} {input.SAM} -o {output.filtBAM}

        printf '\033[1;36m{params.sample}: sorting BAM...\\n\033[0m'
        samtools sort -@ {params.CPUs} {output.filtBAM} -o {output.filtBAM_sorted}
        samtools index -@ {params.CPUs} -b {output.filtBAM_sorted} {output.filtBAM_sorted_index}

        printf '\033[1;36m{params.sample}: Getting flagstat from unfiltered BAM...\\n\033[0m'
        samtools flagstat {output.filtBAM_sorted} -@ {params.CPUs} > {output.flagstat_on_unfiltered_BAM}

        printf '\033[1;36m{params.sample}: Removing MT reads from BAM...\\n\033[0m'
        samtools idxstats {output.filtBAM_sorted} | cut -f 1 | grep -v ^chrM | grep -v ^M | xargs samtools view -@ {params.CPUs} -b {output.filtBAM_sorted} > {output.filtBAM_sorted_woMT}

        printf '\033[1;36m{params.sample}: BAM indexing...\\n\033[0m'
        samtools index -@ {params.CPUs} -b {output.filtBAM_sorted_woMT} {output.filtBAM_sorted_woMT_index}

        printf '\033[1;36m{params.sample}: Getting flagstat from BAM without MT-DNA...\\n\033[0m'
        samtools flagstat {output.filtBAM_sorted_woMT} -@ {params.CPUs} > {output.flagstat_on_filtered_woMT_BAM}

        rm {input.SAM}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# BAM duplicates removal and relative flagstat | BAM reads shifting
rule E_bam_deduplication:
    input:
        BAM = ancient(os.path.join("02_BAM/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT.bam"])))
    output:
        dedup_BAM = os.path.join("03_BAM_{DUP}/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}.bam"])),
        dedup_BAM_index = os.path.join("03_BAM_{DUP}/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}.bai"])),
        dedup_BAM_metrics = os.path.join("03_BAM_{DUP}/metrics", "{SAMPLES}_metrics_woMT_{DUP}_bam.txt"),
        dedup_BAM_flagstat = os.path.join("03_BAM_{DUP}/flagstat/", "{SAMPLES}_flagstat_filtered_bam_woMT_{DUP}.txt")
    params:
        build_BAMdedup_metrics = os.path.dirname("03_BAM_{DUP}/metrics/"),
        build_BAMdedup_unshifted = os.path.dirname("03_BAM_{DUP}/unshifted_bams/"),
        build_BAMdedup_flagstat = os.path.dirname("03_BAM_{DUP}/flagstat/"),
        build_BAMdedup_multiQC = os.path.dirname("03_BAM_{DUP}/fastQC/multiQC_{DUP}_bams/"),
        build_BAMdedup_fragmentPlots = os.path.dirname("03_BAM_{DUP}/fragmentSizeDistribution_plots/"),
        sample = "{SAMPLES}",
        max_records = config["PICARD_max_records_in_ram"],
        max_file_handles = config["PICARD_max_file_handles_for_read_ends_map"],
        rm_dup = rm_duplicates,
        CPUs = config["SAMtools_threads"]
    threads:
        config["SAMtools_threads"]
    shell:
        """
        mkdir -p {params.build_BAMdedup_metrics}
        mkdir -p {params.build_BAMdedup_unshifted}
        mkdir -p {params.build_BAMdedup_flagstat}
        mkdir -p {params.build_BAMdedup_multiQC}
        mkdir -p {params.build_BAMdedup_fragmentPlots}

        printf '\033[1;36m{params.sample}: Removing/Marking BAM duplicates...\\n\033[0m'

        picard MarkDuplicates \
        --CREATE_INDEX true \
        --INPUT {input.BAM} \
        --OUTPUT {output.dedup_BAM} \
        --METRICS_FILE {output.dedup_BAM_metrics} \
        --ASSUME_SORT_ORDER coordinate \
        --REMOVE_DUPLICATES {params.rm_dup} \
        --VALIDATION_STRINGENCY STRICT \
        --MAX_RECORDS_IN_RAM {params.max_records} \
        --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP {params.max_file_handles}

        samtools flagstat {output.dedup_BAM} -@ {params.CPUs} > {output.dedup_BAM_flagstat}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# BAM reads shifting
rule F_bam_shifting:
    input:
        dedup_BAM = ancient(os.path.join("03_BAM_{DUP}/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}.bam"]))),
        dedup_BAM_index = ancient(os.path.join("03_BAM_{DUP}/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}.bai"])))
    output:
        dedup_BAM_shifted_toSort = temp(os.path.join("03_BAM_{DUP}/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}_shifted.ToSort.bam"]))),
        dedup_BAM_shifted_sorted = os.path.join("03_BAM_{DUP}/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_sorted.bam"])),
        dedup_BAM_shifted_sorted_index = os.path.join("03_BAM_{DUP}/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_sorted.bam.bai"])),
        dedup_BAM_shifted_sorted_flagstat = os.path.join("03_BAM_{DUP}/flagstat/", "{SAMPLES}_flagstat_woMT_{DUP}_shifted_sorted.txt")
    params:
        sample = "{SAMPLES}",
        minFragmentLength = str(config["minFragmentLength"]),
        maxFragmentLength = str(config["maxFragmentLength"]),
        CPUs = config["SAMtools_threads"]
    threads:
        config["SAMtools_threads"]
    shell:
        """
        printf '\033[1;36m{params.sample}: Shifting reads in BAM...\\n\033[0m'
        alignmentSieve -p {params.CPUs} --ATACshift --bam {input.dedup_BAM} --outFile {output.dedup_BAM_shifted_toSort} --minFragmentLength {params.minFragmentLength} --maxFragmentLength {params.maxFragmentLength}

        printf '\033[1;36m{params.sample}: Sorting shifted BAM...\\n\033[0m'
        samtools sort -@ {params.CPUs} {output.dedup_BAM_shifted_toSort} -o {output.dedup_BAM_shifted_sorted}
        samtools index -@ {params.CPUs} -b {output.dedup_BAM_shifted_sorted} {output.dedup_BAM_shifted_sorted_index}

        printf '\033[1;36m{params.sample}: Getting flagstat from shifted BAM...\\n\033[0m'
        samtools flagstat {output.dedup_BAM_shifted_sorted} -@ {params.CPUs} > {output.dedup_BAM_shifted_sorted_flagstat}

        echo '------------------------------------------------------------------------'
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# FastQC on BAMs
rule G_fastQC_BAMs:
    input:
        dedup_BAM = ancient(os.path.join("03_BAM_{DUP}/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}.bam"]))),
        dedup_BAM_index = ancient(os.path.join("03_BAM_{DUP}/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}.bai"])))
    output:
        html = os.path.join("03_BAM_{DUP}/fastQC/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}_fastqc.html"])),
        zip = os.path.join("03_BAM_{DUP}/fastQC/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}_fastqc.zip"]))
    params:
        fastQC_BAMs_outdir = os.path.join(config["output_directory"], "03_BAM_{DUP}/fastQC/"),
        sample = "{SAMPLES}",
        CPUs = config["fastQC_threads"]
    threads:
        config["fastQC_threads"]
    shell:
        """
        printf '\033[1;36m{params.sample}: Performing fastQC on deduplicated bam...\\n\033[0m'
        fastqc -t {params.CPUs} --outdir {params.fastQC_BAMs_outdir} {input.dedup_BAM}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform multiQC for BAMs
rule H_multiQC_BAMs:
    input:
        BAM_fastqc_zip = ancient(expand(os.path.join("03_BAM_{dup}/fastQC/", "{sample}_mapQ{MAPQ}_sorted_woMT_{dup}_fastqc.zip"), sample=SAMPLENAMES, dup=DUP, MAPQ=str(config["mapQ_cutoff"])))
    output:
        multiqcReportBAM = "03_BAM_{DUP}/fastQC/multiQC_{DUP}_bams/multiQC_report_BAMs_{DUP}.html"
    params:
        fastQC_BAM_reports = os.path.join("03_BAM_{DUP}/fastQC/", "*.zip"),
        multiQC_BAM_outdir = os.path.join(config["output_directory"], "03_BAM_{DUP}/fastQC/multiQC_{DUP}_bams/")
    shell:
        """
        printf '\033[1;36mGenerating multiQC report from deduplicated bam quality test...\\n\033[0m'
        multiqc -f --outdir {params.multiQC_BAM_outdir} -n multiQC_report_BAMs_{DUP}.html {params.fastQC_BAM_reports}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform fragment size distribution plot
rule I1_fragment_size_distribution:
    input:
        BAM = ancient(os.path.join("03_BAM_{DUP}/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}.bam"]))),
        BAM_index = ancient(os.path.join("03_BAM_{DUP}/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}.bai"]))),
        multiqcReportBAM = ancient("03_BAM_{DUP}/fastQC/multiQC_{DUP}_bams/multiQC_report_BAMs_{DUP}.html")
    output:
        plot = os.path.join("03_BAM_{DUP}/fragmentSizeDistribution_plots/", "{SAMPLES}_fragment_size_distribution.pdf")
    params:
        sample = "{SAMPLES}",
        plotFormat = config["plot_format"],
        binSize = str(config["window_length"]),
        blacklist = config["blacklist_file"],
        CPUs = config["bamPEFragmentSize_threads"],
        max_fragment_length = config["max_fragment_length"]
    threads:
        config["bamPEFragmentSize_threads"]
    shell:
        """
        printf '\033[1;36m{params.sample}: Plotting the fragment size distribution...\\n\033[0m'

        bamPEFragmentSize \
        -p {params.CPUs} \
        -b {input.BAM} \
        --plotFileFormat {params.plotFormat} \
        --plotTitle {params.sample} \
        --samplesLabel {params.sample} \
        --binSize {params.binSize} \
        --maxFragmentLength {params.max_fragment_length} \
        --blackListFileName {params.blacklist} \
        -o {output.plot}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Perform fragment size distribution plot
rule I2_fragment_size_distribution_report:
    input:
        plots = ancient(expand(os.path.join("03_BAM_{dup}/fragmentSizeDistribution_plots/", "{sample}_fragment_size_distribution.pdf"), sample = SAMPLENAMES, dup=DUP))
    output:
        report_pdf ="03_BAM_{DUP}/fragmentSizeDistribution_plots/ALL.SAMPLES_fragmentSizeDistribution_plots.pdf"
    threads:
        1
    shell:
        """
        printf '\033[1;36mMerging fragmentSizeDistribution reports in a unique PDF...\\n\033[0m'
        pdfcombine 03_BAM_{DUP}/fragmentSizeDistribution_plots/*_fragment_size_distribution.pdf -o {output.report_pdf} -sf
        """


# ----------------------------------------------------------------------------------------
# bigWig generation from BAM
rule J1_bigWig_normalization_woCNVcorrection:
    input:
        dedup_BAM_shifted_sorted = ancient(os.path.join("03_BAM_{DUP}/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_sorted.bam"]))),
        dedup_BAM_shifted_sorted_index = ancient(os.path.join("03_BAM_{DUP}/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_sorted.bam.bai"])))
    output:
        norm_bw = os.path.join("04_Normalization/normalized_bigWigs/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_normalized_bs", str(config["bigWig_binSize"]), ".bw"]))
    params:
        sample = "{SAMPLES}",
        binSize = config["bigWig_binSize"],
        normalization_method = config["normalization_method"],
        effective_genomeSize = config["effective_genomeSize"],
        ignore_for_normalization = config["ignore_for_normalization"],
        blacklist = config["blacklist_file"],
        CPUs = config["bamCoverage_threads"]
    threads:
        config["bamCoverage_threads"]
    shell:
        """
        printf "\033[1;36m{params.sample}: generation of the normalized bigWig file...\\n\033[0m"

        bamCoverage -p {params.CPUs} \
        --bam {input.dedup_BAM_shifted_sorted} \
        --binSize {params.binSize} \
        --normalizeUsing {params.normalization_method} \
        --effectiveGenomeSize {params.effective_genomeSize} \
        --ignoreForNormalization {params.ignore_for_normalization} \
        --extendReads \
        --blackListFileName {params.blacklist} \
        -of "bigwig" \
        -o {output.norm_bw}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# SEACR peakCalling on uncorrected bams
rule J3_SEACR_peakCalling:
    input:
        dedup_BAM_shifted_sorted = ancient(os.path.join("03_BAM_{DUP}/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_sorted.bam"]))),
        dedup_BAM_shifted_sorted_index = ancient(os.path.join("03_BAM_{DUP}/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_sorted.bam.bai"]))),
        norm_bw = ancient(os.path.join("04_Normalization/normalized_bigWigs/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_normalized_bs", str(config["bigWig_binSize"]), ".bw"])))
    output:
        norm_bdg = temp(os.path.join("04_Normalization/normalized_bigWigs/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_normalized_bs", str(config["bigWig_binSize"]), ".bg"]))),
        seacr_peaks = os.path.join(PEAKSDIR, ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_seacr_top", str(config["AUC_fration"]), ".stringent.bed"]))
    params:
        dir = ''.join([PEAKSDIR, "log/"]),
        basename_seacr_peaks = os.path.join(PEAKSDIR, ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_seacr_top", str(config["AUC_fration"])])),
        AUC_fration = str(config["AUC_fration"]),
        sample = "{SAMPLES}"
    log:
        out = os.path.join(PEAKSDIR, ''.join(["log/{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_woMT_{DUP}_shifted_seacr_top", str(config["AUC_fration"]), ".log"]))
    shell:
        """
        mkdir -p {params.dir}

        printf '\033[1;36m{params.sample}: Converting normalized bigWig to bedGraph (for SEACR peak calling)...\\n\033[0m'
        bigWigToBedGraph {input.norm_bw} {output.norm_bdg}

        printf '\033[1;36m{params.sample}: Calling peaks by SEACR...\\n\033[0m'
        CONDAPATH=$(conda info | grep 'active env location : ' | sed 's/active env location : //')
        $CONDAPATH/bin/SEACR_1.3.sh {output.norm_bdg} {params.AUC_fration} non stringent {params.basename_seacr_peaks} 2> {log.out}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Computation of the counts summary table
rule K_counts_summary:
    input:
        R1 = ancient(expand(os.path.join(config["runs_directory"], "".join(["{sample}", config['runs_suffix'][0], config['fastq_extension']])), sample = SAMPLENAMES)),
        R2 = ancient(expand(os.path.join(config["runs_directory"], "".join(["{sample}", config['runs_suffix'][1], config['fastq_extension']])), sample = SAMPLENAMES)),
        flagstat_on_unfiltered_BAM = ancient(expand(os.path.join("02_BAM/flagstat/", "{sample}_flagstat_UNfiltered_bam.txt"), sample = SAMPLENAMES)),
        flagstat_on_filtered_woMT_BAM = ancient(expand(os.path.join("02_BAM/flagstat/", "{sample}_flagstat_filtered_bam_woMT.txt"), sample = SAMPLENAMES)),
        dedup_BAM_flagstat = ancient(expand(os.path.join("03_BAM_{dup}/flagstat/", "{sample}_flagstat_filtered_bam_woMT_{dup}.txt"), sample = SAMPLENAMES, dup=DUP)),
        dedup_BAM_shifted_sorted_flagstat = ancient(expand(os.path.join("03_BAM_{dup}/flagstat/", "{sample}_flagstat_woMT_{dup}_shifted_sorted.txt"), sample = SAMPLENAMES, dup=DUP)),
        #scaling_factors = ancient("04_Normalization/scalingFactor/scalingFactor_results.txt"),
        peaks_file = ancient(expand(os.path.join(PEAKSDIR, "{sample}_mapQ{MAPQ}_woMT_{dup}_shifted_seacr_top{AUC}.stringent.bed"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]), dup=DUP, AUC=str(config["AUC_fration"])))
    output:
        summary_file = os.path.join(SUMMARYDIR, "Counts/counts_summary.txt"),
        summary_file_temp = temp(os.path.join(SUMMARYDIR, "Counts/summary_file.temp"))
    params:
        build_summary_directory = os.path.dirname(SUMMARYDIR),
        R1_suffix = config['runs_suffix'][0],
        R2_suffix = config['runs_suffix'][1],
        sample_list = str(' '.join(SAMPLENAMES)),
        multiQC_report = "01_fastQC_raw/multiQC_raw/multiQC_report_fastqRaw_data/multiqc_general_stats.txt",
        peaks_dir = PEAKSDIR,
        FRiP_threshold = config["FRiP_threshold"]
    threads: 1
    shell:
        """
        mkdir -p {params.build_summary_directory}/Counts/subread_featureCounts_output/

        printf '\033[1;36mGeneration of a general counts summary table...\\n\033[0m'
        printf Sample'\\t'Reads_R1'\\t'Reads_R2'\\t'Reads_total'\\t'unfiltered_BAM'\\t'Percentage_MT'\\t'dedup_BAM'\\t'duplicated_reads'\\t'shifted_BAM'\\t'loss_post_shifting'\\t'n.peaks'\\t'FRiP.perc'\\t'FRiP.quality'\\n' > {output.summary_file}

        for NAME in {params.sample_list}
        do
            printf '\033[1;36m     - %s: adding stats to summary table...\\n\033[0m' $NAME

            mkdir -p {params.build_summary_directory}/Counts/subread_featureCounts_output/${{NAME}}/

            R1=$(grep ${{NAME}}{params.R1_suffix} {params.multiQC_report} | cut -f 6 | sed 's/\.[0-9]\+//')
            R2=$(grep ${{NAME}}{params.R2_suffix} {params.multiQC_report} | cut -f 6 | sed 's/\.[0-9]\+//')
            TOTAL=$((R1 + R2))

            unfilteredBAM=$(grep mapped 02_BAM/flagstat/${{NAME}}_flagstat_UNfiltered_bam.txt | head -n 1 | cut -f 1 -d ' ')
            woMT_BAM=$(grep mapped 02_BAM/flagstat/${{NAME}}_flagstat_filtered_bam_woMT.txt | head -n 1 | cut -f 1 -d ' ')
            percMT=$(echo "scale=1; (100 - (($woMT_BAM/$unfilteredBAM) * 100))" | bc)

            dedupBAM=$(grep mapped 03_BAM_{DUP}/flagstat/${{NAME}}_flagstat_filtered_bam_woMT_{DUP}.txt | head -n 1 | cut -f 1 -d ' ')
            dedupREADS=$((woMT_BAM - dedupBAM))

            shiftedBAM=$(grep mapped 03_BAM_{DUP}/flagstat/${{NAME}}_flagstat_woMT_{DUP}_shifted_sorted.txt | head -n 1 | cut -f 1 -d ' ')
            lossReads=$((dedupBAM - shiftedBAM))

            peaks=$(wc -l {params.peaks_dir}${{NAME}}*.*bed | cut -f 1 -d ' ')

            awk 'BEGIN{{FS=OFS="\\t"; print "GeneID\\tChr\\tStart\\tEnd\\tStrand"}}{{print $4, $1, $2+1, $3, "."}}' {params.peaks_dir}${{NAME}}*.*bed > {params.peaks_dir}${{NAME}}.saf
            FEATURECOUNTSLOG={params.build_summary_directory}/Counts/subread_featureCounts_output/${{NAME}}/${{NAME}}.readCountInPeaks.log
            featureCounts -p -a {params.peaks_dir}${{NAME}}.saf -F SAF -o {params.build_summary_directory}/Counts/subread_featureCounts_output/${{NAME}}/${{NAME}}.readCountInPeaks 03_BAM_{DUP}/${{NAME}}_mapQ*_woMT_{DUP}_shifted_sorted.bam 2> ${{FEATURECOUNTSLOG}}
            rm {params.peaks_dir}${{NAME}}.saf
            frip=$(grep 'Successfully assigned alignments' ${{FEATURECOUNTSLOG}} | sed -e 's/.*(//' | sed 's/%.*$//')
            fripScore=$(echo $frip | sed 's/\\..*$//')
            fripLabel=$(if [ $fripScore -ge {params.FRiP_threshold} ]; then echo 'good'; else echo 'bad'; fi)


            printf ${{NAME}}'\\t'$R1'\\t'$R2'\\t'$TOTAL'\\t'$unfilteredBAM'\\t'$percMT'\\t'$dedupBAM'\\t'$dedupREADS'\\t'$shiftedBAM'\\t'$lossReads'\\t'$peaks'\\t'$frip'\\t'$fripLabel'\\n' >> {output.summary_file}
        done

        uniq -u {output.summary_file} > {output.summary_file_temp}
        (head -n 1 {output.summary_file_temp} && tail -n +2 {output.summary_file_temp} | sort -k 1) > {output.summary_file}
        """


# ----------------------------------------------------------------------------------------
# Generation of samples PCA and Heatmap
rule L1_multiBigwigSummary:
    input:
        norm_bw = ancient(expand(os.path.join("04_Normalization/normalized_bigWigs/", "{sample}_mapQ{MAPQ}_woMT_{dup}_shifted_normalized_bs{binSize}.bw"), sample=SAMPLENAMES,  dup=DUP, MAPQ=str(config["mapQ_cutoff"]), binSize=str(BINS)))
    output:
        matrix = os.path.join(SUMMARYDIR, "Sample_comparisons/multiBigWigSummary_matrix_allSamples.npz")
    params:
        make_directory = os.path.dirname(''.join([SUMMARYDIR, "/Sample_comparisons/"])),
        labels = ' '.join(SAMPLENAMES),
        window = config["binning_window_size"],
        blacklist = config["blacklist_file"],
        CPUs = config["multiBigwigSummary_threads"]
    threads:
        config["multiBigwigSummary_threads"]
    shell:
        """
        printf '\033[1;36mComparing the whole signal among samples...\\n\033[0m'

        mkdir -p {params.make_directory}

        multiBigwigSummary bins -p {params.CPUs} -b {input.norm_bw} --labels {params.labels} --binSize {params.window} --blackListFileName {params.blacklist} -o {output.matrix}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Generation of samples PCA and Heatmap
rule L2_PCA_and_samples_correlation:
    input:
        matrix = os.path.join(SUMMARYDIR, "Sample_comparisons/multiBigWigSummary_matrix_allSamples.npz")
    output:
        PCA = os.path.join(SUMMARYDIR, "Sample_comparisons/PCA_on_BigWigs_wholeGenome.pdf"),
        hetamap_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_spearmanMethod.pdf"),
        hetamap_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_heatmap_on_BigWigs_wholeGenome_pearsonMethod.pdf"),
        scatterplot_spearman = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_spearmanMethod.pdf"),
        scatterplot_pearson = os.path.join(SUMMARYDIR, "Sample_comparisons/Sample_correlation/Correlation_scatterplot_on_BigWigs_wholeGenome_pearsonMethod.pdf")
    params:
        make_directory = os.path.dirname(''.join([SUMMARYDIR, "Sample_comparisons/Sample_correlation/"])),
        heatmap_color = config["heatmap_color"]
    shell:
        """
        printf '\033[1;36mPlotting the correlation and variability of the whole signal among samples...\\n\033[0m'

        mkdir -p {params.make_directory}

        printf '\033[1;36m    - plotting PCA...\\n\033[0m'
        plotPCA -in {input.matrix} -o {output.PCA} -T 'PCA on BigWigs (whole genome)' --plotFileFormat 'pdf'


        printf '\033[1;36m    - plotting Spearman correlation heatmap...\\n\033[0m'
        plotCorrelation -in {input.matrix} \
        --corMethod spearman \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Spearman correlation of BigWigs" \
        --whatToPlot heatmap \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.hetamap_spearman}

        printf '\033[1;36m    - plotting Pearson correlation heatmap...\\n\033[0m'
        plotCorrelation -in {input.matrix} \
        --corMethod pearson \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson correlation of BigWigs" \
        --whatToPlot heatmap \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.hetamap_pearson}



        printf '\033[1;36m    - plotting Spearman correlation scatterplot...\\n\033[0m'
        plotCorrelation -in {input.matrix} \
        --corMethod spearman \
        --log1p \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Spearman correlation of BigWigs - ln values" \
        --whatToPlot scatterplot \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.scatterplot_spearman}

        printf '\033[1;36m    - plotting Pearson correlation scatterplot...\\n\033[0m'
        plotCorrelation -in {input.matrix} \
        --corMethod pearson \
        --log1p \
        --skipZeros \
        --removeOutliers \
        --plotTitle "Pearson correlation of BigWigs - ln values" \
        --whatToPlot scatterplot \
        --colorMap {params.heatmap_color} \
        --plotNumbers \
        --plotFileFormat 'pdf' \
        -o {output.scatterplot_pearson}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Generation of Lorenz curves
rule L3_Lorenz_curve:
    input:
        dedup_BAM_shifted_sorted = ancient(expand(os.path.join("03_BAM_{dup}/", "{sample}_mapQ{MAPQ}_woMT_{dup}_shifted_sorted.bam"), sample=SAMPLENAMES, dup=DUP, MAPQ=str(config["mapQ_cutoff"])))
    output:
        lorenz_plot = os.path.join(SUMMARYDIR, "Lorenz_curve_deeptools.plotFingreprint_allSamples.pdf")
    params:
        heatmap_color = config["heatmap_color"],
        all_bams = ' '.join(expand(os.path.join("03_BAM_{dup}/", "{sample}_mapQ{MAPQ}_woMT_{dup}_shifted_sorted.bam"), sample=SAMPLENAMES, dup=DUP, MAPQ=str(config["mapQ_cutoff"]))),
        labels = ' '.join(SAMPLENAMES),
        blacklist = config["blacklist_file"],
        binSize = config["plotFingerprint_binSize"],
        sampledRegions = config["plotFingerprint_sampledRegions"],
        extra_params = config["plotFingerprint_extra_parameters"],
        CPUs = config["multiBigwigSummary_threads"]
    threads:
        config["plotFingerprint_threads"]
    shell:
        """
        printf '\033[1;36mPlotting the Lorenz curves-Fingerprint for all samples...\\n\033[0m'

        plotFingerprint \
        --bamfiles {params.all_bams} \
        --plotFile {output.lorenz_plot} \
        --labels {params.labels} \
        --blackListFileName {params.blacklist} \
        --binSize {params.binSize} \
        --numberOfSamples {params.sampledRegions}\
        -p {params.CPUs} {params.extra_params}
        """
# ----------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------
# Absolute peaks file and relative matrix score generation for SEACR
rule M_all_peaks_file_and_score_matrix:
    input:
        norm_bw = ancient(expand(os.path.join("04_Normalization/normalized_bigWigs/", "{sample}_mapQ{MAPQ}_woMT_{dup}_shifted_normalized_bs{binSize}.bw"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]), dup=DUP, binSize=str(BINS))),
        peaks_file = (expand(os.path.join(PEAKSDIR, "{sample}_mapQ{MAPQ}_woMT_{dup}_shifted_seacr_top{AUC}.stringent.bed"), sample=SAMPLENAMES, MAPQ=str(config["mapQ_cutoff"]), dup=DUP, AUC=str(config["AUC_fration"])))
    output:
        concatenation_bed = temp(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/temp_all_samples_peaks_concatenation.bed")),
        concatenation_bed_sorted = temp(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/temp_all_samples_peaks_concatenation_sorted.bed")),
        concatenation_bed_collapsed = temp(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/temp_all_samples_peaks_concatenation_collapsed.bed")),
        concatenation_bed_collapsed_sorted = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/all_samples_peaks_concatenation_collapsed_sorted.bed"),
        score_matrix_peaks = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_SEACR.npz"),
        score_matrix_peaks_table = os.path.join(SUMMARYDIR, ''.join(["Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_table_", PEAKCALLER, ".tsv"]))
    params:
        make_directory = os.path.dirname(''.join([SUMMARYDIR, "Sample_comparisons/Peak_comparison/"])),
        peaks_dir = PEAKSDIR,
        labels = ' '.join(SAMPLENAMES),
        blacklist = config["blacklist_file"],
        CPUs = config["multiBigwigSummary_threads"]
    threads:
        config["multiBigwigSummary_threads"]
    shell:
        """
        printf '\033[1;36mGenerating a file result of the merge of all the SEACR peaks...\\n\033[0m'

        mkdir -p {params.make_directory}

        cat {params.peaks_dir}*.*bed >> {output.concatenation_bed}
        sort -V -k1,1 -k2,2 -k5,5 {output.concatenation_bed} | cut -f 1,2,3 > {output.concatenation_bed_sorted}

        bedtools merge -i {output.concatenation_bed_sorted} | uniq > {output.concatenation_bed_collapsed}
        sort -V -k1,1 -k2,2 -k5,5 {output.concatenation_bed_collapsed} > {output.concatenation_bed_collapsed_sorted}


        printf '\033[1;36mComputing the score matrix for all the SEACR peaks per each sample...\\n\033[0m'

        multiBigwigSummary BED-file \
        -p {params.CPUs} \
        -b {input.norm_bw} \
        -o {output.score_matrix_peaks} \
        --BED {output.concatenation_bed_collapsed_sorted} \
        --blackListFileName {params.blacklist} \
        --outRawCounts {output.score_matrix_peaks_table} \
        --labels {params.labels}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Compute peaks z-scores and plot heatmap for SEACR
rule N_peaks_zScores_and_heatmap:
    input:
        score_matrix_peaks_table = ancient(os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/peaks_score_matrix_all_samples_table_{PEAKCALLER}.tsv"))
    output:
        rawScores_hetamap = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Heatmaps/Heatmap_on_log1p.rawScores_for_{PEAKCALLER}.peaks_union_population.pdf")
    params:
        make_directory = os.path.dirname(''.join([SUMMARYDIR, "Sample_comparisons/Peak_comparison/Heatmaps/"])),
        heatmap_color = config["heatmap_color"],
        zScore_heatmap_color = config["zScore_heatmap_color"],
        heatmap_basename_rawScores = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Heatmaps/Heatmap_on_log1p.rawScores_for_{PEAKCALLER}.peaks_union_population"),
        heatmap_basename_zScore = os.path.join(SUMMARYDIR, "Sample_comparisons/Peak_comparison/Heatmaps/Heatmap_on_zScores_for_{PEAKCALLER}.peaks_union_population"),
        n_samples = len(SAMPLENAMES),
        peak_caller = PEAKCALLER
    threads:
        config["multiBigwigSummary_threads"]
    run:
        # Messege
        shell("printf '\033[1;36mLoading {PEAKCALLER} peak raw score table...\\n\033[0m'")
        shell("mkdir -p {params.make_directory}")

        # Import multiBigWig summary table
        import pandas as pd
        matrix = pd.read_csv(str(input.score_matrix_peaks_table),  sep='\s+', engine='python')

        # Use peak coordinates as ID ofr each row
        matrix["peak_ID"] = matrix[matrix.columns[:3]].apply(lambda x: '_'.join(x.dropna().astype(str)),axis=1)
        matrix = matrix[matrix.columns[3:]]
        matrix = matrix.set_index('peak_ID')

        # Conversion in log1p
        import numpy as np
        matrix_log1p = np.log1p(matrix)

        # Required to avoid the use of X11
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt

        # Required to avoid errors in scipy iterations
        import sys
        sys.setrecursionlimit(100000)


        # Messege
        shell("printf '\033[1;36mPlotting the {PEAKCALLER} peak raw score hetamaps...\\n\033[0m'")

        # Generation of the rawScore heatmap and clustering
        from bioinfokit import analys, visuz
        visuz.gene_exp.hmap(df=matrix_log1p,
                            cmap=params.heatmap_color,
                            rowclus=True,
                            colclus=(params.n_samples > 1),
                            figtype="pdf",
                            ylabel=False,
                            figname=str(params.heatmap_basename_rawScores),
                            dim=(params.n_samples, 9), # W * H
                            tickfont=(6, 4))


        ## Generation of the zScore heatmap and clustering
        #from bioinfokit import analys, visuz
        #visuz.gene_exp.hmap(df=matrix,
                            #cmap=params.zScore_heatmap_color,
                            #rowclus=True,
                            #colclus=True,
                            #zscore=0,
                            #figtype ="pdf",
                            #ylabel=False,
                            #figname=str(params.heatmap_basename_zScore),
                            #dim=(6, 12),
                            #tickfont=(6, 4))

        # ---------------- Manual way to compute zScore heatmap Z=(rowScore - rowMean)/rowSD -----------------
        if params.n_samples > 1:
            # Messege
            shell("printf '\033[1;36mComputing the zScores for {PEAKCALLER} peaks...\\n\033[0m'")

            import pandas as pd
            matrix = pd.read_csv(str(input.score_matrix_peaks_table),  sep='\s+', engine='python')
            stat_tb = pd.DataFrame({'rowMeans': matrix[matrix.columns[3:]].mean(axis=1),
                                    'SD':  matrix[matrix.columns[3:]].std(axis=1)})

            scores = []
            for i in list(range(3,len(matrix.columns))):
                scores.append((matrix[matrix.columns[i]] - stat_tb["rowMeans"]) / stat_tb["SD"])


            zScores = pd.DataFrame(scores).transpose()
            zScores.columns = list(matrix.columns)[3:len(matrix.columns)]
            zScores = zScores.fillna(0)
            zScores['peak_ID'] = matrix[matrix.columns[:3]].apply(lambda x: '_'.join(x.dropna().astype(str)),axis=1)
            zScores = zScores.set_index('peak_ID')


            # Messege
            shell("printf '\033[1;36mPlotting the {PEAKCALLER} peak zScores hetamaps...\\n\033[0m'")

            visuz.gene_exp.hmap(df=zScores,
                                cmap=params.zScore_heatmap_color,
                                rowclus=True,
                                colclus=True,
                                figtype ="pdf",
                                ylabel=False,
                                figname=str(params.heatmap_basename_zScore),
                                dim=(params.n_samples, 9),
                                tickfont=(6, 4))
        #------- ------- ------- ------- ------- ------- ------- ------- ------- -------
# ----------------------------------------------------------------------------------------



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>  VARIANT CALLING  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# ----------------------------------------------------------------------------------------
# Create reference genome dictionary
rule O_Make_reference_genome_dictionary:
    input:
        genome = config["genome_fasta"]
    output:
        genome_dict = ''.join([re.sub("[.]([a-z]|[A-Z])*$", "",config["genome_fasta"]),'.dict'])
    shell:
        """
        picard CreateSequenceDictionary REFERENCE={input.genome} OUTPUT={output.genome_dict}
        """
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# run base score recalibration (BSQR) of the bams
rule P_GATK_bam_base_quality_score_recalibration:
    input:
        genome_dict = ancient(''.join([re.sub("[.]([a-z]|[A-Z])*$", "",config["genome_fasta"]),'.dict'])),
        dedup_BAM = ancient(os.path.join("03_BAM_{DUP}/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}.bam"]))),
        dedup_BAM_index = ancient(os.path.join("03_BAM_{DUP}/unshifted_bams/", ''.join(["{SAMPLES}_mapQ", str(config["mapQ_cutoff"]), "_sorted_woMT_{DUP}.bai"])))
    params:
        sample = "{SAMPLES}",
        gatk_directory = GATKDIR,
        genome = config["genome_fasta"],
        dbsnp = config["dbsnp_file"],
        CPUs = config["SAMtools_threads"]
    threads:
        config["SAMtools_threads"]
    output:
        bsqr_table = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_bsqr.table"),
        dedup_BAM_bsqr = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_bsqr.bam"),
        dedup_BAM_bsqr_index = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_bsqr.bam.bai")
    shell:
        """
        printf '\033[1;36m{params.sample}: Base Quality Score Recalibration of the deduplicated unshifted bam...\\n\033[0m'
        mkdir -p {params.gatk_directory}{params.sample}

        gatk --java-options '-Xmx4G' BaseRecalibrator \
        --input {input.dedup_BAM} \
        --known-sites {params.dbsnp} \
        --output {output.bsqr_table} \
        --reference {params.genome}

        gatk --java-options '-Xmx4G' ApplyBQSR \
        -R {params.genome} \
        -I {input.dedup_BAM} \
        --bqsr-recal-file {output.bsqr_table} \
        -O {output.dedup_BAM_bsqr}

        printf '\033[1;36m{params.sample}: Indexing recalibrated bam...\\n\033[0m'
        samtools index -@ {params.CPUs} -b {output.dedup_BAM_bsqr} {output.dedup_BAM_bsqr_index}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# run gatk haplotype caller
rule Q_GATK_haplotype_calling:
    input:
        concatenation_bed_collapsed_sorted = ancient(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_collapsed_sorted.bed")),
        dedup_BAM_bsqr = ancient(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_bsqr.bam")),
        dedup_BAM_bsqr_index = ancient(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_bsqr.bam.bai"))
    params:
        genome = config["genome_fasta"],
        sample = "{SAMPLES}",
        to_copy_bed = os.path.join(GATKDIR, "all_samples_peaks_concatenation_collapsed_sorted.bed")
    output:
        gvcf = temp(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk.g.vcf.gz"))
    shell:
        """
        cp {input.concatenation_bed_collapsed_sorted} {params.to_copy_bed}

        printf '\033[1;36m{params.sample}: GATK Haplotype calling...\\n\033[0m'

        gatk --java-options '-Xmx4g' HaplotypeCaller \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -I {input.dedup_BAM_bsqr} \
        -O {output.gvcf} \
        -ERC GVCF \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# correct the genotypes that come out of haplotype caller
rule R_GATK_haplotype_calling_correction:
    input:
        concatenation_bed_collapsed_sorted = ancient(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_collapsed_sorted.bed")),
        gvcf = ancient(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk.g.vcf.gz"))
    params:
        sample = "{SAMPLES}",
        genome = config["genome_fasta"]
    output:
        vcf = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk.vcf.gz")
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK Haplotype call correction...\\n\033[0m'

        gatk --java-options '-Xmx4g' GenotypeGVCFs \
        --include-non-variant-sites \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.gvcf} \
        -O {output.vcf}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Select SNPs
rule S_GATK_call_SNPs:
    input:
        concatenation_bed_collapsed_sorted = ancient(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_collapsed_sorted.bed")),
        vcf = ancient(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk.vcf.gz"))
    params:
        sample = "{SAMPLES}",
        genome = config["genome_fasta"]
    output:
        snp = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk-snp.vcf")
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK SNP calling...\\n\033[0m'

        gatk --java-options "-Xmx4g" SelectVariants \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.vcf} \
        --select-type SNP \
        --select-type NO_VARIATION \
        --select-type-to-exclude INDEL \
        --select-type-to-exclude MIXED \
        --select-type-to-exclude SYMBOLIC \
        --select-type-to-exclude MNP \
        -O {output.snp}
        """
# ----------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------
# Call indels
rule T_GATK_call_indels:
    input:
        concatenation_bed_collapsed_sorted = ancient(os.path.join(SUMMARYDIR, "temp_all_samples_peaks_concatenation_collapsed_sorted.bed")),
        vcf = ancient(os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk.vcf.gz"))
    params:
        sample = "{SAMPLES}",
        genome = config["genome_fasta"]
    output:
        indels = os.path.join(GATKDIR, "{SAMPLES}/{SAMPLES}_gatk-indel.vcf")
    shell:
        """
        printf '\033[1;36m{params.sample}: GATK Indel calling...\\n\033[0m'

        gatk --java-options "-Xmx4g" SelectVariants \
        -L {input.concatenation_bed_collapsed_sorted} \
        -R {params.genome} \
        -V {input.vcf} \
        --select-type INDEL \
        --select-type NO_VARIATION \
        --select-type-to-exclude SNP \
        --select-type-to-exclude MIXED \
        --select-type-to-exclude SYMBOLIC \
        --select-type-to-exclude MNP \
        -O {output.indels}
        """
# ----------------------------------------------------------------------------------------





#    gatk-indel:
# cols:
#    info: [QD,FS,MQ,AC,ExcessHet]
#    format: [DP,GQ]


#snp_filter: ['gatk-snp~DP>10']
#indel_filter: ['gatk-indel~DP>10']
# ----------------------------------------------------------------------------------------
