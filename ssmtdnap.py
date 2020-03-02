"""
example:
    #  check run options
    python ssmtdnap.py


    # initial config generation
    python ssmtdnap.py \\
        --project_root ~/projects/KAZ_WG/KAZ_WG_mtDNA \\
        --fastq_dirs_list ~/icebox/fastq_gz/KAZ_WG/ \\
        --sample_delimiter . \\
        --fastq_extension .fastq.gz \\
        --R1_fastq_extension .R1.fastq.gz \\
        --R2_fastq_extension .R2.fastq.gz \\
        --script_dir_name scripts

    # optional: generate config with --add_tokens option
    # to add tokens to each step, for rerun from last failed step


    # precise config tuning

    # scripts generation
    python ssmtdnap.py -j ~/projects/KAZ_WG/KAZ_WG_hg19/scripts/default_settings.json

    # submit scripts to workload manager

    cd ~/projects/KAZ_WG/KAZ_WG_hg19/scripts/; for i in $( ls *.ss.sh ); do qsub $i; done
"""
__VERSION__ = "0.0.4"

import os
import argparse
import json

from collections import defaultdict


__NOT_READY__ = "NOT_READY"
__READY__ = "READY"
__ALMOST_READY__ = "ALMOST_READY"


def main():
    settings = parse_arguments_to_settings()
    if settings["ready"] == __ALMOST_READY__:
        save_settings(settings)
    elif settings["ready"] == __READY__:
        run_pipeline(settings)
    else:
        print(__doc__)


def parse_arguments_to_settings():
    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--settings_json", default=None, required=False)
    parser.add_argument("--project_root", default=None, required=False)
    parser.add_argument("--fastq_dirs_list", default=[], required=False, nargs="+")
    parser.add_argument("--sample_delimiter", default="_", required=False)
    parser.add_argument("--fastq_extension", default=".fastq.gz", required=False)
    parser.add_argument("--R1_fastq_extension", default=".R1.fastq.gz", required=False)
    parser.add_argument("--R2_fastq_extension", default=".R2.fastq.gz", required=False)
    parser.add_argument("--script_dir_name", default="scripts", required=False)
    parser.add_argument("--run_annovar", action="store_true")
    parser.add_argument("--add_tokens", action="store_true")
    parser.add_argument("--debug", action="store_true")
    #
    args = parser.parse_args()
    if args.settings_json:
        settings = json.load(open(args.settings_json))   # config exist, import it
        __samples_dict__ = load_fastq_samples(settings)  # find all fastq files
        __samples_list__ = settings["samples_list"]  # get list of target samples from config
        settings["samples_dict"] = {  # filter target samples from all fastq
            list_key: dict_value
            for list_key in __samples_list__
            for dict_key, dict_value in __samples_dict__.items()
            if list_key == dict_key or list_key + "_m" == dict_key
        }
        settings["ready"] = __READY__
    elif args.project_root:
        settings = {
            "settings_json": args.settings_json,
            "project_root": args.project_root,
            "fastq_dirs_list": args.fastq_dirs_list,
            "sample_delimiter": args.sample_delimiter,
            "fastq_extension": args.fastq_extension,
            "R1_fastq_extension": args.R1_fastq_extension,
            "R2_fastq_extension": args.R2_fastq_extension,
            "script_dir_name": args.script_dir_name,
            "add_tokens": args.add_tokens,
            "debug": args.debug,
            "ready":__ALMOST_READY__,
        }
    else:
        settings = {
            "ready":__NOT_READY__,
        }
    return settings


def load_fastq_samples(settings):
    fastq_dirs_list = settings["fastq_dirs_list"]
    sample_delimiter = settings["sample_delimiter"]
    fastq_extension = settings["fastq_extension"]
    R1_fastq_extension = settings["R1_fastq_extension"]
    R2_fastq_extension = settings["R2_fastq_extension"]
    #
    res = defaultdict(lambda: defaultdict(str))
    for fastq in get_files_generator(fastq_dirs_list, fastq_extension):
        sample = os.path.basename(fastq).split(sample_delimiter)[0]
        if fastq.endswith(R1_fastq_extension):
            res[sample]["read1"] = fastq
        elif fastq.endswith(R2_fastq_extension):
            res[sample]["read2"] = fastq
    res = {
        key: value
        for key, value in res.items()
        if key + "_m" not in res
    }
    return res


def get_files_generator(dirs_list, extension=""):
    for path in dirs_list:
        for data_file in os.listdir(path):
            if data_file:
                data_path = os.path.join(path, data_file)
                if os.path.isfile(data_path) and data_path.endswith(extension):
                    yield data_path
                elif os.path.isdir(data_path):
                    yield from get_files_generator([data_path], extension)


def save_settings(settings):
    settings["project_script_dir"] = os.path.join(
        settings["project_root"],
        settings["script_dir_name"],
    )
    settings.update(get_default_settings(settings))
    mkdir(settings["project_script_dir"])
    json_file = os.path.join(
        settings["project_script_dir"],
        "default_settings.json",
    )
    with open(json_file, "w") as f:
        default_settings_str = json.dumps(settings, indent=4, sort_keys=True)
        f.write(default_settings_str)
    # debug print
    print(f"# ls  {settings['project_script_dir']}")
    print(f"# more {json_file}")
    print(f"# nano {json_file}")
    print(f"# python ssmtdnap.py -j {json_file}")


def get_default_settings(d):
    default_settings_dict = {
        "number_of_threads": "4",
        "fastq_dirs_list": d["fastq_dirs_list"],
        "sample_delimiter": d["sample_delimiter"],
        "fastq_extension": d["fastq_extension"],
        "R1_fastq_extension": d["R1_fastq_extension"],
        "R2_fastq_extension": d["R2_fastq_extension"],
        "samples_list": sorted(load_fastq_samples(d)) if  d["fastq_dirs_list"] else [],
        "project_root": d["project_root"],
        "tools": {
            "fastqc": "",
            "bwa": "/home/Pipeline/bwa/bwa-0.7.12/bwa",
            "samtools": "/home/Pipeline/samtools/samtools-1.2/samtools",
            "bcftools": "/home/Pipeline/bcftools/bcftools-1.2/bcftools",
            "java": "/home/Pipeline/java/jdk1.8.0_101/bin/java",
            "picard": "/home/Pipeline/java/jdk1.8.0_101/bin/java -jar /home/Pipeline/picard/picard_1_130/picard.jar",
            "gatk": "/home/Pipeline/java/jdk1.8.0_101/bin/java -jar /home/Pipeline/gatk/GATK_3_8_1/GenomeAnalysisTK.jar",
            "vcf_concat": "/home/Pipeline/vcftools/vcftools_0.1.13/bin/vcf-concat",
            "vcf_sort": "/home/Pipeline/vcftools/vcftools_0.1.13/bin/vcf-sort",
            "vcf_merge:": "/home/Pipeline/vcftools/vcftools_0.1.13/bin/vcf-merge",
            "bgzip": "/home/Pipeline/tabix/tabix-0.2.6/bgzip",
            "tabix": "/home/Pipeline/tabix/tabix-0.2.6/tabix",
            "#annovar": "",
            "#snpedia": "",
        },
        "databases": {
            "ref": "/home/PublicData/broadinstitute/2.8/hg19/ucsc.hg19.ref/ucsc.hg19.fasta",
            "ref_ucsc_hg19": "/home/PublicData/broadinstitute/2.8/hg19/ucsc.hg19.ref/ucsc.hg19.fasta",
            "ref_mtDNA": "/home/PublicData/h.sapiens_mtDNA/HS_mtDNA.fa",
            "gold_indel": "/home/PublicData/broadinstitute/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf",
            "oneKG_indel": "/home/PublicData/broadinstitute/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf",
            "oneKG_snp": "/home/PublicData/broadinstitute/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf",
            "dbsnp": "/home/PublicData/broadinstitute/2.8/hg19/dbsnp_138.hg19.vcf",
            "#dbsnp": "/home/PublicData/dbsnp/human_9606_b150_GRCh37p13/All_20170710.liftover.2.17.11.vcf",
            "onmi_snp": "/home/PublicData/broadinstitute/2.8/hg19/1000G_omni2.5.hg19.sites.vcf",
            "hapmap_snp": "/home/PublicData/broadinstitute/2.8/hg19/hapmap_3.3.hg19.sites.vcf",
            "target_region": "/home/PublicData/Agilent_v4_71m_reduced.bed",
            "#ensembl_ref_dir": "/home/PublicData/ensembl_GRCh37_75/",
            "#ensembl_ref_fa": "/home/PublicData/ensembl_GRCh37_75/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa",
            "#ensembl_ref_gtf": "/home/PublicData/ensembl_GRCh37_75/Homo_sapiens.GRCh37.75.gtf"
        },
    }
    return default_settings_dict


def mkdir(dir_name):
    os.makedirs(dir_name, exist_ok=True)


def run_pipeline(settings):
    for sample in sorted(settings["samples_dict"]):
        sample_settings = settings
        sample_settings["sample"] = sample
        sample_settings = get_settings_for_SSAP(sample_settings)
        cmd_list = get_cmd_list_for_SSAP(sample_settings)
        write_cmd_list_to_file(sample_settings, cmd_list)
    else:
        run_ms_merge(settings)
    # debug print
    print(f"# ls {settings['project_script_dir']}")
    print(f"# cd {settings['project_script_dir']}; for i in $( ls *.ss.sh ); do echo $i; done" )
    print(f"# cd {settings['project_script_dir']}; for i in $( ls *.ss.sh ); do qsub $i; done" )


def run_ms_merge(settings):
    project_root = settings["project_root"]
    project_name = os.path.basename(project_root)
    num = len(settings["samples_dict"])
    sample = f"{project_name}_{num}"
    sample_dir = os.path.join(project_root, sample)
    samples_string = " ".join(
        os.path.join(project_root, s, s + ".call.vcf.gz")
        for s in sorted(settings["samples_dict"])
    )
    ms_mtDNA_merged_vcf = os.path.join(sample_dir, sample + ".merged.vcf")
    cmd_list = [
        f"mkdir -p {sample_dir}",
        f"vcf-merge {samples_string} > {ms_mtDNA_merged_vcf}",
    ]
    script_file = os.path.join(settings["project_script_dir"], sample + ".ms.sh")
    with open(script_file, "w") as f:
        f.write("#!/bin/bash\n\n")
        for cmd in cmd_list:
            new_line = cmd + "\n\n"
            f.write(new_line)
    # debug print
    print(f"# cd {settings['project_script_dir']}; for i in $( ls *.ms.sh ); do qsub $i; done" )


def get_settings_for_SSAP(sample_dict):
    sample = sample_dict["sample"]
    sample_dir = os.path.join(sample_dict["project_root"], sample)
    sample_path_prefix = os.path.join(sample_dir, sample)
    _dict = {
        "sample": sample,
        "sample_dir": sample_dir,

        "project_script_dir" : sample_dict["project_script_dir"],
        "number_of_threads": sample_dict["number_of_threads"],

        "bwa": sample_dict["tools"]["bwa"],
        "samtools": sample_dict["tools"]["samtools"],
        "bcftools": sample_dict["tools"]["bcftools"],
        "bgzip": sample_dict["tools"]["bgzip"],
        "tabix": sample_dict["tools"]["tabix"],

        "ref": sample_dict["databases"]["ref_mtDNA"],

        "read1": sample_dict["samples_dict"][sample]["read1"],
        "read2": sample_dict["samples_dict"][sample]["read2"],

        "sam": sample_path_prefix + ".mem.sam",
        "bam": sample_path_prefix + ".view.bam",
        "sorted_bam": sample_path_prefix + ".sorted.bam",
        "sorted_tmp": sample_path_prefix + ".sorted.tmp",

        "call_vcf": sample_path_prefix + ".call.vcf",
        "stats_txt": sample_path_prefix + ".stats.txt",
        "call_vcf_gz": sample_path_prefix + ".call.vcf.gz",

        "add_tokens": sample_dict["add_tokens"],
        "debug": sample_dict["debug"],
    }
    return _dict

###############################################################################
def get_cmd_list_for_SSAP(sample_settings):
    if sample_settings["debug"]:
        print("\n\nget_cmd_list_for_SSAP", sample_settings)
    cmd_list = [
        get_mkdir_cmd(sample_settings),
        get_cmd_bwa_mem_sam(sample_settings),
        get_cmd_samtools_view_bam(sample_settings),
        get_cmd_samtools_sort_bam(sample_settings),
        get_cmd_samtools_mpileup_bcftools_call(sample_settings),
        get_cmd_bcftools_filter_bcftools_stats(sample_settings),
        get_cmd_bgzip_vcf(sample_settings),
        get_cmd_tabix_gz(sample_settings),
    ]
    return cmd_list


###############################################################################
def get_mkdir_cmd(d):
    return "mkdir -p {sample_dir}".format(**d)


###############################################################################
def reduce_spaces_and_newlines(s):
    s = s.replace("\n", " ")
    s = " ".join([i for i in s.split(" ") if i])
    return s

def get_cmd(d):
    d["token"] = "{sample_dir}/token.{sample}.{token_suffix}".format(**d)
    d["flags"] = " && ".join([" [ -f {:s} ] ".format(i) for i in d["files_list"]]) + " && [ ! -f {token} ] ".format(**d)
    if d["add_tokens"]:
        cmd = """
            {flags} &&
            dt1=`date +%y%m%d_%H%M%S` && echo $dt1 {token} &&
            {main_cmd} &&
            du {out_file} > {out_file}.$dt1.du &&
            md5sum {out_file} > {out_file}.$dt1.md5 &&
            dt2=`date +%y%m%d_%H%M%S` &&
            echo $dt1 $dt2 > {token} ||
            echo "TOKEN SKIPPED {token}"
            """.format(**d)
    else:
        cmd = "{main_cmd}".format(**d)
    return reduce_spaces_and_newlines(cmd)


def clear_after_competion(d):
    d["token_suffix"] = "vcftools_sort_vcf"
    d["token"] = "{sample_dir}/token.{sample}.{token_suffix}".format(**d)
    if d["add_tokens"]:
        cmd = """ [ -f {token} ] && [ -f {vcftools_sorted_vcf} ] &&
            rm -f
            {sam} {bam} {sorted_bam} {ARRG_bam}
            {IR_bam}
            {gatk_AR_SNP_vcf} {gatk_AR_INDEL_vcf}
            {gatk_SV_SNP_raw_vcf}  {gatk_SV_INDEL_raw_vcf}
            {vcftools_concat_vcf}
            """.format(**d)
            # {gatk_SV_SNP_fil_vcf} {gatk_SV_INDEL_fil_vcf} # for vcftools concatenate - sometimes library doesnt loaded well
    else:
        cmd = ""
    return reduce_spaces_and_newlines(cmd)

###############################################################################
def get_cmd_bwa_mem_sam(d):
    d["token_suffix"] = "bwa_read1_read2_2_sam"
    d["files_list"] = [d["read1"], d["read2"]]
    d["out_file"] = d["sam"]
    d["main_cmd"] = bash_bwa_mem_sam(d)
    return get_cmd(d)

def bash_bwa_mem_sam(d):
    return """{bwa} mem -M -t {number_of_threads} {ref} {read1} {read2} > {sam}""".format(**d)


###############################################################################
def get_cmd_samtools_view_bam(d):
    d["files_list"] = [d["sam"]]
    d["token_suffix"] = "samtools_sam2bam"
    d["out_file"] = d["bam"]
    d["main_cmd"] = bash_samtools_view_bam(d)
    return get_cmd(d)


def bash_samtools_view_bam(d):
    return """{samtools} view -bT {ref} {sam} > {bam}""".format(**d)


###############################################################################
def get_cmd_samtools_sort_bam(d):
    d["files_list"] = [d["bam"]]
    d["token_suffix"] = "samtools_bam_2_sorted_bam"
    d["out_file"] = d["sorted_bam"]
    d["main_cmd"] = bash_samtools_sort_bam(d)
    return get_cmd(d)


def bash_samtools_sort_bam(d):
    return """{samtools} sort -l 9 -O bam -T {sorted_tmp} {bam} > {sorted_bam}""".format(**d)


###############################################################################
def get_cmd_samtools_mpileup_bcftools_call(d):
    d["files_list"] = [d["sorted_bam"]]
    d["token_suffix"] = "get_cmd_samtools_mpileup_bcftools_call"
    d["out_file"] = d["call_vcf"]
    d["main_cmd"] = bash_samtools_mpileup_bcftools_call(d)
    return get_cmd(d)


def bash_samtools_mpileup_bcftools_call(d):
    return """ {samtools} mpileup -uf {ref} {sorted_bam} | {bcftools} call -m -Ov > {call_vcf} """.format(**d)


###############################################################################
def get_cmd_bcftools_filter_bcftools_stats(d):
    d["files_list"] = [d["call_vcf"]]
    d["token_suffix"] = "get_cmd_bcftools_filter_bcftools_stats"
    d["out_file"] = d["stats_txt"]
    d["main_cmd"] = bash_bcftools_filter_bcftools_stats(d)
    return get_cmd(d)


def bash_bcftools_filter_bcftools_stats(d):
    return """ {bcftools} filter -i"%QUAL>20" {call_vcf} | {bcftools} stats > {stats_txt}  """.format(**d)


###############################################################################
def get_cmd_bgzip_vcf(d):
    d["files_list"] = [d["call_vcf"]]
    d["token_suffix"] = "get_cmd_bgzip_vcf"
    d["out_file"] = d["call_vcf_gz"]
    d["main_cmd"] = bash_bgzip_vcf(d)
    return get_cmd(d)


def bash_bgzip_vcf(d):
    return """ {bgzip} -c {call_vcf} > {call_vcf_gz} """.format(**d)


###############################################################################
def get_cmd_tabix_gz(d):
    d["files_list"] = [d["call_vcf_gz"]]
    d["token_suffix"] = "get_cmd_tabix_gz"
    d["out_file"] = d["call_vcf_gz"] + ".tbi"
    d["main_cmd"] = bash_tabix_gz(d)
    return get_cmd(d)


def bash_tabix_gz(d):
    return """ {tabix} {call_vcf_gz} """.format(**d)


###############################################################################
def get_cmd_list_for_mtDNA_MSAP_sample_dict(_in_dict_):
    cmd_list = [
        get_mkdir_cmd(_in_dict_),
        "vcf-merge {samples_string} > {MS_mtDNA_merged_vcf}".format(**_in_dict_),
    ]
    return cmd_list


###############################################################################
def write_cmd_list_to_file(sample_settings, cmd_list):
    script_file = os.path.join(
        sample_settings["project_script_dir"],
        sample_settings["sample"] + ".ss.sh",
    )
    with open(script_file, "w") as f:
        f.write("#!/bin/bash\n\n")
        for cmd in cmd_list:
            new_line = cmd + "\n\n"
            f.write(new_line)


###############################################################################
if __name__ == "__main__":
    main()
