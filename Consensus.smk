import sys
import os

SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, "scripts"))
import tools

JULIA_COMMAND = f"julia --startup-file=no --project={SNAKEDIR}"

######################################################
# GLOBAL CONSTANTS
######################################################
def abspath(x):
    return os.path.abspath(os.path.expanduser(x))

if "readdir" not in config:
    raise KeyError("You must supply read path: '--config readdir=path/to/reads'")

READDIR = abspath(config["readdir"])

# For some reason it all fucks up if the read directory is a child of the running
# directory, so we check that here.
_dir = READDIR
while _dir != os.path.dirname(_dir):
    if _dir == os.getcwd():
        raise ValueError("Error: Read path cannot be child directory of running directory")
    _dir = os.path.dirname(_dir)

if "platform" not in config or config["platform"] not in ["illumina", "nanopore"]:
    raise KeyError("You must supply platform: '--config platform=[\"illumina/nanopore\"]")

IS_NANOPORE = config["platform"] == "nanopore"
IS_ILLUMINA = config["platform"] == "illumina"
assert IS_NANOPORE ^ IS_ILLUMINA # only one must be true at a time

if IS_NANOPORE:
    if "pore" not in config:
        raise KeyError("On nanopore platform, you must supply pore: --config pore=9/10")
    if str(config["pore"]) not in ["9", "10"]:
        raise ValueError(f"Pore must be 9 or 10, not '{config['pore']}'")
    PORE = str(config["pore"])

# Pass in the reference directory
if "ref" not in config:
    raise KeyError("You must supply reference directory: '--config ref=/path/to/ref'")

REFDIR = abspath(config["ref"])
if not os.path.isdir(REFDIR):
    raise NotADirectoryError(REFDIR)

if IS_ILLUMINA:
    READS = tools.get_read_pairs(READDIR)
else:
    READS = tools.get_nanopore_reads(READDIR)

SAMPLENAMES = sorted(READS.keys())

REFOUTDIR = os.path.join(REFDIR, "refout")

# We have to create this directories either in a rule or outside the DAG to not
# mess up the DAG.
if not os.path.isdir(REFOUTDIR):
    os.mkdir(REFOUTDIR)

SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

######################################################
# Start of pipeline
######################################################

def done_input(wildcards):
    # Add report and the commit
    inputs = ["report.txt"]

    for samplename in SAMPLENAMES:
        inputs.append(f"consensus/{samplename}/consensus.fna")

    return inputs

rule all:
    input: done_input
    output: "commit.txt"
    params: SNAKEDIR
    shell: "git -C {params} rev-parse --short HEAD > {output} && cp {params}/copy_readme.md README.md"

#################################
# REFERENCE-ONLY PART OF PIPELINE
#################################
rule index_ref:
    input: REFDIR + "/refs.fna"
    output:
        comp=REFOUTDIR + "/refs.comp.b",
        name=REFOUTDIR + "/refs.name",
        length=REFOUTDIR + "/refs.length.b",
        seq=REFOUTDIR + "/refs.seq.b"
    params:
        outpath=REFOUTDIR + "/refs",
    log: "tmp/log/kma_ref.log"
    
    shell: "kma index -nbp -k 12 -Sparse - -i {input} -o {params.outpath} 2> {log}"

############################
# CONSENSUS PART OF PIPELINE
############################
rule gzip:
    input: "{base}"
    output: "{base}.gz"
    shell: "gzip -k {input}"

if IS_ILLUMINA:
    ruleorder: second_kma_map > first_kma_map > gzip

    rule fastp:
        input:
            fw=lambda wildcards: READS[wildcards.samplename][0],
            rv=lambda wildcards: READS[wildcards.samplename][1],
        output:
            fw=temp('tmp/trim/{samplename}/fw.fq'),
            rv=temp('tmp/trim/{samplename}/rv.fq'),
            html='tmp/trim/{samplename}/report.html',
            json='tmp/trim/{samplename}/report.json'
        log: "tmp/log/fastp/{samplename}.log"
        threads: 2
        shell:
            'fastp -i {input.fw} -I {input.rv} '
            '-o {output.fw} -O {output.rv} --html {output.html} --json {output.json} '
            '--disable_adapter_trimming --trim_poly_g --cut_tail --cut_front --low_complexity_filter '
            '--complexity_threshold 50 --thread {threads} 2> {log}'

    rule map_best_template:
        input:
            fw=rules.fastp.output.fw,
            rv=rules.fastp.output.rv,
            index=rules.index_ref.output
        output: "tmp/aln/{samplename}/sparse.spa"
        params:
            db=REFOUTDIR + "/refs",
            outbase="tmp/aln/{samplename}/sparse"
        threads: 2
        log: "tmp/log/aln/{samplename}.initial.log"
        shell:
            # Here, not sure if I should sort by template cov (-ss c) or not.
            # Pro: We mostly care about having a fully-covered template, not a partially
            # covered with high depth at the covered areas:
            # Con: Majority vote will win away, so it'll just fuck up if we pick a
            # uniformly, but low covered reference anyway
            "kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
            "-ss c -t {threads} -Sparse 2> {log}"        

elif IS_NANOPORE:
    rule fastp:
        input: lambda wildcards: READS[wildcards.samplename]
        output:
            reads=temp('tmp/trim/{samplename}/reads.fq'),
            html='tmp/trim/{samplename}/report.html',
            json='tmp/trim/{samplename}/report.json'
        log: "tmp/log/fastp/{samplename}.log"
        threads: 2
        shell:
            'fastp -i {input} -o {output.reads} --html {output.html} '
            '--json {output.json} --disable_adapter_trimming  --disable_trim_poly_g '
            '--cut_window_size 10 --cut_mean_quality 10 --cut_tail --cut_front --low_complexity_filter  '
            '--complexity_threshold 50 --length_limit 2400 --length_required 100 '
            '--average_qual 12 --thread {threads} 2> {log}'

    rule map_best_template:
        input:
            reads=rules.fastp.output.reads,
            index=rules.index_ref.output
        output: "tmp/aln/{samplename}/sparse.spa"
        params:
            db=REFOUTDIR + "/refs", # same as index_reffile param
            outbase="tmp/aln/{samplename}/sparse"
        threads: 1
        log: "tmp/log/aln/{samplename}.initial.log"
        shell:
            # See above comment in rule with same name
            "kma -i {input.reads} -o {params.outbase} -t_db {params.db} "
            "-ss c -t {threads} -Sparse 2> {log}"

### Both platforms
rule collect_best_templates:
    input: expand("tmp/aln/{samplename}/sparse.spa", samplename=SAMPLENAMES)
    output: temp(expand("tmp/aln/{samplename}/cat.fna", samplename=SAMPLENAMES))
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/gather_spa.jl",
        refpath=REFDIR
    shell: "{params.juliacmd} {params.scriptpath} tmp/aln {params.refpath}"

rule first_kma_index:
    input: "tmp/aln/{samplename}/cat.fna"
    output:
        comp=temp("tmp/aln/{samplename}/cat.comp.b"),
        name=temp("tmp/aln/{samplename}/cat.name"),
        length=temp("tmp/aln/{samplename}/cat.length.b"),
        seq=temp("tmp/aln/{samplename}/cat.seq.b")
    params:
        t_db="tmp/aln/{samplename}/cat"
    log: "tmp/log/aln/kma1_index_{samplename}.log"
    # The pipeline is very sensitive to the value of k here.
    # Too low means the mapping is excruciatingly slow,
    # too high results in poor mapping quality.
    shell: "kma index -nbp -k 10 -i {input} -o {params.t_db} 2> {log}"

if IS_ILLUMINA:
    rule first_kma_map:
        input:
            fw=rules.fastp.output.fw,
            rv=rules.fastp.output.rv,
            index=rules.first_kma_index.output,
        output:
            res="tmp/aln/{samplename}/kma1.res",
            fsa="tmp/aln/{samplename}/kma1.fsa",
            mat="tmp/aln/{samplename}/kma1.mat.gz",
        params:
            db="tmp/aln/{samplename}/cat",
            outbase="tmp/aln/{samplename}/kma1",
        log: "tmp/log/aln/kma1_map_{samplename}.log"
        threads: 2
        run:
            shell("kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
            "-t {threads} -1t1 -gapopen -5 -nf -matrix 2> {log}")

elif IS_NANOPORE:
    rule first_kma_map:
        input:
            reads=rules.fastp.output.reads,
            index=rules.first_kma_index.output,
        output:
            res="tmp/aln/{samplename}/kma1.res",
            fsa="tmp/aln/{samplename}/kma1.fsa",
            mat="tmp/aln/{samplename}/kma1.mat.gz",
        params:
            db="tmp/aln/{samplename}/cat",
            outbase="tmp/aln/{samplename}/kma1",
        log: "tmp/log/aln/kma1_map_{samplename}.log"
        threads: 2
        run:
            shell("kma -i {input.reads} -o {params.outbase} -t_db {params.db} "
            "-t {threads} -1t1 -bcNano -nf -matrix 2> {log}")

# Both platforms
rule remove_primers:
    input:
        con=rules.first_kma_map.output.fsa,
        primers=f"{REFDIR}/primers.fna"
    output: temp("tmp/aln/{samplename}/cat.trimmed.fna")
    log: "tmp/log/consensus/remove_primers_{samplename}.txt"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/trim_consensus.jl",
        minmatches=4,
        fuzzylen=8,
    shell: """if [ -s {input.primers} ]; then
    {params.juliacmd} {params.scriptpath} {input.primers} \
{input.con} {output} {params.minmatches} {params.fuzzylen} > {log}
else
    cp {input.con} {output}
fi
"""

if IS_ILLUMINA:
    # We now re-map to the created consensus sequence in order to accurately
    # estimate depths and coverage, and get a more reliable assembly seq.
    rule second_kma_index:
        input: rules.remove_primers.output
        output:
            comp=temp("tmp/aln/{samplename}/cat.trimmed.comp.b"),
            name=temp("tmp/aln/{samplename}/cat.trimmed.name"),
            length=temp("tmp/aln/{samplename}/cat.trimmed.length.b"),
            seq=temp("tmp/aln/{samplename}/cat.trimmed.seq.b")
        params:
            t_db="tmp/aln/{samplename}/cat.trimmed"
        log: "tmp/log/aln/kma2_index_{samplename}.log"
        shell: "kma index -nbp -i {input} -o {params.t_db} 2> {log}"

    # And now we KMA map to that index again
    rule second_kma_map:
        input:
            fw=rules.fastp.output.fw,
            rv=rules.fastp.output.rv,
            index=rules.second_kma_index.output,
        output:
            res="tmp/aln/{samplename}/kma2.res",
            fsa="tmp/aln/{samplename}/kma2.fsa",
            mat="tmp/aln/{samplename}/kma2.mat.gz",
        params:
            db="tmp/aln/{samplename}/cat.trimmed",
            outbase="tmp/aln/{samplename}/kma2",
        log: "tmp/log/aln/kma2_map_{samplename}.log"
        threads: 2
        run:
            shell("kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db} "
            "-t {threads} -1t1 -gapopen -5 -nf -matrix 2> {log}")

    rule create_report:
        input:
            matrix=expand("tmp/aln/{samplename}/kma1.mat.gz", samplename=SAMPLENAMES),
            assembly=expand("tmp/aln/{samplename}/kma2.fsa", samplename=SAMPLENAMES),
            res=expand("tmp/aln/{samplename}/kma2.res", samplename=SAMPLENAMES)
        output:
            consensus=expand("consensus/{samplename}/{type}.{nuc}",
                samplename=SAMPLENAMES, type=["consensus", "curated"], nuc=["fna", "faa"]
            ),
            report="report.txt",
        params:
            juliacmd=JULIA_COMMAND,
            scriptpath=f"{SNAKEDIR}/scripts/report.jl",
            refdir=REFDIR
        log: "tmp/log/report.txt"
        threads: workflow.cores
        run:
            shell(f"{params.juliacmd} -t {threads} {params.scriptpath} illumina . {params.refdir} > {log}")

elif IS_NANOPORE:
    rule medaka:
        input: 
            reads=rules.fastp.output.reads,
            draft=rules.remove_primers.output
        output: directory("tmp/aln/{samplename}/medaka")
        log: "tmp/log/aln/medaka_{samplename}.log"
        threads: 2
        params:
            model=lambda wc: "r941_min_high_g360" if PORE == 9 else "r103_min_high_g360"
        shell:
            "medaka_consensus -i {input.reads} -d {input.draft} -o {output} -t {threads} "
            "-m {params.model} 2> {log}"

    rule clean_medaka:
        input: rules.medaka.output
        output: "tmp/aln/{samplename}/moved.txt"
        shell: "rm {input}/*.bam {input}/*.bai {input}/*.hdf && touch {output}"

    rule create_report:
        input:
            assembly=expand("tmp/aln/{samplename}/moved.txt", samplename=SAMPLENAMES),
        output:
            consensus=expand("consensus/{samplename}/{type}.{nuc}",
                samplename=SAMPLENAMES, type=["consensus", "curated"], nuc=["fna", "faa"]
            ),
            report="report.txt",
        params:
            juliacmd=JULIA_COMMAND,
            scriptpath=f"{SNAKEDIR}/scripts/report.jl",
            refdir=REFSEQDIR
        log: "tmp/log/report.txt"
        threads: workflow.cores
        run:
            shell(f"{params.juliacmd} -t {threads} {params.scriptpath} nanopore . {params.refdir} > {log}")
