import sys
import os

SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, "scripts"))
import tools

JULIA_COMMAND = f"JULIA_LOAD_PATH='{SNAKEDIR}' julia --startup-file=no"

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
    inputs = ["report_consensus.txt"]

    for samplename in SAMPLENAMES:
        inputs.append(f"sequences/{samplename}/all.fna")
        inputs.append(f"depths/{samplename}_template.pdf")

    return inputs

rule all:
    input: done_input
    output: "commit_consensus.txt"
    params: SNAKEDIR
    shell: "git -C {params:q} rev-parse --short HEAD > {output} && cp {params:q}/copy_readme.md README_CONSENSUS.md"

#################################
# REFERENCE-ONLY PART OF PIPELINE
#################################
rule create_ref_fna_jls:
    input: REFDIR + "/refs.json"
    output:
        fna=temp(REFOUTDIR + "/refs.fna"),
        jls=REFOUTDIR + "/ref_segment_map.jls",
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/make_jlrefs.jl"
    shell: '{params.juliacmd} {params.scriptpath} {input:q} {output.fna:q} {output.jls:q}'

rule index_ref:
    input: REFOUTDIR + "/refs.fna"
    output:
        comp=REFOUTDIR + "/refs.comp.b",
        name=REFOUTDIR + "/refs.name",
        length=REFOUTDIR + "/refs.length.b",
        seq=REFOUTDIR + "/refs.seq.b"
    params:
        outpath=REFOUTDIR + "/refs",
    log: "tmp/log/kma_ref.log"
    # Set k from 16 to 12 for more sensitivity, but slower. May not be needed.
    shell: "kma index -nbp -k 12 -i {input:q} -o {params.outpath:q} 2> {log}"

############################
# CONSENSUS PART OF PIPELINE
############################
if IS_ILLUMINA:
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
            'fastp -i {input.fw:q} -I {input.rv:q} '
            '-o {output.fw:q} -O {output.rv:q} --html {output.html:q} --json {output.json:q} '
            '--disable_adapter_trimming --trim_poly_g --cut_tail --cut_front --low_complexity_filter '
            '--complexity_threshold 50 --thread {threads} 2> {log:q}'

    rule map_best_template:
        input:
            fw=rules.fastp.output.fw,
            rv=rules.fastp.output.rv,
            index=rules.index_ref.output
        output: "tmp/aln/{samplename}/initial.res"
        params:
            db=REFOUTDIR + "/refs",
            outbase="tmp/aln/{samplename}/initial"
        threads: 2
        log: "tmp/log/aln/{samplename}.initial.log"
        shell:
            "kma -ipe {input.fw} {input.rv} -o {params.outbase} -t_db {params.db:q} "
            "-t {threads} -mrs 0.3 -ConClave 2 -nc -nf 2> {log}"        

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
            'fastp -i {input:q} -o {output.reads} --html {output.html} '
            '--json {output.json} --disable_adapter_trimming  --disable_trim_poly_g '
            '--cut_window_size 10 --cut_mean_quality 10 --cut_tail --cut_front --low_complexity_filter  '
            '--complexity_threshold 50 --length_limit 2400 --length_required 100 '
            '--average_qual 12 --thread {threads} 2> {log}'

    rule map_best_template:
        input:
            reads=rules.fastp.output.reads,
            index=rules.index_ref.output
        output: "tmp/aln/{samplename}/initial.res"
        params:
            db=REFOUTDIR + "/refs", # same as index_reffile param
            outbase="tmp/aln/{samplename}/initial"
        threads: 2
        log: "tmp/log/aln/{samplename}.initial.log"
        shell:
            "kma -i {input.reads} -o {params.outbase:q} -t_db {params.db:q} "
            "-t {threads} -mrs 0.3 -ConClave 2 -nc -nf 2> {log}"

### Both platforms
# This rule downloads and installs all Julia packages needed
rule instantiate:
    output: touch(REFOUTDIR + "/cons_instantiated")
    params: JULIA_COMMAND
    shell: "{params} -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'"

rule collect_best_templates:
    input:
        spas=expand("tmp/aln/{samplename}/initial.res", samplename=SAMPLENAMES),
        inst=REFOUTDIR + "/cons_instantiated"
    output: temp(expand("tmp/aln/{samplename}/cat.fna", samplename=SAMPLENAMES))
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/gather_res.jl",
        refpath=REFDIR
    shell: "{params.juliacmd} {params.scriptpath:q} tmp/aln {params.refpath:q}"

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
    shell: "kma index -nbp -k 10 -i {input:q} -o {params.t_db} 2> {log}"

if IS_ILLUMINA:
    rule first_kma_map:
        input:
            fw=rules.fastp.output.fw,
            rv=rules.fastp.output.rv,
            index=rules.first_kma_index.output,
        output:
            res="tmp/aln/{samplename}/kma_1.res",
            fsa="tmp/aln/{samplename}/kma_1.fsa",
            mat="tmp/aln/{samplename}/kma_1.mat.gz",
        params:
            db="tmp/aln/{samplename}/cat",
            outbase="tmp/aln/{samplename}/kma_1",
        log: "tmp/log/aln/kma1_map_{samplename}.log"
        threads: 2
        shell:
            "kma -ipe {input.fw:q} {input.rv:q} -o {params.outbase} -t_db {params.db} "
            "-t {threads} -mrs 0.3 -1t1 -ConClave 2 -gapopen -5 -nf -matrix 2> {log}"

elif IS_NANOPORE:
    rule first_kma_map:
        input:
            reads=rules.fastp.output.reads,
            index=rules.first_kma_index.output,
        output:
            res="tmp/aln/{samplename}/kma_1.res",
            fsa="tmp/aln/{samplename}/kma_1.fsa",
            mat="tmp/aln/{samplename}/kma_1.mat.gz",
        params:
            db="tmp/aln/{samplename}/cat",
            outbase="tmp/aln/{samplename}/kma_1",
        log: "tmp/log/aln/kma1_map_{samplename}.log"
        threads: 2
        shell:
            "kma -i {input.reads:q} -o {params.outbase} -t_db {params.db} "
            "-mp 20 -bc 0.7 -t {threads} -mrs 0.3 -1t1 -ConClave 2 -bcNano -nf -matrix 2> {log}"

def iterative_reads(wc):
    if IS_ILLUMINA:
        return [f'tmp/trim/{wc.samplename}/fw.fq', f'tmp/trim/{wc.samplename}/rv.fq']
    else:
        return f'tmp/trim/{wc.samplename}/reads.fq'

rule iterative_assembly:
    input:
        reads=rules.fastp.output,
        asm=rules.first_kma_map.output.fsa,
        res=rules.first_kma_map.output.res,
        jls=rules.create_ref_fna_jls.output.jls
    output:
        asm="tmp/aln/{samplename}/kma_final.fsa",
        res="tmp/aln/{samplename}/kma_final.res",
        mat="tmp/aln/{samplename}/kma_final.mat.gz",
        conv="tmp/aln/{samplename}/convergence.tsv"
    log: "tmp/log/aln/{samplename}/iterative.log"
    threads: 2
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/iter_asm.jl",
        samplename=lambda wc: wc.samplename,
        segment_map=f"{REFOUTDIR}/ref_segment_map.jls",
        outdir=lambda wc: f"tmp/aln/{wc.samplename}",
        logdir=lambda wc: f"tmp/log/aln/{wc.samplename}",
        k=10,
        threshold=0.995 if IS_ILLUMINA else 0.990,
        reads=iterative_reads
    shell:
        "{params.juliacmd} -t {threads} {params.scriptpath:q} "
        "{params.samplename} {params.segment_map:q} {input.asm} {input.res} {params.outdir} "
        "{params.logdir} {params.k} {params.threshold} {params.reads} 2> {log}"

# Both platforms
# We need to map to multiple templates per segments to catch superinfections.
# but I haven't nailed down the heuristics for in gather_res.jl for detecting
# superinfections. If they are not truly present, the multiple assemblies will
# quickly converge. We trim primers here and remove converged assemblies
rule remove_primers:
    input:
        con=rules.iterative_assembly.output.asm,
        primers=f"{REFDIR}/primers.fna"
    output: temp("tmp/aln/{samplename}/assembly.fna")
    log: "tmp/log/consensus/remove_primers_{samplename}.txt"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/trim_consensus.jl",
        minmatches=4,
        fuzzylen=8
    shell: 
        "{params.juliacmd} {params.scriptpath:q} {input.primers:q} "
        "{input.con} {output} {params.minmatches} {params.fuzzylen} > {log}"

rule create_report:
    input:
        matrix=expand("tmp/aln/{samplename}/kma_final.mat.gz", samplename=SAMPLENAMES),
        assembly=expand("tmp/aln/{samplename}/assembly.fna", samplename=SAMPLENAMES),
        res=expand("tmp/aln/{samplename}/kma_final.res", samplename=SAMPLENAMES)
    output:
        consensus=expand("sequences/{samplename}/all.fna", samplename=SAMPLENAMES),
        depth=expand("depths/{samplename}_template.pdf", samplename=SAMPLENAMES),
        report="report_consensus.txt"
    params:
        juliacmd=JULIA_COMMAND,
        scriptpath=f"{SNAKEDIR}/scripts/report.jl",
        refdir=REFDIR,
        platform='illumina' if IS_ILLUMINA else 'nanopore'
    log: "tmp/log/report_consensus.txt"
    threads: workflow.cores
    shell: 
        "{params.juliacmd} -t {threads} {params.scriptpath:q} "
        "{params.platform} . {params.refdir:q} > {log}"
