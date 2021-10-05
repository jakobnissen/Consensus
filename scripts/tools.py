import os
import re

def get_read_pairs(directory):
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"Could not locate read directory: \"{os.path.abspath(directory)}\"")

    try:
        filenames = sorted(next(os.walk(directory))[2])
    except StopIteration:
        raise FileNotFoundError("Must specify at least one read pair. "
                                "Are you sure you got the directory correct? "
                                f"{os.path.abspath(directory)}")

    # See Illumina's "Naming convention" documentation.
    # Currently, this logic does not take multiple lanes into account: We assume
    # only lane 1 is used.
    PATTERN = r"(.+)_S(\d+)_L001_R([12])_001\.fastq\.gz"
    pattern = re.compile(PATTERN)

    reads = dict()
    for filename in filenames:
        match = pattern.match(filename)
        if match is None:
            raise ValueError(
                ("Filename {} does not match pattern ".format(filename) + PATTERN +
                " Are you sure it is an Illumina FASTQ file?")
                )

        samplename, samplenumber, readnumber = match.groups()

        # Samplenumber 0 are the undetermined reads, always. These cannot meaningfully
        # be worked with.
        if samplenumber == 0:
            continue

        if samplename not in reads:
            reads[samplename] = []
        reads[samplename].append(os.path.join(directory, filename))

    singletons = [v for v in reads.values() if len(v) != 2]
    if len(singletons) > 0:
        print("Non-paired read files")
        for group in singletons:
            for file in group:
                print(file)
            print("")

        raise ValueError("Some FASTQ files not found in pairs")

    return reads

def get_nanopore_reads(directory):
    if not os.path.isdir(directory):
        raise FileNotFoundError(f"{directory} is not a directory.")

    filenames = sorted(next(os.walk(directory))[2])
    if len(filenames) == 0:
        raise ValueError(f"No files found in {directory}")

    reads = dict()
    for filename in filenames:
        for ending in [".fq.gz", ".fq", ".fastq.gz", ".fastq"]:
            if filename.endswith(ending):
                reads[filename[:-len(ending)]] = os.path.join(directory, filename)
                break
        else: # no break
            raise ValueError(f"File {filename} is not named like a FASTQ file.")
    
    return reads
