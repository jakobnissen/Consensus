# Contents of this directory
This directory contains files from a run of the SSI Influenza consensus pipeline.
Read more on the Consensus page of the InfluenzaDocumentation page.
Contact Jakob Nybo Nissen for details.
There may be additional files unlisted here, most likely from the phylogeny pipeline.
See its equivalent README file.

### report_consensus.txt
This auto-generated report summarizes all consensus sequences. Its contents include
* Short summary read statistics for each sample
* Depth and coverage of each constructed influenza segment
* Identity of segment to the chosen reference it assembled against
* Any problems with the alignment to reference, assembly, and proteins in each segment.

Each segment may be marked "FAIL" if it fails automated quality controls.
Segments that are not marked FAIL can with reasonable certaincy be considered alright.

### consensus_versions.txt
This file displays the version of the consensus pipeline and the Julia language used to
generate the files in this directory.

### depths
In this directory, there are depth plots for every constructed segment. Each sample has two
files, `SAMPLE_template.pdf` and `SAMPLE_assembly.pdf`. Both show the depth (y-axis) along
each segment (x-axis, scaled by segment length).

The `template` file measures depth relative to the original reference it assembled against.
The `assembly` measures depth relative to itself. These may differ if there are indels or
significant nucleotide differences between reference and assembly

The raw values displayed in the plots are available in the .tsv files. Here, the "order"
field is relevant if there are multiple different types of the same segment in the same
sample. Then, the segments are ordered (approximately) by depth such that the segment
with order 1 is the most abundant.

### sequences
Each sample has 6 files in this directory: `all`, `secondary` and `primary`, each of these
with two endings, `.fna` and `.faa`. `.fna` signifies DNA, `.faa` amino acids.
`primary` only contains all segments that passed quality control and which have been designated
primary segments by the pipeline. If there are multiple instances of the same segment in
the same sample, one of each segment is determined to be the primary segment.
`secondary` are all the remaining passed segments not marked as primary.

The `all` file contains everything in `primary` and `secondary`, as well as all segments
that did NOT pass quality control.

### tmp
This directory contains files that are considered internal to the pipeline.
That means the content is subject to change and not documented.
Most importantly, the file `tmp/internal.jls.gz` contains a lot of relevant data
for each assembly, and can be loaded by a computer program and manipulated.
