# Motivation: Pia wanted me to append dates to the FASTA files from swine samples.
# If sample date is not available, instead append recievedate.
# Instead of modifying the pipeline, it's easier to rename the files after the fact.
# function read_datefile(io::IO, samples::Set{String})