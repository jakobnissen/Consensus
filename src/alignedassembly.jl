function load_aligned_assemblies(
    paths::Vector{String},
    samples::Vector{Sample},
    jsonpath::AbstractString,
)::Vector{Vector{AlignedAssembly}}
    asms = map(path -> load_assembly(path), paths)
    refs = find_references(asms, jsonpath)
    zip(asms, refs, samples) |> Map() do (asmv, refv, sample)
        map(collect(zip(asmv, refv))) do (asm, ref)
            alnasm = AlignedAssembly(asm, ref)
            add_alnasm_errors!(alnasm)
            alnasm
        end
    end |> Folds.collect
end

function load_assembly(path::AbstractString)::Vector{Assembly}
    map(open(collect, FASTA.Reader, path)) do record
        temp_asm = Assembly(record, nothing)
        (accession, segment) = split_segment(temp_asm.name)
        Assembly(accession, temp_asm.seq, some(segment), temp_asm.insignificant)
    end
end

function find_references(
    asms::Vector{Vector{Assembly}},
    jsonpath::AbstractString
)::Vector{Vector{Reference}}

    accessions = Set{String}()
    for asmv in asms, asm in asmv
        push!(accessions, asm.name)
    end

    byaccession = Dict{String, Reference}()
    added_accessions = Set{String}()
    for reference in Influenza.load_references(jsonpath)
        if reference.name in accessions
            push!(added_accessions, reference.name)
            byaccession[reference.name] = reference
        end
    end
    missing_acc = setdiff(accessions, added_accessions)
    if !isempty(missing_acc)
        error("Accession $(first(missing_acc)) missing from $jsonpath")
    end

    return map(asms) do asmv
        map(asmv) do asm
            byaccession[asm.name]
        end
    end
end

function add_alnasm_errors!(alnasm::AlignedAssembly)
    # Low identity to reference
    alnasm.identity < 0.8 && push!(alnasm.errors, Influenza.ErrorLowIdentity(alnasm.identity))
    return nothing
end

# Add error to aligned assembly if the convergence report show some segments
# have not converged
function convergence_check(alnasms::Vector{AlignedAssembly}, convpath::String)
    byid = Dict{Tuple{String, Segment}, AlignedAssembly}()
    for alnasm in alnasms
        @assert !haskey(byid, alnasm.reference.name) # for safety
        byid[(alnasm.reference.name, alnasm.reference.segment)] = alnasm
    end
    open(convpath) do io
        header = readline(io)
        @assert strip(header) == "template\tsegment\tisconverged\tid"
        for row in eachline(io)
            identifier, segment, isconverged_, id_ = split(row, '\t')
            alnasm = byid[(identifier, parse(Segment, segment))]
            if !(parse(Bool, isconverged_))
                id = parse(Float64, id_)
                push!(alnasm.errors, Influenza.ErrorAssemblyNotConverged(id))
            end
        end
    end
end

function write_sequences(
    seq_dir_name::AbstractString,
    sample::Sample,
    alnasms::Vector{AlignedAssembly},
    passed::Vector{Bool},
    primary::Vector{Bool}
)
    isdir(seq_dir_name) || mkpath(seq_dir_name)
    file_contents = Dict{String, Vector{FASTA.Record}}()
    for i in eachindex(alnasms, passed, primary)
        alnasm = alnasms[i]
        asm = alnasm.assembly
        segment = alnasm.reference.segment
        samplename = nameof(sample)

        # All DNA
        dnastr = dna_str(asm.seq, asm.insignificant)
        primstr = primary[i] ? "primary" : "secondary"
        rec = FASTA.Record("$(samplename)_$(segment)_$(primstr)", dnastr)
        push!(get!(valtype(file_contents), file_contents, "all.fna"), rec)

        # Primary / secondary DNA
        if passed[i]
            filename = primary[i] ? "primary.fna" : "secondary.fna"
            rec = FASTA.Record("$(samplename)_$(segment)", dnastr)
            push!(get!(valtype(file_contents), file_contents, filename), rec)
        end

        # Protein
        m_aaseqs = Influenza.translate_proteins(alnasm)
        for (protein, m_aaseq) in zip(alnasm.proteins, m_aaseqs)
            orfs = @unwrap_or protein.orfs continue
            aaseq = @unwrap_or m_aaseq continue
            id = "$(samplename)_$(protein.variant)"
            orfstr = join(["$(first(i))-$(last(i))" for i in orfs], ',')
            rec = FASTA.Record("$(id)_$(primstr)_$(orfstr)", aaseq)
            push!(get!(valtype(file_contents), file_contents, "all.faa"), rec)
            if passed[i]
                filename = primary[i] ? "primary.faa" : "secondary.faa"
                rec = FASTA.Record("$(id)_$(orfstr)", aaseq)
                push!(get!(valtype(file_contents), file_contents, filename), rec)
            end
        end

        # Write files
        for (filename, recs) in file_contents
            open(FASTA.Writer, joinpath(seq_dir_name, filename)) do writer
                for rec in recs
                    write(writer, rec)
                end
            end
        end
    end
end

"This creates a FASTA record with lowercase letters on uncertain positions"
function dna_str(seq::LongDNASeq, insig::Option{<:AbstractVector{Bool}})
    buf = IOBuffer()
    print(buf, seq)
    vec = take!(buf)
    bvec = unwrap_or(insig, nothing)
    if bvec !== nothing
        for i in eachindex(vec, bvec)
            # This sets ASCII letters to lowercase
            vec[i] += (0x20 * bvec[i])
        end
    end
    return String(vec)
end
