function load_aligned_assemblies(
    paths::Vector{String},
    jsonpath::AbstractString,
)::Vector{Vector{AlignedAssembly}}
    asms = map(path -> load_assembly(path), paths)
    refs = find_references(asms, jsonpath)
    result = Vector{Vector{AlignedAssembly}}(undef, length(paths))
    Threads.@threads for i in eachindex(result, asms, refs, paths)
        v = map(collect(zip(asms[i], refs[i]))) do (asm, ref)
            alnasm = AlignedAssembly(asm, ref)
            add_alnasm_errors!(alnasm)
            alnasm
        end
        result[i] = v
    end
    return result
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
    order::Vector{UInt8}
)
    isdir(seq_dir_name) || mkpath(seq_dir_name)
    file_contents = Dict{String, Vector{FASTA.Record}}()
    for i in eachindex(alnasms, passed, order)
        alnasm = alnasms[i]
        asm = alnasm.assembly
        segment = alnasm.reference.segment
        samplename = nameof(sample)
        ord = order[i]
        primary = ord == 1

        # All DNA
        dnastr = dna_str(asm.seq, asm.insignificant)
        sec_header = "$(samplename)_$(segment)_$(ord)"
        prim_header = "$(samplename)_$(segment)"
        rec = FASTA.Record(sec_header, dnastr)
        push!(get!(valtype(file_contents), file_contents, "all.fna"), rec)

        # Primary / secondary DNA
        if passed[i]
            filename = primary ? "primary.fna" : "secondary.fna"
            header = primary ? prim_header : sec_header
            rec = FASTA.Record(header, dnastr)
            push!(get!(valtype(file_contents), file_contents, filename), rec)
        end

        # Protein
        m_aaseqs = Influenza.translate_proteins(alnasm)
        for (protein, m_aaseq) in zip(alnasm.proteins, m_aaseqs)
            orfs = @unwrap_or protein.orfs continue
            aaseq = @unwrap_or m_aaseq continue
            id = "$(samplename)_$(protein.variant)"
            orfstr = join(["$(first(i))-$(last(i))" for i in orfs], ',')
            sec_header = "$(id)_$(ord)_$(orfstr)"
            prim_header = "$(id)_$(orfstr)"
            rec = FASTA.Record(sec_header, aaseq)
            push!(get!(valtype(file_contents), file_contents, "all.faa"), rec)
            if passed[i]
                filename = primary ? "primary.faa" : "secondary.faa"
                header = primary ? prim_header : sec_header
                rec = FASTA.Record(header, aaseq)
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

