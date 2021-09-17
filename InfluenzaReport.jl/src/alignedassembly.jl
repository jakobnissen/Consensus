function load_aligned_assemblies(
    asm_paths::Vector{String},
    jsonpath::AbstractString,
    is_kma::Bool,
)::Vector{SegmentTuple{Option{AlignedAssembly}}}
    asms = map(path -> load_assembly(path, is_kma), asm_paths)
    refs = find_references(asms, jsonpath)

    zip(asms, refs) |> Map() do (asm, ref)
        ntuple = SegmentTuple(zip(asm, ref))
        map(ntuple) do (maybe_asm, maybe_ref)
            (is_error(maybe_asm) || is_error(maybe_ref)) && return none(AlignedAssembly)
            alnasm = AlignedAssembly(unwrap(maybe_asm), unwrap(maybe_ref))
            add_alnasm_errors!(alnasm)
            some(alnasm)
        end
    end |> Folds.collect
end

function load_assembly(path::AbstractString, kma::Bool)::SegmentTuple{Option{Assembly}}
    records = open(collect, FASTA.Reader, path)
    result = fill(none(Assembly), N_SEGMENTS)
    for record in records
        temp_asm = Assembly(record, nothing)
        (accession, segment) = let
            if kma
                # Format should be HEADER_SEGMENT
                s = try_split_segment(temp_asm.name)
                if s === nothing
                    error("In $path, found header \"$(temp_asm.name)\", expected pattern \"HEADER_SEGMENT\"")
                else
                    s
                end
            else
                s = try_parse_medaka_header(temp_asm.name)
                if s === nothing
                    error("In $path, found header \"$(temp_asm.name)\", expected pattern \"HEADER_SEGMENT_segment0 HEADER_SEGMENT[stuff]\"")
                else
                    s
                end
            end
        end
        segment_index = reinterpret(UInt8, segment) + 0x01
        is_error(result[segment_index]) || error("Segment $segment present twice in $path")
        result[segment_index] = some(Assembly(accession, temp_asm.seq, some(segment), temp_asm.insignificant))
    end
    return SegmentTuple(result)
end

function try_parse_medaka_header(s::Union{String, SubString{String}})
    ps = findall("_segment0", s)
    # if _segment0 is contained in the header N times, it will appear 2N+1 times total.
    isodd(length(isempty(ps))) || return nothing
    underscore_pos = first(ps[cld(length(ps), 2)])
    stripped_header = SubString(s, 1:prevind(s, underscore_pos))
    return try_split_segment(stripped_header)
end

function try_split_segment(s::Union{String, SubString{String}})
    p = findlast(isequal(UInt8('_')), codeunits(s))
    p === nothing && return nothing
    seg = tryparse(Segment, SubString(s, p+1:lastindex(s)))
    seg === nothing && return nothing
    return (SubString(s, 1, prevind(s, p)), seg)
end

function find_references(
    asms::Vector{SegmentTuple{Option{Assembly}}},
    jsonpath::AbstractString
)::Vector{SegmentTuple{Option{Reference}}}

    accessions = Set{String}()
    for asmt in asms, masm in asmt
        acc = (@unwrap_or masm continue).name
        push!(accessions, acc)
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

    return map(asms) do asmt
        map(asmt) do masm
            if is_error(masm)
                none(Reference)
            else
                some(byaccession[unwrap(masm).name])
            end
        end
    end
end

function add_alnasm_errors!(alnasm::AlignedAssembly)
    # Low identity to reference
    alnasm.identity < 0.8 && push!(alnasm.errors, Influenza.ErrorLowIdentity(alnasm.identity))
    return nothing
end

"This creates a FASTA record with lowercase letters on uncertain positions"
function asm_dna_record(asm::Assembly, header::String)
    record = FASTA.Record(header, asm.seq)
    insignificant = @unwrap_or asm.insignificant (return record)

    data = record.data
    @assert length(record.sequence) == length(insignificant)
    for (isbad, seqpos) in zip(insignificant, record.sequence)
        if isbad
            data[seqpos] += 0x20 # this sets it to lowercase for ASCII
        end
    end
    return record
end

# Add error to aligned assembly if kma2.res does not show they
# have converged
function kma2_identity_check(
    alnasms::Vector{SegmentTuple{Option{AlignedAssembly}}},
    respaths::Vector{String}
)::Nothing
    # TODO: Parallelize?
    for (m_alnasm, respath) in zip(alnasms, respaths)
        seen = falses(N_SEGMENTS)
        res = open(respath) do io
            KMATools.parse_res(io, respath)
        end
        for row in res
            segment::Segment = let
                p = findlast('_', row.template)
                s = p === nothing ? nothing : tryparse(Segment, strip(row.template[p+1:end]))
                s !== nothing ? s : error("Could not parse segment in file $respath, line $(row.template)")
            end
            index = UInt8(segment) + 0x01
            if seen[index]
                error("Duplicate segment $segment in file $respath")
            end
            seen[index] = true
            alnasm = @unwrap_or m_alnasm[index] continue

            # Identity can be low either due to low template id or low query id
            # these may differ due to indels. We take the minimum of these.
            # Strictly speaking, this means an N nt insersion in template plus an
            # N nt deletion would only count as N differences, but it's OK here.
            id = min(row.tid, row.qid)
            if id < 0.995
                push!(alnasm.errors, Influenza.ErrorAssemblyNotConverged(id))
            end
        end
    end
end

function write_files(
    cons_dirname::AbstractString,
    tmp_dirname::AbstractString,
    samplenames::Vector{String},
    alnasms::Vector{SegmentTuple{Option{AlignedAssembly}}},
    passes::Vector{SegmentTuple{Bool}}
)
    open(joinpath(tmp_dirname, "aln_asms.jl"), "w") do io
        serialize_alnasms(io, alnasms, passes)
    end
    write_consensus(cons_dirname, samplenames, alnasms, passes)
end

function write_consensus(
    dirname::AbstractString,
    samplenames::Vector{String},
    alnasms::Vector{SegmentTuple{Option{AlignedAssembly}}},
    passes::Vector{SegmentTuple{Bool}},
)
    isdir(dirname) || mkdir(dirname)
    for (samplename, alnasm_tup, pass_tup) in zip(samplenames, alnasms, passes)
        subdir = joinpath(dirname, samplename)
        isdir(subdir) || mkdir(subdir)
        
        cons_dna_writer = open(FASTA.Writer, joinpath(subdir, "consensus.fna"))
        cons_aa_writer =  open(FASTA.Writer, joinpath(subdir, "consensus.faa"))
        cura_dna_writer = open(FASTA.Writer, joinpath(subdir, "curated.fna"))
        cura_aa_writer =  open(FASTA.Writer, joinpath(subdir, "curated.faa"))

        for (_i, (m_alnasm, is_passed)) in enumerate(zip(alnasm_tup, pass_tup))
            alnasm = @unwrap_or m_alnasm continue
            segment = Segment(_i - 0x01)

            # Write DNA
            asm = alnasm.assembly
            dna_record = asm_dna_record(asm, samplename * '_' * string(segment))
            write(cons_dna_writer, dna_record)
            is_passed && write(cura_dna_writer, dna_record)

            # Write proteins
            m_aaseqs = translate_proteins(alnasm)
            for (protein, m_aaseq) in zip(alnasm.proteins, m_aaseqs)
                orfs = @unwrap_or protein.orfs continue
                aaseq = @unwrap_or m_aaseq continue
                header = (
                    samplename * '_' * string(protein.variant) * '_' *
                    join(["$(first(i))-$(last(i))" for i in orfs], ',')
                )
                record = FASTA.Record(header, aaseq)
                write(cons_aa_writer, record)
                is_passed && write(cura_aa_writer, record)
            end
        end

        close(cons_dna_writer)
        close(cons_aa_writer)
        close(cura_dna_writer)
        close(cura_aa_writer)
    end
end


