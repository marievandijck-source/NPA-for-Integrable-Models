# --------------------------------------------------
# Entry of the Moment Matrix
# --------------------------------------------------

struct MomentEntry
    coefficient::ComplexF64
    word::Word
end

# --------------------------------------------------
# Moment Matrix
# --------------------------------------------------

function moment_matrix_structure(basis::Vector{Word})
    m = length(basis)
    structure = Matrix{MomentEntry}(undef, m, m)

    for i in eachindex(basis)
        for j in eachindex(basis)
            coefficient, word =
                adjoint_product(basis[i], basis[j])

            structure[i, j] =
                MomentEntry(coefficient, word)
        end
    end

    return structure
end

# --------------------------------------------------
# Collect and sort words from the Moment Matrix
# --------------------------------------------------

function word_sort_key(w::Word)
    return (
        length(w.ops),
        Tuple(canonical_index(c) for c in w.ops),
    )
end 

function moment_words(
    structure::AbstractMatrix{MomentEntry},
)
    words = Set{Word}()

    for entry in structure
        push!(words, entry.word)
    end

    return sort!(
        collect(words);
        by=word_sort_key,
    )
end