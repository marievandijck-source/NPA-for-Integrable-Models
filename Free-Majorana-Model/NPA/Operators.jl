using Combinatorics: combinations

# --------------------------------------------------
# Majorana operators
# --------------------------------------------------

@enum MajoranaType γ1 = 1 γ2 = 2

struct Majorana
    mode::Int
    type::MajoranaType

    function Majorana(mode::Int, type::MajoranaType)
        mode >= 1 || throw(ArgumentError("mode must be at least 1"))
        return new(mode, type)
    end
end

canonical_index(c::Majorana) =
    2 * (c.mode - 1) + Int(c.type)

# --------------------------------------------------
# Words
# --------------------------------------------------

struct Word
    ops::Tuple{Vararg{Majorana}}
end

Word(ops::AbstractVector{Majorana}) = Word(Tuple(ops))
Word(ops::Majorana...) = Word(ops)

# --------------------------------------------------
# Operators
# --------------------------------------------------

struct Operator{T}
    coeffs::Dict{Word,T}
end

Base.:(==)(A::Operator, B::Operator) =
    A.coeffs == B.coeffs

# --------------------------------------------------
# Basis generation
# --------------------------------------------------

function generate_majoranas(n_modes::Int)
    n_modes >= 1 ||
        throw(ArgumentError("n_modes must be at least 1"))

    return [
        Majorana(mode, type)
        for mode in 1:n_modes
        for type in (γ1, γ2) 
    ]
end

function word_basis(n_modes::Int, max_degree::Int)
    max_degree >= 0 ||
        throw(ArgumentError("max_degree must be nonnegative"))

    majoranas = generate_majoranas(n_modes)
    highest_degree = min(max_degree, length(majoranas))

    basis = Word[Word()]

    for degree in 1:highest_degree
        for combination in combinations(majoranas, degree)
            push!(basis, Word(combination))
        end
    end

    return basis
end