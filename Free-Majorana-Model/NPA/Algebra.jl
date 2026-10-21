# --------------------------------------------------
# Implement Clifford-Algebra
# --------------------------------------------------

# --------------------------------------------------
# Canonicalize 
# --------------------------------------------------
Base.isless(a::Majorana,b::Majorana) =
    (a.mode,a.type) < (b.mode,b.type)

function permutation_sign(ops)
    inversions = 0

    for i in 1:length(ops)-1
        for j in i+1:length(ops)
            if canonical_index(ops[i]) > canonical_index(ops[j])
                inversions += 1
            end
        end
    end

    return iseven(inversions) ? 1 : -1
end

function canonical_order(ops)
    sign = permutation_sign(ops)
    sorted_ops = sort(collect(ops), by=canonical_index)

    return sign, sorted_ops
end

function canonicalize(word::Word)
    ops = collect(word.ops)

    sign = permutation_sign(ops)
    sort!(ops, by=canonical_index)

    reduced = Majorana[]

    for op in ops
        if !isempty(reduced) && reduced[end] == op
            pop!(reduced)  
        else
            push!(reduced, op)
        end
    end

    return sign, Word(Tuple(reduced))
end

# --------------------------------------------------
# Multiply
# --------------------------------------------------

Base.:*(w1::Word, w2::Word) =
    canonicalize(Word((w1.ops..., w2.ops...)))

function Base.:*(A::Operator, B::Operator)
    result = Dict{Word,Complex64}()

    for (wordA, coeffA) in A.coeffs
        for (wordB, coeffB) in B.coeffs
            sign, product_word = wordA * wordB

            new_coeff = ComplexF64(sign * coeffA * coeffB)

            result[product_word] =
                get(result, product_word, 0.0 + 0.0im) + new_coeff
        end
    end

    filter!(
        pair -> !isapprox(pair.second, 0.0; atol=1e-12),
        result,
    )

    return Operator(result)
end

# --------------------------------------------------
# Dagger
# --------------------------------------------------
Base.adjoint(c::Majorana) = c 

function Base.adjoint(w::Word)
    sign, wdag = canonicalize(Word(reverse(w.ops)))
    return sign, wdag
end

function Base.adjoint(A::Operator)
    result = Dict{Word,ComplexF64}()

    for (word, coeff) in A.coeffs
        sign, canonical_word = adjoint(word)
        new_coeff = sign * conj(coeff)

        result[canonical_word] =
            get(result, canonical_word, 0.0 + 0.0im) +
            new_coeff
    end

    filter!(
        pair -> !isapprox(pair.second, 0.0; atol=1e-12),
        result,
    )

    return Operator(result)
end

# --------------------------------------------------
# Commutation rules
# --------------------------------------------------

function Base.:+(A::Operator, B::Operator)
    result = Dict{Word,ComplexF64}()

    for (word, coeff) in A.coeffs
        result[word] =
            get(result, word, 0.0 + 0.0im) + coeff
    end

    for (word, coeff) in B.coeffs
        result[word] =
            get(result, word, 0.0 + 0.0im) + coeff
    end

    filter!(
        pair -> !isapprox(pair.second, 0.0; atol=1e-12),
        result,
    )

    return Operator(result)
end


function Base.:-(A::Operator, B::Operator)
    result = Dict{Word,ComplexF64}()

    for (word, coeff) in A.coeffs
        result[word] =
            get(result, word, 0.0 + 0.0im) + coeff
    end

    for (word, coeff) in B.coeffs
        result[word] =
            get(result, word, 0.0 + 0.0im) - coeff
    end

    filter!(
        pair -> !isapprox(pair.second, 0.0; atol=1e-12),
        result,
    )

    return Operator(result)
end

function commutator(A::Operator, B::Operator)
    return A*B - B*A
end

function anticommutator(A::Operator, B::Operator)
    return A*B + B*A
end

# --------------------------------------------------
# Adjoint product 
# --------------------------------------------------

function adjoint_product(u::Word, v::Word)
    sign_u, udag = adjoint(u)
    sign_product, product_word = udag * v

    coefficient = ComplexF64(sign_u * sign_product)

    return coefficient, product_word
end 