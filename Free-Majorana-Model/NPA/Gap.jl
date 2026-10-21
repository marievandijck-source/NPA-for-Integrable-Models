# Test voor een random state
# Test voor andere powers
# Test voor general H 
# Make figure


# include("Operators.jl")
# import .Operators: Majorana, Word, Operator, γ1, γ2
# include("Algebra.jl")

# using Pkg; Pkg.add("PauliStrings")
using PauliStrings
using IterTools
using Combinatorics: combinations
using Plots
using LaTeXStrings

function paired_word_basis(N::Int)
    basis = Word[]

    # Empty subset -> identity
    # push!(basis, Word())

    for r in 1:N
        for S in combinations(1:N, r)

            ops = Majorana[]

            for k in S
                push!(ops, Majorana(k, γ1))
                push!(ops, Majorana(k, γ2))
            end

            push!(basis, Word(ops))
        end
    end

    return basis
end

function build_H_equal_coeffs(N, s, c)
    Hcoeffs = Dict{Word,ComplexF64}()

    id = Word()

    for r in 1:N
        for S in combinations(1:N, r)

            ops = Majorana[]

            for k in S
                push!(ops, Majorana(k, γ1))
                push!(ops, Majorana(k, γ2))
            end

            wS = Word(ops)

            # x_S = i^|S| * (∏ s_k) * ∏ d_k0 d_k1
            prefactor = (1im)^r * prod(s[k] for k in S)

            # + c * identity
            Hcoeffs[id] =
                get(Hcoeffs, id, 0.0 + 0im) + c

            # - c * x_S
            Hcoeffs[wS] =
                get(Hcoeffs, wS, 0.0 + 0im) - c * prefactor
        end
    end
    
    return Operator(Hcoeffs)
end

function is_consecutive_pbc(S, N)
    r = length(S)

    if r == 1
        return true
    end

    gaps = [diff(S); S[1] + N - S[end]]

    return count(gaps .== 1) == r - 1
end
    

function energy(z,n, coefficient; neighbours::Bool=false)
    N = length(z)
    E = 0.0

    for r in 1:n
        for S in combinations(1:N, r)
            if neighbours && !is_consecutive_pbc(S, N)
                continue
            end

            xS = prod(z[k] for k in S)
            
            if coefficient isa Number
                c = coefficient
            else 
                c = coefficient[r]
            end

            E += c * (1 - xS) # we here removed a square which results in a factor 1/2
        end
    end

    return E
end

function spectrum(N, n, coefficient; neighbours::Bool=false)
    states = product(fill((-1, 1), N)...)

    spec = []

    for z in states
        E = energy(z, n, coefficient; neighbours=neighbours)

        push!(spec, (
            energy = E,
            state = collect(z)
        ))
    end

    sort!(spec, by = x -> x.energy)

    return spec
end

N = 10
basis = paired_word_basis(N)

#coefficient = N/(2^N-1)
coefficient = 1/(2^N-1)

s = ones(Int, N) #Choose target GS

H = build_H_equal_coeffs(N, s, coefficient)

states = collect(product(fill((-1, 1), N)...))

spec = spectrum(N, N, coefficient)

sort!(spec, by = x -> x.energy)

print(spec[1:2])
#print(spec[end-4:end])
# Analytical 
E = 2*2^(N-1)/(2^N-1)
println("Gap = $E")

# Construct H up to second order
n=2
c2 = 2/(N^2+N)

spec2 = spectrum(N,n, c2)

sort!(spec2, by = x -> x.energy)

print(spec2[1:2])

#Analytical 

E2 = 2*N*2/(N^2+N)
print(E2)


# Construct H with structure of NN


cNN = [2/(N*(N+1)), (N-1)/(N^2+N)]

specNN = spectrum(N,n, cNN; neighbours=true)

sort!(specNN, by = x -> x.energy)

print(specNN[1:2])

ENN = 4/(N+1)

println(ENN)

# Constant Hamiltonian of first order

c1 = 1/N

n=1

spec1 = spectrum(N,n, c1; neighbours=true)

sort!(spec1, by = x -> x.energy)

print(spec1[1:2])

E1 = 2/(N)

println(E1)

# All connected sets 

cC = 1/(N*(N+1))

n=N

specC = spectrum(N,n, cC; neighbours=true)

sort!(specC, by = x -> x.energy)

print(specC[1:2])

EC = 2/(N)

println(EC)


#=
# --------------------------------------------------
# Compute numerical gaps
# --------------------------------------------------

Ns = 3:16

gap_all    = Float64[]
gap_order2 = Float64[]
gap_NN     = Float64[]
gap_order1 = Float64[]

for N in Ns
print(N)
    # 1. All interaction orders
    c_all = 1 / (2^N - 1)
    spec = spectrum(N, N, c_all)
    push!(gap_all, get_gap(spec))

    # 2. All-to-all interactions up to order 2
    c2 = 2 / (N^2 + N)
    spec2 = spectrum(N, 2, c2)
    push!(gap_order2, get_gap(spec2))

    # 3. Nearest-neighbour interactions up to order 2
    cNN = [
        2 / (N * (N + 1)),
        (N - 1) / (N^2 + N)
    ]

    specNN = spectrum(N, 2, cNN; neighbours=true)
    push!(gap_NN, get_gap(specNN))

    # 4. Only first-order interactions
    c1 = 1 / N
    spec1 = spectrum(N, 1, c1)
    push!(gap_order1, get_gap(spec1))
end


# --------------------------------------------------
# Analytical functions
# --------------------------------------------------

Δ_all(N)    = 2.0^N / (2.0^N - 1)
Δ_order2(N) = 4 / (N + 1)
Δ_NN(N)     = 4 / (N + 1)
Δ_order1(N) = 2 / N

Ncontinuous = range(
    minimum(Ns),
    maximum(Ns),
    length=500
)


# --------------------------------------------------
# Plot
#
# Analytical + numerical are plotted consecutively,
# so they are grouped nicely in the legend.
# --------------------------------------------------

p = plot(
    xlabel="System size N",
    ylabel="Energy gap Δ",
    title="Energy gap versus system size",
    legend=:outerright,
    size=(1100, 550),
    framestyle=:box
)


# ==================================================
# 1. ALL ORDERS — BLACK
# ==================================================
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
plot!(
    p,
    Ncontinuous,
    Δ_all.(Ncontinuous),
    label="All orders — analytical:  Δ = 2ᴺ/(2ᴺ − 1)",
    linewidth=2,
    linecolor="black"
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
scatter!(
    p,
    Ns,
    gap_all,
    label="All orders — numerical",
    markercolor="black",
    markerstrokecolor="black",
    markersize=5,
    marker=:circle
)

# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)

# ==================================================
# 2. ORDER 2, ALL-TO-ALL — LIGHT BLUE
# ==================================================

plot!(
    p,
    Ncontinuous,
    Δ_order2.(Ncontinuous),
    label="Order ≤ 2, all-to-all — analytical:  Δ = 4/(N + 1)",
    linewidth=2,
    linecolor="lightblue"
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
scatter!(
    p,
    Ns,
    gap_order2,
    label="Order ≤ 2, all-to-all — numerical",
    markersize=6,
    markercolor="lightblue",
    markerstrokecolor="lightblue",
    marker=:square
)

# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)

# ==================================================
# 3. ORDER 2, NN — DARK BLUE
# ==================================================

plot!(
    p,
    Ncontinuous,
    Δ_NN.(Ncontinuous),
    label="Order ≤ 2, NN — analytical:  Δ = 4/(N + 1)",
    linewidth=2,
    linecolor="darkblue",
    linestyle=:dash
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
scatter!(
    p,
    Ns,
    gap_NN,
    label="Order ≤ 2, NN — numerical",
    markersize=5,
    markercolor="darkblue",
    markerstrokecolor="darkblue",
    marker=:circle
)

# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
# ==================================================
# 4. FIRST ORDER — BLUE
# ==================================================

plot!(
    p,
    Ncontinuous,
    Δ_order1.(Ncontinuous),
    label="First order — analytical:  Δ = 2/N",
    linewidth=2,
    linecolor="blue"
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
scatter!(
    p,
    Ns,
    gap_order1,
    label="First order — numerical",
    markersize=5,
    markercolor="blue",
    markerstrokecolor="blue",
    marker=:circle
)
# Empty legend row
plot!(p, [NaN], [NaN],
    label=" ",
    linealpha=0,
    markerstrokecolor=:transparent
)
display(p)

=#
#See if a higher degree also makes no differnece between the order and the connectedness of the set 
#Find the relation between the degree and the way the gap scales → introduce the degree as a function of N