using JuMP
using LinearAlgebra
using Test
using SCS  
using Pkg
import MathOptInterface as MOI


include("Operators.jl")
include("Algebra.jl")
include("Moment-Matrix.jl")
include("SDP.jl")

# --------------------------------------------------
# Define the problem
# --------------------------------------------------

n_modes = 4
n = 2 * n_modes
degree = 3

# --------------------------------------------------
# Build c^T H c 
# --------------------------------------------------

basis = word_basis(n_modes, degree) #word_basis(n_modes::Int, max_degree::Int)

H = zeros(Int, n, n)
for i in 1:n, j in 1:n
    if i < j
        H[i, j] = 1
    elseif i > j
        H[i, j] = -1
    end
end

c = generate_majoranas(n_modes)   # order matches H's index order

coeffs = Dict{Word,ComplexF64}()
for i in 1:n, j in (i+1):n
    Hij = H[i, j]
    if Hij != 0
        w = Word(c[i], c[j])
        coeffs[w] = get(coeffs, w, 0.0 + 0.0im) + im * Hij   # <-- factor of i, adjust prefactor to match your convention
    end
end

Hop = Operator(coeffs)

# --------------------------------------------------
# Build the SDP
# --------------------------------------------------
problem = build_moment_sdp(
    basis,
    Hop,
    SCS.Optimizer,
)

set_optimizer_attribute(problem.model, "eps_abs", 1e-10)

result = solve_moment_sdp!(
    problem;
    verbose=false,
)

println("Status: ", result.termination_status)
println("Energy: ", result.objective)