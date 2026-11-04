# -----------------------------------------------------------
# Probability that a subset of size k contains
# an odd number of the n excited sites
# -----------------------------------------------------------

function odd_probability(N::Int, n::Int, k::Int)

    jmin = max(0, k - (N - n))
    jmax = min(n, k)

    numerator = big(0)

    for j in jmin:jmax
        if isodd(j)
            numerator +=
                binomial(big(n), j) *
                binomial(big(N - n), k - j)
        end
    end

    denominator = binomial(big(N), k)

    return Float64(
        BigFloat(numerator) / BigFloat(denominator)
    )
end


# -----------------------------------------------------------
# Optimize the gap
#
# Normalization:
#
#     sum_k c[k] = 1
# -----------------------------------------------------------

function optimize_gap(
    N::Int,
    kmax::Int;
    term_energy = 2.0
)

    @assert 1 <= kmax <= N

    # --------------------------------------------------------
    # P[n,k] = fraction of k-body subsets violated
    #          by an n-excitation state
    # --------------------------------------------------------

    P = zeros(Float64, N, kmax)

    for n in 1:N
        for k in 1:kmax
            P[n, k] = odd_probability(N, n, k)
        end
    end

    # --------------------------------------------------------
    # Linear program
    # --------------------------------------------------------

    model = Model(HiGHS.Optimizer)

    set_silent(model)

    @variable(model, c[1:kmax] >= 0)
    @variable(model, Δ >= 0)


    # --------------------------------------------------------
    # Normalization
    # --------------------------------------------------------

    @constraint(
        model,
        sum(c[k] for k in 1:kmax) == 1
    )


    # --------------------------------------------------------
    # All excited sectors lie above the gap
    # --------------------------------------------------------

    for n in 1:N

        @constraint(
            model,
            term_energy *
            sum(
                P[n, k] * c[k]
                for k in 1:kmax
            )
            >= Δ
        )

    end

    # --------------------------------------------------------
    # Maximize minimum excitation energy
    # --------------------------------------------------------

    @objective(
        model,
        Max,
        Δ
    )

    optimize!(model)

    # --------------------------------------------------------
    # Results
    # --------------------------------------------------------

    c_opt = value.(c)

    gap = value(Δ)

    energies = [
        term_energy *
        sum(
            P[n, k] * c_opt[k]
            for k in 1:kmax
        )
        for n in 1:N
    ]

    return (
        gap = gap,
        coefficients = c_opt,
        energies = energies,
        probabilities = P
    )
end

function fit_loglog(Ns, gaps)

    logN = log.(Float64.(Ns))
    loggap = log.(gaps)

    # log(Δ) = a + b log(N)
    X = hcat(
        ones(length(logN)),
        logN
    )

    coeff = X \ loggap

    intercept = coeff[1]
    slope = coeff[2]

    prefactor = exp(intercept)

    # Predicted values
    loggap_fit = X * coeff
    gap_fit = prefactor .* Float64.(Ns).^slope

    # R² in log-space
    SS_res = sum(
        (loggap .- loggap_fit).^2
    )

    SS_tot = sum(
        (loggap .- mean(loggap)).^2
    )

    R2 = if SS_tot > 0
        1 - SS_res / SS_tot
    else
        NaN
    end

    return (
        slope = slope,
        intercept = intercept,
        prefactor = prefactor,
        R2 = R2,
        fit = gap_fit
    )
end


function print_gap_result(
    label::String,
    N::Int,
    kmax::Int,
    x:: Real,
    result;
    tolerance=1e-9,
)
    coefficients = result.coefficients
    active_orders = findall(c -> abs(c) > tolerance, coefficients)

    println()
    println("="^70)
    println("N = $N")
    println("Scaling: $label")
    println("x: $x")
    println("kmax = $kmax")

    @printf("Optimal gap Δ = %.10f\n", result.gap)

    limiting_sector = argmin(result.energies)

    @printf(
        "Limiting excitation sector: n = %d, energie = %.10f\n",
        limiting_sector,
        result.energies[limiting_sector],
    )

    println("cₖ larger than $tolerance:")

    if isempty(active_orders)
        println("  No cₖ larger than $tolerance")
    else
        println("  k       cₖ (total)       cₖ/binomial(N,k)")
        println("  ---------------------------------------------------------")

        for k in active_orders
            ck = coefficients[k]

            coefficient_per_subset = Float64(
                BigFloat(ck) / BigFloat(binomial(big(N), k))
            )

            @printf(
                "  %-4d    %18.10e       %18.10e\n",
                k,
                ck,
                coefficient_per_subset,
            )
        end
    end

    println("="^70)
end

function push_coefficients!(
    df, N::Int, x, kmax::Int, result;
    family::String = "", tolerance = 1e-9
)
    for k in findall(c -> abs(c) > tolerance, result.coefficients)
        ck = result.coefficients[k]
        push!(df, (
            N = N, family = family, x = x, kmax = kmax, k = k,
            c_k = ck,
            c_k_per_subset = Float64(BigFloat(ck) / BigFloat(binomial(big(N), k)))
        ))
    end
end