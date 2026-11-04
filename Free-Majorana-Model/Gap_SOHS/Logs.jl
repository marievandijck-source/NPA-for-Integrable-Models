using Plots
using JuMP
using HiGHS
using Measures
using CSV
using DataFrames
using Statistics
using Printf

include("Functions.jl")

# ------------------------------------------------------------
# Plotting gaps for kmax = logaritmic function
# ------------------------------------------------------------

# -----------------------------------------------------------
# Parameters
# -----------------------------------------------------------

Ns = 50:20:100

# ------------------------------------------------------------
# Choose x values
#
# kmax = logartimic functions with powers x
# ------------------------------------------------------------

x_min  = 1
x_max  = 4
x_step = 1

xs = x_min:x_step:x_max

colors = palette(:auto, length(xs))

# -----------------------------------------------------------
# Store gaps
# -----------------------------------------------------------

gaps_log = Dict(x => Float64[] for x in xs)

gaps_loglog = Dict(x => Float64[] for x in xs)

gaps_Nlog = Dict(x => Float64[] for x in xs)

# -----------------------------------------------------------
# DataFrame for all numerical results
# -----------------------------------------------------------

data = DataFrame(
    N = Int[],
    family = String[], # Distinction between different logaritmic functions
    x = Int[],
    kmax = Int[],
    gap = Float64[]
)

coeff_data = DataFrame(
    N = Int[], family = String[], x = Int[], kmax = Int[],
    k = Int[], c_k = Float64[], c_k_per_subset = Float64[]
)

# -----------------------------------------------------------
# Compute gaps
# -----------------------------------------------------------

for N in Ns

    println()
    println("----------------------------")
    println("N = $N")
    println("----------------------------")
    flush(stdout)
    for x in xs

        # -----------------------------------------------------------
        # kmax = log(N)^x
        # -----------------------------------------------------------

        kmax_log = round(Int, log(N)^x)

        kmax_log = clamp(kmax_log, 1, N)

        result_log = optimize_gap(N, kmax_log)

        print_gap_result("log(N)^x", N, kmax_log, x, result_log)

        push!(
            gaps_log[x],
            result_log.gap
        )

        push!(
            data,
            (
                N = N,
                family = "log(N)^x",
                x = x,
                kmax = kmax_log,
                gap = result_log.gap
            )
        )

        push_coefficients!(coeff_data, N, x, kmax_log, result_log; family = "log(N)^x")

        # -----------------------------------------------------------
        # kmax = log(log(N))^x
        # -----------------------------------------------------------

        kmax_loglog = round(Int, log(log(N))^x)

        kmax_loglog = clamp(kmax_loglog, 1, N)

        result_loglog = optimize_gap(N, kmax_loglog)

        print_gap_result("log(log(N))^x", N, kmax_loglog, x, result_loglog)

        push!(
            gaps_loglog[x],
            result_loglog.gap
        )

        push!(
            data,
            (
                N = N,
                family = "log(log(N))^x",
                x = x,
                kmax = kmax_loglog,
                gap = result_loglog.gap
            )
        )

        push_coefficients!(coeff_data, N, x, kmax_loglog, result_loglog; family = "log(log(N))^x")

        # -----------------------------------------------------------
        # kmax = N / log(N)^(1/x)
        # -----------------------------------------------------------

        kmax_Nlog = round(Int, N / (log(N)^(1.0 / x)))

        kmax_Nlog = clamp(kmax_Nlog, 1, N)

        result_Nlog = optimize_gap(N, kmax_Nlog)

        print_gap_result("N / log(N)^(1/x)", N, kmax_Nlog, x, result_Nlog)

        push!(gaps_Nlog[x], result_Nlog.gap)

        push!(
            data,
            (
                N = N,
                family = "N/log(N)^(1/x)",
                x = x,
                kmax = kmax_Nlog,
                gap = result_Nlog.gap
            )
        )

        push_coefficients!(coeff_data, N, x, kmax_Nlog, result_Nlog; family = "N/log(N)^(1/x)")

    end
    flush(stdout)
end

# -----------------------------------------------------------
# Create results folder
# -----------------------------------------------------------

mkpath("results")

# -----------------------------------------------------------
# Save all numerical data
# -----------------------------------------------------------

CSV.write("results/gap_scaling_all_data.csv", data)

CSV.write("results/gap_scaling_coefficients.csv", coeff_data)

# -----------------------------------------------------------
# Linear fits on log-log scale
# -----------------------------------------------------------

fit_results = DataFrame(
    family = String[],
    x = Int[],
    slope = Float64[],
    intercept = Float64[],
    prefactor = Float64[],
    R2 = Float64[]
)

# -----------------------------------------------------------
# Store fitted curves
# -----------------------------------------------------------

fits_log = Dict{Int, Vector{Float64}}()
fits_loglog = Dict{Int, Vector{Float64}}()
fits_Nlog = Dict{Int, Vector{Float64}}()

Nfloat = Float64.(collect(Ns))
logN = log.(Nfloat)

# -----------------------------------------------------------
# Perform fits for all three families using fit_loglog
# -----------------------------------------------------------

for x in xs

    # --------------------------------------------------------
    # 1. kmax = log(N)^x
    # --------------------------------------------------------

    fit1 = fit_loglog(
        collect(Ns),
        gaps_log[x]
    )

    fits_log[x] = fit1.fit

    push!(
        fit_results,
        (
            family = "log(N)^x",
            x = x,
            slope = fit1.slope,
            intercept = fit1.intercept,
            prefactor = fit1.prefactor,
            R2 = fit1.R2
        )
    )


    # --------------------------------------------------------
    # 2. kmax = log(log(N))^x
    # --------------------------------------------------------

    fit2 = fit_loglog(
        collect(Ns),
        gaps_loglog[x]
    )

    fits_loglog[x] = fit2.fit

    push!(
        fit_results,
        (
            family = "log(log(N))^x",
            x = x,
            slope = fit2.slope,
            intercept = fit2.intercept,
            prefactor = fit2.prefactor,
            R2 = fit2.R2
        )
    )


    # --------------------------------------------------------
    # 3. kmax = N / log(N)^(1/x)
    # --------------------------------------------------------

    fit3 = fit_loglog(
        collect(Ns),
        gaps_Nlog[x]
    )

    fits_Nlog[x] = fit3.fit

    push!(
        fit_results,
        (
            family = "N/log(N)^(1/x)",
            x = x,
            slope = fit3.slope,
            intercept = fit3.intercept,
            prefactor = fit3.prefactor,
            R2 = fit3.R2
        )
    )
end


# -----------------------------------------------------------
# Save fit parameters
# -----------------------------------------------------------
mkpath("results")

CSV.write(
    "results/loglog_fit_parameters_logs.csv",
    fit_results
)

println()
println("Fit parameters:")
println(fit_results)


# -----------------------------------------------------------
# Plot 1:
# kmax = log(N)^x
# -----------------------------------------------------------

p1 = plot(
    xlabel = "System size N",
    ylabel = "Energy gap Δ",
    title = "kₘₐₓ = log(N)^x",

    xscale = :log10,
    yscale = :log10,

    framestyle = :box,
    legend = :topright,

    guidefontsize = 12,
    tickfontsize = 10,
    titlefontsize = 13,
    legendfontsize = 9,

    left_margin   = 10mm,
    right_margin  = 4mm,
    bottom_margin = 10mm,
    top_margin    = 7mm
)

for (i, x) in enumerate(xs)
    scatter!(
        p1,
        Ns,
        gaps_log[x],
        label = "x = $x",
        markersize = 4,
        marker = :circle,
        color = colors[i]
    )

    plot!(
        p1,
        Ns,
        fits_log[x],
        label = false,
        linewidth = 1.5,
        linestyle = :dash,
        color = colors[i]
    )
end

# -----------------------------------------------------------
# Plot 2:
# kmax = log(log(N))^x
# -----------------------------------------------------------

p2 = plot(
    xlabel = "System size N",
    ylabel = "Energy gap Δ",
    title = "kₘₐₓ = log(log(N))^x",

    xscale = :log10,
    yscale = :log10,

    framestyle = :box,
    legend = :topright,

    guidefontsize = 12,
    tickfontsize = 10,
    titlefontsize = 13,
    legendfontsize = 9,

    left_margin   = 10mm,
    right_margin  = 4mm,
    bottom_margin = 10mm,
    top_margin    = 7mm
)

for (i, x) in enumerate(xs)
    scatter!(
        p2,
        Ns,
        gaps_loglog[x],
        label = "x = $x",
        markersize = 4,
        marker = :circle,
        color = colors[i]
    )

       plot!(
        p2,
        Ns,
        fits_loglog[x],
        label = false,
        linewidth = 1.5,
        linestyle = :dash,
        color = colors[i]
    )
end


# -----------------------------------------------------------
# Plot 3:
# kmax = N / log(N)^(1/x)
# -----------------------------------------------------------

p3 = plot(
    xlabel = "System size N",
    ylabel = "Energy gap Δ",
    title = "kₘₐₓ = N / log(N)^(1/x)",

    xscale = :log10,
    yscale = :log10,

    framestyle = :box,
    legend = :topright,

    guidefontsize = 12,
    tickfontsize = 10,
    titlefontsize = 13,
    legendfontsize = 9,

    left_margin   = 10mm,
    right_margin  = 5mm,
    bottom_margin = 10mm,
    top_margin    = 7mm
)

for (i, x) in enumerate(xs)
    scatter!(
        p3,
        Ns,
        gaps_Nlog[x],
        label = "x = $x",
        markersize = 4,
        marker = :circle,
        color = colors[i]
    )

     plot!(
        p3,
        Ns,
        fits_Nlog[x],
        label = false,
        linewidth = 1.5,
        linestyle = :dash,
        color = colors[i]
    )
end


# -----------------------------------------------------------
# Combined figure
# -----------------------------------------------------------

p_all = plot(
    p1,
    p2,
    p3,

    layout = (1, 3),

    # More vertical and horizontal space
    size = (1800, 700),

    # Outer margin of complete figure
    left_margin   = 5mm,
    right_margin  = 5mm,
    bottom_margin = 15mm,
    top_margin    = 5mm
)

display(p_all)

# -----------------------------------------------------------
# Save combined figure
# -----------------------------------------------------------

savefig(
    p_all,
    "results/gap_scaling_three_families.png"
)

savefig(
    p_all,
    "results/gap_scaling_three_families.pdf"
)

println("Results saved in results/")