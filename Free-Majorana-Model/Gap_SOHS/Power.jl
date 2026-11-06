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
# Plotting gaps for kmax = N^x
# ------------------------------------------------------------

# -----------------------------------------------------------
# Parameters
# -----------------------------------------------------------

Ns = 50:20:200

# ------------------------------------------------------------
# Choose x values
#
# kmax = N^x
# ------------------------------------------------------------

x_min  = 0.1
x_max  = 1
x_step = 0.1

xs = x_min:x_step:x_max

colors = palette(:auto, length(xs))

# ------------------------------------------------------------
# Store gaps
# ------------------------------------------------------------

gaps = Dict(x => Float64[] for x in xs)

# ------------------------------------------------------------
# DataFrame for all numerical results
# ------------------------------------------------------------

data = DataFrame(
    N = Int[],
    family = String[],
    x = Float64[],
    kmax = Int[],
    gap = Float64[]
)

coeff_data = DataFrame(
    N = Int[], family = String[], x = Float64[], kmax = Int[],
    k = Int[], c_k = Float64[], c_k_per_subset = Float64[]
)

# -----------------------------------------------------------
# Compute gaps
# -----------------------------------------------------------

for N in Ns

    flush(stdout)

    for x in xs

        # -----------------------------------------------------------
        # kmax = N^x
        # -----------------------------------------------------------

        kmax = round(Int, N^x)

        kmax = clamp(kmax, 1, N)

        result = optimize_gap(N, kmax)

        print_gap_result("N^x", N, kmax, x, result)

        push!(gaps[x], result.gap)

        push!(
            data,
            (
                N = N,
                family = "N^x",
                x = x,
                kmax = kmax,
                gap = result.gap
            ))
        
        push_coefficients!(coeff_data, N, x, kmax, result; family = "N^x")
        
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

CSV.write("results/gap_scaling_N_power_x_all_data.csv", data)

CSV.write("results/gap_scaling_N_power_x_coefficients.csv", coeff_data)

# -----------------------------------------------------------
# Linear scale plot
# -----------------------------------------------------------

p1 = plot(
    xlabel = "System size N", ylabel = "Energy gap Δ",
    title = "Linear scale",
    legend = false, framestyle = :box,
    ylims = (0.0, 1.05),
    guidefontsize = 12, tickfontsize = 10, titlefontsize = 14,
    left_margin = 8mm, right_margin = 8mm, bottom_margin = 10mm, top_margin = 5mm
)

for (i, x) in enumerate(xs)
    scatter!(p1, Ns, gaps[x], label = "kₘₐₓ = N^$x",
             markersize = 4, marker = :circle, color = colors[i])
end

# -----------------------------------------------------------
# Linear fits on log-log scale
# -----------------------------------------------------------

fit_results = DataFrame(
    family = String[], x = Float64[], slope = Float64[],
    intercept = Float64[], prefactor = Float64[], R2 = Float64[]
)

fits = Dict{Float64, Vector{Float64}}()

for x in xs
    fit = fit_loglog(collect(Ns), gaps[x])
    fits[x] = fit.fit
    push!(fit_results, (
        family = "N^x", x = x, slope = fit.slope,
        intercept = fit.intercept, prefactor = fit.prefactor, R2 = fit.R2
    ))
end

CSV.write("results/loglog_fit_parameters_N_power_x.csv", fit_results)

println()
println("Fit parameters:")
println(fit_results)

# -----------------------------------------------------------
# Log-log figure
# -----------------------------------------------------------

p2 = plot(
    xlabel = "System size N",
    ylabel = "Energy gap Δ",
    title = "Logarithmic scale",

    legend = :outerright,
    framestyle = :box,

    xscale = :log10,
    yscale = :log10,

    guidefontsize = 12,
    tickfontsize = 10,
    titlefontsize = 14,
    legendfontsize = 9,

    left_margin   = 10mm,
    right_margin  = 10mm,
    bottom_margin = 10mm,
    top_margin    = 5mm
)

for (i,x) in enumerate(xs)
    scatter!(
        p2,
        Ns,
        gaps[x],
        label = "kₘₐₓ = N^$x",
        markersize = 4,
        marker = :circle,
        color = colors[i]
    )
    # Power-law fit
    plot!(
        p2,
        Ns,
        fits[x],
        label = false,
        linewidth = 1.5,
        linestyle = :dash,
        color = colors[i]
    )
end


# -----------------------------------------------------------
# Combined figure
# -----------------------------------------------------------

p = plot(
    p1,
    p2,
    layout = grid(
        1, 2,
        widths = [0.44, 0.56]
    ),
    size = (1550, 550),
    margin = 5mm
)

display(p)


savefig(p, "results/gap_scaling_N_power_x.png")
savefig(p, "results/gap_scaling_N_power_x.pdf")