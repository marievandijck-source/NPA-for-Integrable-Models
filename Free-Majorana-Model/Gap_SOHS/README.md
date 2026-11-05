# Gap Scaling

Julia code to study how the optimal energy gap Δ scales with the system size N when the interaction order is limited to `kmax`.

For every system size N and every chosen `kmax(N)`, a linear program (LP) finds the weights `c_k` of the k-body terms (k = 1, …, kmax) that maximise the gap.

## Scripts

| File | Description |
|------|-------------|
| `Functions.jl` | Shared functions: `odd_probability`, `optimize_gap`, `fit_loglog`, `print_gap_result`, `push_coefficients!` |
| `Dev.jl` | Scaling for `kmax = N / x` |
| `Power.jl` | Scaling for `kmax = N^x` |
| `Logs.jl` | Scaling for three logarithmic families: `log(N)^x`, `log(log(N))^x` and `N / log(N)^(1/x)` |
|`run_gap.sh`| File to run code on HPC (doduo UGent). All files run in one job.|

`Dev.jl`, `Power.jl` and `Logs.jl` are standalone scripts. Each one loads `Functions.jl` with `include("Functions.jl")`, so keep all files in the same folder.

## Method

`See Free_Majorana_Model_Notes.pdf`

**Optimisation (LP, solved with HiGHS).**

$$
\Delta = \max_{c}\ \min_{n}\ E_n \quad \text{subject to} \quad c_k \ge 0,\ \sum_k c_k = 1
$$

## Dependencies
`Plots`, `JuMP`, `HiGHS`, `Measures`, `CSV`, `DataFrames`, `Statistics`, `Printf` (and `Combinatorics`, which `Dev.jl` imports)

### Parameters

The settings are at the top of each script:

- `Ns`: system sizes, e.g. `Ns = 50:20:200`
- `x_min`, `x_max`, `x_step`: the range of the exponent / divisor `x`

## Output

Every script writes three CSV files and a figure (PNG and PDF) to `results/`.

| Script | Gaps | Coefficients | Fit parameters | Figure |
|--------|------|--------------|----------------|--------|
| `Dev.jl` | `gap_scaling_N_over_x_all_data.csv` | `gap_scaling_N_over_x_coefficients.csv` | `loglog_fit_parameters_N_over_x.csv` | `gap_scaling_N_over_x.png/pdf` |
| `Power.jl` | `gap_scaling_N_power_x_all_data.csv` | `gap_scaling_N_power_x_coefficients.csv` | `loglog_fit_parameters_N_power_x.csv` | `gap_scaling_N_power_x.png/pdf` |
| `Logs.jl` | `gap_scaling_all_data.csv` | `gap_scaling_coefficients.csv` | `loglog_fit_parameters_logs.csv` | `gap_scaling_three_families.png/pdf` |

### Columns

**Gaps (`*_all_data.csv`)**: one row per (N, x) combination

| Column | Meaning |
|--------|---------|
| `N` | system size |
| `family` | scaling family used for `kmax` |
| `x` | exponent / divisor |
| `kmax` | maximal interaction order used |
| `gap` | optimal gap Δ |

**Coefficients (`*_coefficients.csv`)**: one row per (N, x, k) with a non-zero weight (`|c_k| > 1e-9`)

| Column | Meaning |
|--------|---------|
| `N`, `family`, `x`, `kmax` | as above |
| `k` | interaction order |
| `c_k` | optimal total weight of all k-body terms |
| `c_k_per_subset` | `c_k / binomial(N, k)`, the weight per individual k-body subset |

**Fit parameters (`loglog_fit_parameters_*.csv`)**: one row per (family, x)

| Column | Meaning |
|--------|---------|
| `slope` | exponent `b` in Δ ≈ A · N^b |
| `intercept` | `a = log A` |
| `prefactor` | `A = exp(a)` |
| `R2` | goodness of fit in log space |

## Scaling families

| Script | Family | `kmax` |
|--------|--------|--------|
| `Dev.jl` | `N/x` | `round(N / x)` |
| `Power.jl` | `N^x` | `round(N^x)` |
| `Logs.jl` | `log(N)^x` | `round(log(N)^x)` |
| `Logs.jl` | `log(log(N))^x` | `round(log(log(N))^x)` |
| `Logs.jl` | `N/log(N)^(1/x)` | `round(N / log(N)^(1/x))` |

In all cases `kmax` is clamped to the range `[1, N]`.

## Numerical notes

- The binomial coefficients are computed exactly with `BigInt` and converted to `Float64` only after the division, so `P(n,k)` is accurate to machine precision.
- The accuracy of the gap itself is limited by the LP solver tolerance (HiGHS default about 1e-7). For very small gaps this can matter for the fit; tighten it with `set_attribute(model, "primal_feasibility_tolerance", 1e-10)` if needed.
- All gaps must be strictly positive, because the fit takes `log(gap)`.
