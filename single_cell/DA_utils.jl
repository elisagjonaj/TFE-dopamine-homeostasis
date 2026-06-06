#=
================================================================================
DA_utils.jl - Utility Functions for DA Neuron Model Analysis
================================================================================

This file is taken directly from Arthur Fyon's codebase and was not modified 
as part of this work. It contains utility functions for analyzing dopaminergic 
neuron model simulations, including:

    - Numerical algorithms (bisection)
    - Spike detection and firing pattern analysis
    - Interspike interval (ISI) extraction
    - Burstiness quantification
    - Visualization tools for scatter matrices and PCA heatmaps

Dependencies:
    Statistics, Plots, StatsPlots, LaTeXStrings, Printf

Author: Arthur Fyon
=#

using Statistics, Plots, StatsPlots, LaTeXStrings, Printf

# ==============================================================================
# Numerical Algorithms
# ==============================================================================

"""
    bisection(f, a, b, tol=1e-6, max_iter=100)

Find the root of function `f` in the interval [a, b] using the bisection method.

# Arguments
- `f`        : Function to find root of
- `a`        : Left endpoint of interval
- `b`        : Right endpoint of interval
- `tol`      : Tolerance for convergence (default: 1e-6)
- `max_iter` : Maximum number of iterations (default: 100)

# Returns
- `mid`  : Approximate root
- `iter` : Number of iterations performed

# Throws
- `ErrorException` if f(a) and f(b) have the same sign
"""
function bisection(f, a, b, tol=1e-6, max_iter=100)
    # Validate initial interval
    if f(a) * f(b) > 0
        error("The function must change sign on the interval [a, b]")
    end

    # Initialize
    mid = (a + b) / 2.0
    iter = 0

    # Bisection iteration
    while abs(f(mid)) > tol && iter < max_iter
        mid = (a + b) / 2.0

        # Update interval based on sign
        if f(a) * f(mid) < 0
            b = mid
        else
            a = mid
        end

        iter += 1
    end

    return mid, iter
end

# ==============================================================================
# Firing Pattern Analysis Functions
# ==============================================================================

"""
    extract_frequency(V, t)

Extract the mean spiking frequency from a voltage trace.

# Arguments
- `V` : Membrane potential time series (mV)
- `t` : Time vector (ms)

# Returns
- Spiking frequency in Hz, or `NaN` if neuron is silent (< 2 spikes)

# Algorithm
Uses threshold crossing detection with:
- Spike onset threshold: +10 mV
- Spike offset threshold: 0 mV
"""
function extract_frequency(V, t)
    # Thresholds for spike detection
    spike_up_threshold = 10.0
    spike_down_threshold = 0.0

    # Detect spike times
    spike_detected = 0
    spike_times = Float64[]
    
    for i in 1:length(V)
        if V[i] > spike_up_threshold && spike_detected == 0
            push!(spike_times, t[i])
            spike_detected = 1
        end
        if V[i] < spike_down_threshold && spike_detected == 1
            spike_detected = 0
        end
    end

    # Check for sufficient spikes
    if length(spike_times) < 2
        return NaN
    end

    # Compute interspike intervals
    ISI = diff(spike_times)

    if length(ISI) < 2
        return NaN
    end

    # Compute mean frequency
    T = mean(ISI) / 1000  # Convert to seconds
    f = 1 / T              # Frequency in Hz

    return f
end

"""
    extract_ISI(V, t)

Extract all interspike intervals from a voltage trace.

# Arguments
- `V` : Membrane potential time series (mV)
- `t` : Time vector (ms)

# Returns
- Vector of interspike intervals (ms), or `NaN` if insufficient spikes

# Note
Uses lower thresholds than `extract_frequency`:
- Spike onset threshold: -20 mV
- Spike offset threshold: -30 mV
"""
function extract_ISI(V, t)
    # Thresholds for spike detection (lower than extract_frequency)
    spike_up_threshold = -20.0
    spike_down_threshold = -30.0

    # Detect spike times
    spike_detected = 0
    spike_times = Float64[]
    
    for i in 1:length(V)
        if V[i] > spike_up_threshold && spike_detected == 0
            push!(spike_times, t[i])
            spike_detected = 1
        end
        if V[i] < spike_down_threshold && spike_detected == 1
            spike_detected = 0
        end
    end

    # Check for sufficient spikes
    if length(spike_times) < 2
        return NaN
    end

    # Compute interspike intervals
    ISI = diff(spike_times)

    if length(ISI) < 2
        return NaN
    end

    return ISI
end

"""
    extract_burstiness(V, t)

Extract bursting characteristics from a voltage trace.

# Arguments
- `V` : Membrane potential time series (mV)
- `t` : Time vector (ms)

# Returns
Tuple of four values:
- `burstiness`        : Burst intensity metric (spikes/burst × intraburst_f / interburst_T)
- `nb_spike_per_burst`: Mean number of spikes per burst
- `intraburst_f`      : Mean frequency within bursts (Hz)
- `interburst_f`      : Burst frequency (Hz)

Returns `(NaN, NaN, NaN, NaN)` if the neuron is not bursting.

# Algorithm
1. Detects all spikes using threshold crossing
2. Computes ISIs and separates them into intraburst/interburst using the midpoint
3. Identifies burst boundaries and computes statistics
"""
function extract_burstiness(V, t)
    # Thresholds for spike detection
    spike_up_threshold = 10.0
    spike_down_threshold = 0.0

    # Detect spike times
    spike_detected = 0
    spike_times = Float64[]
    
    for i in 1:length(V)
        if V[i] > spike_up_threshold && spike_detected == 0
            push!(spike_times, t[i])
            spike_detected = 1
        end
        if V[i] < spike_down_threshold && spike_detected == 1
            spike_detected = 0
        end
    end

    # Check for at least one spike
    if length(spike_times) < 1
        return NaN, NaN, NaN, NaN
    end

    # Compute interspike intervals
    ISI = diff(spike_times)

    # Threshold to separate intraburst from interburst ISIs
    max_ISI = maximum(ISI)
    min_ISI = minimum(ISI)
    half_ISI = (max_ISI + min_ISI) / 2

    # If ISI too uniform, neuron is tonic spiking (not bursting)
    if max_ISI - min_ISI < 10
        return NaN, NaN, NaN, NaN
    end

    # Find indices marking the first spike of each burst (long preceding ISI)
    first_spike_burst = findall(x -> x > half_ISI, ISI)

    # Compute interburst frequency
    Ts = ISI[first_spike_burst]
    interburst_T = mean(Ts) / 1000  # in seconds
    interburst_f = 1 / interburst_T  # in Hz

    # Compute number of spikes per burst
    nb_spike_burst = diff(first_spike_burst)

    # Check for sufficient bursts
    if length(nb_spike_burst) < 2
        return NaN, NaN, NaN, NaN
    end
    
    nb_spike_per_burst = round(mean(nb_spike_burst))

    # Validate bursting behavior
    if nb_spike_per_burst < 1.5 || nb_spike_per_burst > 500
        return NaN, NaN, NaN, NaN
    end

    # Compute intraburst frequency
    intra_spike_burst = findall(x -> x < half_ISI, ISI)
    Ts_intraburst = ISI[intra_spike_burst]
    T_intraburst = mean(Ts_intraburst) / 1000  # in seconds
    intraburst_f = 1 / T_intraburst             # in Hz

    # Compute burstiness metric
    burstiness = (nb_spike_per_burst * intraburst_f) / interburst_T

    return burstiness, nb_spike_per_burst, intraburst_f, interburst_f
end

# ==============================================================================
# Visualization Functions
# ==============================================================================

"""
    heatmap_dir(dir_val, nb_channels)

Create a heatmap visualization of principal component loadings.

# Arguments
- `dir_val`     : Eigen decomposition result (with `.values` and `.vectors` fields)
- `nb_channels` : Number of channels/dimensions

# Returns
- Plots.jl combined plot object showing PC loadings as heatmaps

# Description
Displays the absolute loading weights for each principal component as a series
of heatmaps, with variance explained shown as percentage labels.
"""
function heatmap_dir(dir_val, nb_channels)
    total_var = sum(dir_val.values)
    eig_val_decreasing = reverse(dir_val.values) ./ total_var
    
    # First principal component
    bin = 1
    val = eig_val_decreasing[1] * 100
    str_val = @sprintf "%d" val
    p1 = heatmap(1:1, 1:nb_channels, 
                 reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ 
                        norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]),
                 colorbar=false, clim=(0, 1), tickfontsize=18, guidefontsize=15,
                 yticks=(1:nb_channels, reverse([L"\bar{g}_\mathrm{Na}", L"\bar{g}_\mathrm{Kd}", 
                        L"\bar{g}_\mathrm{CaL}", L"\bar{g}_\mathrm{CaN}", L"\bar{g}_\mathrm{ERG}", 
                        L"g_\mathrm{leak}"])))

    # Second principal component
    bin = 2
    val = eig_val_decreasing[2] * 100
    str_val = @sprintf "%d" val
    p2 = heatmap(1:1, 1:nb_channels, 
                 reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ 
                        norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Third principal component
    bin = 3
    val = eig_val_decreasing[3] * 100
    str_val = @sprintf "%d" val
    p3 = heatmap(1:1, 1:nb_channels, 
                 reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ 
                        norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Fourth principal component
    bin = 4
    val = eig_val_decreasing[4] * 100
    str_val = @sprintf "%d" val
    p4 = heatmap(1:1, 1:nb_channels, 
                 reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ 
                        norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Fifth principal component
    bin = 5
    val = eig_val_decreasing[5] * 100
    str_val = @sprintf "%d" val
    p5 = heatmap(1:1, 1:nb_channels, 
                 reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ 
                        norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Sixth principal component
    bin = 6
    val = eig_val_decreasing[6] * 100
    str_val = @sprintf "%d" val
    p6 = heatmap(1:1, 1:nb_channels, 
                 reverse(reshape(abs.(dir_val.vectors[:, nb_channels - bin + 1]) ./ 
                        norm(dir_val.vectors[:, nb_channels - bin + 1]), nb_channels, 1)),
                 grid=false, xlabel=L"%$str_val\%",
                 axis=false, ticks=false, c=cgrad([:gray93, :orangered3]), guidefontsize=15,
                 colorbar=false, clim=(0, 1))

    # Colorbar
    colors = -1.0:0.002:1.0
    p90 = heatmap(1:1, colors, reshape(colors, length(colors), 1), size=(200, 500),
                  grid=false, axis=false, xticks=false, colorbar=false,
                  c=cgrad([:gray93, :orangered3]), clim=(0, 1), ymirror=true,
                  yticks=(0:0.2:1, [L"0\%", L"20\%", L"40\%", L"60\%", L"80\%", L"100\%"]), 
                  ylims=(0, 1), yaxis=true)

    # Combine plots
    CC = plot(p1, p2, p3, p4, p5, p6, layout=(1, nb_channels), size=(600, 500))

    return CC
end

"""
    scatter_matrix3x3(g_all, maxs, color_p, m_shape, names; kwargs...)

Create a 3×3 scatter plot matrix for visualizing conductance correlations.

# Arguments
- `g_all`    : Matrix of conductance values (N × 4)
- `maxs`     : Vector of axis limits for each dimension
- `color_p`  : Point color
- `m_shape`  : Marker shape
- `names`    : Vector of dimension names (LaTeX strings)

# Keyword Arguments
- `flag`     : Overlay type (0=correlation line, 1=PC1, 2=PC1+PC2)
- `dir_val`  : Eigen decomposition for PC overlay
- `mean_vec` : Mean vector for PC overlay
- `s1`       : Scale factor for PC1
- `s2`       : Scale factor for PC2
- `flag2`    : If 1, overlay second dataset
- `g_all2`   : Second dataset for overlay
- `color_p2` : Color for second dataset
- `m_shape2` : Marker shape for second dataset

# Returns
- Plots.jl combined plot object
"""
function scatter_matrix3x3(g_all, maxs, color_p, m_shape, names; 
                           flag=0, dir_val=Nothing, mean_vec=Nothing, 
                           s1=Nothing, flag2=0, g_all2=Nothing, 
                           color_p2=:gray70, m_shape2=Nothing, s2=Nothing)
    
    cors = NaN * ones(3, 3)
    msw_main = m_shape == :cross ? 1 : 0

    # Panel (1,2): dimension 1 vs 2
    p12 = scatter(g_all[:, 1], g_all[:, 2], label="", markerstrokewidth=msw_main, 
                  color=color_p, top_margin=12Plots.mm,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 1], g_all2[:, 2], label="", color=color_p2, markerstrokewidth=0.0,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    annotate!(maxs[1]/2, maxs[2]*1.3, Plots.text(names[1], :black, :center, 18))
    xlims!((0, maxs[1]))
    ylims!((0, maxs[2]))

    if flag == 0
        line_12 = fit(g_all[:, 1], g_all[:, 2], 1)
        s0, sn = minimum(g_all[:, 1]), maximum(g_all[:, 1])
        plot!([s0, sn], [line_12(s0), line_12(sn)], linewidth=2, label="", linecolor="black")
        cors[1, 1] = cor(g_all[:, 1], g_all[:, 2])
    elseif flag == 1
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, 
               mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, 
               mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, 
               mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, 
               mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values, 
               mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values],
              [mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values, 
               mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end

    # Panel (1,3): dimension 1 vs 3
    p13 = scatter(g_all[:, 1], g_all[:, 3], label="", markerstrokewidth=msw_main, color=color_p,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 1], g_all2[:, 3], label="", color=color_p2, markerstrokewidth=0.0,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    xlims!((0, maxs[1]))
    ylims!((0, maxs[3]))

    if flag == 0
        line_13 = fit(g_all[:, 1], g_all[:, 3], 1)
        s0, sn = minimum(g_all[:, 1]), maximum(g_all[:, 1])
        plot!([s0, sn], [line_13(s0), line_13(sn)], linewidth=2, label="", linecolor="black")
        cors[2, 1] = cor(g_all[:, 1], g_all[:, 3])
    elseif flag == 1
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, 
               mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, 
               mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, 
               mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, 
               mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values, 
               mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values],
              [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values, 
               mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end

    # Panel (1,4): dimension 1 vs 4
    p14 = scatter(g_all[:, 1], g_all[:, 4], label="", markerstrokewidth=msw_main, color=color_p,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 1], g_all2[:, 4], label="", color=color_p2, markerstrokewidth=0.0,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    xlims!((0, maxs[1]))
    ylims!((0, maxs[4]))

    if flag == 0
        line_14 = fit(g_all[:, 1], g_all[:, 4], 1)
        s0, sn = minimum(g_all[:, 1]), maximum(g_all[:, 1])
        plot!([s0, sn], [line_14(s0), line_14(sn)], linewidth=2, label="", linecolor="black")
        cors[3, 1] = cor(g_all[:, 1], g_all[:, 4])
    elseif flag == 1
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, 
               mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, 
               mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, 
               mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, 
               mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[1] - s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values, 
               mean_vec[1] + s2*dir_val.vectors[:, nb_channels-1][1]*dir_val.values],
              [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values, 
               mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end

    # Panel (2,3): dimension 2 vs 3
    p23 = scatter(g_all[:, 2], g_all[:, 3], label="", markerstrokewidth=msw_main, color=color_p,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 2], g_all2[:, 3], label="", color=color_p2, markerstrokewidth=0.0,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    xlims!((0, maxs[2]))
    ylims!((0, maxs[3]))

    if flag == 0
        line_23 = fit(g_all[:, 2], g_all[:, 3], 1)
        s0, sn = minimum(g_all[:, 2]), maximum(g_all[:, 2])
        plot!([s0, sn], [line_23(s0), line_23(sn)], linewidth=2, label="", linecolor="black")
        cors[2, 2] = cor(g_all[:, 2], g_all[:, 3])
    elseif flag == 1
        plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, 
               mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, 
               mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, 
               mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, 
               mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values, 
               mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values],
              [mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values, 
               mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end

    # Panel (2,4): dimension 2 vs 4
    p24 = scatter(g_all[:, 2], g_all[:, 4], label="", markerstrokewidth=msw_main, color=color_p,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 2], g_all2[:, 4], label="", color=color_p2, markerstrokewidth=0.0,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    xlims!((0, maxs[2]))
    ylims!((0, maxs[4]))

    if flag == 0
        line_24 = fit(g_all[:, 2], g_all[:, 4], 1)
        s0, sn = minimum(g_all[:, 2]), maximum(g_all[:, 2])
        plot!([s0, sn], [line_24(s0), line_24(sn)], linewidth=2, label="", linecolor="black")
        cors[3, 2] = cor(g_all[:, 2], g_all[:, 4])
    elseif flag == 1
        plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, 
               mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, 
               mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, 
               mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, 
               mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[2] - s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values, 
               mean_vec[2] + s2*dir_val.vectors[:, nb_channels-1][2]*dir_val.values],
              [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values, 
               mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end

    # Panel (3,4): dimension 3 vs 4
    p34 = scatter(g_all[:, 3], g_all[:, 4], label="", markerstrokewidth=msw_main, color=color_p, 
                  right_margin=15Plots.mm,
                  grid=false, ticks=false, tickfontsize=10, markershape=m_shape, guidefontsize=18)
    if flag2 == 1
        scatter!(g_all2[:, 3], g_all2[:, 4], label="", color=color_p2, markerstrokewidth=0.0,
                 grid=false, ticks=false, tickfontsize=10, markershape=m_shape2, guidefontsize=18)
    end
    annotate!(maxs[3]*1.3, maxs[4]/2, Plots.text(names[4], :black, :center, 18))
    xlims!((0, maxs[3]))
    ylims!((0, maxs[4]))

    if flag == 0
        line_34 = fit(g_all[:, 3], g_all[:, 4], 1)
        s0, sn = minimum(g_all[:, 3]), maximum(g_all[:, 3])
        plot!([s0, sn], [line_34(s0), line_34(sn)], linewidth=2, label="", linecolor="black")
        cors[3, 3] = cor(g_all[:, 3], g_all[:, 4])
        display(cors)
    elseif flag == 1
        plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, 
               mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, 
               mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    elseif flag == 2
        plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, 
               mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, 
               mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
        plot!([mean_vec[3] - s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values, 
               mean_vec[3] + s2*dir_val.vectors[:, nb_channels-1][3]*dir_val.values],
              [mean_vec[4] - s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values, 
               mean_vec[4] + s2*dir_val.vectors[:, nb_channels-1][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:dash)
    end

    # Label panels
    p21 = plot(axis=false, ticks=false, labels=false)
    xlims!((-1, 1))
    ylims!((-1, 1))
    annotate!(0, 0, Plots.text(names[2], :black, :center, 18))

    p32 = plot(axis=false, ticks=false, labels=false)
    xlims!((-1, 1))
    ylims!((-1, 1))
    annotate!(0, 0, Plots.text(names[3], :black, :center, 18))

    # Combine into scatter matrix
    CC = plot(p12, p21, p13, p23, p32, p14, p24, p34, size=(500, 500),
              layout=@layout([° ° _; ° ° °; ° ° °]), margin=3Plots.mm)

    return CC
end

"""
    scatter_matrix3x3_zcolor(g_all, maxs, Rin, m_shape, names; kwargs...)

Create a 3×3 scatter plot matrix with color-coded points by input resistance.

# Arguments
- `g_all`   : Matrix of conductance values (N × 4)
- `maxs`    : Vector of axis limits for each dimension
- `Rin`     : Input resistance values for color coding
- `m_shape` : Marker shape
- `names`   : Vector of dimension names (LaTeX strings)

# Keyword Arguments
- `flag`     : Overlay type (0=correlation line, 1=PC1)
- `dir_val`  : Eigen decomposition for PC overlay
- `mean_vec` : Mean vector for PC overlay
- `s1`       : Scale factor for PC1

# Returns
- Plots.jl combined plot object with thermal colormap for Rin
"""
function scatter_matrix3x3_zcolor(g_all, maxs, Rin, m_shape, names; 
                                   flag=0, dir_val=Nothing, mean_vec=Nothing, s1=Nothing)
    
    cors = NaN * ones(3, 3)
    msw_main = m_shape == :cross ? 1 : 0

    # Panel (1,2): dimension 1 vs 2
    p12 = scatter(g_all[:, 1], g_all[:, 2], label="", markerstrokewidth=msw_main, zcolor=Rin, 
                  top_margin=12Plots.mm, legend=false,
                  c=cgrad(:thermal, rev=false), clims=(13, 45), grid=false, ticks=false, 
                  tickfontsize=10, markershape=m_shape, guidefontsize=18)
    annotate!(maxs[1]/2, maxs[2]*1.3, Plots.text(names[1], :black, :center, 18))
    xlims!((0, maxs[1]))
    ylims!((0, maxs[2]))

    if flag == 0
        line_12 = fit(g_all[:, 1], g_all[:, 2], 1)
        s0, sn = minimum(g_all[:, 1]), maximum(g_all[:, 1])
        plot!([s0, sn], [line_12(s0), line_12(sn)], linewidth=2, label="", linecolor="black")
        cors[1, 1] = cor(g_all[:, 1], g_all[:, 2])
    elseif flag == 1
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, 
               mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, 
               mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    end

    # Panel (1,3): dimension 1 vs 3
    p13 = scatter(g_all[:, 1], g_all[:, 3], label="", markerstrokewidth=msw_main, zcolor=Rin, 
                  legend=false,
                  c=cgrad(:thermal, rev=false), clims=(13, 45), grid=false, ticks=false, 
                  tickfontsize=10, markershape=m_shape, guidefontsize=18)
    xlims!((0, maxs[1]))
    ylims!((0, maxs[3]))

    if flag == 0
        line_13 = fit(g_all[:, 1], g_all[:, 3], 1)
        s0, sn = minimum(g_all[:, 1]), maximum(g_all[:, 1])
        plot!([s0, sn], [line_13(s0), line_13(sn)], linewidth=2, label="", linecolor="black")
        cors[2, 1] = cor(g_all[:, 1], g_all[:, 3])
    elseif flag == 1
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, 
               mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, 
               mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    end

    # Panel (1,4): dimension 1 vs 4
    p14 = scatter(g_all[:, 1], g_all[:, 4], label="", markerstrokewidth=msw_main, zcolor=Rin, 
                  legend=false,
                  c=cgrad(:thermal, rev=false), clims=(13, 45), grid=false, ticks=false, 
                  tickfontsize=10, markershape=m_shape, guidefontsize=18)
    xlims!((0, maxs[1]))
    ylims!((0, maxs[4]))

    if flag == 0
        line_14 = fit(g_all[:, 1], g_all[:, 4], 1)
        s0, sn = minimum(g_all[:, 1]), maximum(g_all[:, 1])
        plot!([s0, sn], [line_14(s0), line_14(sn)], linewidth=2, label="", linecolor="black")
        cors[3, 1] = cor(g_all[:, 1], g_all[:, 4])
    elseif flag == 1
        plot!([mean_vec[1] - s1*dir_val.vectors[:, nb_channels][1]*dir_val.values, 
               mean_vec[1] + s1*dir_val.vectors[:, nb_channels][1]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, 
               mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    end

    # Panel (2,3): dimension 2 vs 3
    p23 = scatter(g_all[:, 2], g_all[:, 3], label="", markerstrokewidth=msw_main, zcolor=Rin, 
                  legend=false,
                  c=cgrad(:thermal, rev=false), clims=(13, 45), grid=false, ticks=false, 
                  tickfontsize=10, markershape=m_shape, guidefontsize=18)
    xlims!((0, maxs[2]))
    ylims!((0, maxs[3]))

    if flag == 0
        line_23 = fit(g_all[:, 2], g_all[:, 3], 1)
        s0, sn = minimum(g_all[:, 2]), maximum(g_all[:, 2])
        plot!([s0, sn], [line_23(s0), line_23(sn)], linewidth=2, label="", linecolor="black")
        cors[2, 2] = cor(g_all[:, 2], g_all[:, 3])
    elseif flag == 1
        plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, 
               mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              [mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, 
               mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    end

    # Panel (2,4): dimension 2 vs 4
    p24 = scatter(g_all[:, 2], g_all[:, 4], label="", markerstrokewidth=msw_main, zcolor=Rin, 
                  legend=false,
                  c=cgrad(:thermal, rev=false), clims=(13, 45), grid=false, ticks=false, 
                  tickfontsize=10, markershape=m_shape, guidefontsize=18)
    xlims!((0, maxs[2]))
    ylims!((0, maxs[4]))

    if flag == 0
        line_24 = fit(g_all[:, 2], g_all[:, 4], 1)
        s0, sn = minimum(g_all[:, 2]), maximum(g_all[:, 2])
        plot!([s0, sn], [line_24(s0), line_24(sn)], linewidth=2, label="", linecolor="black")
        cors[3, 2] = cor(g_all[:, 2], g_all[:, 4])
    elseif flag == 1
        plot!([mean_vec[2] - s1*dir_val.vectors[:, nb_channels][2]*dir_val.values, 
               mean_vec[2] + s1*dir_val.vectors[:, nb_channels][2]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, 
               mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    end

    # Panel (3,4): dimension 3 vs 4
    p34 = scatter(g_all[:, 3], g_all[:, 4], label="", markerstrokewidth=msw_main, zcolor=Rin, 
                  right_margin=15Plots.mm,
                  c=cgrad(:thermal, rev=false), clims=(13, 45), grid=false, ticks=false, 
                  tickfontsize=10, markershape=m_shape, guidefontsize=18, legend=false)
    annotate!(maxs[3]*1.3, maxs[4]/2, Plots.text(names[4], :black, :center, 18))
    xlims!((0, maxs[3]))
    ylims!((0, maxs[4]))

    if flag == 0
        line_34 = fit(g_all[:, 3], g_all[:, 4], 1)
        s0, sn = minimum(g_all[:, 3]), maximum(g_all[:, 3])
        plot!([s0, sn], [line_34(s0), line_34(sn)], linewidth=2, label="", linecolor="black")
        cors[3, 3] = cor(g_all[:, 3], g_all[:, 4])
        display(cors)
    elseif flag == 1
        plot!([mean_vec[3] - s1*dir_val.vectors[:, nb_channels][3]*dir_val.values, 
               mean_vec[3] + s1*dir_val.vectors[:, nb_channels][3]*dir_val.values],
              [mean_vec[4] - s1*dir_val.vectors[:, nb_channels][4]*dir_val.values, 
               mean_vec[4] + s1*dir_val.vectors[:, nb_channels][4]*dir_val.values],
              arrow=false, color=:black, linewidth=2, label="", linestyle=:solid)
    end

    # Label panels
    p21 = plot(axis=false, ticks=false, labels=false)
    xlims!((-1, 1))
    ylims!((-1, 1))
    annotate!(0, 0, Plots.text(names[2], :black, :center, 18))

    p32 = plot(axis=false, ticks=false, labels=false)
    xlims!((-1, 1))
    ylims!((-1, 1))
    annotate!(0, 0, Plots.text(names[3], :black, :center, 18))

    # Combine into scatter matrix
    CC = plot(p12, p21, p13, p23, p32, p14, p24, p34, size=(500, 500),
              layout=@layout([° ° _; ° ° °; ° ° °]), margin=3Plots.mm)

    return CC
end
