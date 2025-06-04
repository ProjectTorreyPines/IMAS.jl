"""
    robust_outlier_removal(data::AbstractMatrix{T}; 
        spatial_window::Tuple{Int,Int}=(3,3),
        threshold::Float64=3.0,
        method::Symbol=:mad,
        boundary_mode::Symbol=:reflect,
        min_neighbors::Int=4,
        max_iterations::Int=1,
        return_mask::Bool=false) where T<:Real

Advanced outlier removal with flexible stencil sizes and robust statistics.

# Arguments

  - `data`: 2D matrix of measurements (space × time)
  - `spatial_window`: Tuple (height, width) defining the neighborhood size
  - `threshold`: Number of robust standard deviations for outlier detection
  - `method`: Statistical method (`:mad`, `:iqr`, `:huber`, `:trimmed_mean`)
  - `boundary_mode`: How to handle edges (`:reflect`, `:symmetric`, `:constant`, `:nearest`)
  - `min_neighbors`: Minimum number of valid neighbors required for replacement
  - `max_iterations`: Number of iterative passes (useful for cluster outliers)

# Returns

  - Cleaned 2D matrix with outliers replaced
  - Mask indicating which points were modified (optional second return)
"""
function robust_outlier_removal(data::AbstractMatrix{T};
    spatial_window::Tuple{Int,Int}=(3, 3),
    threshold::Float64=3.0,
    method::Symbol=:mad,
    boundary_mode::Symbol=:reflect,
    min_neighbors::Int=4,
    max_iterations::Int=1,
    return_mask::Bool=false) where {T<:Real}

    cleaned = copy(data)
    m, n = size(data)
    h, w = spatial_window

    # Ensure odd window sizes for symmetry
    h = h + (h % 2 == 0 ? 1 : 0)
    w = w + (w % 2 == 0 ? 1 : 0)

    h_half = h ÷ 2
    w_half = w ÷ 2

    outlier_mask = falses(m, n)

    println("Outlier removal: $(spatial_window) window, threshold=$threshold, method=$method")

    for iteration in 1:max_iterations
        outliers_found = 0

        for i in 1:m, j in 1:n
            # Extract neighborhood with boundary handling
            neighborhood = extract_neighborhood(cleaned, i, j, h_half, w_half, boundary_mode)

            # Skip if too few neighbors
            valid_neighbors = sum(.!isnan.(neighborhood))
            if valid_neighbors < min_neighbors
                continue
            end

            # Compute robust center and scale
            center, scale = compute_robust_statistics(neighborhood, method)

            # Check if current point is an outlier
            current_val = cleaned[i, j]
            if !isnan(current_val) && abs(current_val - center) > threshold * scale
                # Replace with robust estimate
                cleaned[i, j] = center
                outlier_mask[i, j] = true
                outliers_found += 1
            end
        end

        if outliers_found == 0
            println("Converged after $iteration iteration(s)")
            break
        else
            println("Iteration $iteration: $outliers_found outliers corrected")
        end
    end

    if return_mask
        return cleaned, outlier_mask
    else
        return cleaned
    end
end

"""
    extract_neighborhood(data::AbstractMatrix{T}, i::Int, j::Int, h_half::Int, w_half::Int, boundary_mode::Symbol) where T<:Real

Extract neighborhood around point (i,j) with specified boundary handling.
"""
function extract_neighborhood(data::AbstractMatrix{T}, i::Int, j::Int, h_half::Int, w_half::Int, boundary_mode::Symbol) where {T<:Real}
    m, n = size(data)

    # # Define neighborhood indices
    # i_start, i_end = i - h_half, i + h_half
    # j_start, j_end = j - w_half, j + w_half

    neighborhood = Matrix{Float64}(undef, 2 * h_half + 1, 2 * w_half + 1)

    for di in -h_half:h_half, dj in -w_half:w_half
        ii, jj = i + di, j + dj

        # Handle boundaries
        if boundary_mode == :reflect
            ii = reflect_index(ii, m)
            jj = reflect_index(jj, n)
        elseif boundary_mode == :symmetric
            ii = symmetric_index(ii, m)
            jj = symmetric_index(jj, n)
        elseif boundary_mode == :nearest
            ii = clamp(ii, 1, m)
            jj = clamp(jj, 1, n)
        elseif boundary_mode == :constant
            if ii < 1 || ii > m || jj < 1 || jj > n
                neighborhood[di+h_half+1, dj+w_half+1] = NaN
                continue
            end
        end

        if 1 <= ii <= m && 1 <= jj <= n
            neighborhood[di+h_half+1, dj+w_half+1] = data[ii, jj]
        else
            neighborhood[di+h_half+1, dj+w_half+1] = NaN
        end
    end

    return neighborhood
end

function reflect_index(idx::Int, max_idx::Int)
    if idx < 1
        return 2 - idx
    elseif idx > max_idx
        return 2 * max_idx - idx
    else
        return idx
    end
end

function symmetric_index(idx::Int, max_idx::Int)
    if idx < 1
        return 1 - idx
    elseif idx > max_idx
        return 2 * max_idx - idx + 1
    else
        return idx
    end
end

"""
    compute_robust_statistics(neighborhood::AbstractMatrix{Float64}, method::Symbol)

Compute robust center and scale estimates for outlier detection.

Method can be :mad, :iqr, :huber, or :trimmed_mean
"""
function compute_robust_statistics(neighborhood::AbstractMatrix{Float64}, method::Symbol)
    # Remove NaN values
    valid_vals = neighborhood[.!isnan.(neighborhood)]

    if length(valid_vals) == 0
        return NaN, NaN
    end

    if method == :mad
        # Median Absolute Deviation
        center = Statistics.median(valid_vals)
        scale = Statistics.median(abs.(valid_vals .- center)) * 1.4826  # Scale factor for normal distribution
        scale = max(scale, 1e-8)  # Avoid division by zero

    elseif method == :iqr
        # Interquartile Range
        center = Statistics.median(valid_vals)
        q1, q3 = Statistics.quantile(valid_vals, [0.25, 0.75])
        scale = (q3 - q1) / 1.349  # Scale factor for normal distribution
        scale = max(scale, 1e-8)

    elseif method == :huber
        # Huber robust estimator (iterative)
        center = Statistics.median(valid_vals)
        scale = Statistics.mad(valid_vals, center) * 1.4826
        for _ in 1:5  # Few iterations usually enough
            weights = huber_weights.(abs.(valid_vals .- center) ./ scale)
            center = sum(weights .* valid_vals) / sum(weights)
            residuals = abs.(valid_vals .- center)
            scale = sqrt(sum(weights .* residuals .^ 2) / sum(weights)) * 1.134
            scale = max(scale, 1e-8)
        end

    elseif method == :trimmed_mean
        # Trimmed mean (remove 10% from each tail)
        trim_frac = 0.1
        sorted_vals = sort(valid_vals)
        n_trim = max(1, round(Int, trim_frac * length(sorted_vals)))
        trimmed = sorted_vals[n_trim+1:end-n_trim]
        center = Statistics.mean(trimmed)
        scale = Statistics.std(trimmed)
        scale = max(scale, 1e-8)

    else
        error("Unknown method: $method. Use :mad, :iqr, :huber, or :trimmed_mean")
    end

    return center, scale
end

"""
    Huber weight function for robust estimation
"""
function huber_weights(x::Float64, k::Float64=1.345)
    if abs(x) <= k
        return 1.0
    else
        return k / abs(x)
    end
end

"""
    adaptive_outlier_removal!(data::AbstractMatrix{T};
        min_window::Tuple{Int,Int}=(3, 3),
        max_window::Tuple{Int,Int}=(7, 7),
        base_window::Tuple{Int,Int}=min_window,
        threshold::Float64=3.0,
        adaptivity::Symbol=:gradient) where {T<:Real}

Adaptive outlier removal that adjusts window size based on local data characteristics.
"""
function adaptive_outlier_removal!(data::AbstractMatrix{T};
    min_window::Tuple{Int,Int}=(3, 3),
    max_window::Tuple{Int,Int}=(7, 7),
    base_window::Tuple{Int,Int}=min_window,
    threshold::Float64=3.0,
    adaptivity::Symbol=:gradient) where {T<:Real}

    m, n = size(data)
    cleaned = IMASdd.ThreadSafeDicts.ThreadSafeDict{Tuple{Int,Int},T}()

    #Threads.@threads 
    for i in 1:m, j in 1:n
        # Determine local window size based on adaptivity criterion
        if adaptivity == :gradient
            # Smaller windows in high-gradient regions (preserve edges)
            local_grad = compute_local_gradient(data, i, j)
            window_scale = 1.0 / (1.0 + local_grad)
        elseif adaptivity == :variance
            # Smaller windows in high-variance regions
            local_var = compute_local_variance(data, i, j, base_window)
            window_scale = 1.0 / (1.0 + local_var / Statistics.mean(data)^2)
        else
            window_scale = 1.0
        end
        if isnan(window_scale)
            window_scale = 1.0
        end

        # Scale window size
        h_scaled = round(Int, base_window[1] * window_scale)
        w_scaled = round(Int, base_window[2] * window_scale)

        # Clamp to min/max window sizes
        h_scaled = clamp(h_scaled, min_window[1], max_window[1])
        w_scaled = clamp(w_scaled, min_window[2], max_window[2])

        # Ensure odd sizes
        h_scaled += (h_scaled % 2 == 0 ? 1 : 0)
        w_scaled += (w_scaled % 2 == 0 ? 1 : 0)

        # Apply outlier removal with adaptive window
        h_half = h_scaled ÷ 2
        w_half = w_scaled ÷ 2

        neighborhood = extract_neighborhood(data, i, j, h_half, w_half, :reflect)
        center, scale = compute_robust_statistics(neighborhood, :mad)

        current_val = data[i, j]
        if Float64(current_val) == 0.0 # we assume exact zero is a NaN
            cleaned[i, j] = NaN
        elseif !isnan(current_val) && abs(current_val - center) > threshold * scale
            cleaned[i, j] = center
        end
    end

    for ((i, j), value) in cleaned
        data[i, j] = value
    end

    return data
end

function compute_local_gradient(data::AbstractMatrix{T}, i::Int, j::Int) where {T<:Real}
    m, n = size(data)

    # Compute finite differences
    grad_i = 0.0
    grad_j = 0.0

    if i > 1 && i < m
        grad_i = abs(data[i+1, j] - data[i-1, j]) / 2
    end
    if j > 1 && j < n
        grad_j = abs(data[i, j+1] - data[i, j-1]) / 2
    end

    return sqrt(grad_i^2 + grad_j^2)
end

function compute_local_variance(data::AbstractMatrix{T}, i::Int, j::Int, window::Tuple{Int,Int}) where {T<:Real}
    h_half = window[1] ÷ 2
    w_half = window[2] ÷ 2

    neighborhood = extract_neighborhood(data, i, j, h_half, w_half, :reflect)
    valid_vals = neighborhood[.!isnan.(neighborhood)]

    return length(valid_vals) > 1 ? Statistics.var(valid_vals) : 0.0
end

function adaptive_outlier_removal!(tg::Vector{Vector{IMASnodeRepr{T}}}) where {T<:Real}
    for time_group in tg
        leaf1 = first(time_group)
        for field in keys(leaf1.ids)
            if field != :time && hasdata(leaf1.ids, field) && typeof(getproperty(leaf1.ids, field)) <: Vector
                data = hcat((getproperty(leaf.ids, field) for leaf in time_group if hasdata(leaf.ids, field))...)
                clean_data = adaptive_outlier_removal!(data; adaptivity=:gradient)
                for (k, leaf) in enumerate(time_group)
                    if hasdata(leaf.ids, field)
                        setproperty!(leaf.ids, field, clean_data[:, k])
                    end
                end
            end
        end
    end
    return tg
end
