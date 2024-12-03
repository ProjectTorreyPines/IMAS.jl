document[Symbol("Physics tf")] = Symbol[]

"""
    tf_ripple(r, R_tf::Real, N_tf::Integer)

Evaluate fraction of toroidal magnetic field ripple at `r` [m]
generated from `N_tf` toroidal field coils with outer leg at `R_tf` [m]
"""
function tf_ripple(r, R_tf::Real, N_tf::Integer)
    eta = (r ./ R_tf) .^ N_tf
    return eta ./ (1.0 .- eta)
end

@compat public tf_ripple
push!(document[Symbol("Physics tf")], :tf_ripple)

"""
    R_tf_ripple(r, ripple::Real, N_tf::Integer)

Evaluate location of toroidal field coils outer leg `R_tf`` [m] at which `N_tf`toroidal field coils generate a given fraction of toroidal magnetic field ripple at`r` [m]
"""
function R_tf_ripple(r, ripple::Real, N_tf::Integer)
    return r .* (ripple ./ (ripple .+ 1.0)) .^ (-1 / N_tf)
end

@compat public R_tf_ripple
push!(document[Symbol("Physics tf")], :R_tf_ripple)

"""
    top_view_outline(tf::IMAS.build__tf, n::Int=1; cutouts::Bool=false)

returns x,y outline (top view) of the n'th TF coil
"""
function top_view_outline(tf::IMAS.build__tf, n::Int=1; cutouts::Bool=false)
    layers = parent(tf).layer
    TF = get_build_layers(layers; type=IMAS._tf_)

    N = mod(abs(n) - 1, tf.coils_n) + 1

    ϕROT = collect(range(0.0, 2pi, tf.coils_n + 1))[1:end-1]
    ϕwedge = 2pi / tf.coils_n / 2

    if cutouts
        ϕstart = atan(TF[1].end_radius * sin(ϕwedge), TF[2].start_radius)
        ϕend2 = atan(TF[1].end_radius * sin(ϕwedge), TF[2].end_radius * 2.0)

        @assert n != 0
        if n > 0
            x = [
                0.1 * cos(ϕwedge),
                TF[1].end_radius * cos(ϕwedge),
                TF[2].end_radius * 2.0
            ]

            y = [
                0.1 * sin(ϕwedge),
                TF[1].end_radius * sin(ϕwedge),
                TF[2].end_radius * 2.0 * sin(ϕend2)
            ]
        else
            x1 = [
                TF[1].end_radius * cos(ϕwedge),
                TF[2].end_radius * 2.0
            ]

            y1 = [
                TF[1].end_radius * sin(ϕwedge),
                TF[2].end_radius * 2.0 * sin(ϕend2)
            ]

            ϕrot = ϕROT[2]
            x2 = x1 .* cos(ϕrot) .+ y1 .* sin(ϕrot)
            y2 = x1 .* sin(ϕrot) .- y1 .* cos(ϕrot)

            x = [x1; reverse(x2)]
            y = [y1; reverse(y2)]
            x[end] = x[1]
            y[end] = y[1]

            ϕrot = ϕROT[N]
            xrot = x .* cos(ϕrot) .- y .* sin(ϕrot)
            yrot = x .* sin(ϕrot) .+ y .* cos(ϕrot)
            return (x=xrot, y=yrot)
        end

    else
        ϕstart = atan(TF[1].end_radius * sin(ϕwedge), TF[2].start_radius)
        ϕend = atan(TF[1].end_radius * sin(ϕwedge), TF[2].end_radius)

        x = [
            TF[1].start_radius * cos(ϕwedge),
            TF[1].end_radius * cos(ϕwedge),
            TF[2].start_radius,
            TF[2].end_radius
        ]

        y = [
            TF[1].start_radius * sin(ϕwedge),
            TF[1].end_radius * sin(ϕwedge),
            TF[2].start_radius * sin(ϕstart),
            TF[2].end_radius * sin(ϕend)
        ]
    end

    x = [x; reverse(x); x[1]]
    y = [y; reverse(-y); y[1]]

    ϕrot = ϕROT[N]
    xrot = x .* cos(ϕrot) .- y .* sin(ϕrot)
    yrot = x .* sin(ϕrot) .+ y .* cos(ϕrot)
    return (x=xrot, y=yrot)
end

@compat public top_view_outline
push!(document[Symbol("Physics tf")], :top_view_outline)
