"""
    tf_ripple(r, R_tf::Real, N_tf::Integer)

Evaluate fraction of toroidal magnetic field ripple at `r` [m]
generated from `N_tf` toroidal field coils with outer leg at `R_tf` [m]
"""
function tf_ripple(r, R_tf::Real, N_tf::Integer)
    eta = (r ./ R_tf) .^ N_tf
    return eta ./ (1.0 .- eta)
end

"""
    R_tf_ripple(r, ripple::Real, N_tf::Integer)

Evaluate location of toroidal field coils outer leg `R_tf`` [m] at which `N_tf`toroidal field coils generate a given fraction of toroidal magnetic field ripple at`r` [m]
"""
function R_tf_ripple(r, ripple::Real, N_tf::Integer)
    return r .* (ripple ./ (ripple .+ 1.0)) .^ (-1 / N_tf)
end

"""
    top_outline(tf::IMAS.build__tf)

returns x,y outline (top view) of the tf coil
"""
function top_outline(tf::IMAS.build__tf, n::Int=1; cutouts::Bool=false)
    layers = parent(tf).layer
    TF = get_build_layers(layers; type=IMAS._tf_)

    ϕwedge = 2pi / tf.coils_n / 2

    if cutouts
        ϕstart = atan(TF[1].end_radius * sin(ϕwedge), TF[2].start_radius)
        ϕend2 = atan(TF[1].end_radius * sin(ϕwedge), TF[2].end_radius * 2.0)

        x = [
            TF[1].start_radius / 2.0 * cos(ϕwedge),
            TF[1].end_radius * cos(ϕwedge),
            TF[2].end_radius * 2.0
        ]
    
        y = [
            TF[1].start_radius / 2.0 * sin(ϕwedge),
            TF[1].end_radius * sin(ϕwedge),
            TF[2].end_radius * 2.0 * sin(ϕend2)
        ]

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

    if n != 1
        ϕrot = collect(range(0.0, 2pi, tf.coils_n + 1))[1:end-1][n]
        xrot = x .* cos(ϕrot) .- y .* sin(ϕrot)
        yrot = x .* sin(ϕrot) .+ y .* cos(ϕrot)
        return (x=xrot, y=yrot)
    else
        return (x=x, y=y)
    end
end