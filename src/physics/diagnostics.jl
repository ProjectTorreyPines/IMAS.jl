"""
    lenght_line_of_sight(icls::IMAS.interferometer__channel___line_of_sight)

Returns lenght of a given interferometer channel line of sight
"""
function lenght_line_of_sight(icls::IMAS.interferometer__channel___line_of_sight)
    x1 = icls.first_point.r * cos(icls.first_point.phi)
    y1 = icls.first_point.r * sin(icls.first_point.phi)
    z1 = icls.first_point.z

    x2 = icls.second_point.r * cos(icls.second_point.phi)
    y2 = icls.second_point.r * sin(icls.second_point.phi)
    z2 = icls.second_point.z

    d = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)

    if !isempty(icls.third_point)
        x3 = icls.third_point.r * cos(icls.third_point.phi)
        y3 = icls.third_point.r * sin(icls.third_point.phi)
        z3 = icls.third_point.z

        d += sqrt((x3 - x2)^2 + (y3 - y2)^2 + (z3 - z2)^2)
    end

    return d
end