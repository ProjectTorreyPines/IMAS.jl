"""
    controller(ct::IMAS.controllers, name::String)

Pick a controller based on its name
"""
function controller(ct::IMAS.controllers, name::String)
    index = findfirst(c -> c.name == name, ct.linear_controller)
    if index === nothing
        return nothing
    end
    return ct.linear_controller[index]
end

"""
    pid_controller(controller::IMAS.controllers__linear_controller, kP::T, kI::T, kD::T) where {T<:Real}

Initializes SISO PID controller
"""
function pid_controller(controller::IMAS.controllers__linear_controller, kP::T, kI::T, kD::T) where {T<:Real}
    controller.input_names = ["error"]
    controller.output_names = ["control"]
    controller.pid.p.time = [-Inf]
    controller.pid.p.data = [kP;;;]
    controller.pid.i.time = [-Inf]
    controller.pid.i.data = [kI;;;]
    controller.pid.d.time = [-Inf]
    controller.pid.d.data = [kD;;;]
    controller.inputs.time = T[]
    controller.inputs.data = zeros(T, (1, 0))
    controller.outputs.time = T[]
    controller.outputs.data = zeros(T, (1, 0))
    return controller
end

"""
    (controller::controllers__linear_controller{T})(setpoint::T, value::T, time0::Float64) where {T<:Real}

Operates the linear controller
"""
function (controller::controllers__linear_controller{T})(setpoint::T, value::T, time0::Float64) where {T<:Real}

    P = controller.pid.p.data[1]
    I = controller.pid.i.data[1]
    D = controller.pid.d.data[1]

    dd = top_dd(controller)

    if stream_has_controller(dd, controller.name)
        control = stream_controller(dd, controller.name, 1.0; time=time0, setpoint, value, P, I, D)

    else
        error = setpoint - value
        push!(controller.inputs.time, time0)
        controller.inputs.data = hcat(controller.inputs.data, [error])

        if size(controller.inputs.data)[2] < 2
            integral = zero(T)
            derivative = zero(T)
        else
            integral = sum(
                (controller.inputs.data[1, k+1] + controller.inputs.data[1, k]) / 2.0 * (controller.inputs.time[k+1] - controller.inputs.time[k]) for
                k in 1:length(controller.inputs.time)-1
            )
            derivative = (controller.inputs.data[1, end] - controller.inputs.data[1, end-1]) / (controller.inputs.time[end] - controller.inputs.time[end-1])
        end
        control = (P * error) + (I * integral) + (D * derivative)
    end

    push!(controller.outputs.time, time0)
    controller.outputs.data = hcat(controller.outputs.data, [control])

    return control
end

function stream_controller(dd::IMAS.dd, controller_name::String, timeout::Float64; kw...)
    stream_fuse2ctrl = "$(controller_name)__fuse2ctrl"
    stream_ctrl2fuse = "$(controller_name)__ctrl2fuse"
    IMAS.stream_push!(dd, stream_fuse2ctrl; kw...)
    return IMAS.stream_pop!(dd, stream_ctrl2fuse; timeout)[:control]
end

function stream_has_controller(::Nothing, controller_name::String)
    return false
end

function stream_has_controller(dd::IMAS.dd, controller_name::String)
    return is_streaming(dd) && stream_has_service_provider(dd, controller_name)
end
