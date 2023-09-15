"""
    controllers__linear_controller(kP::T, kI::T, kD::T) where {T<:Real}

Constructor that initializes SISO PID controller
"""
function controllers__linear_controller(kP::T, kI::T, kD::T) where {T<:Real}
    controller = controllers__linear_controller{T}()
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
    (controller::controllers__linear_controller{T})(set_point::T, value::T, time::Float64) where {T<:Real}

Operates the linear controller
"""
function (controller::controllers__linear_controller{T})(set_point::T, value::T, time::Float64) where {T<:Real}
    error = set_point - value
    push!(controller.inputs.time, time)
    controller.inputs.data = hcat(controller.inputs.data, [error])

    if size(controller.inputs.data)[2] < 2
        integral = 0.0
        derivative = 0.0
    else
        integral = sum((controller.inputs.data[1, k+1] + controller.inputs.data[1, k]) / 2.0 * (controller.inputs.time[k+1] - controller.inputs.time[k]) for k in 1:length(controller.inputs.time)-1)
        derivative = (controller.inputs.data[1, end] - controller.inputs.data[1, end-1]) / (controller.inputs.time[end] - controller.inputs.time[end-1])
    end
    control = (controller.pid.p.data[1] * error) + (controller.pid.i.data[1] * integral) + (controller.pid.d.data[1] * derivative)

    push!(controller.outputs.time, time)
    controller.outputs.data = hcat(controller.outputs.data, [control])

    return control
end