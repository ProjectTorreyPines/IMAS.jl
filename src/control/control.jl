document[Symbol("Control")] = Symbol[]

import Dates
import ExpiringCaches
import FuseExchangeProtocol as FXP

"""
    controller(ct::IMAS.controllers, name::String)

Pick a controller based on its name, returns `nothing` if not found

NOTE: for now only looks under `dd.controllers.linear_controller`
"""
function controller(ct::IMAS.controllers, name::String)
    index = findfirst(c -> c.name == name, ct.linear_controller)
    if index === nothing
        return nothing
    end
    return ct.linear_controller[index]
end

@compat public controller
push!(document[Symbol("Control")], :controller)

"""
    pid_controller(controller::IMAS.controllers__linear_controller, P::T, I::T, D::T) where {T<:Real}

Initializes PID controllers
"""
function pid_controller(controller::IMAS.controllers__linear_controller, P::T, I::T, D::T) where {T<:Real}
    controller.pid.p.time = [-Inf]
    controller.pid.p.data = [P;;;]
    controller.pid.i.time = [-Inf]
    controller.pid.i.data = [I;;;]
    controller.pid.d.time = [-Inf]
    controller.pid.d.data = [D;;;]
    controller.inputs.time = T[]
    controller.inputs.data = zeros(T, (1, 0))
    controller.outputs.time = T[]
    controller.outputs.data = zeros(T, (1, 0))
    return controller
end

"""
    (controller::controllers__linear_controller{T})(control_name::Any, control_algorithm::Val{:PID}, setpoint::T, value::T, time0::Float64) where {T<:Real}

Operates PID controllers
"""
function (controller::controllers__linear_controller{T})(control_name::Any, control_algorithm::Val{:PID}, setpoint::T, value::T, time0::Float64) where {T<:Real}
    P = controller.pid.p.data[1]
    I = controller.pid.i.data[1]
    D = controller.pid.d.data[1]

    error = setpoint - value
    push!(controller.inputs.time, time0)
    controller.inputs.data = hcat(controller.inputs.data, [error])

    if fxp_has_controller(controller)
        control = fxp_controller(controller, 1.0; time=time0, setpoint, value, P, I, D)
    else
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

@compat public controllers__linear_controller
push!(document[Symbol("Control")], :controllers__linear_controller)

function fxp_request_service(controller::controllers__linear_controller)::Bool
    if controller.name in keys(FXP_CONTROLLERS)
        fxp_serviced = true
    else
        fxp_serviced = fxp_request_service(top_dd(controller), "$(controller.name)_control")
    end
    if fxp_serviced
        controller.description = "fxp_serviced"
    else
        controller.description = "local_controller"
    end
    return fxp_serviced
end

function fxp_controller(controller::controllers__linear_controller, timeout::Float64; kw...)
    dd = top_dd(controller)
    client = fxp_client(dd)
    if controller.name in keys(FXP_CONTROLLERS)
        service_name = FXP_CONTROLLERS[controller.name]["service_name"]
        session_id = FXP_CONTROLLERS[controller.name]["session_id"]
    else
        service_name = "$(controller.name)_control"
        session_id = fxp_identifier(dd)
    end
    FXP.json_push(client, session_id, service_name, :requestor; kw...)
    return FXP.json_pop(client, session_id, service_name, :requestor; timeout)[:control]
end

ExpiringCaches.@cacheable Dates.Second(1) function fxp_has_controller(controller::controllers__linear_controller)::Bool
    if ismissing(controller, :description) || controller.description != "fxp_serviced"
        return false
    end
    if controller.name in keys(FXP_CONTROLLERS)
        return true
    end
    dd = top_dd(controller)
    if dd === nothing
        return false
    end
    return fxp_has_service_provider(dd, "$(controller.name)_control")
end
