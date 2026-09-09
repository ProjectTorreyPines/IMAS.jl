using IMAS
using Test

"""
Stand-in for a satellite package (e.g. IFEdd) that defines its own concrete
container on top of IMAS's abstract `DD`. Mirrors the shape `GenerateDD` emits,
trimmed to a single reused IMAS IDS.

Lets IMAS verify that its `DD`-generic code paths work for containers other than
`IMAS.dd`, without depending on any satellite package.
"""
module FakeSatellite

import IMAS

const TSD = IMAS.IMASdd.ThreadSafeDicts

mutable struct FilledFields____my_own_dd <: IMAS.FilledFields
    var"requirements"::Bool
end

mutable struct my_own_dd{T} <: IMAS.DD{T}
    var"requirements"::IMAS.requirements{T}
    global_time::Float64
    _aux::TSD.ThreadSafeDict{Symbol,Any}
    _name::Symbol
    _filled::FilledFields____my_own_dd
    _frozen::Bool
    _threads_lock::ReentrantLock
    _in_expression::TSD.ThreadSafeDict{Int,Vector{Symbol}}
    _parent::WeakRef
end

function my_own_dd{T}(; frozen::Bool=false) where {T}
    ids = my_own_dd{T}(
        IMAS.requirements{T}(; frozen),
        0.0,
        TSD.ThreadSafeDict{Symbol,Any}(),
        Symbol(""),
        FilledFields____my_own_dd(false),
        frozen,
        ReentrantLock(),
        TSD.ThreadSafeDict{Int,Vector{Symbol}}(),
        WeakRef(nothing))
    setfield!(ids.requirements, :_parent, WeakRef(ids))
    return ids
end

my_own_dd(; frozen::Bool=false) = my_own_dd{Float64}(; frozen)

# satellites register their own paths into the shared registry; keyed by our own
# type, so this cannot collide with IMAS's entries
merge!(IMAS._all_info, Dict(
    (my_own_dd, :requirements) => IMAS.Info(String[], "-", "STRUCTURE", "Reused IMAS requirements IDS", true, String[])
))

end # module FakeSatellite

# Exercises both `IMAS.dd` and a foreign `<: IMAS.DD` container through the same
# assertions. `nested` matters because a bare `deepcopy(dd)` is served by IMASdd's
# `Base.deepcopy(::DD)`, while a dd reached through a container goes through
# `IMAS.deepcopy_internal` in src/control/fxp.jl — a different code path.
function check_deepcopy(make)
    original = make()
    original.global_time = 3.0
    original.requirements.cost = 7.0

    copied = deepcopy(original)
    @test typeof(copied) === typeof(original)
    @test copied !== original
    @test copied.global_time == 3.0
    @test copied.requirements.cost == 7.0

    # a real copy, not an alias
    copied.requirements.cost = 9.0
    @test original.requirements.cost == 7.0

    # the method the container paths funnel into
    direct = Base.deepcopy_internal(original, IdDict())
    @test typeof(direct) === typeof(original)
    @test direct.requirements.cost == 7.0

    for (label, wrapped, unwrap) in (
        ("Vector", [original], r -> r[1]),
        ("NamedTuple", (dd=original,), r -> r.dd),
        ("Dict", Dict(:k => original), r -> r[:k]))
        inner = unwrap(deepcopy(wrapped))
        @test typeof(inner) === typeof(original)
        @test inner.requirements.cost == 7.0
    end

    # `_frozen` travels with the copy
    @test getfield(deepcopy([make(; frozen=true)])[1], :_frozen)

    return nothing
end

@testset "satellite DD" begin

    @testset "type contract" begin
        @test isabstracttype(IMAS.DD)
        @test IMAS.dd{Float64} <: IMAS.DD{Float64}
        @test FakeSatellite.my_own_dd{Float64} <: IMAS.DD{Float64}

        # a satellite container is a sibling of IMAS.dd, not a subtype of it
        @test !(FakeSatellite.my_own_dd{Float64} <: IMAS.dd)
        @test !(IMAS.dd{Float64} <: FakeSatellite.my_own_dd)
    end

    @testset "deepcopy — IMAS.dd (baseline)" begin
        check_deepcopy(IMAS.dd)
    end

    @testset "deepcopy — user-defined DD subtype" begin
        check_deepcopy(FakeSatellite.my_own_dd)
    end

    # `deepcopy_internal` exists so that a deep-copied dd keeps talking to the same
    # FXP client rather than cloning it
    @testset "deepcopy_internal shares _aux[:fxp]" begin
        for make in (IMAS.dd, FakeSatellite.my_own_dd)
            original = make()
            client = Ref(0)   # stands in for a Jedis client
            getfield(original, :_aux)[:fxp] = client
            copied = Base.deepcopy_internal(original, IdDict())
            @test getfield(copied, :_aux)[:fxp] === client
        end
    end

end
