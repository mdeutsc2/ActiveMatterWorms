using StaticArrays
using ProgressMeter
using ArgParse
using CellListMap
import CellListMap.wrap_relative_to


mutable struct Filaments
    x::Vector{SVector{2,Float64}}
    v::Vector{SVector{2,Float64}}
    id::Int
end

function update_forces!(x, y, i, j, d2, forces, cutoff)
    if (i > 20)
        r = y - x
        dudr = 10^6 * 4 * r * (d2 - cutoff^2)
        forces[i] += dudr
        forces[j] -= dudr
    else
        forces[i] += zeros(SVector{2,Float64})
        forces[j] -= zeros(SVector{2,Float64})
    end
    return forces
end

function init_worms(N::Int)
    Vec2D = SVector{2,Float64}
    #positions = rand(Vec2D,N)
    positions = zeros(Vec2D,N)
    for i in 1:N
        phi = randn()*0.9
        theta = 2*pi*randn()*0.9
        positions[i] = Vec2D(cos(phi)*sin(theta),sin(phi)*sin(theta))
    end
    return positions
end

function init_system(;N::Int=200)
    positions = init_worms(N)
    unitcell = [1.0, 1.0]
    cutoff = 0.1
    system = ParticleSystem(
        positions=positions,
        cutoff=cutoff,
        unitcell=unitcell,
        output=similar(positions),
        output_name=:forces,
    )
    return system
end
function simulate(system=init_system(); nsteps::Int=100, isave=1)
    # initial velocities
    velocities = [ randn(eltype(system.positions)) for _ in 1:length(system.positions) ]
    for i in 1:20
        velocities[i] = zeros(SVector{2,Float64})
    end
    dt = 1e-3
    trajectory = typeof(system.positions)[]
    for step in 1:nsteps
	    println(step)
        # compute forces at this step
        map_pairwise!(
            (x,y,i,j,d2,forces) -> update_forces!(x,y,i,j,d2,forces,system.cutoff),
            system
        )
        # Update positions and velocities
        for i in eachindex(system.positions, system.forces)
            if (i > 20)
                f = system.forces[i]
            else
                f = zeros(SVector{2,Float64})
            end
            x = system.positions[i]
            v = velocities[i]
            x = x + v * dt + (f / 2) * dt^2
            v = v + f * dt
            # wrapping to origin for obtaining a pretty animation
            x = wrap_relative_to(x, SVector(0.0, 0.0), system.unitcell)
            # !!! IMPORTANT: Update arrays of positions and velocities
            system.positions[i] = x
            velocities[i] = v
        end
        # Save step for printing
        if step % isave == 0
            push!(trajectory, copy(system.positions))
        end
    end
    return trajectory
end

using Plots
function animate(trajectory)
    c = zeros(40)
    c[20:40] .= 1
    anim = @animate for step in trajectory
        scatter(
            Tuple.(step),
            label=nothing,
            lims=(-0.5, 0.5),
            aspect_ratio=1,
            framestyle=:box,
            zcolor=c,
        )
    end
    gif(anim, "simulation.gif", fps=10)
end

system = init_system(N=40);
@info Threads.nthreads()
system.parallel = true
trajectory = simulate(system,nsteps=50);
animate(trajectory)
