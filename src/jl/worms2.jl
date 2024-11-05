using StaticArrays
using ProgressMeter
using BenchmarkTools
using KernelAbstractions
using ArgParse
using SpatialHashTables
using SpatialHashTables: dist²,numthreads,domainsize

const DT::Float64 = 1e-3
const DT2O2::Float64 = DT*DT/2
const KSPRING::Float64 = 57.146436
const K2SPRING::Float64 = 50.0*KSPRING #100
const K3SPRING::Float64 = 50.0*KSPRING
const KBEND::Float64 = 40.0
const LENGTH0::Float64 = 0.8 # equilibrium length
const LENGTH2::Float64 = 2.0*LENGTH0
const LENGTH3::Float64 = 3.0*LENGTH0
const RCUT::Float64 = 2.5
const RCUTSMALL::Float64 = 2.0^(1.0/6.0)
const SIGMA::Float64 = 2.0
const A::Float64 = 0.14
const GAMMA::Float64 = 6.0 # frictional constant for dissipative force (~1/damp)
const NUM_POINTS::Int = 5000 # number of boundary points
const FLUID_OFFSET::Float64 = RCUTSMALL*SIGMA # z-offset of fluid
#const THETANOW::Float64 = 5.0*pi
Vec3D = SVector{3,Float64}


@inline function is_filament(i,params)
    if i <= params["np"]*params["nworms"]
        return true
    end
    return false
end

@inline dist_sq(xy) = sum(z -> z^2, xy)

@inline function force_fnc(Xi, Xj, Xij, d²)
    return Xij / sqrt(d²)
end 

@inline function force_fnc2(Xi, Xj, Xij, d²)
    return Xij / sqrt(d²)
end

@kernel function forces_multithreaded_kernel!(F, @Const(X), @Const(V), r, grid, nworms,np, nbounds,ww_epsilon,fdep)
    i = @index(Global)
    N = nworms*np
    Xi = X[i]
    Fi = zero(eltype(F))
    r² = r*r #cutoff
    for j in neighbours(grid, Xi, r)
        Xj = X[j]

        Xij = Xi - Xj
        r2 = dist_sq(Xij)
        if (i <= N) && (j <= N) # both are filaments
            iworm = 1 + (i - 1)/np
            ip = i - np*(iworm - 1)
            jworm = 1 + (j - 1)/np
            jp = j - np*(jworm - 1)
            #if not on same filament
            if (iworm != jworm) && (abs(ip-jp) > 2)
                if zero(r) < r2 < RCUTSMALL^2
                    Fi += (-48.0*ww_epsilon*r2^(-7.0) + 24.0*ww_epsilon*r2^(-4.0) + fdep/sqrt(r2))*Xij
                end
            end

            # if zero(r) < r2 < r^2
            #     # inter-worm forces
            #     Fi += force_fnc(Xi, Xj, Xij, d²)

            # elseif  (i < N) && (N < j <= N+nbounds) # i = filament, j = boundary
            #     Fi += force_fnc2(Xi, Xj, Xij, d²)
            # elseif  (j < N) && (N < i <= N+nbounds) # i = boundary, j = filament
            #     Fi += force_fnc2(Xi, Xj, Xij, d²)
            # end
        end 
    end
    F[i] = Fi    
end


function forces_multithreaded!(F, X, V, params, r, grid)
    kernel = forces_multithreaded_kernel!(grid.backend, numthreads(grid))
    kernel(F, X, V, r, grid,params["nworms"], params["np"],NUM_POINTS, params["ww_epsilon"],params["fdep"],ndrange = length(X))
    synchronize(grid.backend)
    return nothing
end

function init_worms(params,N::Int)
    #positions = rand(Vec2D,N)
    positions = zeros(Vec3D,N)
    ireverse = zeros(params["nworms"])
    save = zeros(Vec3D,params["np"])
    bound = zeros(Vec3D,NUM_POINTS)
    hx = params["rwall"]*2.0 + 1.0
    hxo2 = hx/2
    hyo2 = hx/2
    if (params["boundary"] == 1)
        eqidistantThetaValues = zeros(NUM_POINTS)
        deltaTheta = 2*pi/NUM_POINTS
        for i in 1:NUM_POINTS
            eqidistantThetaValues[i] = i * deltaTheta
        end
        for i in 1:NUM_POINTS
            r = params["rwall"]
            bound[i] = Vec3D(r*cos(eqidistantThetaValues[i])+hxo2,
                            r*sin(eqidistantThetaValues[i])+hxo2,
                            0.0)
            positions[i+params["np"]*params["nworms"]] = bound[i]
        end
        thetanow = 5.0*pi
        for iw in 1:params["nworms"]
            #wormid = (iw-1)*np+ip;
            #iw = 1 + ((i - 1)/np):int; // find which worm j is in
            #ip = i - np*(iw - 1); // which particle in the worm is j?
            if (rand() <= 0.5)
                ireverse[iw] = 1.0
            end
            dth = 0.0
            for ip in 1:params["np"]
                r = A*thetanow
                dth = LENGTH0/r
                thetanow += dth
                id = (iw-1)*params["np"]+ip
                #println(id," ",iw," ",ip)
                positions[id] = Vec3D(hxo2 + r*cos(thetanow),
                                      hyo2 + r*sin(thetanow),
                                      rand()*params["L"])
            end
            thetanow += 2.0*dth
        end
        for iw in 1:params["nworms"]
            if ireverse[iw] == 1
                for ip in 1:params["np"]
                    id = (iw-1)*params["np"]+ip
                    save[ip] = positions[id]
                end
                for ip in 1:params["np"]
                    id = (iw-1)*params["np"]+ip
                    positions[id] = save[params["np"]+1-ip]
                end
            end
        end
    end
    return positions
end

function write_xyz(params,dtype,filename,istep,d,positions)
    # dtype, filename, istep, d == dimensions, positions (3d array of positions)
    f = open(filename, "w+")
    #write(f, string(d[1]," ",d[2]," ",d[2],"\n"))
    #write(f,string("# ",length(positions)),"\n")
    write(f,string(length(positions)),"\n")
    write(f,string("# ",istep,"\n"))
    hx = params["rwall"]*2.0 + 1.0
    hxo2 = hx/2
    hyo2 = hx/2
    ic = 0
    for iw in 1:params["nworms"]
        id1 = (iw-1)*params["np"]+1
        idnp = (iw-1)*params["np"]+params["np"]
        dx = positions[id1][1] - hxo2
        dy = positions[id1][2] - hyo2
        xang = atan(dy,dx)
        rx = -sin(xang)
        ry = cos(xang)
        dot = (positions[id1][1] - positions[idnp][1])*rx + (positions[id1][2] - positions[idnp][2])*ry;
        if (dot >= 0.0)
            for ip in 1:params["np"]
                i = (iw-1)*params["np"]+ip
                write(f,string("A ",positions[i][1]," ",positions[i][2]," ",positions[i][3],"\n"))
                ic+=1
            end
        else
            for ip in 1:params["np"]
                i = (iw-1)*params["np"]+ip
                write(f,string("B ",positions[i][1]," ",positions[i][2]," ",positions[i][3],"\n"))
                ic+=1
            end
        end
    end
    for i in params["np"]*params["nworms"]+1:params["np"]*params["nworms"]+NUM_POINTS
        write(f,string("I ",positions[i][1]," ",positions[i][2]," ",positions[i][3],"\n"))
        ic+=1
    end
    close(f)
    println(ic,"/",length(positions)," written")
    return nothing
end

function update_positions!(positions,velocities,F,Fold,params)
    Threads.@threads for i in eachindex(positions)
        if is_filament(i,params)
            positions[i] += velocities[i]*DT + F[i]*DT*DT*0.5
            Fold[i] = F[i]
            F[i] = [0.0,0.0,0.0]
        end
    end
end

function intraworm_forces!(positions,F,params)
    Threads.@threads for iw in 1:params["nworms"]
        for ip in 1:params["np"]-1
            id = (iw-1)*params["np"]+ip
            idp1 = (iw-1)*params["np"]+ip+1
            dx = positions[idp1][1] - positions[id][1]
            dy = positions[idp1][2] - positions[id][2]
            dz = positions[idp1][3] - positions[id][3]
            r = sqrt(dx*dx + dy*dy + dz*dz)
            ff = -KSPRING*(r - LENGTH0)/r
            F[idp1] += Vec3D(ff*dx,ff*dy,ff*dz)
            F[id] -= Vec3D(ff*dx,ff*dy,ff*dz)
        end
    end
    Threads.@threads for iw in 1:params["nworms"]
        for ip in 1:params["np"]-2
            id = (iw-1)*params["np"]+ip
            idp2 = (iw-1)*params["np"]+ip+2
            dx = positions[idp2][1] - positions[id][1]
            dy = positions[idp2][2] - positions[id][2]
            dz = positions[idp2][3] - positions[id][3]
            r = sqrt(dx*dx + dy*dy + dz*dz)
            ff = -K2SPRING*(r - LENGTH2)/r
            F[idp2] += Vec3D(ff*dx,ff*dy,ff*dz)
            F[id] -= Vec3D(ff*dx,ff*dy,ff*dz)
        end
    end
    Threads.@threads for iw in 1:params["nworms"]
        for ip in 1:params["np"]-3
            id = (iw-1)*params["np"]+ip
            idp3 = (iw-1)*params["np"]+ip+3
            dx = positions[idp3][1] - positions[id][1]
            dy = positions[idp3][2] - positions[id][2]
            dz = positions[idp3][3] - positions[id][3]
            r = sqrt(dx*dx + dy*dy + dz*dz)
            ff = -K3SPRING*(r - LENGTH3)/r
            F[idp3] += Vec3D(ff*dx,ff*dy,ff*dz)
            F[id] -= Vec3D(ff*dx,ff*dy,ff*dz)
        end
    end
end

function update_velocities!(positions,velocities,F,Fold,params)
    Threads.@threads for i in eachindex(positions)
        if is_filament(i,params)
            velocities[i] += DT*0.5*(F[i] + Fold[i])
        end
    end
end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--np"
            help = "number of particles/worm"
            arg_type = Int
            default = 80
        "--nworms"
            help = "number of worms"
            arg_type = Int
            default = 800
        "--nsteps"
            help = "number of steps"
            arg_type = Int
            default = 12000000
        "--fdogic"
            help = "active force"
            arg_type = Float64
            default = 0.06
        "--walldrive"
            help = "walldrive"
            arg_type = Bool
            default = false
        "--fdogicwall"
            help = "fdogicwall"
            arg_type = Float64
            default = 0.001
        "--fdep"
            help = "fdep"
            arg_type = Float64
            default = 0.25
        "--dogic_fdep"
            help = "dogic_fdep:  extra attractive force when dogic shearing is present (originall 0.7)"
            arg_type = Float64
            default = 0.25
        "--fdepwall"
            help = "wall depletion force"
            arg_type = Float64
            default = 6.0
        "--rwall"
            help = "wall radius (for disk, adjusted for cycloid)"
            arg_type = Float64
            default = 164.0 #125.0*rcutsmall*sqrt(2.0)
        "--save_interval"
            help = "save interval"
            arg_type = Int
            default = 4000
        "--boundary"
            help = "boundary type"
            arg_type = Int
            default = 2
        "--kbt"
            help = "temperature"
            arg_type = Float64
            default = 1.5
        "--fluid_rho"
            help = "underlying fluid layer density"
            arg_type = Float64
            default = 0.1
        "--L"
            help = "thickness of filament cell"
            arg_type = Float64
            default = 3.2
        "--ww_epsilon"
            help = "worm-worm WCA interaction epsilon"
            arg_type = Float64
            default = 0.5
        "--sw_epsilon"
            help = "solvent worm LJ interaction epsilon"
            arg_type = Float64
            default = 2.0
        "arg1"
            help = "a positional argument"
            required = false
    end
    parsed_args = parse_args(s)
    println("Parameters:")
    for (arg,val) in parsed_args
        println("\t$arg  =>  $val")
    end
    return parsed_args
end

function main()
    @info Threads.nthreads()
    params = parse_commandline()

    dimensions = (1,1,1)::NTuple{3,Int}
    N = params["np"]*params["nworms"] + NUM_POINTS

    positions = init_worms(params,N)
    write_xyz(params,Float32,"test_"*string(000)*".xyz",0,dimensions,positions)
    F = zeros(SVector{3, Float64}, N) # forces
    Fold = zeros(SVector{3, Float64}, N) 
    velocities = randn(SVector{3, Float64}, N)

    gridsize = (10, 10, 10)::NTuple{3,Int}
    grid = BoundedGrid(RCUT,gridsize,positions)
    @info "Starting"
    for istep in 1:10_000 #main loop on nsteps
        e_time = @elapsed begin
        update_positions!(positions,velocities,F,Fold,params)
        intraworm_forces!(positions,F,params)
        forces_multithreaded!(F,positions,velocities,params,RCUT,grid)
        update_velocities!(positions, velocities, F, Fold,params)
        updatecells!(grid,positions)
        end
        println(istep," ",e_time," s")
        if (istep % 100 == 0) 
            println(istep,sum(F))
            write_xyz(params,Float32,"test_"*string(istep)*".xyz",0,dimensions,positions)
        end
    end
    @info "Done!"
end
main()