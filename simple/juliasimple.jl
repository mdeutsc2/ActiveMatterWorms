using CUDA;


mutable struct Bin
    id::Int32
    atoms::Array{Int32,1}
    neighbors::Array{Int32,1}
    ncount::Int32
end

function init_bins!(bins)
    for ibin in 1:numBins
        for jbint in 1:numBins
            binid = (jbin-1)*numBins+ibin
            ibnnab = ibin+1
            jbinnab = jbin
            if (ibinnab <= numBins)
                binidnbor = (jbinnab-1)*numBins)+ibinnab
                bins[binid].neighbors[1] = binidnbor
            end
            
end
function update_pos!(pos,vel,frc)
    Threads.@threads for i in 1:numParticles
        pos[i,1] += vel[i,1]*dt + (frc[i,1]/2)*(dt^2)
        pos[i,2] += vel[i,2]*dt + (frc[i,2]/2)*(dt^2)

        if pos[i,1] > L
            pos[i,1] = pos[i,1]-L
        elseif pos[i,1] < 0.0
            pos[i,1] = pos[i,1] + L
        end
        if pos[i,2] > L
            pos[i,2] = pos[i,2]-L
        elseif pos[i,2] < 0.0
            pos[i,2] = pos[i,2] + L
        end

        vel[i,1] += 0.5*dt*frc[i,1]
        vel[i,2] += 0.5*dt*frc[i,2]
        frc[i,1] = 0.0
        frc[i,2] = 0.0
    end
end

function update_vel!(vel,frc,KE)
    Threads.@threads for i in 1:numParticles
        vel[i,1] += 0.5*dt*frc[i,1]
        vel[i,2] += 0.5*dt*frc[i,2]
        KE[i] = 0.5*(vel[i,1]*vel[i,1] + vel[i,2]*vel[i,2])
    end
end

function calc_forces!(pos,frc)
    n = size(pos)[1]
    for i in 1:n-1
        for j in i:n
            dx = pos[j,1] - pos[i,1]
            dy = pos[j,2] - pos[i,2]
            r2 = (dx*dx + dy*dy)
            if (r2 <= r2cut)
                ffor = -48.0*r2^(-7.0) + 24.0*r2^(-4.0)
                ffx = ffor*dx
                ffy = ffor*dy
                frc[i,1] += ffx
                frc[i,2] += ffy
                frc[j,1] -= ffx
                frc[j,2] -= ffy
            end
        end
    end
end

numParticles = 625
L = 100.0
dt = 0.001
nsteps = 1000
save_interval = 5000
thermo = false
kbt = 0.5

sigma = 1.0
rcutsmall = 2.0^(1.0/6.0)
rcut = 2.5*sigma
r2cut = rcut*rcut
gamma = 3.0
print_interval = 100

pos = zeros((numParticles,2))
vel = zeros((numParticles,2))
frc = zeros((numParticles,2))
numBins = Int(ceil(L/rcut))
bins = Array{Bin,1}(undef,numBins)
for binid in 1:numBins
    tmp_atoms = zeros(Int32,64)
    tmp_nbors = zeros(Int32,4)
    tmp_atoms .= -1
    tmp_nbors .= -1
    bins[binid] = Bin(binid,tmp_atoms,tmp_nbors,0)
end

KE = zeros(numParticles)
if (numParticles != sqrt(numParticles)*sqrt(numParticles))
    println("non-square numParticles")
end

row_length = sqrt(numParticles)
spacing = 1.0
center = L/2 - spacing*(row_length/2)
for i in 1:numParticles
    row = i%row_length
    col = ((i-row)/row_length)-1
    pos[i,1] = center + spacing*rcutsmall*row
    pos[i,2] = center + spacing*rcutsmall*col
    vel[i,1] = 0.1*randn()
    vel[i,2] = 0.1*randn()
    KE[i] = 0.5*(vel[i,1]*vel[i,1] + vel[i,2]*vel[i,2])
end

vxave = sum(vel[:,1])/numParticles
vyave = sum(vel[:,1])/numParticles
vel[:,1] = vel[:,1] .- vxave
vel[:,2] = vel[:,2] .- vyave

calc_forces!(pos,frc)
time_elapsed = 0.0
for istep in 1:nsteps
    global time_elapsed += @elapsed begin
        update_pos!(pos,vel,frc)
        calc_forces!(pos,frc)
        update_vel!(vel,frc,KE)
    end
    if (istep%print_interval == 0)
        println("Step: ",istep,"\t",istep/time_elapsed," iters/s")
    end
end
