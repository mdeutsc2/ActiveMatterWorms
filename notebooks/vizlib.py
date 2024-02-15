#vizlib.py
import numpy as np
import tqdm
import matplotlib.pyplot as plt
from scipy import interpolate

def load_data(filename):
    Ns = []
    with open(filename) as f:
        N = int(f.readline().strip())
        print(N," particles per frame")
        f.seek(0)
        line_count = 0
        for line in f:
            line_count += 1
            if (str(N) in line.strip()) and ("A" not in line.strip()) and ("B" not in line.strip()):
                Ns.append([line.strip(),line_count+1,line_count+N+1]) # number of particles, data start, data end
    print("Loading ",len(Ns)," frames of data...")
    pos_data = np.zeros((len(Ns),N,3))
    vel_data = np.zeros((len(Ns),N,3))
    with open(filename,'r') as f:
        #frame_data = np.zeros((N,6))
        raw_data = f.readlines()
        for iline in range(len(raw_data)):
            line = raw_data[iline]
            raw_data[iline] = line.strip()
        for frame in tqdm.tqdm(range(len(Ns))):
            data_line_start = Ns[frame][1]
            data_line_end = Ns[frame][2]
            raw_frame_data_list = raw_data[data_line_start:data_line_end] # strips "A/B"
            for (i,line) in enumerate(raw_frame_data_list):
                line_data = [float(val) for val in line[1:].split()]
                pos_data[frame,i,:] = line_data[0:3]
                vel_data[frame,i,:] = line_data[3:6]
    return N,len(Ns),pos_data,vel_data

def get_particle_vecs_frame(N,worm_length,iframe,pos_data,disp=True):
    """ calculating the vectors of each of the individual particles"""
    vector_field = np.zeros((N,6))
    num_worms = int(N/worm_length)
    if disp:
        print("num_worms:",num_worms)

    worm_view = np.array_split(pos_data[iframe,:,:],num_worms)
    count = 0
    for iw in range(num_worms):
        for ip in range(len(worm_view[iw])):
            ip1 = ip + 1
            if ip == len(worm_view[iw])-1:
                ip1 = ip-1
            vec = np.zeros(3)
            pos = np.zeros(3)
            pos[0] = worm_view[iw][ip][0]
            pos[1] = worm_view[iw][ip][1]
            pos[2] = worm_view[iw][ip][2]
            vec[0] = worm_view[iw][ip][0] - worm_view[iw][ip1][0] #X
            vec[1] = worm_view[iw][ip][1] - worm_view[iw][ip1][1] #Y
            vec[2] = worm_view[iw][ip][2] - worm_view[iw][ip1][2] #Z
            n = np.sqrt(vec[0] ** 2 + vec[1] ** 2 + vec[2]**2)
            vec = vec/n
            vector_field[count,0:3] = pos
            vector_field[count,3:6] = vec
            count += 1
            #print(ip,worm_view[iw][ip],vec)
    return vector_field

def interp_particle_vecs(vector_field,num_interpolation_points,interpolation_type,disp=True):
    #extend of the particles in +x,-x,+y,-y
    max_x = np.max(vector_field[:,0])
    min_x = np.min(vector_field[:,0])
    max_y = np.max(vector_field[:,1])
    min_y = np.min(vector_field[:,1])
    if disp:
        print("X:",min_x,max_x,max_x-min_x)
        print("Y:",min_y,max_y,max_y-min_y)

    # create structured grid
    #num_interpolation_points = 100
    grid_x,grid_y = np.meshgrid(np.linspace(min_x,max_x,num=num_interpolation_points),
                                np.linspace(min_y,max_y,num=num_interpolation_points))
    # Create points for interpolation
    points = np.column_stack((vector_field[:,0], vector_field[:,1]))
    # Interpolate the x and y components separately using bilinear interpolation
    grid_u = interpolate.griddata(points, vector_field[:,3], (grid_x, grid_y), method=interpolation_type)
    grid_v = interpolate.griddata(points, vector_field[:,4], (grid_x, grid_y), method=interpolation_type)
    # merged everything back together to make single object
    structured_vector_field = np.column_stack((grid_x.flatten(), grid_y.flatten(), grid_u.flatten(), grid_v.flatten()))
    # normalize the vector parts
    n = np.sqrt(structured_vector_field[:,2] ** 2 + structured_vector_field[:,3] ** 2)
    structured_vector_field[:,2] = structured_vector_field[:,2]/n
    structured_vector_field[:,3] = structured_vector_field[:,3]/n

    # replace nans with 0.0
    np.nan_to_num(structured_vector_field[:,2], copy=False, nan=0.0)
    np.nan_to_num(structured_vector_field[:,3], copy=False, nan=0.0)
    return structured_vector_field

def get_energy_simple(vector_grid,num_interpolation_points,disp=True):
    if disp:
        print(vector_grid.shape)
    vecs = vector_grid.reshape(num_interpolation_points,num_interpolation_points,4)
    max_x = np.max(vector_grid[-1,0])
    min_x = np.min(vector_grid[0,0])
    max_y = np.max(vector_grid[-1,1])
    min_y = np.min(vector_grid[1,1])
    u = vecs[:,:,2]
    v = vecs[:,:,3]
    grid_x,grid_y = np.meshgrid(np.linspace(min_x,max_x,num=num_interpolation_points),
                                np.linspace(min_y,max_y,num=num_interpolation_points))
    E = np.zeros((num_interpolation_points,num_interpolation_points))
    for i in np.arange(num_interpolation_points):
        for j in np.arange(num_interpolation_points):
            jnab = j
            inab = i + 1
            if inab >= num_interpolation_points-1:
                inab = i
            dp = u[i,j]*u[inab,jnab] + v[i,j]*v[inab,jnab]
            E[i,j] += 1-(dp*dp)

            jnab = j
            inab = i - 1
            if inab <= 0:
                inab = i
            dp = u[i,j]*u[inab,jnab] + v[i,j]*v[inab,jnab]
            E[i,j] += 1-(dp*dp)

            inab = i
            jnab = j + 1
            if jnab >= num_interpolation_points-1:
                jnab = j
            dp = u[i,j]*u[inab,jnab] + v[i,j]*v[inab,jnab]
            E[i,j] += 1-(dp*dp)

            inab = i
            jnab = j - 1
            if jnab <= 0:
                jnab = j
            dp = u[i,j]*u[inab,jnab] + v[i,j]*v[inab,jnab]
            E[i,j] += 1-(dp*dp)
    E_norm = (E-np.min(E))/(np.max(E)-np.min(E))
    return E_norm