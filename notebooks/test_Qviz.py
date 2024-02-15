import numpy as np
import tqdm
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from scipy.spatial.distance import cdist
import matplotlib.ticker as mticker
import os

filename = "/media/matt/SCRATCH/ActiveMatterWorms/data_02_02_2024/filament_data.npz"
data = np.load(filename)
pos_data = data['a']
vel_data = data['b']
N = data['a'].shape[1]
nframes = data['a'].shape[0]
print(nframes)

def get_vector_frame(iframe,worm_length,pos_data):
    #worm_length = 80
    #iframe = -1
    num_worms = int(N/worm_length)
    print("num_worms:",num_worms)
    vector_field = np.zeros((N,6))
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
    return vector_field

def get_Q(vector_field,num_interpolation_points):
    max_x = np.max(vector_field[:,0])
    min_x = np.min(vector_field[:,0])
    max_y = np.max(vector_field[:,1])
    min_y = np.min(vector_field[:,1])
    print("X:",min_x,max_x,max_x-min_x)
    print("Y:",min_y,max_y,max_y-min_y)

    x_grid,y_grid = np.meshgrid(np.linspace(min_x,max_x,num_interpolation_points),
                                np.linspace(min_y,max_y,num_interpolation_points))

    Qxx = np.zeros((num_interpolation_points,num_interpolation_points))
    Qxy = np.zeros((num_interpolation_points,num_interpolation_points))
    for i in np.arange(0,x_grid.shape[0]-1):
        for j in np.arange(0,y_grid.shape[1]-1):
            # getting the cell boundary points
            mask = (vector_field[:,0] >= x_grid[i,j]) & (vector_field[:,0] < x_grid[i,j+1]) & \
                 (vector_field[:,1] >= y_grid[i,j]) & (vector_field[:,1] < y_grid[i+1,j])
            # selecting particles that are in that cell
            if np.any(mask):
                x_cell_vals = vector_field[mask,3]
                y_cell_vals = vector_field[mask,4]
            else:
                x_cell_vals = 0.0
                y_cell_vals = 0.0

            Qxx[i,j] = 0.5*(3*np.mean(x_cell_vals*x_cell_vals)-1)
            Qxy[i,j] = 0.5*(3*np.mean(x_cell_vals*y_cell_vals))
    np.nan_to_num(Qxx, copy=False, nan=0.0)
    np.nan_to_num(Qxy, copy=False, nan=0.0)
    return Qxx,Qxy,x_grid,y_grid

def get_director(Qxx,Qxy,x_grid,y_grid,num_interpolation_points):
    nx = np.zeros((num_interpolation_points,num_interpolation_points))
    ny = np.zeros((num_interpolation_points,num_interpolation_points))
    for i in np.arange(0,x_grid.shape[0]-1):
        for j in np.arange(0,y_grid.shape[1]-1):
            #S[x,y] = 0.5*np.sqrt(self.trace_Q_squared(self.Q[x,y,0],self.Q[x,y,1]))
            a = Qxx[i,j]
            b = Qxy[i,j]
            eigenvalues,eigenvectors = np.linalg.eig(np.array([[a,b],[b,-a]]))
            nx[i,j] = eigenvectors[0,np.argmax(eigenvalues)]
            ny[i,j] = eigenvectors[1,np.argmax(eigenvalues)]
    return nx,ny

def calc_defects(Qxx,Qxy,x_grid,y_grid,normalize=True):
    dx = x_grid[0,1] - x_grid[0,0]
    dy = y_grid[1,0] - y_grid[0,0]

    dq11_tensor_dx = np.zeros_like(x_grid)
    dq12_tensor_dx = np.zeros_like(x_grid)
    dq11_tensor_dy = np.zeros_like(y_grid)
    dq12_tensor_dy = np.zeros_like(y_grid)

    dq11_tensor_dx[:, 1:-1] = (Qxx[:, 2:] - Qxx[:, :-2]) / (2 * dx)
    dq11_tensor_dx[:, 0] = (Qxx[:, 1] - Qxx[:, 0]) / dx
    dq11_tensor_dx[:, -1] = (Qxx[:, -1] - Qxx[:, -2]) / dx

    dq12_tensor_dx[:, 1:-1] = (Qxy[:, 2:] - Qxy[:, :-2]) / (2 * dx)
    dq12_tensor_dx[:, 0] = (Qxy[:, 1] - Qxy[:, 0]) / dx
    dq12_tensor_dx[:, -1] = (Qxy[:, -1] - Qxy[:, -2]) / dx

    dq11_tensor_dy[1:-1, :] = (Qxx[2:, :] - Qxx[:-2, :]) / (2 * dy)
    dq11_tensor_dy[0, :] = (Qxx[1, :] - Qxx[0, :]) / dy
    dq11_tensor_dy[-1, :] = (Qxx[-1, :] - Qxx[-2, :]) / dy

    dq12_tensor_dy[1:-1, :] = (Qxy[2:, :] - Qxy[:-2, :]) / (2 * dy)
    dq12_tensor_dy[0, :] = (Qxy[1, :] - Qxy[0, :]) / dy
    dq12_tensor_dy[-1, :] = (Qxy[-1, :] - Qxy[-2, :]) / dy

    dij = dq11_tensor_dx * dq12_tensor_dy - dq11_tensor_dy*dq12_tensor_dx
    print(dij.shape)
    #normalizing dij
    if normalize:
        dij = (dij-np.min(dij))/(np.max(dij)-np.min(dij)) - 0.5
    return dij

def plot_frame(iframe,pos_data):
    worm_length = 80
    vector_field = get_vector_frame(iframe,worm_length,pos_data)
    num_interpolation_points = 40
    Qxx,Qxy,x_grid,y_grid = get_Q(vector_field,num_interpolation_points)
    nx,ny = get_director(Qxx,Qxy,x_grid,y_grid,num_interpolation_points)
    dij = calc_defects(Qxx,Qxy,x_grid,y_grid)

    max_x = np.max(vector_field[:,0])
    min_x = np.min(vector_field[:,0])
    max_y = np.max(vector_field[:,1])
    min_y = np.min(vector_field[:,1])
    #drawing the cardioid
    rwall=164
    c_a = 1.5*(rwall/2)
    t = np.linspace(0,2*np.pi,1000)
    cx = c_a * (1 - np.cos(t)) * np.cos(t) + rwall + c_a
    cy = c_a * (1 - np.cos(t)) * np.sin(t) + rwall
    for i in np.arange(0,x_grid.shape[0]):
        for j in np.arange(0,y_grid.shape[1]):
            dx = x_grid[i,j] - rwall - c_a
            dy = y_grid[i,j] - rwall
            r = np.sqrt(dx*dx + dy*dy)
            c_r = 1 - np.cos(np.arctan2(dy,dx))
            if r > c_a*c_r:
                nx[i,j] = np.nan
                ny[i,j] = np.nan


    fig,ax2 = plt.subplots(1,1,figsize=(10,10))
    quiver = ax2.quiver(x_grid, #x
                       y_grid, #y
                       nx, #u
                       ny, #v
                       pivot='mid',
                       headlength=0,
                       headwidth=0,
                       headaxislength=0,
                       color="tab:blue",
                       scale_units='xy',
                       scale=0.25
    )
    #ax2.plot(cx,cy,color='k')
    pcm = ax2.imshow(dij,cmap="seismic",interpolation="bicubic",alpha=0.5,origin="lower",extent=(min_x,max_x,min_y,max_y))
    #pcm = ax2.imshow(dij,cmap="hot",alpha=0.5,origin="lower",extent=(min_x,max_x,min_y,max_y))
    cbar = fig.colorbar(pcm, shrink=0.77, ax=ax2,ticks=[-0.5, -0.25, 0, 0.25, 0.5],format=mticker.FixedFormatter(['-0.5', '', '0', '', '0.5']))
    labels = cbar.ax.get_yticklabels()
    labels[0].set_verticalalignment('top')
    labels[-1].set_verticalalignment('bottom')
    ax2.set_box_aspect(1)
    ax2.set_adjustable("datalim")
    ax2.set_title(iframe)
    plt.tight_layout()
    plt.savefig("frame"+str(iframe).zfill(5)+".png")
    #plt.show()
    
try:  
    os.mkdir("./anim_frames")  
except OSError as error:  
    print(error)   
os.chdir("./anim_frames")
frames = np.arange(nframes)
for iframe in tqdm.tqdm(frames[1200:1400]):
    plot_frame(iframe,pos_data)
