import math,os,argparse
import numpy as np
from ovito.io import import_file
from ovito.modifiers import *
from ovito.vis import Viewport
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support


def render_all_positions(filename):
    ''' renders all the worm and solvent position data'''
    pipeline = import_file(filename,columns = ['Particle Type','Position.X','Position.Y','Position.Z','Velocity.X','Velocity.Y','Velocity.Z'])
    print("render_all_positions")
    print("Running ",filename,pipeline.source)
    # modifiers
    def compute_particle_colors(frame,data):
        # ovito for some idiotic reason wants rgb in 0-1 inclusive so...
        typeA = np.array([255,102,102])/255.0# 1
        typeB = np.array([163,215,248])/255.0 # 2
        typeE = np.array([102,255,51])/255.0 # 5
        typeI = np.array([148,0,148])/255.0 # 3
        typeS = np.array([179,179,0])/255.0 # 4
        rgb = np.empty((data.particles.count,3))

        typeA_loc = np.where(data.particles['Particle Type'][:] == 1)
        typeB_loc = np.where(data.particles['Particle Type'][:] == 2)
        typeE_loc = np.where(data.particles['Particle Type'][:] == 5)
        typeI_loc = np.where(data.particles['Particle Type'][:] == 3)
        typeS_loc = np.where(data.particles['Particle Type'][:] == 4)
        # setting particle type A (1)
        rgb[typeA_loc,:] = typeA
        rgb[typeB_loc,:] = typeB
        rgb[typeE_loc,:] = typeE
        rgb[typeI_loc,:] = typeI
        rgb[typeS_loc,:] = typeS

        data.particles_.create_property('Color',data=rgb)

    pipeline.modifiers.append(ComputePropertyModifier(
        output_property = 'Radius',
        expressions = ['0.5']
    ))
    pipeline.modifiers.append(compute_particle_colors)


    data = pipeline.compute()
    print("# particles",data.particles.count)
    print(list(data.particles.keys()))

    print("Particle Types")
    for ptc_type in data.particles['Particle Type'].types:
        count = np.where(data.particles['Particle Type'][:] == ptc_type.id)[0].shape[0]
        print(ptc_type.id,ptc_type.name,ptc_type.color,count)

    pipeline.add_to_scene()

    vp = Viewport()
    vp.type = Viewport.Type.Top
    vp.camera_pos = (74.5,74.5,0)
    vp.camera_up = (0,1,0)
    vp.camera_dir = (0,0,-1)
    vp.zoom_all(size=(1000,1000))
    if alt_name == None:
        out_name = filename.split('.')[0] + "_all.mp4"
    else:
        out_name = alt_name
    print("writing to: ",out_name)
    vp.render_anim(filename = out_name,
                   size=(1000,1000),
                   fps=60,
                   background = (1.0,1.0,1.0),
                   stop_on_error=True)

def render_solvent_positions(filename,alt_name = None):
    ''' renders ONLY solvent position data'''
    pipeline = import_file(filename,columns = ['Particle Type','Position.X','Position.Y','Position.Z','Velocity.X','Velocity.Y','Velocity.Z'])
    print("render_solvent_positions")
    print("Running ",filename,pipeline.source)
    # modifiers
    def compute_particle_colors(frame,data):
        # ovito for some idiotic reason wants rgb in 0-1 inclusive so...
        typeA = np.array([255,102,102])/255.0# 1
        typeB = np.array([163,215,248])/255.0 # 2
        typeE = np.array([102,255,51])/255.0 # 5
        typeI = np.array([148,0,148])/255.0 # 3
        typeS = np.array([179,179,0])/255.0 # 4
        rgb = np.empty((data.particles.count,3))

        typeA_loc = np.where(data.particles['Particle Type'][:] == 1)
        typeB_loc = np.where(data.particles['Particle Type'][:] == 2)
        typeE_loc = np.where(data.particles['Particle Type'][:] == 5)
        typeI_loc = np.where(data.particles['Particle Type'][:] == 3)
        typeS_loc = np.where(data.particles['Particle Type'][:] == 4)
        # setting particle type A (1)
        rgb[typeA_loc,:] = typeA
        rgb[typeB_loc,:] = typeB
        rgb[typeE_loc,:] = typeE
        rgb[typeI_loc,:] = typeI
        rgb[typeS_loc,:] = typeS

        data.particles_.create_property('Color',data=rgb)

    # def vanish_worms_radius(frame,data):
    #     radius = np.empty((data.particles.count))
    #     radius = 0.5
    #     typeA_loc = np.where(data.particles['Particle Type'][:] == 1)
    #     typeB_loc = np.where(data.particles['Particle Type'][:] == 2)
    #     radius[typeA_loc] = 0.00001
    #     radius[typeB_loc] = 0.00001
    #     data.particles_.create_property('Radius',data=radius)
    def setup_particle_types(frame,data):
        types = data.particles_.particle_types_
        types.type_by_id_(1).name = "A"
        types.type_by_id_(1).radius = 0.5
        types.type_by_id_(1).color = np.array([255,102,102])/255.0

        types.type_by_id_(2).name = "B"
        types.type_by_id_(2).radius = 0.5
        types.type_by_id_(2).color = np.array([163,215,248])/255.0

        types.type_by_id_(3).name = "I"
        types.type_by_id_(3).radius = 0.5
        types.type_by_id_(3).color = np.array([148,0,148])/255.0

        types.type_by_id_(4).name = "S"
        types.type_by_id_(4).radius = 0.5
        types.type_by_id_(4).color = np.array([179,179,0])/255.0

        types.type_by_id_(5).name = "E"
        types.type_by_id_(5).radius = 1.0
        types.type_by_id_(5).color = np.array([102,255,51])/255.0
    pipeline.modifiers.append(setup_particle_types)

    # pipeline.modifiers.append(compute_particle_colors)
    # pipeline.modifiers.append(vanish_worms_radius)

    data = pipeline.compute()
    print("# particles",data.particles.count)
    print(list(data.particles.keys()))

    print("Particle Types")
    for ptc_type in data.particles['Particle Type'].types:
        count = np.where(data.particles['Particle Type'][:] == ptc_type.id)[0].shape[0]
        print(ptc_type.id,ptc_type.name,ptc_type.color,count)

    pipeline.add_to_scene()

    vp = Viewport()
    vp.type = Viewport.Type.Top
    vp.camera_pos = (74.5,74.5,0)
    vp.camera_up = (0,1,0)
    vp.camera_dir = (0,0,-1)
    vp.zoom_all(size=(1000,1000))
    if alt_name == None:
        out_name = filename.split('.')[0] + "_solvent.mp4"
    else:
        out_name = alt_name
    print("writing to: ",out_name)
    vp.render_anim(filename = out_name,
                   size=(1000,1000),
                   fps=60,
                   background = (1.0,1.0,1.0),
                   stop_on_error=True)

def main(args):
    print(args.filename)
    render_solvent_positions(args.filename,"test_anim.mp4")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('filename')
    args = parser.parse_args()
    main(args)
