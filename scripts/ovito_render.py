import math,os,argparse,multiprocessing,tqdm,time
import numpy as np
from ovito.io import import_file
from ovito.modifiers import *
from ovito.vis import Viewport

os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support

class Render:
    def all_positions(filename,debug=False,alt_name = None):
        ''' renders all the worm and solvent position data'''
        pipeline = import_file(filename,columns = ['Particle Type','Position.X','Position.Y','Position.Z','Velocity.X','Velocity.Y','Velocity.Z'])
        if debug:
            print("render_all_positions")
            print("Running ",filename,pipeline.source)
        # modifiers
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


        data = pipeline.compute()
        if debug:
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
        if debug:
            print("writing to: ",out_name)
        vp.render_anim(filename = out_name,
                       size=(1000,1000),
                       fps=60,
                       background = (1.0,1.0,1.0),
                       stop_on_error=True)
        return out_name

    def solvent_positions(filename,debug=False,alt_name = None):
        ''' renders ONLY solvent position data'''
        pipeline = import_file(filename,columns = ['Particle Type','Position.X','Position.Y','Position.Z','Velocity.X','Velocity.Y','Velocity.Z'])
        if debug:
            print("render_solvent_positions")
            print("Running ",filename,pipeline.source)
        # modifiers
        def setup_particle_types(frame,data):
            types = data.particles_.particle_types_
            types.type_by_id_(1).name = "A"
            types.type_by_id_(1).radius = 0.001
            types.type_by_id_(1).color = np.array([255,102,102])/255.0

            types.type_by_id_(2).name = "B"
            types.type_by_id_(2).radius = 0.001
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

        data = pipeline.compute()
        if debug:
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
        if debug:
            print("writing to: ",out_name)
        vp.render_anim(filename = out_name,
                       size=(1000,1000),
                       fps=60,
                       background = (1.0,1.0,1.0),
                       stop_on_error=True)
        return out_name

    def worm_positions(filename,debug=False,alt_name = None):
        ''' renders ONLY worm position data'''
        pipeline = import_file(filename,columns = ['Particle Type','Position.X','Position.Y','Position.Z','Velocity.X','Velocity.Y','Velocity.Z'])
        if debug:
            print("render_worm_positions")
            print("Running ",filename,pipeline.source)
        # modifiers
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
            types.type_by_id_(4).radius = 0.001
            types.type_by_id_(4).color = np.array([179,179,0])/255.0

            types.type_by_id_(5).name = "E"
            types.type_by_id_(5).radius = 1.0
            types.type_by_id_(5).color = np.array([102,255,51])/255.0
        pipeline.modifiers.append(setup_particle_types)

        data = pipeline.compute()
        if debug:
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
            out_name = filename.split('.')[0] + "_worms.mp4"
        else:
            out_name = alt_name
        if debug:
            print("writing to: ",out_name)
        vp.render_anim(filename = out_name,
                       size=(1000,1000),
                       fps=60,
                       background = (1.0,1.0,1.0),
                       stop_on_error=True)
        return out_name


def reencode_video(f,args):
    # create temporary directory and move all old mp4s to directory
    fold = f.split('.')[0] + "-264.mp4"
    os.rename(f,fold)
    os.system("ffmpeg -i {0} -vcodec libx265 -crf 18 -c copy -tag:v hvc1 {1}".format(fold,f)) # needs hvc1 media tag for macos playback
    # os.system("ffmpeg -i {0} -vcodec libx265 -crf 18 {1}".format(fold,f))
    os.remove(fold)

def process_file(f,args):
    start_time = time.time()
    debug = args.verbose
    outfiles = []
    if args.all == True:
        args.solvent = True
        args.worms = True
        args.velocity = True
    if args.worms:
        # Render.worm_positions(f,debug,"test_anim_worms.mp4")
        fname = Render.worm_positions(f,debug)
        outfiles.append(fname)
    if args.solvent:
        # Render.solvent_positions(f,debug,"test_anim_solvent.mp4")
        fname = Render.solvent_positions(f,debug)
        outfiles.append(fname)
    if args.all:
        # Render.all_positions(f,debug,"test_anim_all.mp4")
        fname = Render.all_positions(f,debug)
        outfiles.append(fname)
    if args.velocity:
        print("velocity not done yet")
        # Render.velocity_field(f)
    end_time = time.time()
    total_time = end_time - start_time
    if args.ncpus > 1:
        process_name = multiprocessing.current_process().name
        print(f"Processed {f} in {total_time} s on {process_name}",)
    else:
        print(f"Processed {f} in {total_time} s")
    return outfiles

def main(args):
    print(args.filename)
    if args.ncpus < 0:
        args.ncpus = multiprocessing.cpu_count()

    print("# cpus",args.ncpus)
    print("Re-encode? ",args.reencode)

    # SERIAL
    if args.ncpus > 1:
        # PARALLEL
        pool = multiprocessing.Pool(processes=ncpus)
        #pool.starmap(process_file,tqdm.tqdm([(f,args) for f in args.filename],total=len(args.filename)))
        pool.starmap(process_file,[(f,args) for f in args.filename])

        pool.close()
        pool.join()
    else:
        for f in args.filename:
            outfiles = process_file(f,args)
            if args.reencode:
                for f in outfiles:
                    reencode_video(f,args)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
             prog='ovito_render.py',
             description='renders various animations and visualizations of .xyzv data from amatter simuation',
             epilog='Matt Deutsch, Kent State University (2023)'
    )

    parser.add_argument('filename',nargs="+",help="<Required> data filename or list of filenames")
    parser.add_argument("-a","--all",help="render all",action="store_false")
    parser.add_argument("-s","--solvent",help="only render solvent data",action="store_true")
    parser.add_argument("-w","--worms",help="only render worm data",action="store_true")
    parser.add_argument("-v","--velocity",help="only render velocity data",action="store_true")
    parser.add_argument("-n","--ncpus",help="number of cpus for parallel data processing, default=0 (serial), -1=all cpus",default=0)
    parser.add_argument("--verbose",help="verbose logging",action="store_true")
    parser.add_argument("--reencode",help="re-encodes videos into H.265 format with minimal loss in quality",action="store_true")
    args = parser.parse_args()
    main(args)


#ffmpeg -i amatter_28_11_2023-2_worms.mp4 -vcodec libx265 -crf 18 amatter_28_11_2023-2_worms_265.mp4
