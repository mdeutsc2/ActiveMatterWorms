module Output {
    use Amatter3d;
    use IO;
    use IO.FormattedIO;
    use CTypes;
    use Math only atan2,sin,cos;
    import C;
    // IO FUNCTIONS
proc write_xyz(istep:int) {
    var dx :real, dy:real, xang:real, rx:real,ry:real,dot:real;
    var ic:int;
    var filename:string = "amatter%{010u}.xyz".format(istep);
    //var filename = "amatter" + (istep:string) + ".xyz";
    try {
    var xyzfile = open(filename, ioMode.cw);
    var myFileWriter = xyzfile.writer(locking=false);
    // number of active particles + 4 edge-defining particles + boundary + solvent (optional)
    if (fluid_cpl) {
        myFileWriter.writeln(nworms * np + 4 + numPoints + numSol);
    } else {
        myFileWriter.writeln(nworms*np + 4 + numPoints);
    }
    myFileWriter.writeln("# 0");
    ic = 2;
    for iw in 1..nworms {
        dx = worms[iw,1].x - hxo2;
        dy = worms[iw,1].y - hyo2;
        xang = atan2(dy,dx);
        rx = -sin(xang);
        ry = cos(xang);
        dot = (worms[iw,1].x - worms[iw,np].x)*rx + (worms[iw,1].y - worms[iw,np].y)*ry;
        if (dot >= 0.0) {
            for i in 1..np {
                myFileWriter.writeln("A ",worms[iw,i].x," ",worms[iw,i].y," ", worms[iw,i].z);
                //myFileWriter.writeln("A ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        } else {
            for i in 1..np {
                myFileWriter.writeln("B ",worms[iw,i].x," ",worms[iw,i].y," ", worms[iw,i].z);
                //myFileWriter.writeln("B ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        }
    }
    for i in 1..numPoints {
        myFileWriter.writeln("I ",bound[i].x," ",bound[i].y," ",0.0);
        ic += 1;
    }
    if (fluid_cpl) {
        for i in 1..numSol {
            myFileWriter.writeln("S ",solvent[i].x," ",solvent[i].y," ",0.0);
            ic += 1;
        }
    }
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0);
    //ic += 4;
    // myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    //writeln(filename,"\t",ic," lines written");
    //xyzfile.close();
    } catch e: Error {
        writeln(e);
    }
}

proc write_xyzv(istep:int) {
    var dx :real, dy:real, xang:real, rx:real,ry:real,dot:real;
    var ic:int;
    var filename:string;
    try! {
          filename = "amatter%{010u}.xyz".format(istep);
    }
    //var filename = "amatter" + (istep:string) + ".xyz";
    try {
    var xyzfile = open(filename, ioMode.cw);
    var myFileWriter = xyzfile.writer(locking=false);
    // number of active particles + 4 edge-defining particles + boundary + solvent (optional)
    if (fluid_cpl) {
        myFileWriter.writeln(nworms * np + 4 + numPoints + numSol);
    } else {
        myFileWriter.writeln(nworms*np + 4 + numPoints);
    }
    myFileWriter.writeln("# 0");
    ic = 2;
    for iw in 1..nworms {
        dx = worms[iw,1].x - hxo2;
        dy = worms[iw,1].y - hyo2;
        xang = atan2(dy,dx);
        rx = -sin(xang);
        ry = cos(xang);
        dot = (worms[iw,1].x - worms[iw,np].x)*rx + (worms[iw,1].y - worms[iw,np].y)*ry;
        if (dot >= 0.0) {
            for i in 1..np {
                myFileWriter.writeln("A ",worms[iw,i].x," ",worms[iw,i].y," ", worms[iw,i].z," ",worms[iw,i].vx," ",worms[iw,i].vy," ",worms[iw,i].vz);
                //myFileWriter.writeln("A ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        } else {
            for i in 1..np {
                myFileWriter.writeln("B ",worms[iw,i].x," ",worms[iw,i].y," ",worms[iw,i].z," ",worms[iw,i].vx," ",worms[iw,i].vy," ",worms[iw,i].vz);
                //myFileWriter.writeln("B ",x[iw,i]," ",y[iw,i]," ", 0.0," ",vx[iw,i]," ",vy[iw,i]," ",0.0);
                ic += 1;
            }
        }
    }
    for i in 1..numPoints {
        myFileWriter.writeln("I ",bound[i].x," ",bound[i].y," ",0.0," ",0.0," ",0.0," ",0.0);
        ic += 1;
    }
    if (fluid_cpl) {
        for i in 1..numSol {
            myFileWriter.writeln("S ",solvent[i].x," ",solvent[i].y," ",0.0," ",solvent[i].vx," ",solvent[i].vy," ",0.0);
            ic += 1;
        }
    }
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    //ic += 4;
    // myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 - rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 - rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    // myFileWriter.writeln("E ",hxo2 + rwall," ",hyo2 + rwall," ",0.0," ",0.0," ",0.0," ",0.0);
    write_log(logfile,filename+"\t"+ic:string+" lines written",true);
    //xyzfile.close();
    } catch e: Error {
        writeln(e);
    }
}

proc init_macro(filename:string) {
    if C.appendFileExists(filename.c_str()) != 0 {
        var error_str:string = filename+" already exists";
        halt(error_str);
    } else {
        try {
            var datfile = open(filename, ioMode.cw);
            var myFileWriter = datfile.writer(locking=false);
            myFileWriter.writeln("step \t KEworm \t KEworm-vave \t KEsol");
            datfile.fsync();
            write_log(logfile,filename+" header written");
        } catch e: Error {
            writeln(e);
        }
    }
}

proc write_macro(filename: string,istep: int) {
    var ret_int:int;
    var out_str:string = istep:string+"\t"+KEworm_total[istep]:string+"\t"+KEworm_local_total[istep]:string+"\t"+KEsol_total[istep]:string+"\t"+AMworm_total[istep]:string+"\t"+AMsol_total[istep]:string;
    ret_int = C.appendToFile(filename.c_str(),out_str.c_str());
    if ret_int != 0 {
        var err_str:string = "error in appending to "+filename;
        halt(err_str);
    }
}

proc init_log(filename:string) {
    if C.appendFileExists(filename.c_str()) != 0 {
        var error_str:string = filename+" already exists";
        halt(error_str);
    }
}

proc write_log(filename: string, out_str: string, silent: bool = false) {
    if silent == false {
      writeln(out_str);
    }
    var ret_int: int;
    ret_int = C.appendToFile(filename.c_str(),out_str.c_str());
    if ret_int != 0 {
        var err_str:string = "error in appending to "+filename;
        halt(err_str);
    }
}

proc write_params() {
    var filename:string = "params.dat";
    try {
        var paramsfile = open(filename, ioMode.cw);
        var myFileWriter = paramsfile.writer(locking=false);
        myFileWriter.writeln("np\t",np.type:string,"\t",np);
        myFileWriter.writeln("nworms\t",nworms.type:string,"\t",nworms);
        myFileWriter.writeln("nsteps\t",nsteps.type:string,"\t",nsteps);
        myFileWriter.writeln("fdogic\t",fdogic.type:string,"\t",fdogic);
        myFileWriter.writeln("walldrive\t",walldrive.type:string,"\t",walldrive);
        myFileWriter.writeln("fdogicwall\t",fdogicwall.type:string,"\t",fdogicwall);
        myFileWriter.writeln("fdep\t",fdep.type:string,"\t",fdep);
        myFileWriter.writeln("dogic_fdep\t",dogic_fdep.type:string,"\t",dogic_fdep);
        myFileWriter.writeln("fdepwall\t",fdepwall.type:string,"\t",fdepwall);
        myFileWriter.writeln("dt\t",dt.type:string,"\t",dt);
        myFileWriter.writeln("kspring\t",kspring.type:string,"\t",kspring);
        myFileWriter.writeln("k2spring\t",k2spring.type:string,"\t",k2spring);
        myFileWriter.writeln("k3psring\t",k3spring.type:string,"\t",k3spring);
        myFileWriter.writeln("kbend\t",kbend.type:string,"\t",kbend);
        myFileWriter.writeln("length0\t",length0.type:string,"\t",length0);
        myFileWriter.writeln("rcut\t",rcut.type:string,"\t",rcut);
        myFileWriter.writeln("save_interval\t",save_interval.type:string,"\t",save_interval);
        myFileWriter.writeln("boundary\t",bd.type:string,"\t",bd.t);
        myFileWriter.writeln("fluid_cpl\t",fluid_cpl.type:string,"\t",fluid_cpl);
        myFileWriter.writeln("debug\t",debug.type:string,"\t",debug);
        myFileWriter.writeln("thermo\t",thermo.type:string,"\t",thermo);
        myFileWriter.writeln("thermow\t",thermo.type:string,"\t",thermow);
        myFileWriter.writeln("kbt\t",kbt.type:string,"\t",kbt);
        myFileWriter.writeln("sigma\t",sigma.type:string,"\t",sigma);
        myFileWriter.writeln("gamma\t",gamma.type:string,"\t",gamma);
        myFileWriter.writeln("fdrag\t",fdrag.type:string,"\t",fdrag);
        myFileWriter.writeln("numPoints\t",numPoints.type:string,"\t",numPoints);
        myFileWriter.writeln("fluid_rho\t",fluid_rho.type:string,"\t",fluid_rho);
        myFileWriter.writeln("random_init\t",random_init.type:string,"\t",random_init);
        myFileWriter.writeln("numSol\t",numSol.type:string,"\t",numSol);
        myFileWriter.writeln("fluid_offset\t",fluid_offset.type:string,"\t",fluid_offset);
        myFileWriter.writeln("worm_particl_mass\t",worm_particle_mass.type:string,"\t",worm_particle_mass);
        myFileWriter.writeln("sw_epsilon\t",sw_epsilon.type:string,"\t",sw_epsilon);
        myFileWriter.writeln("ww_epsilon\t",ww_epsilon.type:string,"\t",ww_epsilon);
        myFileWriter.writeln("L\t",L.type:string,"\t",L);
        paramsfile.fsync();
        //paramsfile.close();
        write_log(logfile,"params.dat written");
    } catch e: Error {
        writeln(e);
    }
}

proc restart_write(istep:int) {
   /*
    File structures
    all lines starting with '#' are comment lines and contain data about the rest of the simulation
    First # line is the time step of the restart
    First, all worms will be written out, starting with 2 comment lines, the first containing nworms, 2nd containing np
    The data will then be as following with nworms*np lines
    iworm:int ip:int x:real y:real z:real vx:real vy:real vz:real
    The next comment line will denote the beginning of the boundary section and contains numPoints
    The data will then be as follows with numPoints lines
    i:int, x:real, y:real
    The next comment line will denote the beginning of the fluid section (if applicable) and contains numSol
    The data will the be as follows with numSol lines
    i:int x:real, y:real vx:real vy:real
    */
    var filename:string;
    if istep == 0 {
      filename = "amatter_init.restart";
    } else {
      filename = "amatter.restart";
    }
    try {
        var restart_file = open(filename, ioMode.cw);
        var restartWriter = restart_file.writer(locking=false);
        // total number of particles = number of active particles + boundary + solvent (optional)
        restartWriter.writeln("# ",istep);
        restartWriter.writeln("#nworms:",nworms);
        restartWriter.writeln("#np:",np);
        for iw in 1..nworms {
            for i in 1..np {
                restartWriter.writeln(iw," ",i," ",worms[iw,i].x," ",worms[iw,i].y," ", worms[iw,i].z," ",worms[iw,i].vx," ",worms[iw,i].vy," ",worms[iw,i].vz);
            }
        }
        restartWriter.writeln("#numPoints:",numPoints);
        for i in 1..numPoints {
            restartWriter.writeln(i," ",bound[i].x," ",bound[i].y);
        }
        if fluid_cpl {
            restartWriter.writeln("#numSol:",numSol);
            for i in 1..numSol {
                restartWriter.writeln(i," ",solvent[i].x," ",solvent[i].y," ",solvent[i].vx," ",solvent[i].vy);
            }
        }
        //restart_file.fsync();
        //restart_file.close();
    } catch e: Error {
        writeln(e);
    }
}

proc restart_read(filename) {
    var istep_read =  0;
    try {
    var restart_file = open(filename, ioMode.r).reader(locking=false);
    var read_sucess :bool;
    var linestring :string;
    read_sucess = restart_file.readLine(linestring,stripNewline=true);
    writeln("Reading restart from step ", linestring);
    istep_read = linestring.partition(" ")[2]:int;

    read_sucess = restart_file.readLine(linestring,stripNewline=true);
    var nworms_read = linestring.partition(":")[2]:int;
    read_sucess = restart_file.readLine(linestring,stripNewline=true);
    var np_read = linestring.partition(":")[2]:int;
    if worms.shape != (nworms_read,np_read) {
        writeln("nworms: ",nworms," ",nworms_read);
        writeln("np: ",np," ",np_read);
        halt("shape of input from restart wrong!");
    }
    for iline in 1..(nworms_read*np_read) {
        read_sucess = restart_file.readLine(linestring,stripNewline=true);
        var istring : [1..8] string;
        for (i,isep) in zip(1..8,linestring.split(" ",-1)) {
            istring[i] = isep;
        }
        var iw = istring[1]:int;
        var ip = istring[2]:int;
        worms[iw,ip].x = istring[3]:real;
        worms[iw,ip].y = istring[4]:real;
        worms[iw,ip].z = istring[5]:real;
        worms[iw,ip].vx = istring[6]:real;
        worms[iw,ip].vy = istring[7]:real;
        worms[iw,ip].vz = istring[8]:real;
        worms[iw,ip].fx = 0.0;
        worms[iw,ip].fy = 0.0;
        worms[iw,ip].fz = 0.0;
        worms[iw,ip].fxold = 0.0;
        worms[iw,ip].fyold = 0.0;
        worms[iw,ip].fzold = 0.0;
        worms[iw,ip].ptype = 1;
        worms[iw,ip].m = worm_particle_mass;
    }
    read_sucess = restart_file.readLine(linestring,stripNewline=true);
    var numPoints_read = linestring.partition(":")[2]:int;
    if (numPoints != numPoints_read) {
        writeln("numPoints:",numPoints," ",numPoints_read);
        halt("shape of input from restart wrong");
    }
    for iline in 1..numPoints_read {
        read_sucess = restart_file.readLine(linestring,stripNewline=true);
        var istring : [1..3] string;
        for (i,isep) in zip(1..3,linestring.split(" ",-1)) {
            istring[i] = isep;
        }
        var ib = istring[1]:int;
        bound[ib].x = istring[2]:real;
        bound[ib].y = istring[3]:real;
        bound[ib].z = 0.0;
        bound[ib].ptype = 3;
    }
    //assert(bound==bound2);
    if (fluid_cpl) {
        read_sucess = restart_file.readLine(linestring,stripNewline=true);
        var numSol_read = linestring.partition(":")[2]:int;
        if (numSol != numSol_read) {
            writeln("numSol:",numSol," ",numSol_read);
            halt("shape of input from restart wrong");
        }
        for iline in 1..numSol_read {
            read_sucess = restart_file.readLine(linestring,stripNewline=true);
            var istring : [1..5] string;
            for (i,isep) in zip(1..5,linestring.split(" ",-1)) {
                istring[i] = isep;
            }
            var ip = istring[1]:int;
            solvent[ip].x = istring[2]:real;
            solvent[ip].y = istring[3]:real;
            solvent[ip].z = 0.0;
            solvent[ip].vx = istring[4]:real;
            solvent[ip].vy = istring[5]:real;
            solvent[ip].vz = 0.0;
            solvent[ip].m = 1.0;
            solvent[ip].ptype = 2;
        }
        //assert(solvent==solvent2);
    }
    //restart_file.fsync();
    //restart_file.close();
    } catch e: Error {
        writeln(e);
    }
    return istep_read;
}
}