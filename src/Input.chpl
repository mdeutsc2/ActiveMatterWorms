module Input {
    use CTypes;
    import C;
    use IO;

    proc write_xyzv(istep:int) {
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
    //write_log(logfile,filename+"\t"+ic:string+" lines written",true);
    //xyzfile.close();
    } catch e: Error {
        writeln(e);
    }
}
}