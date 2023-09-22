use Math;
use Random;
use IO;
use IO.FormattedIO;
use Time; // for stopwatch
use List;
use MemDiagnostics;


config const debug = false, 
            L = 100.0, // size of the simulation box in x and y
            nsteps = 2,//0000,
            dt = 0.00025,//0.001,
            save_interval = 50,
            numParticles = 144, // number of particles
            thermo = false,
            kbt = 0.5; //reduced temperature

var hxo2 = L/2,
    sigma = 1.0,
    rcut = 2.5*sigma,
    r2cut = rcut**2,
    rcutsmall = sigma*2.0**(1.0/6.0),
    r2cutsmall = rcutsmall*rcutsmall,
    pi = 4.0*atan(1.0),
    gamma = 1.5,
    print_interval = 10;

var ptc_init_counter = 1;
record Particle {
    var id: int;
    var x,y,z: real;
    var vx,vy,vz: real;
    var vxave,vyave,vzave: real;
    var fx,fy,fz: real;
    var fxold,fyold,fzold: real;
    var m: real; //mass
    var ptype: int; //  ptype=1 (active), ptype=2(solvent), ptype=3 (boundary), ptype=-1 (unassigned) 
    proc init() {
        this.id = ptc_init_counter;
        ptc_init_counter += 1;
        this.x = 0.0;
        this.y = 0.0;
        this.z = 0.0;
        this.vx = 0.0;
        this.vy = 0.0;
        this.vz = 0.0;
        this.vxave = 0.0;
        this.vyave = 0.0;
        this.vzave = 0.0;
        this.fx = 0.0;
        this.fy = 0.0;
        this.fz = 0.0;
        this.fxold = 0.0;
        this.fyold = 0.0;
        this.fzold = 0.0;
        this.ptype = -1;
    }
    proc info() {
        var typestring: string;
        if (this.ptype == 1) {
            typestring = "active";
        } else if (this.ptype == 2) {
            typestring = "solvent";
        } else if (this.ptype == 3) {
            typestring = "boundary";
        } else {
            typestring = "unassigned";
        }
        var s:string = "id# %i \t type: %s \t mass: %s \n pos: %r \t %r \t %r \n vel: %r \t %r \t %r \n force: %r \t %r \t %r".format(this.id,typestring,this.m,this.x,this.y,this.z,this.vx,this.vy,this.vz,this.fx,this.fy,this.fz);
        return s;
    }
    proc p(px: real, py: real, pz: real) {
        this.x = px;
        this.y = py;
        this.z = pz;
    }
    proc p() {
        return (this.x,this.y,this.z);
    }
    proc v(velx: real, vely:real, velz:real) {
        this.vx = velx;
        this.vy = vely;
        this.vz = velz;
    }
    proc v() {
        return (this.vx,this.vy,this.vz);
    }
    proc f(forcex:real,forcey:real,forcez:real) {
        this.fx = forcex;
        this.fy = forcey;
        this.fz = forcez;
    }
    proc f() {
        return (this.fx,this.fy,this.fz);
    }
    proc set(p: Particle) {
        this.x = p.x;
        this.y = p.y;
        this.z = p.z;
        this.vx = p.vx;
        this.vy = p.vy;
        this.vz = p.vz;
        this.vxave = p.vxave;
        this.vyave = p.vyave;
        this.vzave = p.vzave;
        this.fx = p.fx;
        this.fy = p.fy;
        this.fz = p.fz;
        this.fxold = p.fxold;
        this.fyold = p.fyold;
        this.fzold = p.fzold;
    }
}

var bin_init_counter = 1;
record Bin {
    var id: (int,int); //id of each bin
    var atoms: list(int); // list of particle id's in each bin
    var neighbors: [1..4][1..2] int; // indices of each bin's neighboring bin 
    var ncount: int; // count of number of particles in each bin
    var x: (real,real); // precalculate the max_x and min_x values for the space the box occupies
    var y: (real,real); // this is for easier neighbor list creation

    proc init() { // record initializer
        this.id = (0,0);
        this.ncount = 0;
        this.x = (0.0,0.0);
        this.y = (0.0,0.0);
        for i in 1..4 {
            this.neighbors[i] = [0,0];
        }
    }
}
var solvent: [1..numParticles] Particle;
var KE: [1..numParticles] real;
var KE_total: [1..nsteps] real;

var numBins = ceil(L/rcut):int;
const binSpace = {1..numBins, 1..numBins};
var bins : [binSpace] Bin;
var randStream = new RandomStream(real); // creating random number generator

proc main () {
    writeln("numBins: ",numBins);
    startVerboseMem();
    init_bins();

    // init particles
    // init in lattice
    if (numParticles != sqrt(numParticles)*sqrt(numParticles)) {writeln("non-square numParticles");halt();}
    var row_length = sqrt(numParticles):int;
    var row,col,center,spacing :real;
    spacing = 1.0;
    center = hxo2 - spacing*(row_length/2);
    for i in 1..numParticles {
        row = i % row_length;
        col = ((i - row)/row_length)-1;
        solvent[i].x = center + spacing*rcutsmall*row;
        solvent[i].y = center + spacing*rcutsmall*col;
        solvent[i].vx = 0.01;//*gaussRand(0.0,1.0);
        solvent[i].vy = 0.0;//0.1*gaussRand(0.0,1.0);
        solvent[i].ptype = 2;
        //KE[i] = 0.0;
        KE[i] = 0.5*(solvent[i].vx * solvent[i].vx + solvent[i].vy * solvent[i].vy);
    }
    
    // init randomly
    

    update_cells(0);

    var t = 0.0;
    var total_time = 0.0;
    var ct: stopwatch, wt:stopwatch, xt:stopwatch; //calc time, io time, totaltime
    write_xyz(0);
    calc_forces_old();
    //calc_forces();
    var vel_mag = 0.0;
    //setting up stopwatch
    xt.start();
    for istep in 1..nsteps {
        if (istep % print_interval == 0) {
	    xt.stop();
	    total_time = xt.elapsed();
            KE_total[istep] = (+ reduce KE);
            vel_mag = sqrt((max reduce solvent.vx)**2 + (max reduce solvent.vy)**2);
            writeln("Step: ",istep,"\t",
                    istep/total_time,"iter/s\tCalc:",
                    (ct.elapsed()/total_time)*100," %\tElapsed:",
                    total_time," s\t",
                    KE_total[istep],"\t",
                    (max reduce solvent.x),"\t",
                    (max reduce solvent.y),"\t",
                    vel_mag);                    
        xt.start();
        }
	    ct.start();
        // update positions
        update_position();
        update_cells(istep);

        calc_forces_old();
        //calc_forces();

        update_velocities();
        

        ct.stop();
        if (istep % save_interval == 0){
            wt.start();
            write_xyz(istep);
            wt.stop();
        }
    }
    write_macro(nsteps);
    xt.stop();
    writeln("Total Time:",xt.elapsed()," s");
    stopVerboseMem();
    printMemAllocStats();
}


proc update_position() {
    // update positions
    forall i in 1..numParticles {
        solvent[i].x += solvent[i].vx*dt + (solvent[i].fx/2)*(dt**2);
        solvent[i].y += solvent[i].vy*dt + (solvent[i].fy/2)*(dt**2);
        //periodic boundary conditions
        // if solvent[i].x > L {
        //     solvent[i].x = solvent[i].x - L;
        // } else if solvent[i].x < 0.0 {
        //     solvent[i].x = solvent[i].x + L;
        // }
        // if solvent[i].y > L {
        //     solvent[i].y = solvent[i].y - L;
        // } else if solvent[i].y < 0.0 {
        //     solvent[i].y = solvent[i].y + L;
        // }
        
        solvent[i].vx += 0.5*dt*solvent[i].fx;
        solvent[i].vy += 0.5*dt*solvent[i].fy;
        // collide with boundaries
        if isclose(solvent[i].x, L, rtol = 1e-3, atol = 0.0) {
            solvent[i].vx = -1.0*solvent[i].vx;
        }
        // if isclose(solvent[i].x, 0.0, rtol = 1e-3, atol = 0.0) {
        //     solvent[i].vx = -1.0*solvent[i].vx;
        //     writeln("hit left wall",solvent[i].vx);
        // }
        if solvent[i].x <= 0.01 {
            solvent[i].vx = -1.0*solvent[i].vx;
        }
        if isclose(solvent[i].y, L, rtol = 1e-3, atol = 0.0) {
            solvent[i].vy = -1.0*solvent[i].vy;
        }
        // if isclose(solvent[i].y, 0.0, rtol = 1e-3, atol = 0.0) {
        //     solvent[i].vy = -1.0*solvent[i].vy;
        //     writeln("hit bottom wall",solvent[i].vy);
        // }
        if (solvent[i].y <= 0.01) {
            solvent[i].vy = -1.0*solvent[i].vy;
        }
        solvent[i].fx = 0.0;
        solvent[i].fy = 0.0;
        solvent[i].fz = 0.0;
    }
}

proc lj(i,j) {
    var dx,dy,r2,ffor,ffx,ffy,dvx,dvy,r,rhatx,rhaty,omega,fdissx,fdissy,gauss,frand :real;
    dx = solvent[j].x - solvent[i].x;
    dy = solvent[j].y - solvent[i].y;
    r2 = (dx*dx + dy*dy);
    //if (debug) {writeln("icount: ",icount,"\t",i,"\tjcount: ",jcount,"\t",j,"\t",r2,"\t",r2cut);}
    if (r2 <= r2cut) {
        // LJ force
        ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
        ffx = ffor*dx;
        ffy = ffor*dy;
        if (debug) {
            if (i == 1) || (j == 1) {
                writeln("LJ ",i," ",j," ",sqrt(r2)," ",sqrt(r2cut)," ",ffx," ",ffy," ",dx," ",dy);
            }
        }
        solvent[i].fx += ffx;
        solvent[i].fy += ffy;
        solvent[j].fx -= ffx;
        solvent[j].fy -= ffy;
        if (thermo) {
            // DPD thermostat
            //adding dissipative force
            dvx = solvent[j].vx - solvent[i].vx;
            dvy = solvent[j].vy - solvent[i].vy;
            r = sqrt(r2);
            rhatx = dx/r;
            rhaty = dy/r;
            omega = (1.0-r2/r2cut);
            fdissx = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty)*rhatx; //gamma = 1/damp (proportional to friction force)
            fdissy = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty)*rhaty;
            solvent[i].fx -= fdissx;
            solvent[i].fy -= fdissy;
            solvent[j].fx += fdissx;
            solvent[j].fy += fdissy;
            // adding random forces
            gauss = gaussRand(0.0,1.0); // generates normal random numbers (mean, stddev)
            //gauss = randStream.getNext();
            frand = (1.0/sqrt(dt))*sqrt(omega)*gauss*sqrt(2.0*kbt*gamma);
            solvent[i].fx += frand*rhatx;
            solvent[i].fy += frand*rhaty;
            solvent[j].fx -= frand*rhatx;
            solvent[j].fy -= frand*rhaty;
        }
    }
}

proc wca(i,j) {

}

proc calc_forces_old () {
    // if (thermo) {
    //     // needed for idiot thermostat
    //     var temp_KE_total = (+ reduce KE);
    //     alpha = sqrt(kbt*numParticles)/temp_KE_total;
    //     if (temp_KE_total <= 0.00001) {
    //         alpha = 1.0;
    //     } 
    //     writeln(alpha);
    // }
    var Lo2 = L/2.0;
    for i in 1..numParticles {
        var dx,dy,r2,ffor,ffx,ffy,dvx,dvy,rhatx,rhaty,r,frand,gauss,fdissx,fdissy,omega,alpha,lb,dp,col_prob:real;
        // calculating the force
        for j in (i+1)..numParticles {
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
            // if (dx > Lo2) {
            //     dx = dx - L;
            // }
            // if (dx < -Lo2) {
            //     dx = dx + L;
            // }
            // if (dy > Lo2) {
            //     dy = dy - L;
            // }
            // if (dy < -Lo2) {
            //     dy = dy + L;
            // }
            r2 = (dx*dx + dy*dy);
            if (r2 <= r2cut) {
                // LJ force
                ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
                ffx = ffor*dx;
                ffy = ffor*dy;
                solvent[i].fx += ffx;
                solvent[i].fy += ffy;
                solvent[j].fx -= ffx;
                solvent[j].fy -= ffy;
                if (thermo) {
                    // DPD thermostat
                    //adding dissipative force
                    dvx = solvent[j].vx - solvent[i].vx;
                    dvy = solvent[j].vy - solvent[i].vy;
                    r = sqrt(r2);
                    rhatx = dx/r;
                    rhaty = dy/r;
                    omega = (1.0-r2/r2cut);
                    fdissx = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty)*rhatx; //gamma = 1/damp (proportional to friction force)
                    fdissy = -1.0*gamma*omega*(dvx*rhatx + dvy*rhaty)*rhaty;
                    solvent[i].fx -= fdissx;
                    solvent[i].fy -= fdissy;
                    solvent[j].fx += fdissx;
                    solvent[j].fy += fdissy;
                    // adding random forces
                    gauss = gaussRand(0.0,1.0); // generates normal random numbers (mean, stddev)
                    //gauss = randStream.getNext();
                    frand = (1.0/sqrt(dt))*sqrt(omega)*gauss*sqrt(2.0*kbt*gamma);
                    solvent[i].fx += frand*rhatx;
                    solvent[i].fy += frand*rhaty;
                    solvent[j].fx -= frand*rhatx;
                    solvent[j].fy -= frand*rhaty;
                }
            }
        }
    }
}

proc calc_forces () {
    for ibin in binSpace{
        // calculate forces inside this bin
        if bins[ibin].ncount > 1 { // check if atoms are in this bin
            // iterate over all the atoms in this bin
            //if (debug) {writeln("in bin int ",bins[ibin].id,"\t",bins[ibin].atoms,"\t");}
            for icount in 0..bins[ibin].ncount-2 {
                var i = bins[ibin].atoms[icount];
                //if (debug) {writeln("icount: ",icount,"\t",bins[ibin].atoms[icount]);}
                for jcount in icount+1..bins[ibin].ncount-1 {
                    var j = bins[ibin].atoms[jcount];
                    // if (i == 1) || (j == 1) {
                    //     writeln("in current cell ",ibin);
                    //     writeln(solvent[i].info());
                    //     writeln(solvent[j].info());
                    // }
                    lj(i,j);
                    }
            }
            // calculate forces inside neighbor bins
            for nab in bins[ibin].neighbors { //loop over each neighboring bin (inab,jnab)
                var inab = nab[1];
                var jnab = nab[2];
                //if (debug) {writeln("current",ibin," neighbor ",nab,"\t",bins[inab,jnab].atoms,"\t",bins[inab,jnab].id);}
                if ((inab != -1) || (jnab != -1) && (bins[inab,jnab].ncount > 1)) { // checks if neighbor is valid
                for icount in 1..bins[ibin].ncount { // loop over the atoms in this bin
                    var i = bins[ibin].atoms[icount-1];
                    for jcount in 1..bins[inab,jnab].ncount { // loop over all the atoms in the neighboring bin
                        var j = bins[inab,jnab].atoms[jcount-1];
                        // if (i == 1) || (j == 1) {
                        //     writeln("in neighbor cell ",nab);
                        //     writeln(solvent[i].info());
                        //     writeln(solvent[j].info());
                        // }
                        lj(i,j);
                        }
                    }
                }
                }
            }
        }
}

proc update_velocities() {
    //update velocities
    forall i in 1..numParticles {
        //solvent[i].vx += 0.5*dt*(solvent[i].fxold + solvent[i].fx);
        //solvent[i].vy += 0.5*dt*(solvent[i].fyold + solvent[i].fy);
        //solvent[i].fxold = solvent[i].fx;
        //solvent[i].fyold = solvent[i].fy;
        solvent[i].vx += 0.5*dt*solvent[i].fx;
        solvent[i].vy += 0.5*dt*solvent[i].fy;
        // calculating kinetic energy here too
        KE[i] = 0.5*(solvent[i].vx * solvent[i].vx + solvent[i].vy * solvent[i].vy);
    }
}

// CELL LIST FUNCTIONS
proc update_cells(istep:int) {
    for ibin in binSpace{
        if bins[ibin].ncount > 0{
            bins[ibin].ncount = 0;
            bins[ibin].atoms.clear();
        }
    }

    for i in 1..numParticles { // populating particles into bins
        var binposx = 1+floor(solvent[i].x/rcut):int;
        var binposy = 1+floor(solvent[i].y/rcut):int;
        //if (debug) {writeln(i,"\t",solvent[i].x,"\t",solvent[i].y,"\t",binposx,"\t",binposy);}
        if (debug) {
            if (binposx > numBins) {
                dump_particles();
            } else if (binposy > numBins) {
                dump_particles();
            }else if (binposx < 1) {
                dump_particles();
            } else if (binposy < 1) {
                dump_particles();
            }
        }
        if (binposx > numBins) {
            write_xyz(istep);
        } else if (binposy > numBins) {
            write_xyz(istep);
        }else if (binposx < 1) {
            write_xyz(istep);
        } else if (binposy < 1) {
            write_xyz(istep);
        }
        bins[binposx,binposy].ncount += 1;
        bins[binposx,binposy].atoms.pushBack(solvent[i].id);
        //writeln(p.id,"\t",posx,"\t",posy,"\t",bins[posx,posy].id," ",bins[posx,posy].ncount);
    }
    if (debug) {
        writeln("recalculating bins");
        for ibin in binSpace {
            if bins[ibin].ncount > 0 {
            writeln(bins[ibin].id,"\t",bins[ibin].atoms);
            }
        }
    }
    /*
        //bins[24,24].atoms.clear();
        writeln("*",bins[24,24].atoms,bins[24,24].atoms.isEmpty());
        for a in bins[24,24].atoms {
            writeln(solvent[a].info());
        }
        var ic = 0;
        writeln(bins[24,20].ncount);
        var inds = bins[24,20].atoms.indices;
        writeln(bins[24,20].atoms);
        writeln(bins[24,20].atoms[inds]);
        var inds2 = inds[1..];
        writeln(inds2);
        writeln(bins[24,20].atoms[inds2]); 
        for atom1 in bins[24,20].atoms[inds] { 
            for atom2 in bins[24,20].atoms[inds2] {
                if solvent[atom1].id != solvent[atom2].id { //make sure we 
                    writeln(atom1,"\t",solvent[atom1].id,"\t",atom2,"\t",solvent[atom2].id);
                    ic += 1;
                }
            }
        }
        writeln(ic);
    */
}

proc init_bins() {
    // writeln((4/numBins):int +1); //row
    // writeln(4%numBins); // col
    for ibin in binSpace {
        bins[ibin].x[0] = (ibin[0]-1)*rcut;
        bins[ibin].x[1] = ibin[0]*rcut;
        bins[ibin].y[0] = (ibin[1]-1)*rcut;
        bins[ibin].y[1] = ibin[1]*rcut;
        bins[ibin].id[0] = ibin[0];
        bins[ibin].id[1] = ibin[1];
        var i = ibin[0];
        var j = ibin[1];
        //  4   3   2
        //     *   1
        //  
        /*  1 = i+1,j
            2 = i+1,j+1
            3 = i,j+1
            4 = i-1,j+1
        */
        var ip1,im1,jp1,jm1:int;
        ip1 = i + 1;
        im1 = i - 1;
        jp1 = j + 1;
        jm1 = j - 1;
        if ip1 > numBins {
            ip1 = -1; // no neighbors in this direction
        }
        if jp1 > numBins {
            jp1 = -1;
        }
        if im1 == 0 {
            im1 = -1;
        }
        if jm1 == 0{
            jm1 = -1;
        }
        bins[ibin].neighbors[1] = [ip1,j];
        bins[ibin].neighbors[2] = [ip1,jp1];
        bins[ibin].neighbors[3] = [i,jp1];
        bins[ibin].neighbors[4] = [im1,jp1];
    }
}

// I/O FUNCTIONS
proc write_xyz(istep:int) {
    var ic:int;
    var filename:string = "simple%{07u}.xyz".format(istep);
    //var filename = "amatter" + (istep:string) + ".xyz";
    try {
        var xyzfile = open(filename, ioMode.cw);
        var myFileWriter = xyzfile.writer();
        myFileWriter.writeln(numParticles + 4);
        myFileWriter.writeln("# 0");
        ic = 2;
        for i in 1..numParticles {
            myFileWriter.writeln("S ",solvent[i].x," ",solvent[i].y," ",0.0);
            ic += 1;
        }
        myFileWriter.writeln("E ",0.0," ",0.0," ",0.0);
        myFileWriter.writeln("E ",L," ",0.0," ",0.0);
        myFileWriter.writeln("E ",0.0," ",L," ",0.0);
        myFileWriter.writeln("E ",L," ",L," ",0.0);
        ic += 4;
        writeln(filename,"\t",ic," lines written");
    } catch e: Error {
        writeln(e);
    }
}

proc write_macro(nsteps: int) {
    var filename:string = "energies.dat";
    //var vel_mag = sqrt((max reduce solvent.vx)**2 + (max reduce solvent.vy)**2);
    try{
        var datfile = open(filename, ioMode.cw);
        var myFileWriter = datfile.writer();
        myFileWriter.writeln("step \t KE \t maxX \t maxY \t vel_mag");
        for istep in 1..nsteps {
            myFileWriter.writeln(istep,"\t",
                                KE_total[istep],"\t");
                                //(max reduce solvent.x),"\t",
                                //(max reduce solvent.y),"\t",
                                //vel_mag);
        }
        datfile.fsync();
        writeln("energies.dat written");
    } catch e: Error {
        writeln(e);
    }
}

proc dump_particles() {
    var filename:string = "particles.dump";
    try {
        var dfile = open(filename, ioMode.cw);
        var myFileWriter = dfile.writer();
        for i in 1..numParticles {
            myFileWriter.writeln(solvent[i].info());
        }
    } catch e: Error {
        writeln(e);
    }
}

// HELPER FUNCTIONS
proc gaussRand(mean: real, stddev: real): real {
    var u1 = randStream.getNext();
    var u2 = randStream.getNext();
    var z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2);  // Box-Muller transform
    var gauss = mean + stddev * z0;
    if (gauss > 2.5) {
        gauss = 2.5;
    } else if (gauss < -2.5) {
        gauss = -2.5;
    }
    return gauss;
}