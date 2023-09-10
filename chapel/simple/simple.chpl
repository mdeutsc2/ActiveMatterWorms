use Math;
use Random;
use IO;
use IO.FormattedIO;
use Time; // for stopwatch

config const L = 100.0, // size of the simulation box in x and y
            nsteps = 10000,
            dt = 0.01,
            save_interval = 250,
            numParticles = 144, // number of particles
            thermo = true,
            kbt = 1.0; //reduced temperature

var hxo2 = L/2,
    sigma = 1.0,
    rcut = 2.5*sigma,
    r2cut = rcut**2,
    rcutsmall = sigma*2.0**(1.0/6.0),
    r2cutsmall = rcutsmall*rcutsmall,
    pi = 4.0*atan(1.0),
    gamma = 1;

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

var solvent: [1..numParticles] Particle;
var KE: [1..numParticles] real;
var KE_total: [1..nsteps] real;

var randStream = new RandomStream(real); // creating random number generator

proc main () {
    // init particles
    // init in lattice
    if (numParticles != sqrt(numParticles)*sqrt(numParticles)) {writeln("non-square numParticles");halt();}
    var row_length = sqrt(numParticles):int;
    var row,col,center,spacing :real;
    spacing = 2.0;
    center = hxo2 - spacing*(row_length/2);
    for i in 1..numParticles {
        row = i % row_length;
        col = ((i - row)/row_length)-1;
        solvent[i].x = center + spacing*rcutsmall*row;
        solvent[i].y = center + spacing*rcutsmall*col;
        KE[i] = 0.0;
    }
    // init randomly

    write_xyz(0);
    calc_forces();
    for istep in 1..nsteps {
        // update positions
        update_position();

        calc_forces();

        update_velocities();
        
        KE_total[istep] = (+ reduce KE);
        writeln(istep,"\t",KE_total[istep],"\t");
        if (istep % save_interval == 0){
            write_xyz(istep);
        }
    }
    write_macro(nsteps);
}


proc update_position() {
    // update positions
    forall i in 1..numParticles {
        solvent[i].x += solvent[i].vx*dt + (solvent[i].fx/2)*(dt**2);
        solvent[i].y += solvent[i].vy*dt + (solvent[i].fy/2)*(dt**2);
        //periodic boundary conditions
        if solvent[i].x > L {
            solvent[i].x = solvent[i].x - L;
        } else if solvent[i].x < 0.0 {
            solvent[i].x = solvent[i].x + L;
        }
        if solvent[i].y > L {
            solvent[i].y = solvent[i].y - L;
        } else if solvent[i].y < 0.0 {
            solvent[i].y = solvent[i].y + L;
        }

        solvent[i].fx = 0.0;
        solvent[i].fy = 0.0;
        solvent[i].fz = 0.0;
    }
}

proc calc_forces () {
    var dx,dy,r2,ffor,ffx,ffy,dvx,dvy,rhatx,rhaty,r,frand,gauss,fdissx,fdissy,omega:real;
    for i in 1..numParticles {
        // calculating the force
        for j in (i+1)..numParticles {
            dx = solvent[j].x - solvent[i].x;
            dy = solvent[j].y - solvent[i].y;
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
                    /*
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
                    */
                }
            }
        }
        // idiot thermostat
        // rescaling velocities
        //var total_kinetic = (+ reduce KEsol)
        var alpha = sqrt(kbt*numParticles)/(+ reduce KE);
        solvent[i].vx = alpha*solvent[i].vx;
        solvent[i].vy = alpha*solvent[i].vy;
        /*
        // Lowe-Anderson Thermostat
        */
    }
}

proc update_velocities() {
    //update velocities
    forall i in 1..numParticles {
        solvent[i].vx += 0.5*dt*(solvent[i].fxold + solvent[i].fx);
        solvent[i].vy += 0.5*dt*(solvent[i].fyold + solvent[i].fy);
        solvent[i].fxold = solvent[i].fx;
        solvent[i].fyold = solvent[i].fy;
        // calculating kinetic energy here too
        KE[i] = 0.5*(solvent[i].vx * solvent[i].vx + solvent[i].vy * solvent[i].vy);
    }
}

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
    try{
        var datfile = open(filename, ioMode.cw);
        var myFileWriter = datfile.writer();
        myFileWriter.writeln("step \t KEworm \t KEsol \n");
        for istep in 1..nsteps {
            myFileWriter.writeln(istep,"\t",KE_total[istep]);
        }
        datfile.fsync();
        writeln("energies.dat written");
    } catch e: Error {
        writeln(e);
    }
}

proc gaussRand(mean: real, stddev: real): real {
    var u1 = randStream.getNext();
    var u2 = randStream.getNext();
    var z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2);  // Box-Muller transform
    var gauss = mean + stddev * z0;
    return gauss;   
}