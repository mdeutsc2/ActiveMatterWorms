module Structs {
    use List;
    use IO;
    use IO.FormattedIO;

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
        // proc info() {
        //     var typestring: string;
        //     if (this.ptype == 1) {
        //         typestring = "active";
        //     } else if (this.ptype == 2) {
        //         typestring = "solvent";
        //     } else if (this.ptype == 3) {
        //         typestring = "boundary";
        //     } else {
        //         typestring = "unassigned";
        //     }
        //     var s:string = "id# %i \t type: %s \t mass: %s \n pos: %r \t %r \t %r \n vel: %r \t %r \t %r \n force: %r \t %r \t %r".format(this.id,typestring,this.m,this.x,this.y,this.z,this.vx,this.vy,this.vz,this.fx,this.fy,this.fz);
        //     return s;
        // }
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
        //var id: (int,int); //id of each bin
        var id: int; //id of each bin
        var atoms: list(int); // list of particle id's in each bin
        var neighbors: [1..4] int; // indices of each bin's neighboring bin
        var ncount: int; // count of number of particles in each bin
        var x: (real,real); // precalculate the max_x and min_x values for the space the box occupies
        var y: (real,real); // this is for easier neighbor list creation

        proc init() { // record initializer
            this.id = bin_init_counter;
            bin_init_counter += 1;
            this.ncount = 0;
            this.x = (0.0,0.0);
            this.y = (0.0,0.0);
            for i in 1..4 {
                this.neighbors[i] = -1;
            }
        }
    }
}