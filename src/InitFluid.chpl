module InitFluid {
    use Amatter3d;
    use Math;
    use Random;
    use BoundaryTypes;
    import Structs;
    use Time; // for stopwatch

    proc estimate_area(bd:Boundary):real {
        var area:real;
        if (bd.t == BD_TYPE.CIRCLE) {
            area = pi*rwall*rwall;
        } else if (bd.t == BD_TYPE.CARDIOID) {
            var ca = 1.5*(rwall/2);
            area = (6 * pi * ca ** 2);
        } else if (bd.t == BD_TYPE.EPICYCLOID) {
            area = 12*pi*((rwall/4)+1)**2; //FIXME to use MC method of area estimation
        } else if (bd.t == BD_TYPE.EPITROCHOID) {
            // area =  (# of points inside / total # of points) * area of bounding box
            var d = 0.8;//0.99; // change smoothing parameter here
            var k = 2.0; // Change the value of k here
            var r = 41.0;// outer circle
            var R = k*r; // inner circle 
            var num_mc = 1000000;
            var num_inside = 0;
            for i in 1..num_mc {
                // generate random point within bounding square
                var x = 2*rwall*randStream.next() - rwall;
                var y = 2*rwall*randStream.next() - rwall;
                // Check if the point is inside the epitrochoid
                var dx = x - hxo2;
                var dy = y - hyo2;
                var theta = atan2(dy,dx);
                var ex = (R + r) * cos(theta) - (d *r) *cos((R + r) / r * theta);
                var ey = (R + r) * sin(theta) - (d *r)* sin((R + r) / r * theta);
                if (sqrt(dx*dx + dy*dy) < sqrt(ex*ex + ey*ey)) {
                    num_inside += 1;
                }
            }
            writeln("Epitrochoid area:");
            writeln(num_inside/num_mc);
            area = (num_inside/num_mc)*(2.0*rwall)**2.0;
            writeln(area);
            halt();
        } else {
            writeln("estimate_area(): Boundary type unknown/unspecified");
        }
        return area;
    }
    proc init_fluid_count():int {
            var numSol:int;
            if (bd.t == BD_TYPE.CIRCLE) {
                var fluid_a = 2*rwall;
                var nSol_row = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
                numSol = nSol_row*nSol_row; // calcuating the number of particles in the box (exceeds rwall)
                var pos : [1..numSol,1..3] real;
                var fluid_px = fluid_a;
                var fluid_mx = 0;
                var fluid_py = fluid_a;
                var fluid_my = 0;
                var spacing = fluid_a/nSol_row;
                var row_length = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
                var row,col:real;
                for i in 1..numSol {
                        row = i % row_length;
                        col = ((i - row)/row_length)+1;

                        pos[i,1] = fluid_mx + spacing*row + spacing; // x
                        pos[i,2] = fluid_my + spacing*col; // y
                        pos[i,3] = 0.0;
                }
                // mark all solvent particles that are out of bounds
                var rm_count = 0;
                for i in 1..numSol {
                    var dx = pos[i,1] - hxo2;
                    var dy = pos[i,2] - hyo2;
                    var r = sqrt(dx*dx + dy*dy);
                    if r > (0.95*rwall) {
                        pos[i,3] = -1.0;
                        rm_count +=1;
                    }
                }
                // recounting
                writeln(rm_count);
                numSol -= rm_count;
            } else if (bd.t == BD_TYPE.CARDIOID) {
                var fluid_a = 2*rwall;
                var nSol_row = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
                numSol = nSol_row*nSol_row; // calcuating the number of particles in the box (exceeds rwall)
                var pos : [1..numSol,1..3] real;
                var fluid_px = fluid_a;
                var fluid_mx = 0;
                var fluid_py = fluid_a;
                var fluid_my = 0;
                var spacing = fluid_a/nSol_row;
                var row_length = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
                var row,col:real;
                for i in 1..numSol {
                        row = i % row_length;
                        col = ((i - row)/row_length)+1;

                        pos[i,1] = fluid_mx + spacing*col; // x
                        pos[i,2] = fluid_my + spacing*row; // y
                        pos[i,3] = 0.0;
                }
                // mark all solvent particles that are out of bounds
                var rm_count = 0;
                var ca = 1.5*(rwall/2);
                for i in 1..numSol {
                    var dx = pos[i,1] - hxo2;
                    var dy = pos[i,2] - hyo2;
                    var r = sqrt(dx*dx + dy*dy);
                    var cardioidRadius = 1 - cos(atan2(dy, dx));
                    if r > (0.95*ca*cardioidRadius) {
                        pos[i,3] = -1.0;
                        rm_count +=1;
                    }
                }
                // recounting
                //writeln(rm_count);
                numSol -= rm_count;
                writeln("fluid density",numSol/(6 * pi * ca ** 2));
            } else if (bd.t == BD_TYPE.EPICYCLOID) {
                // var fluid_x = (max reduce bound.x) - (min reduce bound.x);
                // var fluid_y = (max reduce bound.y) - (min reduce bound.y);
                // var numSol = ceil(fluid_rho * fluid_x*fluid_y):int;
                // var pos : [1..numSol,1..3] real;
                // var fluid_px = (max reduce bound.x);
                // var fluid_mx = (min reduce bound.x);
                // var fluid_py = (max reduce bound.y);
                // var fluid_my = (min reduce bound.y);
                // var nSol_row = floor(sqrt(floor(fluid_rho*fluid_x*fluid_y)));
                // var spacing = fluid_x/nSol_row;
                // var row,col:real;
                // for i in 1..numSol {
                //     row = i % nSol_row;
                //     col = ((i - row)/nSol_row) + 1;
                //     pos[i,1] = fluid_mx + spacing*col;
                //     pos[i,2] = fluid_my - spacing*row;
                //     pos[i,3] = 0.0;
                // }
                // var rm_count = 0;
                //     for i in 1..numSol {
                //         var dx = pos[i,1] - hxo2;
                //         var dy = pos[i,2] - hyo2;
                //         var r = sqrt(dx*dx + dy*dy);
                //         var cycloid_radius = (0.5*(rwall/4 + 1)**2)*(5-3*cos(2*atan2(dy,dx)));
                //         if r > cycloid_radius {
                //             pos[i,3] = -1.0;
                //             rm_count += 1;
                //         }
                //     }
                //     numSol -= rm_count;
                //     writeln("fluid_density ",numSol/(12*pi*(rwall/4 + 1)**2));
                } else if (bd.t == BD_TYPE.EPITROCHOID) {
                // var fluid_x = (max reduce bound.x) - (min reduce bound.x);
                // var fluid_y = (max reduce bound.y) - (min reduce bound.y);
                // var numSol = ceil(fluid_rho * fluid_x*fluid_y):int;
                // var pos : [1..numSol,1..3] real;
                // var fluid_px = (max reduce bound.x);
                // var fluid_mx = (min reduce bound.x);
                // var fluid_py = (max reduce bound.y);
                // var fluid_my = (min reduce bound.y);
                // var nSol_row = floor(sqrt(floor(fluid_rho*fluid_x*fluid_y)));
                // var spacing = fluid_x/nSol_row;
                // var row,col:real;
                // for i in 1..numSol {
                //     row = i % nSol_row;
                //     col = ((i - row)/nSol_row) + 1;
                //     pos[i,1] = fluid_mx + spacing*col;
                //     pos[i,2] = fluid_my - spacing*row;
                //     pos[i,3] = 0.0;
                // }
                // var rm_count = 0;
                //     for i in 1..numSol {
                //         var dx = pos[i,1] - hxo2;
                //         var dy = pos[i,2] - hyo2;
                //         var r = sqrt(dx*dx + dy*dy);
                //         var cycloid_radius = (0.5*(rwall/4 + 1)**2)*(5-3*cos(2*atan2(dy,dx)));
                //         if r > cycloid_radius {
                //             pos[i,3] = -1.0;
                //             rm_count += 1;
                //         }
                //     }
                //     numSol -= rm_count;
                //     writeln("fluid_density ",numSol/(12*pi*(rwall/4 + 1)**2));
            } else {
                halt("no other boundaries supported yet");
            }
    return numSol;
    }
    proc init_fluid(ref solvent: [] Structs.Particle,ref numSol: int) {
    var random_placement = false;
    // place particles in a square lattice centered in the boundary
    if (bd.t == BD_TYPE.CIRCLE) {
        // circular boundary
        var fluid_a = 2*rwall;
        var nSol_row = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
        var numSol_tmp = nSol_row*nSol_row; // calcuating the number of particles in the box (exceeds rwall)
        var pos : [1..numSol_tmp,1..3] real;
        var fluid_px = fluid_a;
        var fluid_mx = 0;
        var fluid_py = fluid_a;
        var fluid_my = 0;
        var spacing = fluid_a/nSol_row;
        var row_length = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
        var row,col:real;
        for i in 1..numSol_tmp {
                row = i % row_length;
                col = ((i - row)/row_length)+1;

                pos[i,1] = fluid_mx + spacing*row + spacing; // x
                pos[i,2] = fluid_my + spacing*col; // y
                pos[i,3] = 0.0;
        }
        // mark all solvent particles that are out of bounds
        var count = 1;
        for i in 1..numSol_tmp {
            var dx = pos[i,1] - hxo2;
            var dy = pos[i,2] - hyo2;
            var r = sqrt(dx*dx + dy*dy);
            if r < (0.95*rwall) {
                solvent[count].x = pos[i,1];
                solvent[count].y = pos[i,2];
                solvent[count].z = 0.0;
                solvent[count].vx = 0.0;
                solvent[count].vy = 0.0;
                solvent[count].vz = 0.0;
                solvent[count].fx = 0.0;
                solvent[count].fy = 0.0;
                solvent[count].fz = 0.0;
                solvent[count].ptype = 2;
                solvent[count].m = 1.0;
                count +=1;
            }
        }
    } else if (bd.t == BD_TYPE.CARDIOID) {
        // cardioid boundary
        // circular boundary
        var fluid_a = 2*rwall;
        var nSol_row = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
        var numSol_tmp = nSol_row*nSol_row; // calcuating the number of particles in the box (exceeds rwall)
        var pos : [1..numSol_tmp,1..3] real;
        var fluid_px = fluid_a;
        var fluid_mx = 0;
        var fluid_py = fluid_a;
        var fluid_my = 0;
        var spacing = fluid_a/nSol_row;
        write_log(logfile,"spacing "+spacing:string);
        var row_length = floor(sqrt(floor(fluid_rho*fluid_a*fluid_a))):int;
        var row,col:real;
        for i in 1..numSol_tmp {
                row = i % row_length;
                col = ((i - row)/row_length)+1;

                pos[i,1] = fluid_mx + spacing*col; // x
                pos[i,2] = fluid_my + spacing*row; // y
                pos[i,3] = 0.0;
        }
        // mark all solvent particles that are out of bounds
        var count = 1;
        var ca = 1.5*(rwall/2);
        for i in 1..numSol_tmp {
            var dx = pos[i,1] - hxo2;
            var dy = pos[i,2] - hyo2;
            var r = sqrt(dx*dx + dy*dy);
            var cardioidRadius = 1 - cos(atan2(dy, dx));
            if r < (0.95*ca*cardioidRadius) {
                solvent[count].x = pos[i,1]+0.5*ca;
                solvent[count].y = pos[i,2];
                solvent[count].z = 0.0;
                solvent[count].vx = 0.0;
                solvent[count].vy = 0.0;
                solvent[count].vz = 0.0;
                solvent[count].fx = 0.0;
                solvent[count].fy = 0.0;
                solvent[count].fz = 0.0;
                solvent[count].ptype = 2;
                solvent[count].m = 1.0;
                count +=1;
            }
        }
        writeln((max reduce solvent.x)," ",(min reduce solvent.y));
        writeln((max reduce solvent.x)," ",(min reduce solvent.y));
    } else if (boundary == 3) {
        writeln("haven't put in channel fluid");
        halt();
    } else {
        writeln("invalid boundary");
        halt();
    }
    // init macro variable arrays
    for i in 1..numSol {
        KEsol[i] = 0.0;
    }
    if (bd.t == BD_TYPE.CIRCLE){
        writeln("fluid density=",numSol/(pi*rwall**2));
    } else if (bd.t == BD_TYPE.CARDIOID){
        writeln("fluid density=",numSol/(6*pi*(1.5*(rwall/2))**2));
    }
    writeln("fluid init done");
    return solvent;
    }

    proc init_fluid_rsa1(ref solvent: [] Structs.Particle, ref numSol: int){
        var isInCardioid,tooClose:bool;
        var minDist = 1.5;
        var numPlaced = 0; //Counter for the number of particles placed
        var init_timer:stopwatch, pl_timer:stopwatch;
        init_timer.start();
        while numPlaced <= numSol {
            // generate random point within bounding square
            var x = 2*rwall*randStream.next();
            var y = 2*rwall*randStream.next();
            // Check if the point is inside the cardioid
            var dx = x - hxo2;
            var dy = y - hyo2;
            var r = sqrt(dx*dx + dy*dy);
            if r <= 0.95*rwall {
                isInCardioid = true;
            } else {
                isInCardioid = false;
            }
            // check if the point is too close to other particles
            tooClose = false;
            for i in 1..numPlaced {
                var dist = sqrt((solvent[i].x - x)**2 + (solvent[i].y - y)**2);
                if dist < minDist {
                    tooClose = true;
                    break;
                }
            }
            if isInCardioid && !tooClose {
                writeln(numPlaced,"/",numSol," ",numPlaced/(pi*rwall**2)," ",init_timer.elapsed());
                numPlaced += 1;
                solvent[numPlaced].x = x;
                solvent[numPlaced].y = y;
                solvent[numPlaced].z = 0.0;
                solvent[numPlaced].vx = 0.0;
                solvent[numPlaced].vy = 0.0;
                solvent[numPlaced].vz = 0.0;
                solvent[numPlaced].fx = 0.0;
                solvent[numPlaced].fy = 0.0;
                solvent[numPlaced].fz = 0.0;
                solvent[numPlaced].ptype = 2;
                solvent[numPlaced].m = 1.0;
            }
            if init_timer.elapsed() > (720*60) {
                writeln("FLUID INIT FAILED, CHECK minDist AND fluid_rho VALUES");
                halt();
            }
        }
        init_timer.stop();
        writeln(init_timer.elapsed()," s for solvent init");
        return solvent;
    }

    proc init_fluid_rsa2(ref solvent: [] Structs.Particle, ref numSol: int){
    var isInCardioid,tooClose:bool;
    var minDist = 1.5;
    var numPlaced = 0; //Counter for the number of particles placed
    var ca = 1.5*(rwall/2);
    var init_timer:stopwatch, pl_timer:stopwatch;
    init_timer.start();
    while numPlaced <= numSol {
        // generate random point within bounding square
        var x = 2*rwall*randStream.next();
        var y = 2*rwall*randStream.next();
        // Check if the point is inside the cardioid
        var dx = x - hxo2 - 0.95*ca;
        var dy = y - hyo2;
        var r = sqrt(dx*dx + dy*dy);
        var rCardioid = ca * (1 - cos(atan2(dy, dx)));
        if r <= 0.95*rCardioid {
            isInCardioid = true;
        } else {
            isInCardioid = false;
        }
        // check if the point is too close to other particles
        tooClose = false;
        for i in 1..numPlaced {
            var dist = sqrt((solvent[i].x - x)**2 + (solvent[i].y - y)**2);
            if dist < minDist {
                tooClose = true;
                break;
            }
        }
        if isInCardioid && !tooClose {
            writeln(numPlaced,"/",numSol," ",numPlaced/(6 * pi * ca ** 2)," ",init_timer.elapsed());
            numPlaced += 1;
            solvent[numPlaced].x = x;
            solvent[numPlaced].y = y;
            solvent[numPlaced].z = 0.0;
            solvent[numPlaced].vx = 0.0;
            solvent[numPlaced].vy = 0.0;
            solvent[numPlaced].vz = 0.0;
            solvent[numPlaced].fx = 0.0;
            solvent[numPlaced].fy = 0.0;
            solvent[numPlaced].fz = 0.0;
            solvent[numPlaced].ptype = 2;
            solvent[numPlaced].m = 1.0;
        }
        if init_timer.elapsed() > (720*60) {
            writeln("FLUID INIT FAILED, CHECK minDist AND fluid_rho VALUES");
            halt();
        }
    }
    init_timer.stop();
    writeln(init_timer.elapsed()," s for solvent init");
    return solvent;
    }

    proc init_fluid_rsa3(ref solvent: [] Structs.Particle, ref numSol: int) {
        var isInCycloid,tooClose:bool;
        var minDist = 1.5;
        var numPlaced = 0; //Counter for the number of particles placed
        var ca = (rwall/4)+1;
        var init_timer:stopwatch, pl_timer:stopwatch;
        init_timer.start();
        while numPlaced <= numSol {
            // generate random point within bounding square
            var x = 2*rwall*randStream.next();
            var y = 2*rwall*randStream.next();
            // Check if the point is inside the epitrochoid
            var dx = x - hxo2 - 0.95*ca;
            var dy = y - hyo2;
            var rCylc = 1.9*ca*(1+2*sin(t)/2); // somehow a circle? 2*ca?
            //var rCylc = (0.5*(rwall/4 + 1)**2)*(5-3*cos(2*atan2(x,x)));//sqrt((cylc_x-x)**2 + (cylc_y-y)**2);
            var r = sqrt(dx*dx + dy*dy);
            if r <= 0.95*rCylc {
                isInCycloid = true;
            } else {
                isInCycloid = false;
            }
            // check if the point is too close to other particles
            tooClose = false;
            for i in 1..numPlaced {
                var dist = sqrt((solvent[i].x - x)**2 + (solvent[i].y - y)**2);
                if dist < minDist {
                    tooClose = true;
                    break;
                }
            }
            if isInCycloid && !tooClose {
                writeln(numPlaced,"/",numSol," ",numPlaced/(12 * pi * ca ** 2)," ",init_timer.elapsed());
                numPlaced += 1;
                if numPlaced > numSol {
                    break;
                }
                solvent[numPlaced].x = x;
                solvent[numPlaced].y = y;
                solvent[numPlaced].z = 0.0;
                solvent[numPlaced].vx = 0.0;
                solvent[numPlaced].vy = 0.0;
                solvent[numPlaced].vz = 0.0;
                solvent[numPlaced].fx = 0.0;
                solvent[numPlaced].fy = 0.0;
                solvent[numPlaced].fz = 0.0;
                solvent[numPlaced].ptype = 2;
                solvent[numPlaced].m = 1.0;
            }
            if init_timer.elapsed() > (2*720*60) {
                writeln("FLUID INIT FAILED, CHECK minDist AND fluid_rho VALUES");
                halt();
            }
        }
        init_timer.stop();
        writeln(init_timer.elapsed()," s for solvent init");
        return solvent;
    }

    proc init_fluid_rsa4(ref solvent: [] Structs.Particle, ref numSol: int) {
        var isInCycloid,tooClose:bool;
        var minDist = 1.5;
        var numPlaced = 0; //Counter for the number of particles placed
        // epitrochoid
        var d = 0.8;//0.99; // change smoothing parameter here
        var k = 2.0; // Change the value of k here
        var r = 41.0;// outer circle
        var R = k*r; // inner circle
        var xoffset1 = hxo2; //1.25
        var yoffset1 = 0.675*hyo2; // lower lobe
        var yoffset2 = 1.325*hyo2; // upper lobe

        var init_timer:stopwatch, pl_timer:stopwatch;
        init_timer.start();
        while numPlaced <= numSol {
            // generate random point within bounding square
            var x = 2*rwall*randStream.next();
            var y = 2*rwall*randStream.next();
            // Check if the point is inside the epitrochoid
            var dx = x - hxo2;
            var dy = y - hyo2;
            var epi_pt = epitrochoid_coords(atan2(dy,dx),R,r,d);
            var rCylc = sqrt(epi_pt[0]*epi_pt[0] + epi_pt[1]*epi_pt[1]);
            //var rCylc = 1.9*ca*(1+2*sin(t)/2); // somehow a circle? 2*ca?
            //var rCylc = (0.5*(rwall/4 + 1)**2)*(5-3*cos(2*atan2(x,x)));//sqrt((cylc_x-x)**2 + (cylc_y-y)**2);
            var r_pt = sqrt(dx*dx + dy*dy);
            if r_pt <= 0.95*rCylc {
                isInCycloid = true;
            } else {
                isInCycloid = false;
            }
            // check if the point is too close to other particles
            tooClose = false;
            for i in 1..numPlaced {
                var dist = sqrt((solvent[i].x - x)**2 + (solvent[i].y - y)**2);
                if dist < minDist {
                    tooClose = true;
                    break;
                }
            }
            if isInCycloid && !tooClose {
                writeln(numPlaced,"/",numSol," ",init_timer.elapsed());
                numPlaced += 1;
                if numPlaced > numSol {
                    break;
                }
                solvent[numPlaced].x = x;
                solvent[numPlaced].y = y;
                solvent[numPlaced].z = 0.0;
                solvent[numPlaced].vx = 0.0;
                solvent[numPlaced].vy = 0.0;
                solvent[numPlaced].vz = 0.0;
                solvent[numPlaced].fx = 0.0;
                solvent[numPlaced].fy = 0.0;
                solvent[numPlaced].fz = 0.0;
                solvent[numPlaced].ptype = 2;
                solvent[numPlaced].m = 1.0;
            }
            if init_timer.elapsed() > (2*720*60) {
                writeln("FLUID INIT FAILED, CHECK minDist AND fluid_rho VALUES");
                halt();
            }
        }
        init_timer.stop();
        writeln(init_timer.elapsed()," s for solvent init");
        return solvent;
    }
    

}