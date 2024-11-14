module InitWorms {
    use Amatter3d;
    use Math;
    use BoundaryTypes;
    import Structs;
    use Time; // for stopwatch
    use Output;
    
proc epitrochoid_coords(theta,R,r,d) {
    var x = (R + r) * cos(theta) - (d *r) *cos((R + r) / r * theta);
    var y = (R + r) * sin(theta) - (d *r)* sin((R + r) / r * theta);
    return (x, y);
}

proc arc_length(theta1,theta2,R,r,d) { //epitrochoid arc length
    // Arc length estimation function (using numerical integration)
    var n = 25000; // Number of subdivisions for numerical integration
    var dtheta = (theta2 - theta1) / n;
    var length = 0.0;
    for i in 1..n-1 {
        var theta_left = theta1 + (i-1)*dtheta;
        var theta_right = theta1 + i * dtheta;
        var px_left = 0.0;
        var py_left = 0.0; 
        var px_right = 0.0;
        var py_right = 0.0;
        var p_left = epitrochoid_coords(theta_left,R,r,d);
        var p_right = epitrochoid_coords(theta_right,R,r,d);
        var dx = p_right[0]-p_left[0];
        var dy = p_right[1]-p_left[1];
        length += sqrt(dx*dx + dy*dy);
    }
    return length;
}

proc init_worms() {
    //array for locating neighbor cells
    var rand1:real; // temp fix for randSteam
    var r:real, dth:real, xangle:real;
    var thetanow = 5.0*pi :real; // changes the initial radius of annulus
    var rmin = a*thetanow;
    var density = 0.0;
    if (bd.t == BD_TYPE.CIRCLE) {
        density = nworms*np/(pi*rwall**2); //density for a circle
    } else if (bd.t == BD_TYPE.CARDIOID) {
        density = nworms*np/(6*pi*(1.5*(rwall/2))**2); // density for a cardioid
    }
    write_log(logfile,"nworms\t"+nworms:string);
    write_log(logfile,"np\t"+np:string);
    write_log(logfile,"rwall\t"+rwall:string);
    write_log(logfile,"nxcells\t"+nxcell:string+"\t nycells\t"+nycell:string);
    write_log(logfile,"density\t"+density:string);
    write_log(logfile,"dt\t"+dt:string);
    write_log(logfile,"opt numPoints\t"+(ceil(2*pi*rwall)/rcutsmall):string);
    //setting up the worms
    if (bd.t == BD_TYPE.CIRCLE){
        // circular bound
        var equidistantThetaValues: [1..numPoints] real = 0.0;
        var deltaTheta = 2 * pi / numPoints;

        for i in 1..numPoints {
            equidistantThetaValues[i] = i * deltaTheta;
        }

        for i in 1..numPoints {
            bound[i].x = rwall * cos(equidistantThetaValues[i])+hxo2;
            bound[i].y = rwall * sin(equidistantThetaValues[i])+hyo2;
            bound[i].z = 0.0;
            bound[i].ptype = 3;
        }
        for iw in 1..nworms {
            ireverse[iw] = 0;
            rand1 = randStream.next();
            if (rand1 <= 0.5) {
                ireverse[iw] = 1;
            }
            var worm_z_height = randStream.next()*(L);
            for i in 1..np {
                r = a*thetanow;
                dth = length0/r;
                thetanow += dth;
                worms[iw,i].x = hxo2 + r*cos(thetanow);
                worms[iw,i].y = hyo2 + r*sin(thetanow);
                worms[iw,i].z = worm_z_height;//0.5*L + gaussRand(0.0,0.1);//randStream.next()*(L);
                xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
                //TODO give them an initial velocity going around the circle
                worms[iw,i].ptype = 1;
                worms[iw,i].m = 1.0; // setting mass
                // vx[iw,i] = 0.0; // NOTE: this is all set by initializer
                // vy[iw,i] = 0.0;
                worms[iw,i].vxave = 0.0;
                worms[iw,i].vyave = 0.0;
                worms[iw,i].vzave = 0.0;
                // fx[iw,i] = 0.0;
                // fy[iw,i] = 0.0;
                // fxold[iw,i] = 0.0;
                // fyold[iw,i] = 0.0;
            }
            thetanow += 2.0*dth;
        }
    } else if (bd.t == BD_TYPE.CARDIOID) {
        // cardioid boundary
        var equidistantArcLengths: [1..numPoints] real;
        var thetaValues: [1..numPoints] real;
        var ca = 1.5*(rwall/2);
        var totalArcLength = 8 * ca;
        write_log(logfile,"cardioid width: "+(ca*2):string);
        write_log(logfile,"cardioid cirum: "+totalArcLength:string+"\t"+(totalArcLength/r2cutsmall):string);
        for i in 1..numPoints {
            equidistantArcLengths[i] = totalArcLength * (i - 1) / (numPoints - 1);
            thetaValues[i] = 4 * asin(sqrt(equidistantArcLengths[i] / totalArcLength));
        }

        for i in 1..numPoints {
            bound[i].x = ca * (1 - cos(thetaValues[i])) * cos(thetaValues[i]) + hxo2 + ca;
            bound[i].y = ca * (1 - cos(thetaValues[i])) * sin(thetaValues[i]) + hyo2;
            bound[i].z = 0.0;
            bound[i].ptype = 3;
        }
        //now place worm particles
        // high-density placement 
        var n_levels = 4; // number of levels of filaments
        var theta_now_init = thetanow;
        for ilevel in 1..n_levels {
            thetanow = theta_now_init;
            var wormstart = (nworms/n_levels)*(ilevel-1)+1;
            var wormend = wormstart + (nworms/n_levels)-1;
            writeln(ilevel," ",nworms/n_levels);
            for iw in wormstart..wormend {
                ireverse[iw] = 0;
                rand1 = randStream.next();
                if (rand1 <= 0.5) {
                    ireverse[iw] = 1;
                }
                var worm_z_height = ((0.5*randStream.next())*L)/n_levels + (L/n_levels)*(ilevel-1);
                // if worm_z_height >= L {
                //     worm_z_height= L-0.01;
                // }
                // if worm_z_height <= 0.0 {
                //     worm_z_height = 0.01;
                // }
                //writeln(ilevel*L/n_levels);
                for i in 1..np {
                    r = a*thetanow;
                    dth = (length0)/r;
                    thetanow += dth;
                    worms[iw,i].x = hxo2 + r*cos(thetanow);
                    //worms[iw,i].x = r * (1 - cos(thetanow)) * cos(thetanow) + hxo2 + ca;
                    worms[iw,i].y = hyo2 + r*sin(thetanow);
                    //worms[iw,i].y = r * (1 - cos(thetanow)) * sin(thetanow) + hyo2;
                    worms[iw,i].z = worm_z_height;//0.5*L + gaussRand(0.0,0.1);//randStream.next()*(L);
                    xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
                    //TODO give them an initial velocity going around the circle
                    worms[iw,i].ptype = 1;
                    worms[iw,i].m = 1.0; // setting mass
                    // vx[iw,i] = 0.0;
                    // vy[iw,i] = 0.0;
                    worms[iw,i].vxave = 0.0;
                    worms[iw,i].vyave = 0.0;
                    worms[iw,i].vzave = 0.0;
                    // fx[iw,i] = 0.0;
                    // fy[iw,i] = 0.0;
                    // fxold[iw,i] = 0.0;
                    // fyold[iw,i] = 0.0;
                }
                thetanow += 2.0*dth;
            }
        }
    } else if (bd.t == BD_TYPE.EPICYCLOID) {
        // epicycloid boundary
        var k = 2; // Change the value of k here
        var xoffset1 = 1.25*hxo2;
        var yoffset1 = 0.55*hyo2;
        var yoffset2 = 1.45*hyo2;
        writeln("EPICYCLOID1 BOUNDARY k = ",k);
        var equidistantArcLengths: [1..numPoints] real;
        var thetaValues: [1..numPoints] real;
        var ca = (rwall/4)-1;
        writeln("ca=",ca," numpoints=",numPoints);
        var totalArcLength = 4*(k+1)*ca;

        for i in 1..numPoints {
            // note that 1.5*i is necessary for even-ish spacing
            equidistantArcLengths[i] = totalArcLength * (1.5*i - 1) / (numPoints - 1);
            //equidistantArcLengths[i] = totalArcLength * (i) / (numPoints);
            thetaValues[i] = 4 * asin(sqrt(equidistantArcLengths[i] / totalArcLength));
        }
        for i in 1..numPoints {
            bound[i].x = ca * (k + 1) * cos(thetaValues[i]) - ca * cos((k + 1) * thetaValues[i]) + hxo2 + ca;
            bound[i].y = ca * (k + 1) * sin(thetaValues[i]) - ca * sin((k + 1) * thetaValues[i]) + hyo2;
            bound[i].z = 0.0;
            bound[i].ptype = 3;
        }
        var nfailed = 0;
        //  LOWER LOBE high-density placement 
        var n_levels = 4; // number of levels of filaments
        var theta_now_init = thetanow;
        for ilevel in 1..n_levels {
            thetanow = theta_now_init;
            var wormstart = (nworms/n_levels)*(ilevel-1)+1;
            var wormend = wormstart + (nworms/n_levels)-1;
            writeln(ilevel," ",nworms/n_levels);
            for iw in wormstart..wormend {
                ireverse[iw] = 0;
                rand1 = randStream.next();
                if (rand1 <= 0.5) {
                    ireverse[iw] = 1;
                }
                var worm_z_height = ((0.5*randStream.next())*L)/n_levels + (L/n_levels)*(ilevel-1);
                var thetanow_temp = thetanow;
                var place_failed = false;
                for i in 1..np {
                    r = a*thetanow_temp;
                    thetanow_temp += dth;
                    var xtmp = xoffset1 + r*cos(thetanow);
                    var ytmp = yoffset1 + r*sin(thetanow);
                    if (ytmp >= hyo2) {
                        place_failed = true;
                        //break;
                    }
                }
                if place_failed {
                    writeln("Place failed! ",iw);
                    for i in 1..np {
                        r = a*thetanow;
                        dth = (length0)/r;
                        thetanow +=dth;
                    }
                    nfailed += 1;
                    break;
                } else {
                    for i in 1..np {
                        r = a*thetanow;
                        dth = (length0)/r;
                        thetanow += dth;
                        worms[iw,i].x = xoffset1 + r*cos(thetanow);
                        worms[iw,i].y = yoffset1 + r*sin(thetanow);
                        worms[iw,i].z = worm_z_height;//0.5*L + gaussRand(0.0,0.1);//randStream.next()*(L);
                        xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
                        //TODO give them an initial velocity going around the circle
                        worms[iw,i].ptype = 1;
                        worms[iw,i].m = 1.0; // setting mass
                        worms[iw,i].vxave = 0.0;
                        worms[iw,i].vyave = 0.0;
                        worms[iw,i].vzave = 0.0;
                    }
                }
                thetanow += 2.0*dth;
            }
        }
        //  HIGHER LOBE high-density placement 
        // for ilevel in 1..n_levels {
        //     thetanow = theta_now_init;
        //     var wormstart = (nworms/n_levels)*(ilevel-1)+1;
        //     var wormend = wormstart + (nworms/n_levels)-1;
        //     writeln(ilevel," ",nworms/n_levels);
        //     for iw in wormstart..wormend {
        //         ireverse[iw] = 0;
        //         rand1 = randStream.next();
        //         if (rand1 <= 0.5) {
        //             ireverse[iw] = 1;
        //         }
        //         var worm_z_height = ((0.5*randStream.next())*L)/n_levels + (L/n_levels)*(ilevel-1);
        //         var thetanow_temp = thetanow;
        //         var place_failed = false;
        //         for i in 1..np {
        //             r = a*thetanow_temp;
        //             thetanow_temp += dth;
        //             var xtmp = xoffset1 + r*cos(thetanow);
        //             var ytmp = yoffset2 + r*sin(thetanow);
        //             if (ytmp <= hyo2) {
        //                 place_failed = true;
        //                 break;
        //             }
        //         }
        //         if place_failed {
        //             writeln("Place failed!");
        //             for i in 1..np {
        //                 r = a*thetanow;
        //                 dth = (length0)/r;
        //                 thetanow +=dth;
        //             }
        //             nfailed += 1;
        //             break;
        //         } else {
        //             for i in 1..np {
        //                 r = a*thetanow;
        //                 dth = (length0)/r;
        //                 thetanow += dth;
        //                 worms[iw,i].x = xoffset1 + r*cos(thetanow);
        //                 worms[iw,i].y = yoffset2 + r*sin(thetanow);
        //                 worms[iw,i].z = worm_z_height;//0.5*L + gaussRand(0.0,0.1);//randStream.next()*(L);
        //                 xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
        //                 //TODO give them an initial velocity going around the circle
        //                 worms[iw,i].ptype = 1;
        //                 worms[iw,i].m = 1.0; // setting mass
        //                 worms[iw,i].vxave = 0.0;
        //                 worms[iw,i].vyave = 0.0;
        //                 worms[iw,i].vzave = 0.0;
        //             }
        //         }
        //         thetanow += 2.0*dth;
        //     }
        // }
        writeln("Nworms failed ",nfailed);
        //halt();
    } else if (bd.t == BD_TYPE.EPITROCHOID) {
        // epitrochoid
        var d = 0.8;//0.99; // change smoothing parameter here
        var k = 2.0; // Change the value of k here
        var r = 41.0;// outer circle
        var R = k*r; // inner circle 
        var xoffset1 = hxo2; //1.25
        var yoffset1 = 0.55*hyo2;//0.675*hyo2; // lower lobe
        var yoffset2 = 1.45*hyo2;//1.325*hyo2; // upper lobe
        writeln("outer circle ideal size: ",floor((2.0*rwall)/8.0)," ",r);
        writeln("epitrochoid arc length ",arc_length(0.0,2.0*pi,R,r,d));
        writeln("arc length/spacing ", arc_length(0.0,2.0*pi,R,r,d)/bdSpacing);
        writeln("numPoints ",numPoints);
        if (numPoints > arc_length(0.0,2.0*pi,R,r,d)/bdSpacing) {
            writeln("numPoints > given spacing!!!");
        }
        // OLD-METHOD
        // var equidistantArcLengths: [1..numPoints] real;
        // var thetaValues: [1..numPoints] real;
        // var totalArcLength = arc_length(0.0,2.0*pi,R,r,d);
        // write_log(logfile,"epitrochoid width: "+(r+R):string);
        // write_log(logfile,"epitrochoid cirum: "+totalArcLength:string+"\t"+(totalArcLength/r2cutsmall):string);
        // for i in 1..numPoints {
        //     equidistantArcLengths[i] = totalArcLength * (i - 1) / (numPoints - 1);
        //     thetaValues[i] = 4 * asin(sqrt(equidistantArcLengths[i] / totalArcLength));
        // }
        // var max_y = 0.0;
        // var min_y = 0.0;
        // for i in 1..numPoints {
        //     var x = (R + r) * cos(thetaValues[i]) - (d*r) * cos((R + r) / r * thetaValues[i]) + hxo2;
        //     var y = (R + r) * sin(thetaValues[i]) - (d*r) * sin((R + r) / r * thetaValues[i]) + hyo2;
        //     //points[i] = (x, y);
        //     if (y > max_y) {
        //         max_y = y;
        //     }
        //     if (y < min_y) {
        //         min_y = y;
        //     }
        //     bound[i].x = x;
        //     bound[i].y = y;
        //     bound[i].z = 0.0;
        //     bound[i].ptype = 3;
        // }
        // writeln("bound extents: ",max_y," ",min_y);
        var total_length = 0.0;
        var itheta = 0.0;
        var x = (R + r) * cos(itheta) - (d*r) * cos((R + r) / r * itheta) + hxo2;
        var y = (R + r) * sin(itheta) - (d*r) * sin((R + r) / r * itheta) + hyo2;
        var iptc = 1;
        bound[iptc].x = x;
        bound[iptc].y = y;
        bound[iptc].z = 0.0;
        bound[iptc].ptype = 3;
        while total_length <= arc_length(0.0,2.0*pi,R,r,d) {
            var step_size = 0.001;
            var next_theta = itheta + step_size;
            var current_length = arc_length(0.0,next_theta,R,r,d);
            if (current_length - total_length >= bdSpacing) {
                bound[iptc].x = (R + r) * cos(next_theta) - (d*r) * cos((R + r) / r * next_theta) + hxo2;
                bound[iptc].y = (R + r) * sin(next_theta) - (d*r) * sin((R + r) / r * next_theta) + hyo2;
                bound[iptc].z = 0.0;
                bound[iptc].ptype = 3;
                iptc += 1;
                total_length = current_length;
                itheta = next_theta;
            } else {
                itheta = next_theta;
            }
        }

        //  LOWER LOBE high-density placement 
        var n_levels = 8; // number of levels of filaments
        //n_levels = 2*n_levels; // have to alternate number of filaments
        var theta_now_init = thetanow;
        for ilevel in 1..n_levels {
            thetanow = theta_now_init;
            var wormstart = (nworms/n_levels)*(ilevel-1)+1;
            var wormend = wormstart + (nworms/n_levels)-1;
            writeln(ilevel," ",nworms/n_levels);
            for iw in wormstart..wormend {
                ireverse[iw] = 0;
                rand1 = randStream.next();
                if (rand1 <= 0.5) {
                    ireverse[iw] = 1;
                }
                var worm_z_height = ((0.5*randStream.next())*L)/n_levels + (L/n_levels)*(ilevel-1);
                for i in 1..np {
                    r = a*thetanow;
                    dth = (length0)/r;
                    thetanow += dth;
                    worms[iw,i].x = xoffset1 + r*cos(thetanow);
                    if ilevel%(k:int) == 0 {
                        worms[iw,i].y = yoffset1 + r*sin(thetanow);
                        worms[iw,i].z = worm_z_height;
                    } else {
                        worms[iw,i].y = yoffset2 + r*sin(thetanow);
                        worms[iw,i].z = worm_z_height;
                    }
                    //worms[iw,i].z = worm_z_height;//0.5*L + gaussRand(0.0,0.1);//randStream.next()*(L);
                    xangle = atan2(worms[iw,i].y - hyo2, worms[iw,i].x - hxo2);
                    //TODO give them an initial velocity going around the circle
                    worms[iw,i].ptype = 1;
                    worms[iw,i].m = 1.0; // setting mass
                    worms[iw,i].vxave = 0.0;
                    worms[iw,i].vyave = 0.0;
                    worms[iw,i].vzave = 0.0;
                }
                thetanow += 2.0*dth;
            }
        }
    }
    //reverse some of the worms and give them crazy colors
    for iw in 1..nworms {
        if (ireverse[iw] == 1) {
            for i in 1..np {
                savex[i] = worms[iw,i].x;
                savey[i] = worms[iw,i].y;
                savez[i] = worms[iw,i].z;
            }
            for ip in 1..np {
                worms[iw,ip].x = savex[np+1-ip];
                worms[iw,ip].y = savey[np+1-ip];
                worms[iw,ip].z = savez[np+1-ip];
            }
        }
    }
    //init macro variable arrays
    for i in 1..nworms*np{
        KEworm[i] = 0.0;
    }
    write_log(logfile,"init done");
}
}