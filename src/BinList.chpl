module BinList {
    use Amatter3d;
    use Output;
    use Helper;
// CELL LIST FUNCTIONS
proc update_cells(istep:int) {
    // ncount array set to zero
    // particle list set to empty

    // resetting all other bins lists for worms and solvent
    for binid in binSpace{
        if bins[binid].ncount > 0 {
            bins[binid].ncount = 0; // number of boundary atoms should not change
            bins[binid].scount = 0;
            bins[binid].wcount = 0;
            bins[binid].bcount = 0;
            bins[binid].atoms.clear();
            bins[binid].types.clear();
        }
    }
    var ibin,jbin,binid :int;
    for i in 1..numSol { // populating solvent particles into bins
        ibin=ceil(solvent[i].x/rcut):int;
        jbin=ceil(solvent[i].y/rcut):int;
        binid=(jbin-1)*numBins+ibin;
        if (binid > numBins*numBins) {
            write_log(logfile,solvent[i].info());
            sudden_halt(istep);
        } else if (binid < 1) {
            write_log(logfile,solvent[i].info());
            sudden_halt(istep);
        }
        //append i to particle_list for binid
        bins[binid].scount += 1;
        bins[binid].ncount += 1;
        bins[binid].atoms.pushBack(solvent[i].id);
        bins[binid].types.pushBack(solvent[i].ptype);
        //if (debug) {writeln(i,"\t",solvent[i].x,"\t",solvent[i].y,"\t",binposx,"\t",binposy);}
    }
    var wormid : int;
    for iw in 1..nworms { // populating worm particles into bins
        for ip in 1..np {
            ibin=ceil(worms[iw,ip].x/rcut):int;
            jbin=ceil(worms[iw,ip].y/rcut):int;
            binid=(jbin-1)*numBins+ibin;
            wormid = (iw-1)*np+ip;
            if (binid > numBins*numBins) {
                write_log(logfile,worms[iw,ip].info());
                sudden_halt(istep);
            } else if (binid < 1) {
                write_log(logfile,worms[iw,ip].info());
                sudden_halt(istep);
            }
            //append i to particle_list for binid
            bins[binid].wcount += 1;
            bins[binid].ncount += 1;
            bins[binid].atoms.pushBack(wormid);
            bins[binid].types.pushBack(worms[iw,ip].ptype);
            //if (debug) {writeln(i,"\t",solvent[i].x,"\t",solvent[i].y,"\t",binposx,"\t",binposy);}
        }
    }
    for ib in 1..numPoints {
        if (bound[ib].x == 0.0) && (bound[ib].y == 0.0) {
            ibin=ceil((bound[ib].x+0.01)/rcut):int;
            jbin=ceil((bound[ib].y+0.01)/rcut):int;
            binid=(jbin-1)*numBins+ibin;
        } else {
            ibin=ceil(bound[ib].x/rcut):int;
            jbin=ceil(bound[ib].y/rcut):int;
            binid=(jbin-1)*numBins+ibin;
        }
        //append i to particle_list for binid
        bins[binid].bcount += 1;
        bins[binid].ncount += 1;
        bins[binid].atoms.pushBack(bound[ib].id);
        bins[binid].types.pushBack(bound[ib].ptype);
    }
}

proc init_bins() {
    // writeln((4/numBins):int +1); //row
    // writeln(4%numBins); // col
    var binid,ibinnab,jbinnab,binidnbor:int;
    for ibin in 1..numBins {
        for jbin in 1..numBins {
            binid = (jbin-1)*numBins+ibin;
            //  4   3   2
            //     (*)   1
            //
            /*  1 = i+1,j
                2 = i+1,j+1
                3 = i,j+1
                4 = i-1,j+1
            */
            ibinnab=ibin+1;
            jbinnab=jbin;
            if (ibinnab <= numBins) {
                binidnbor=(jbinnab-1)*numBins+ibinnab;
                bins[binid].neighbors[1] = binidnbor;
            } else {
                bins[binid].neighbors[1] = -1;
            }

            ibinnab=ibin+1;
            jbinnab=jbin+1;
            if (ibinnab <= numBins) && (jbinnab <= numBins) {
                binidnbor=(jbinnab-1)*numBins+ibinnab;
                bins[binid].neighbors[2] = binidnbor;
            } else {
                bins[binid].neighbors[2] = -1;
            }

            ibinnab=ibin;
            jbinnab=jbin+1;
            if (jbinnab <= numBins) {
                binidnbor=(jbinnab-1)*numBins+ibinnab;
                bins[binid].neighbors[3] = binidnbor;
            } else {
                bins[binid].neighbors[3] = -1;
            }

            ibinnab=ibin-1;
            jbinnab=jbin+1;
            if (ibinnab > 0) && (jbinnab <= numBins){
                binidnbor=(jbinnab-1)*numBins+ibinnab;
                bins[binid].neighbors[4] = binidnbor;
            } else {
                bins[binid].neighbors[4] = -1;
            }

        }
    }
        // bins[ibin].x[0] = (ibin[0]-1)*rcut;
        // bins[ibin].x[1] = ibin[0]*rcut;
        // bins[ibin].y[0] = (ibin[1]-1)*rcut;
        // bins[ibin].y[1] = ibin[1]*rcut;
}

proc init_binspace() {
    // making lists of binids for different
    var icount = 0;
     for ibin in 1..numBins by 2 {
         for jbin in 1..numBins {
            icount += 1;
            binSpaceiodd[icount] = (jbin-1)*numBins+ibin;
        }
    }
    icount = 0;
    for ibin in 2..numBins by 2 {
        for jbin in 1..numBins {
            icount +=1;
            binSpaceieven[icount] = (jbin-1)*numBins+ibin;
        }
    }
    icount = 0;
    for ibin in 1..numBins {
        for jbin in 1..numBins by 2 {
            icount += 1;
            binSpacejodd[icount] = (jbin-1)*numBins+ibin;
        }
    }
    icount = 0;
    for ibin in 1..numBins {
        for jbin in 2..numBins by 2 {
            icount += 1;
            binSpacejeven[icount] = (jbin-1)*numBins+ibin;
        }
    }
}
}