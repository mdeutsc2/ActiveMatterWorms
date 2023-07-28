
// configuration
config const np = 80,
             nworms = 1125,
             nsteps = 250,
             fdogic = 0.06,
             fdogicwall = 0.0,
             fdep = 1.0,
             fdepwall = 0.0,
             diss = 0.08,
             dt = 0.02,
             kspring = 57.146436,
             kbend = 40.0,
             length0 = 0.8,
             rcut = 2.5,
             save_interval = 10,
             fluid_cpl = false;

// variables
const r2cut = rcut*rcut,
      rcutsmall = 2.0**(1/6),
      rwall = 125.0*rcutsmall*sqrt(2.0),
      pi = 3.1415926,
      twopi = 2*pi,
      pio4 = pi*0.25,
      density = nworms*np/(pi*rwall**2), //density for a circle
      hx = 2.0*rwall + 1.0,
      hy = hx,
      hyo2 = hy/2,
      hxo2 = hx/2,
      nxcell = (hx/rcutsmall) - 1,
      nycell = (hy/rcutsmall) - 1,
      dcells = hx/float(nxcell)
      ncells = nxcell*nycell;


//main


//functions




