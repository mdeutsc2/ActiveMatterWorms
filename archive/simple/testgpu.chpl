use GPU;
use Time;

config const n = 16;
config const L = 100.0;
const dt = 0.001;
const device = here.gpus[0];
const X = 1;
const Y = 2;

writeln(here.gpus);
writeln(device);

record Particle {
  var x,y: real;
  var l: [1..64] int;

  proc init() {
    this.x = 0.0;
    this.y = 0.0;
    this.l = 0;
  }
}

on device {
  /* var A : [1..n] real;
  var B : [1..n] real;
  forall i in 1..n {
    B[i] = 1.0;
  }
  forall i in 1..n {
    A[i] = B[i] * 2.0;
  }
  writeln(A);
  writeln(B);

  var p : [1..n] Particle;
  writeln("p.x\t",p.x);
  writeln("p.y\t",p.y);
  forall i in 1..n {
    p[i].x = p[i].x + i:real;
    p[i].y = A[i];
  }
  writeln("after kernel");
  writeln("p.x\t",p.x);
  writeln("p.y\t",p.y);
  writeln(p[1].l);
  forall i in 1..n {
    B[i] = 0.0;
  }
  A = test(B);
  writeln(A); */


  var forces: [1..n, 1..2] real;
  var velocities: [1..n, 1..2] real;
  var positions: [1..n, 1..2] real;
  var KE : [1..n] real;
  var pDomain = {1..n};
  // init particles in lattice
  if (n != sqrt(n)*sqrt(n)) {writeln("non-square n");halt();}
  var row_length = sqrt(n):int;
  var row,col,center,spacing :real;
  var hxo2 = L/2;
  var rcutsmall = 2.0**(1.0/6.0);
  spacing = 1.0;
  center = hxo2 - spacing*(row_length/2);
  foreach i in 1..n {
      row = i % row_length;
      col = ((i - row)/row_length)-1;
      positions[i,X] = center + spacing*rcutsmall*row;
      positions[i,Y] = center + spacing*rcutsmall*col;
      velocities[i,X] = 0.1;//*gaussRand(0.0,1.0);
      velocities[i,Y] = 0.1;//*gaussRand(0.0,1.0);
      //KE[i] = 0.0;
      KE[i] = 0.5*(velocities[i,X] * velocities[i,X] + velocities[i,Y] * velocities[i,Y]);
  }
  var vxave = 0.0;
  var vyave = 0.0;
  for i in 1..n {
    vxave += velocities[i,X];
    vyave += velocities[i,Y];
  }
  vxave = vxave/n;
  vyave = vyave/n;
  foreach i in 1..n {
      velocities[i,X] = velocities[i,X] - vxave;
      velocities[i,Y] = velocities[i,Y] - vyave;
  }
  writeln("init done");
  calc_forces(positions,forces);
  var xt : stopwatch;
  var total_time = 0.0;
  writeln("starting");
  xt.start();
  for istep in 1..5000 {
    if (istep % 100 == 0) {
      xt.stop();
      total_time = xt.elapsed();
      writeln("Step: ",istep,"\t",
              istep/total_time," iter/s\tElapsed ",
              total_time," s");
      xt.start();
    }
    update_position(positions,velocities,forces);
    calc_forces(positions,forces);
    update_velocities(velocities,forces,KE);
  }
}

proc update_position(positions,velocities,forces) {
  foreach i in 1..n {
    positions[i,X] += velocities[i,X]*dt + (forces[i,X]/2)*(dt**2);
    positions[i,Y] += velocities[i,Y]*dt + (forces[i,Y]/2)*(dt**2);
    if positions[i,X] > L {
      positions[i,X] = positions[i,X] - L;
    } else if positions[i,X] < 0.0 {
      positions[i,X] = positions[i,X] + L;
    }
    if positions[i,Y] > L {
      positions[i,Y] = positions[i,Y] - L;
    } else if positions[i,Y] < 0.0 {
      positions[i,Y] = positions[i,Y] + L;
    }
    velocities[i,X] += 0.5*dt*forces[i,X];
    velocities[i,Y] += 0.5*dt*forces[i,Y];
    forces[i,X] = 0.0;
    forces[i,Y] = 0.0;
  }
}

proc update_velocities(velocities,forces,KE) {
  foreach i in 1..n {
    velocities[i,X] += 0.5*dt*forces[i,X];
    velocities[i,Y] += 0.5*dt*forces[i,Y];
    KE[i] = 0.5*(velocities[i,X] * velocities[i,X] + velocities[i,Y] * velocities[i,Y]);
  }
}

proc calc_forces(positions,forces) {
  forall i in 1..n with (+ reduce forces) {
    var dx,dy,r2,ffor,ffx,ffy:real;
    for j in 1..n {
      if j <= i {
        continue;
      }
      dx = positions[j,1] - positions[i,1];
      dy = positions[j,2] - positions[i,2];
      r2 = (dx*dx + dy*dy);
      if (r2 <= 2.5) {
        ffor = -48.0*r2**(-7.0) + 24.0*r2**(-4.0);
        ffx = ffor*dx;
        ffy = ffor*dy;
        forces[i,1] += ffx;
        forces[i,2] += ffy;
        forces[j,1] -= ffx;
        forces[j,2] -= ffy;
      }
    }
  }
}

proc test(A) {
  forall i in 1..n {
    A[i] += i;
  }
  return A;
}
