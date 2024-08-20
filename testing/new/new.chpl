use Random;

// Constants
const numParticles = 100; // Number of particles
const numSteps = 1000;    // Number of time steps
const dt = 0.01;          // Time step
const boxLength = 10.0;   // Length of simulation box
const cutoff = 2.5;       // Cutoff distance for LJ potential

// Particle structure
record Particle {
  var x: [1..3] real;  // Position
  var v: [1..3] real;  // Velocity
}
var randStream = new RandomStream(real); // creating random number generator

// Initialize particles
proc initParticles(): [1..numParticles] Particle {
  var particles: [1..numParticles] Particle;
  var rand1 = randStream.getNext();

  for p in 1..numParticles {
    for i in 1..3 {
      particles[p].x[i] = randStream.getNext()*boxLength;
      particles[p].v[i] = (randStream.getNext()*2.0)-1.0;//randGen.uniform(-1.0, 1.0);
    }
  }
  return particles;
}

// Lennard-Jones potential
proc ljPotential(r: real): real {
  var r6 = r**6;
  var r12 = r6 * r6;
  return 4.0 * (1.0 / r12 - 1.0 / r6);
}

// Dot product function
proc dot(vec1: [1..3] real, vec2: [1..3] real): real {
  var result: real = 0.0;
  for i in 1..3 {
    result += vec1[i] * vec2[i];
  }
  return result;
}

// Calculate distance between particles
proc distance(p1: Particle, p2: Particle): real {
  var dx = p2.x[1] - p1.x[1];
  var dy = p2.x[2] - p1.x[2];
  var dz = p2.x[3] - p1.x[3];
  return sqrt(dx*dx + dy*dy + dz*dz);
}

// Main simulation loop
proc main() {
  var particles = initParticles();

  // Molecular dynamics loop
  for step in 1..numSteps {
    // Update positions
    coforall p in 1..numParticles {
      for i in 1..3 {
        particles[p].x[i] += particles[p].v[i] * dt;
        // Apply periodic boundary conditions
        particles[p].x[i] -= boxLength * floor(particles[p].x[i] / boxLength);
      }
    }

    // Calculate forces and update velocities
    coforall p1 in 1..numParticles {
      for p2 in (p1 + 1)..numParticles {
        var dist = distance(particles[p1], particles[p2]);
        if dist < cutoff {
          var rVec: [1..3] real;
          for i in 1..3 {
            rVec[i] = particles[p2].x[i] - particles[p1].x[i];
          }
          var r = sqrt(dot(rVec, rVec));
          var f = ljPotential(r) / r;
          for i in 1..3 {
            atomic {
            particles[p1].v[i] += f * rVec[i] * dt;
            particles[p2].v[i] -= f * rVec[i] * dt;
            }
          }
        }
      }
    }

    // Print particle positions every 100 steps
    if step % 100 == 0 {
      writeln("Step: ", step);
      for p in 1..numParticles {
        writeln("Particle ", p, " position: ", particles[p].x);
      }
    }
  }
}

// Run the simulation
main();

