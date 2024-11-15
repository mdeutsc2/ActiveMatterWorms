use Random;
use List;
const numParticles = 1000;  // Total number of particles
const boxSize = 10.0;        // Size of the simulation box
const cellSize = 2.5;        // Size of each cell

// Define a particle structure
record Particle {
var position: [1..3] real;
var id: int;
}


proc getCellIndex(p: Particle, cellSize: real) {
    return [(p.position[1] / cellSize):int + 1,(p.position[2] / cellSize):int + 1,(p.position[3] / cellSize):int + 1];
}

// Function to update the particle's position
proc updateParticlePosition(inout p: Particle, newPos, inout cellList) {
  var oldCellIndex = getCellIndex(p, cellSize);
  p.position[1] = newPos[0];
  var newCellIndex = getCellIndex(p, cellSize);

  // If the particle has moved to a new cell, update the cell list
  if (oldCellIndex[0] != newCellIndex[0] || oldCellIndex[1] != newCellIndex[1] || oldCellIndex[2] != newCellIndex[2]) {
    // Remove the particle from the old cell
    cellList[oldCellIndex[0], oldCellIndex[1], oldCellIndex[2]].remove(p.id);
    
    // Append the particle ID to the new cell
    cellList[newCellIndex[0], newCellIndex[1], newCellIndex[2]].pushBack(p.id);
  }
}

// Main function
proc main() {
// Initialize particles with random positions
var particles: [1..numParticles] Particle;
var cl_nb: [1..numParticles] int;
var bf_nb: [1..numParticles] int;
var randStream = new randomStream(real);
for i in 1..numParticles {
    particles[i].position = (randStream.next()*boxSize, randStream.next()*boxSize, randStream.next()*boxSize);
    particles[i].id = i;
    cl_nb[i] = 0;
    bf_nb[i] = 0;
}

// Calculate number of cells in each dimension
var numCells = ((boxSize/cellSize):int,(boxSize/cellSize):int,(boxSize/cellSize):int);

// Create a cell list (3D array of lists)
var cellList: [1..numCells[0], 1..numCells[1], 1..numCells[2]] list(int);

// Populate cell list
for p in particles {
    var cellIndex = getCellIndex(p,cellSize);
    cellList[cellIndex[0], cellIndex[1], cellIndex[2]].pushBack(p.id);
}

// Find neighbors for each particle
for p in particles {
    var cellIndex = getCellIndex(p,cellSize);
    
    var nnab = 0;
    // Check neighboring cells (including the cell itself)
    for dx in -1..1 {
        for dy in -1..1 {
            for dz in -1..1 {
                var neighborCellIndex = [
                    cellIndex[0] + dx,
                    cellIndex[1] + dy,
                    cellIndex[2] + dz
                ];

                // Check bounds
                if neighborCellIndex[0] >= 1 && neighborCellIndex[0] <= numCells[0] &&
                    neighborCellIndex[1] >= 1 && neighborCellIndex[1] <= numCells[1] &&
                    neighborCellIndex[2] >= 1 && neighborCellIndex[2] <= numCells[2] {
                    // List the neighbors
                    for neighborId in cellList[neighborCellIndex[0], neighborCellIndex[1], neighborCellIndex[2]] {
                        if neighborId != 0 && neighborId != p.id {
                            //writeln("  Particle ", neighborId);
                            var distance = sqrt((p.position[1] - particles[neighborId].position[1])**2 +  
                            (p.position[2] - particles[neighborId].position[2])**2 +
                            (p.position[3] - particles[neighborId].position[3])**2);
                            if distance < 2.5 {
                                nnab += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    cl_nb[p.id] = nnab;
    // writeln("Neighbors for Particle ", p.id, " : ",nnab);
}

for i in 990..1000 {
    var newPos = (randStream.next()*boxSize, randStream.next()*boxSize, randStream.next()*boxSize);
    updateParticlePosition(particles[i], newPos,cellList);
    // var oldCellIndex = getCellIndex(particles[i], cellSize);
    // // p.position[1] = newx;
    // // p.position[2] = newy;
    // // p.position[3] = newz;
    // particles[i].position = newPos;
    // var newCellIndex = getCellIndex(particles[i], cellSize);
    // // writeln(oldCellIndex," ",newCellIndex,oldCellIndex == newCellIndex);
    // // If the particle has moved to a new cell, update the cell list
    // if (oldCellIndex[0] != newCellIndex[0] || oldCellIndex[1] != newCellIndex[1] || oldCellIndex[2] != newCellIndex[2]) {
    //     // Remove the particle from the old cell
    //     cellList[oldCellIndex[0], oldCellIndex[1], oldCellIndex[2]].remove(particles[i].id);
        
    //     // Append the particle ID to the new cell
    //     cellList[newCellIndex[0], newCellIndex[1], newCellIndex[2]].pushBack(particles[i].id);
    // }
}
writeln("----------");
for i in 1..1000 {
    var nnab = 0;
    for j in i+1..1000 {
        var distance = sqrt((particles[i].position[1] - particles[j].position[1])**2 +  
                            (particles[i].position[2] - particles[j].position[2])**2 +
                            (particles[i].position[3] - particles[j].position[3])**2);
        if distance < 2.5 {
            nnab += 1;
        }
    }
    bf_nb[i] = nnab;
}
for p in particles {
    writeln("Neighbors for Particle ", p.id, " : ",cl_nb[p.id]," : ",bf_nb[p.id]," ",cl_nb[p.id] == bf_nb[p.id]);
}
}