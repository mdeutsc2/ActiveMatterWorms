from amat2 import AmatterSim

sim = AmatterSim()
sim.init_worms()
sim.init_arrays()
print(sim.P)

sim.run()
print(sim.P)
