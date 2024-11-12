#/bin/bash
rm energies.dat
rm params.dat
rm amatter.log
rm *.restart

make clean
make
#./runsim.sh --boundary=3 --np=120 --nworms=3200 --kbt=0.6 --fluid_rho=0.05 --dt=0.005 --save_interval=12000 --nsteps=24000000 --L=6.4 --gamma=24.0 --fdrag=0.04 

./amatter.x --boundary=4 --np=120 --nworms=3200 --kbt=0.6 --fluid_rho=0.05 --dt=0.005 --save_interval=12000 --nsteps=24000000 --L=6.4 --gamma=24.0 --fdrag=0.04 