# run forward simulation to obtain observed data
make clean
make install  MODEL=hh_toy_true.c
mpirun -np 8 nice -19 ../bin/ifos3d ./in_and_out/ifos3d_toy_FW.json | tee ./in_and_out/ifos3D.out

cp model/toy.vs_it0 model/toy.vs.true
cp model/toy.vp_it0 model/toy.vp.true
cp model/toy.rho_it0 model/toy.rho.true


# invert observed data using homogeneous starting model
make clean
make install MODEL=hh_toy_start.c
 mpirun -np 8 nice -19 ../bin/ifos3d ./in_and_out/ifos3d_toy.json | tee ./in_and_out/ifos3D_INV.out
