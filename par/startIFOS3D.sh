#lamboot -v lamhosts
mpirun -np 8 nice -19 ../bin/ifos3d ./in_and_out/ifos3d_toy.inp | tee ./in_and_out/ifos3D.out
