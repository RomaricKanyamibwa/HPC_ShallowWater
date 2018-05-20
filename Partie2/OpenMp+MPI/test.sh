
echo "-----------------Sans export-----------------" >> perform.txt

echo "-----------------Test 2048*2048-----------------" >> perform.txt

../../../Version_Seq/bin/shalw   -x 2048 -y 2048 -t 80

mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc  -x 2048 -y 2048 -t 80
mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  -x 2048 -y 2048 -t 80
mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc --non_block_comm  -x 2048 -y 2048 -t 80
mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  --non_block_comm  -x 2048 -y 2048 -t 80

mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc  -x 2048 -y 2048 -t 80
mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  -x 2048 -y 2048 -t 80
mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc --non_block_comm  -x 2048 -y 2048 -t 80
mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  --non_block_comm  -x 2048 -y 2048 -t 80

echo "-----------------Test 4096*4096-----------------" >> perform.txt

../../../Version_Seq/bin/shalw   -x 4096 -y 4096 -t 80

mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc  -x 4096 -y 4096 -t 80
mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  -x 4096 -y 4096 -t 80
mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc --non_block_comm  -x 4096 -y 4096 -t 80
mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  --non_block_comm  -x 4096 -y 4096 -t 80

mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc  -x 4096 -y 4096 -t 80
mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  -x 4096 -y 4096 -t 80
mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc --non_block_comm  -x 4096 -y 4096 -t 80
mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  --non_block_comm  -x 4096 -y 4096 -t 80

echo "-----------------Test 8192*8192-----------------" >> perform.txt

../../../Version_Seq/bin/shalw   -x 8192 -y 8192 -t 20

mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc  -x 8192 -y 8192 -t 20
mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  -x 8192 -y 8192 -t 20
mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc --non_block_comm  -x 8192 -y 8192 -t 20
mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw  --non_block_comm  -x 8192 -y 8192 -t 20

mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc  -x 8192 -y 8192 -t 20
mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  -x 8192 -y 8192 -t 20
mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  --bloc --non_block_comm  -x 8192 -y 8192 -t 20
mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw  --non_block_comm  -x 8192 -y 4096 -t 20


# echo "-----------------Avec export-----------------" >> perform.txt

# mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw --export --export-path /tmp/3600594/ --bloc 
# mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw --export --export-path /tmp/3600594/
# mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw --export --export-path /tmp/3600594/ --bloc --non_block_comm
# mpirun -n 2 -bynode -hostfile hostfile  -bynode ./bin/shalw --export --export-path /tmp/3600594/ --non_block_comm

# mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw --export --export-path /tmp/3600594/ --bloc
# mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw --export --export-path /tmp/3600594/
# mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw --export --export-path /tmp/3600594/ --bloc --non_block_comm
# mpirun -n 4 -bynode -hostfile hostfile  -bynode ./bin/shalw --export --export-path /tmp/3600594/ --non_block_comm
