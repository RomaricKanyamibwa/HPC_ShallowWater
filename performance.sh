
#echo "----------------Avec hostfile----------------" >> perform.txt

echo "-----------------Test 1024*1024-----------------" >> perform.txt

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40

echo "-----------------Non-Blocking Test-----------------" >> perform.txt

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40 --non_block_comm
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40 --non_block_comm
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40 --non_block_comm
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40 --non_block_comm
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40 --non_block_comm

echo "-----------------Test 2048*2048-----------------" >> perform.txt

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40

echo "-----------------Non-Blocking Test-----------------" >> perform.txt

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40 --non_block_comm
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40 --non_block_comm
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40 --non_block_comm
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40 --non_block_comm
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40 --non_block_comm

echo "-----------------Test 4096*4096-----------------" >> perform.txt

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20

echo "-----------------Non-Blocking Test-----------------" >> perform.txt

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20 --non_block_comm
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20 --non_block_comm
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20 --non_block_comm
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20 --non_block_comm
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20 --non_block_comm

echo "-----------------Test 8192*8192-----------------" >> perform.txt

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20 
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20

echo "-----------------Non-Blocking Test-----------------" >> perform.txt

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20 --non_block_comm
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20 --non_block_comm
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20 --non_block_comm
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20 --non_block_comm
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20 --non_block_comm
