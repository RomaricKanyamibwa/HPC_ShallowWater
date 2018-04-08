
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

echo "-----------------MPI_IO Test-----------------" >> perform.txt

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40

rm -f shalw_*

mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40 --non_block_comm 
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40 --mpi_io
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40 --mpi_io_non_block
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40 --non_block_comm --mpi_io
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40 --non_block_comm --mpi_io_non_block

rm -f shalw_*

mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40 --non_block_comm 
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40 --mpi_io
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40 --mpi_io_non_block
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40 --non_block_comm --mpi_io
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 40 --non_block_comm --mpi_io_non_block

rm -f shalw_*

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80

rm -f shalw_*

mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80 --non_block_comm 
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80 --mpi_io
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80 --mpi_io_non_block
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80 --non_block_comm --mpi_io
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80 --non_block_comm --mpi_io_non_block

rm -f shalw_*

mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80 --non_block_comm 
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80 --mpi_io
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80 --mpi_io_non_block
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80 --non_block_comm --mpi_io
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw --export -x 512 -y 512 -t 80 --non_block_comm --mpi_io_non_block

rm -f shalw_*

