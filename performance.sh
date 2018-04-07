
#echo "----------------Avec hostfile----------------\n" >> perform.txt

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 1024 -y 1024 -t 40

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 2048 -y 2048 -t 40

mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20
mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20
mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20
mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20
mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 4096 -y 4096 -t 20

if [ "$#" -eq 1 ]; then
	mpirun -n 1 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20
	mpirun -n 2 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20
	mpirun -n 4 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20
	mpirun -n 8 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20
	mpirun -n 16 -hostfile hostfile -bynode ./bin/shalw  -x 8192 -y 8192 -t 20
fi