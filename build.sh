mkdir /tmp/3600594/
make
proc=1
if [ "$#" -eq 1 ]; then
	proc=$1
fi
#mpirun -n 1 -hostfile hostfile -bynode  ./bin/shalw -x 256 -y 256 -t 40 --export --export-path /tmp/3600594/
mpirun -n $proc ./bin/shalw -x 128 -y 128 -t 1000 --export --export-path /tmp/3600594/
echo $proc
if [ "$#" -eq 2 ]; then
    echo "Exec avec python"
    ./visu.py /tmp/3600594/shalw_128x128_T1000_NP$proc.sav
fi