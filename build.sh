mkdir /tmp/3600594/
make
mpirun -n 1 -hostfile hostfile -bynode  ./bin/shalw --export --export-path /tmp/3600594/
if [ "$#" -eq 2 ]; then
    echo "Exec avec python"
    ./visu.py /tmp/3600594/shalw_256x256_T1000.sav
fi