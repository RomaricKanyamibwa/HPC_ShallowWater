
#echo "----------------Avec hostfile----------------" >> perform.txt

echo "-----------------Test 1024*1024-----------------" >> perform.txt

../../../Version_Seq/bin/shalw  -x 1024 -y 1024 -t 40
../../../Version_Seq/bin/shalw  -x 1024 -y 1024 -t 40 --cache
./bin/shalw  -x 1024 -y 1024 -t 40

echo "-----------------Test 2048*2048-----------------" >> perform.txt

echo "****Decomposition par Bande****" >> perform.txt

../../../Version_Seq/bin/shalw  -x 2048 -y 2048 -t 40
../../../Version_Seq/bin/shalw  -x 2048 -y 2048 -t 40 --cache
./bin/shalw  -x 2048 -y 2048 -t 40 


echo "-----------------Test 4096*4096-----------------" >> perform.txt

echo "****Decomposition par Bande****" >> perform.txt

../../../Version_Seq/bin/shalw  -x 4096 -y 4096 -t 20
../../../Version_Seq/bin/shalw  -x 4096 -y 4096 -t 20 --cache
./bin/shalw  -x 4096 -y 4096 -t 20


echo "-----------------Test 8192*8192-----------------" >> perform.txt

echo "****Decomposition par Bande****" >> perform.txt

../../../Version_Seq/bin/shalw  -x 8192 -y 8192 -t 20 
../../../Version_Seq/bin/shalw  -x 8192 -y 8192 -t 20 --cache
./bin/shalw  -x 8192 -y 8192 -t 20



echo "-----------------IO Test-----------------" >> perform.txt

echo "-----------------t =40 -----------------" >> perform.txt

../../../Version_Seq/bin/shalw --export -x 512 -y 512 -t 40
../../../Version_Seq/bin/shalw --export -x 512 -y 512 -t 40 --cache
./bin/shalw --export -x 512 -y 512 -t 40
rm -f shalw_*

echo "-----------------t =80 -----------------" >> perform.txt

../../../Version_Seq/bin/shalw --export -x 512 -y 512 -t 80
../../../Version_Seq/bin/shalw --export -x 512 -y 512 -t 80 --cache
./bin/shalw --export -x 512 -y 512 -t 80

rm -f shalw_*


