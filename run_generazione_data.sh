./compila.sh -o data data_gen.cxx header.cxx
./data
sed -e  s/Event_ID://g -e s/m://g -e s/q://g data_retta.dat >> data_tmp.dat
rm data_retta.dat
mv data_tmp.dat data_retta.dat 

