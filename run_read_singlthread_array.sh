./compila.sh -o read_singlethread_array read_singlethread_array.cxx header.cxx 
./read_singlethread_array 
sed -e s/Evento//g -e s/Massimo//g -e s/CoeffAng//g -e s/offset//g reco_retta.dat >> reco_tmp.dat
rm reco_retta.dat
mv reco_tmp.dat reco_retta.dat
