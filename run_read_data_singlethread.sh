./compila.sh -o read_data_singlethread read_data_singlethread.cxx header.cxx
./read_data_singlethread
sed -e s/Evento//g -e s/Massimo//g -e s/CoeffAng//g -e s/offset//g reco_retta.dat >> reco_tmp.dat
rm reco_retta.dat
mv reco_tmp.dat reco_retta.dat
