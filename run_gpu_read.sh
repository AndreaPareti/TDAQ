./cuda_compila.sh gpu_read.cu header.cxx
./gpu_read | grep -w Evento >> reco_retta.dat
sed -e s/Evento//g -e s/Massimo//g -e s/CoeffAng//g -e s/offset//g reco_retta.dat >> reco_tmp.dat
rm reco_retta.dat
mv reco_tmp.dat reco_retta.dat
