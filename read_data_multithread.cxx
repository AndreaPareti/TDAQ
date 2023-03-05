#include<iostream>
#include<stdlib.h>
#include<stdint.h>
#include<inttypes.h>
#include<iomanip>
#include<string.h>
#include<TH1F.h>
#include<vector>
#include<TH2F.h>
#include<TH2I.h>
#include<algorithm>
#include<chrono>
#include<thread>
#include<mutex>
#include<fstream>
#include"header.h"
const long int num_bins_a = 1024;	//dump a 1445
const long int num_bins_b = 1024;	//dump a 1445
std::mutex mtx;

void coord(std::vector<float> &vec_x, std::vector<float> &vec_y, int size_xy,  std::vector<event> vec_eventi, int num_pixel, float dimensioni_tracker, float dist_interazione, float dist_piani)
//void coord(float *vec_x, float *vec_y, int size_xy,  std::vector<event> vec_eventi, int num_pixel, float dimensioni_tracker, float dist_interazione, float dist_piani)
{
	for(int i=0; i<size_xy; i++)
		{
			vec_x[i] = ((dimensioni_tracker*vec_eventi[i].channel_ID())/num_pixel)-(dimensioni_tracker/2);
			vec_y[i] = dist_interazione + dist_piani*(vec_eventi[i].layer()-1);	
			//std::cout << "X: " << vec_x[i] << "\ty: " << vec_y[i] << std::endl;
		}			
}


void full_header(std::ifstream &data)
{
	// LEGGO PRIMA PAROLA, RITORNO ERRORE SE QUESTA NON E' "BABACACA"
	uint32_t word;
	data.read((char*) &word, sizeof(word));
	//std::cout << "0x" << std::hex << std::setw(8) << std::setfill('0') << word << std::endl;
	uint32_t w = uint32_t(word);
	if(w==0xbabacaca){std::cout << "Checkword Full Header Corretta" << std::endl;}
	else
	{
		std::cerr << "Checkword di inizio file non corrisponde! " ;
	}
	int size_full_header;
	data.read((char*) &size_full_header, sizeof(size_full_header));
	//size_full_header=int(word2);
	//std::cout << "0x" << std::hex << std::setw(8) << std::setfill('0') << word2 << std::endl;
	std::cout << "Dimensione Full Header: " << size_full_header << std::endl;
	uint32_t* buf = new uint32_t[size_full_header-2];
	for(int i=0; i<(size_full_header-2); i++)
	{
		data.read((char*) &buf[i], 4);
	}
	std::cout << "Data Format Versione: " << "0x" << std::hex << std::setw(8) << std::setfill('0') << buf[0] << std::endl;
	std::cout << "Identificazione Detector: " << "0x" << std::hex << std::setw(8) << std::setfill('0') << buf[1] << std::endl;
}

void fragment(std::ifstream &data, float dimensioni_tracker, int num_pixel, std::vector<event> &vec_eventi)
{
	uint32_t event_ID;
	data.read((char*) &event_ID, sizeof(event_ID));
	if(event_ID==0x00ff00ff){return;}
	event_ID = int(event_ID);
	std::cout << "\nEvento numero: " << std::dec << event_ID << std::endl;
	while(true)
	{
		uint32_t new_word;
		data.read((char*) &new_word, sizeof(new_word));
		// SE LA PAROLA CORRISPONDE AL TRAILER DI FINE EVENTO RITORNO 
		if(new_word==0x00ff00ff){ break;}
		event hit(new_word, event_ID);
		vec_eventi.push_back(hit);
	}
std::cout<< "Dimensione payload: " << vec_eventi.size() << std::endl;
}

void hough_cpu_parallel(matrix &m, int min_bin, int max_bin, std::vector<float> vettore_di_bin, std::vector<float> vec_x, std::vector<float> vec_y)
//void hough_cpu_parallel(matrix &m, int min_bin, int max_bin, float *vettore_di_bin, float *vec_x, float *vec_y)
{
	std::lock_guard<std::mutex> guard(mtx);
	for(int j=min_bin; j<=max_bin; j++)
	{
		for(int k=0; k<m.h; k++)
		{
//			std::lock_guard<std::mutex> guard(mtx);
			if(vec_x[k]==0){;}			
			else m.e[k*m.w + j] = (vec_x[k]-vettore_di_bin[j])/vec_y[k];
		}	 
	}
}


void ricerca_max_cpu_parallel(int array_a[num_bins_a][num_bins_b], std::vector<float> vettore_di_bin, std::vector<float> vettore_coeff_ang, int min_bin, int max_bin, int threshold, std::ofstream &reco_retta, uint32_t event_ID)
{

	//int threshold = 20;
	int array_max[num_bins_a];
	int loc[num_bins_a];
	for(int i=min_bin; i<=max_bin; i++)
	{
		int location = 1;
		for(int c=0; c<num_bins_b; c++)
		{
			if(array_a[i][c] > array_a[i][location])
			{
				location = c;
				array_max[i] = array_a[i][location];
			}
		}
		loc[i] = location; 
		std::lock_guard<std::mutex> guard(mtx);

		//std::cout << array_max[i] << " ";
		if(i>2 and array_max[i-1] > threshold and array_max[i-2] < array_max[i-1] and array_max[i-1] > array_max[i])
		{
			reco_retta << "Evento " << event_ID << "\tMassimo " << array_max[i-1] << "\tCoeffAng " << vettore_coeff_ang[i-1] << "\toffset: " << vettore_di_bin[loc[i-1]]  << std::endl;
		}
		if(i>2 and array_max[i-1]>threshold and i<(num_bins_b-2) and array_max[i-2] == array_max[i-1] and array_max[i-2]>array_max[i-3] and array_max[i-2]>array_max[i])
		{
			reco_retta << "Evento " << event_ID << "\tMassimo " << array_max[i-1] << "\tCoeffAng " << vettore_coeff_ang[i-1] << "\toffset: " << vettore_di_bin[loc[i-1]]  << std::endl;
		
		}
	}
}

void histo_fill_cpu_parallel(TH2F *istogramma,  int array_a[num_bins_a][num_bins_b], int size_xy, matrix m, std::vector<float> vettore_di_bin, int min_bin, int max_bin)
{
	// RIEMPIO ISTOGRAMMA DELLE FREQUENZE CON I COEFICIENTI ANGOLARI TROVATI
	for(int i=min_bin; i<=max_bin; i++)
	{
		for(int j=0; j<size_xy; j++)	
		{
			istogramma->Fill(m.e[i+j*m.w], vettore_di_bin[i]);
		}
		for(int k=0; k<num_bins_a; k++)
		{
			array_a[k][i] = istogramma->GetBinContent(k,i);
		}
	}
}
int main()
{
	// Quanti core (reali o virtuali) abbiamo a disposizione?
    const int ncores = std::thread::hardware_concurrency();
    std::cout << "La parallelizzazione massima su questo PC è " << ncores << std::endl;
	
	// LETTURA FILE DI CONFIGURAZIONE 
	std::string nomi_param;
	std::vector<float> parametri;
	std::ifstream infile("config.txt");
	float par_value;
	while(infile >> nomi_param >> par_value)
	{
		parametri.push_back(par_value);
	}
	infile.close();
	// INSERISCO PARAMETRI DEL RIVELATORE E DI GENERAZIONE SEGNALI
	int num_piani= int(parametri[0]);
	float dist_interazione = parametri[1];
	float dist_piani = parametri[2];
	float dimensioni_tracker = parametri[3];
	float dim_pixel = parametri[4];
	float max_offset = parametri[10];
	int threshold = int(parametri[11]);
	int num_pixel = int((dimensioni_tracker*1000)/dim_pixel);

	float max_coeff_ang = 3;

	// INIZIO LETTURA DI FILE DATA
	// PARTO LEGGENDO FULL HEADER
	std::ifstream data("file.dat", std::ios::binary);
	std::ofstream reco_retta("reco_retta.dat");
	full_header(data);

	// RIEMPIO UN VETTORE I CUI ELEMENTI SONO I BIN DELL'ISTOGRAMMA IN CUI VIENE SUDDIVISO L'OFFSET
	std::vector<float> vettore_di_bin(num_bins_b);
	for(int i=0; i<num_bins_b; i++)
	{
		vettore_di_bin[i] = -max_offset + i*(2*max_offset/num_bins_b);
	}
	std::vector<float> vettore_coeff_ang(num_bins_a);
	//float vettore_coeff_ang[num_bins_a];
	for(int j=0; j<num_bins_a; j++)
	{
		vettore_coeff_ang[j] = -max_coeff_ang + j*(2*max_coeff_ang/num_bins_a);
//		std::cout  << vettore_coeff_ang[j] << " " ;
	}

	// INIZIO LETTURA FRAMMENTI DATA
	while(!data.eof())
	{
		uint32_t word1;
		data.read((char*) &word1, sizeof(word1));
		uint32_t w1 = uint32_t(word1);
		std::vector<event> vec_eventi;
		if(w1==0xdeadcafe)
		{
		// CREAZIONE ISTOGRAMMA 2D: SU UN ASSE METTO IL BIN CORRIPONDENTE ALLA SUDDIVISIONE DELL'OFFSET
		// SUL SECONDO ASSE INSERISCO I COEFFICIENTI ANGOLARI TROVATI
			TH2F* istogramma = new TH2F("h2", "h2; Coefficiente angolare;  Offset" , num_bins_a, -max_coeff_ang, max_coeff_ang, num_bins_b, -(max_offset), (max_offset));
			istogramma->AddDirectory(false);
			// LEGGO IL PAYLOAD E RIEMPIO UN VETTORE EVENT(EVENT_ID, LAYER, CHANNEL_ID)
			fragment(data, dimensioni_tracker, num_pixel, vec_eventi);
			const int size_xy = vec_eventi.size();

			std::cout << "frammento ricostruito correttamente" << std::endl;
			// INIZIALIZZO E RIEMPIO VETTORI DA MANDARE A GPU E CPU
			std::vector<float> vec_x(size_xy);
			std::vector<float> vec_y(size_xy);

			// PER CIASCUN EVENTO ESTRAGGO LE COORDINATE (X, Z)

			coord(vec_x, vec_y, size_xy, vec_eventi, num_pixel, dimensioni_tracker, dist_interazione, dist_piani);
			uint32_t event_ID = vec_eventi[1].event_ID();

			// CREO MATRICE I CUI ELEMENTI SONO I COEFFICIENTI ANGOLARI DELLE RETTE
			// CHE CONNETTONO CIASCUN BIN ALLE COORDINATE (X,Y) DELL'HIT
			matrix m;
			m.w = num_bins_b;
			m.h = size_xy;
			m.e = (float*) malloc(m.w*m.h*sizeof(float));
			
			// ### RIEMPIMENTO MATRICE USANDO CPU ###
			int array_a[num_bins_a][num_bins_b];
			// USANDO CPU E MULTITHREADING -> uso 5 thread e assumo che il numero di bin scelto sia multiplo
			std::vector<std::thread> vT; 
			int nthread = ncores;
			int width = m.w/nthread;
			auto start_hough = std::chrono::steady_clock::now();
			//divido num_bins_b in ncores intervalli
			for(int count=0; count<nthread; count++)
			{
				int min_bin = count*width;
				int max_bin = (count+1)*width - 1;
				//hough_cpu_parallel(m, min_bin, max_bin, vettore_di_bin, vec_x, vec_y);
				vT.push_back( std::thread(hough_cpu_parallel, std::ref(m), min_bin, max_bin, vettore_di_bin, vec_x, vec_y));
			}
			for(auto &t: vT)
				t.join();
			auto end_hough = std::chrono::steady_clock::now();
			auto t_hough = std::chrono::duration_cast<std::chrono::milliseconds> (end_hough - start_hough).count();
			std::cout << "Time for hough matrix calculation: " << t_hough << " ms" << std::endl;
			// PREPRARO ARRAY BIDIMENSIONALE IL CUI CONTENUTO SONO LE FREQUENZE DELL'ISTOGRAMMA
			
//			histo_fill_cpu(istogramma,  array_a, size_xy, m, vettore_di_bin);
			auto start_histo = std::chrono::steady_clock::now();
			std::vector<std::thread> v;
			for(int count=0; count<nthread; count++)
			{
				int min_bin = count*width;
				int max_bin = (count+1)*width - 1;
				v.push_back( std::thread(histo_fill_cpu_parallel, istogramma, std::ref(array_a), size_xy, m, vettore_di_bin, min_bin, max_bin));
			}
			for(auto &boh: v)
				boh.join();
			auto end_histo = std::chrono::steady_clock::now();
			auto t_histo = std::chrono::duration_cast<std::chrono::milliseconds> (end_histo-start_histo).count();
			std::cout << "Histo filling time: " << t_histo << " ms" << std::endl;


			std::cout << "Ora parte ricerca massimi: " << std::endl;

			auto start_max = std::chrono::steady_clock::now();
			std::vector<std::thread> altro_vettore;
			int width_a = num_bins_a/ncores;
			for(int count=0; count<nthread; count++)
			{
				int min_bin = count*width_a;
				int max_bin = (count+1)*width_a;
				altro_vettore.push_back( std::thread(ricerca_max_cpu_parallel, array_a, vettore_di_bin, vettore_coeff_ang, min_bin, max_bin, threshold, std::ref(reco_retta), event_ID));
			}
			for(auto &altro_boh: altro_vettore)
				altro_boh.join();
			auto end_max = std::chrono::steady_clock::now();
			auto t_max = std::chrono::duration_cast<std::chrono::microseconds> (end_max-start_max).count();
			std::cout << "Time to find local maxima: " << t_max << " micros" << std::endl;



			std::cout << "fine Evento " << event_ID << std::endl;
			
			std::fill(vec_x.begin(), vec_x.end(), 0);
			std::fill(vec_y.begin(), vec_y.end(), 0);

			memset(array_a, 0, sizeof(array_a));
			istogramma->Reset();
			reco_retta.close();
		}

	

	}	



	std::cout << "finito!" << std::endl;
	return 0;
}
