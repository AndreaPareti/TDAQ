#include<iostream>
#include<stdlib.h>
#include<stdint.h>
#include<inttypes.h>
#include<string.h>
#include<TH1F.h>
#include<vector>
#include<iomanip>
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

//void coord(std::vector<float> &vec_x, std::vector<float> &vec_y, int size_xy,  std::vector<event> vec_eventi, int num_pixel, float dimensioni_tracker, float dist_interazione, float dist_piani)
void coord(float *vec_x, float *vec_y, int size_xy,  std::vector<event> vec_eventi, int num_pixel, float dimensioni_tracker, float dist_interazione, float dist_piani)
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
	std::cout << "\nEvent ID: " << std::dec << event_ID << std::endl;
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

void hough_cpu(matrix &m, float* vettore_di_bin, float *vec_x, float *vec_y)
{
	for(int k=0; k<m.h; k++) // i indicatore della riga
		for(int j=0; j<m.w; j++) // j indicarore colonna
		{
			if(vec_x[k]==0){;}
			else m.e[k*m.w + j] = (vec_x[k]-vettore_di_bin[j])/vec_y[k];
		}
}



void ricerca_max_cpu_array(unsigned int *histo, float *vettore_di_bin, float*vettore_coeff_ang, unsigned int event_ID, int threshold, std::ofstream &reco_retta)
{
	int array_max[num_bins_a];
	int loc[num_bins_a];
	for(int bin_a=0; bin_a<num_bins_a; bin_a++)
	{
		int location=1;
		for(int c=0; c<num_bins_b; c++)
		{
			if(histo[bin_a*num_bins_b+c] > histo[bin_a*num_bins_b+location])
			{
				location=c;
			}
			loc[bin_a]=location;
			array_max[bin_a]=histo[bin_a*num_bins_b+location];
		}
	}
	for(int bin_a=0; bin_a<num_bins_a; bin_a++)
	{
        if(bin_a>1 and bin_a<num_bins_a and array_max[bin_a]>threshold and array_max[bin_a-1]<array_max[bin_a] and array_max[bin_a]>array_max[bin_a+1])
        {
        	reco_retta << "Evento " << event_ID << "\tMassimo " << array_max[bin_a-1] << "\tCoeffAng" << vettore_coeff_ang[bin_a] << "\toffset " << vettore_di_bin[loc[bin_a]] << std::endl;
        }	
	    else if(bin_a>1 and bin_a < num_bins_a-1 and array_max[bin_a]>threshold and array_max[bin_a-1]==array_max[bin_a] 
        	and array_max[bin_a-2]<array_max[bin_a-1] and array_max[bin_a]>array_max[bin_a+1])
        	{
           reco_retta << "Evento " << event_ID << "\tMassimo " << array_max[bin_a-1] << "\tCoeffAng" << vettore_coeff_ang[bin_a] << "\toffset " << vettore_di_bin[loc[bin_a]] << std::endl;

         }   
	}

}

void ricerca_max_locale_cpu(unsigned int *histo, float *vettore_di_bin, float *vettore_coeff_ang, unsigned int event_ID, int threshold, std::ofstream &reco_retta)
{
	int array_max[num_bins_a];
	int loc[num_bins_a];

	for(int bin_a=0; bin_a<num_bins_a; bin_a++)
	{
		int location=1;
		for(int c=0; c<num_bins_b; c++)
		{
			if(histo[bin_a*num_bins_b+c]>threshold and c>2 and c<num_bins_b-2 and histo[bin_a*num_bins_b+c]>=histo[bin_a*num_bins_b+(c-1)] and histo[bin_a*num_bins_b+c]>histo[bin_a*num_bins_b+(c-2)] and histo[bin_a*num_bins_b+c]>=histo[bin_a*num_bins_b+(c+1)] and histo[bin_a*num_bins_b+c]>histo[bin_a*num_bins_b+(c+2)])
				if(bin_a>2 and bin_a<num_bins_a-2 and histo[bin_a*num_bins_b+c]>=histo[(bin_a-1)*num_bins_b+c] and histo[bin_a*num_bins_b+c]>histo[(bin_a-2)*num_bins_b+c] and histo[(bin_a)*num_bins_b+c]>=histo[(bin_a+1)*num_bins_b+c] and histo[bin_a*num_bins_b+c]>histo[(bin_a+2)*num_bins_b+c])					
					if(histo[bin_a*num_bins_b+c]>histo[(bin_a-1)*num_bins_b+(c-1)] and histo[bin_a*num_bins_b+c]>histo[(bin_a-1)*num_bins_b+(c+1)] and histo[bin_a*num_bins_b+c]>histo[(bin_a+1)*num_bins_b+(c-1)] and histo[bin_a*num_bins_b+c]>histo[(bin_a+1)*num_bins_b+(c+1)]){
						//if(histo[bin_a*num_bins_b+c]>histo[(bin_a-2)*num_bins_b+(c-1)] and histo[bin_a*num_bins_a+c]>histo[(bin_a-2)*num_bins_b+(c+1)] and histo[bin_a*num_bins_b+c]>histo[(bin_a+2)*num_bins_b+(c-1)] and histo[bin_a*num_bins_b+c]>histo[(bin_a+2)*num_bins_b+(c+1)])
							//if(histo[bin_a*num_bins_b+c]>histo[(bin_a-1)*num_bins_b+(c-2)] and histo[bin_a*num_bins_b+c]>histo[(bin_a-1)*num_bins_b+(c+2)] and histo[bin_a*num_bins_b+c]>histo[(bin_a+1)*num_bins_b+(c-2)] and histo[bin_a*num_bins_b+c]>histo[(bin_a+1)*num_bins_b+(c+2)])
								reco_retta << "Evento " << event_ID << "\tMassimo " << histo[bin_a*num_bins_b+c] << "\tCoeffAng " << vettore_coeff_ang[bin_a] << "\toffset " << vettore_di_bin[c] << std::endl;
											
						}

		}
	}
}



void histo_fill_cpu(TH2F *istogramma,  int array_a[num_bins_a][num_bins_b], int size_xy, matrix m, std::vector<float> vettore_di_bin)
{
	// RIEMPIO ISTOGRAMMA DELLE FREQUENZE CON I COEFICIENTI ANGOLARI TROVATI
	for(int i=0; i<num_bins_b; i++)
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

void histo_fill_cpu_array(unsigned int *histo, matrix &m, int size_xy, float* vettore_di_bin, float* vettore_coeff_ang)
{
	for(int bin_offset=0; bin_offset<num_bins_b; bin_offset++)
	{
		for(int hit=0; hit<size_xy; hit++)
		{
			for(int bin_a=0; bin_a<num_bins_a; bin_a++)
			{
				if(m.e[hit*m.w+bin_offset] < vettore_coeff_ang[bin_a])
				{
					histo[bin_a*num_bins_a+bin_offset]++;
					break;
				}
			}
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
	float *vettore_di_bin, *vettore_coeff_ang;
	vettore_di_bin=(float*)malloc(num_bins_b*sizeof(float));
	vettore_coeff_ang=(float*)malloc(num_bins_a*sizeof(float));
	for(int i=0; i<num_bins_b; i++)
	{
		vettore_di_bin[i] = -max_offset + i*(2*max_offset/num_bins_b);
	}
	for(int j=0; j<num_bins_a; j++)
	{
		vettore_coeff_ang[j] = -max_coeff_ang + j*(2*max_coeff_ang/num_bins_a);
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
			// LEGGO IL PAYLOAD E RIEMPIO UN VETTORE EVENT(EVENT_ID, LAYER, CHANNEL_ID)
			fragment(data, dimensioni_tracker, num_pixel, vec_eventi);
			const int size_xy = vec_eventi.size();

			std::cout << "frammento ricostruito correttamente" << std::endl;
			// INIZIALIZZO E RIEMPIO VETTORI DA MANDARE A GPU E CPU
			int size = size_xy*sizeof(int);
			float *vec_x, *vec_y;
			vec_x = (float *)malloc(size);
			vec_y = (float *)malloc(size);


			// PER CIASCUN EVENTO ESTRAGGO LE COORDINATE (X, Z)
			coord(vec_x, vec_y, size_xy, vec_eventi, num_pixel, dimensioni_tracker, dist_interazione, dist_piani);

			// CREO MATRICE I CUI ELEMENTI SONO I COEFFICIENTI ANGOLARI DELLE RETTE
			// CHE CONNETTONO CIASCUN BIN ALLE COORDINATE (X,Y) DELL'HIT
			matrix m;
			m.w = num_bins_b;
			m.h = size_xy;
			m.e = (float*) malloc(m.w*m.h*sizeof(float));

			// Partenza cronometro
			auto start_hough = std::chrono::steady_clock::now();
			// ### RIEMPIMENTO MATRICE USANDO CPU ###
			hough_cpu(m, vettore_di_bin, vec_x, vec_y);
			auto end_hough = std::chrono::steady_clock::now();
			auto t_hough = std::chrono::duration_cast<std::chrono::microseconds> (end_hough-start_hough).count();
			std::cout << "Hough matrix calculation time: " << t_hough << " us" << std::endl;
			// PREPRARO ARRAY BIDIMENSIONALE IL CUI CONTENUTO SONO LE FREQUENZE DELL'ISTOGRAMMA
			unsigned int *histo;
			histo = (unsigned int*) malloc(num_bins_a*num_bins_b*sizeof(unsigned int));

			auto start_fill = std::chrono::steady_clock::now();
			histo_fill_cpu_array(histo, m, size_xy, vettore_di_bin, vettore_coeff_ang);
			auto end_fill = std::chrono::steady_clock::now();
			auto t_fill = std::chrono::duration_cast<std::chrono::milliseconds> (end_fill-start_fill).count();
			std::cout << "Histogram Fill Time: " << t_fill << " ms" << std::endl;

			uint32_t event_ID = vec_eventi[1].event_ID();
			std::cout << "Ora parte ricerca massimi: " << std::endl;
			

			auto start_max = std::chrono::steady_clock::now();
			//ricerca_max_cpu_array(histo, vettore_di_bin, vettore_coeff_ang, event_ID, threshold, reco_retta);
			ricerca_max_locale_cpu(histo, vettore_di_bin, vettore_coeff_ang, event_ID, threshold, reco_retta);

			auto end_max = std::chrono::steady_clock::now();
			auto time_max = std::chrono::duration_cast<std::chrono::microseconds> (end_max-start_max).count();
			std::cout << "Time to fine local maxima: " << time_max << " us" << std::endl;

			std::cout << "Numero di Hit: " << size_xy << std::endl;
			//std::cout << "fine  " << event_ID << std::endl;
			
			memset(vec_x, 0, sizeof(vec_x));
			memset(vec_y, 0, sizeof(vec_y));
			memset(histo, 0, sizeof(histo));		
		}

	

	}	



	std::cout << "finito!" << std::endl;
	return 0;
}
