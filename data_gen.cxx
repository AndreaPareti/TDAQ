#include<iostream>
#include<stdlib.h>
#include<stdint.h>
#include<inttypes.h>
#include<TMath.h>
#include<TRandom3.h>
#include<iomanip>
#include<TH1F.h>
#include<TApplication.h>
#include<TCanvas.h>
#include<vector>
#include<TH2F.h>
#include<algorithm>
#include<chrono>
#include<thread>
#include<fstream>
#include"header.h"

void eventi(int media_sig, int &num_sig)
{
	TRandom3* rnd_sig = new TRandom3();
	rnd_sig->SetSeed();
	num_sig = rnd_sig->Gaus(media_sig, 1);
	if(num_sig<1){num_sig=1;}
	if(num_sig<0){num_sig=0;}
}

void retta(float max_offset, float dimensione, float dist_interazione, float &q, float &m, int num_piani, float dist_piani)
{
	//generazione offset punto di produzione delle particelle,
	//centrato attorno a (0,0) e distrib. gaussiana
	TRandom3* rnd_q = new TRandom3();
	rnd_q->SetSeed();
	q=rnd_q->Uniform(-max_offset/2, max_offset/2);
	if(q<(-max_offset)){q=-max_offset;}
	else if(q>max_offset){q=max_offset;}

	//generazione coefficiente angolare della traccia
	//distribuzione uniforme tale che avvenga almeno una hit sul primo piano
	TRandom3* rnd_m = new TRandom3();
	rnd_m->SetSeed();
	//m = rnd_m->Uniform((-dimensione/(3*2*dist_interazione)), (dimensione/(3*2*dist_interazione)));
	//m = rnd_m->Uniform((-dimensione/(dist_interazione)), (dimensione/(dist_interazione)));
	//m=rnd_m->Gaus(0, (dimensione/2*dist_interazione));	
	//m = rnd_m->Gaus(0, 1);	
	//m = rnd_m->Uniform(-2.8, 2.8);

	// estremi piano centrale
	float coord_piano_centrale = dist_interazione + (num_piani/2)*dist_piani;
	float estremo_piano = sqrt( (coord_piano_centrale*coord_piano_centrale) + (dimensione*dimensione) );
	float limite = atan(estremo_piano/coord_piano_centrale);
	m = rnd_m->Uniform(-limite*2, limite*2);



}

void smearing(uint32_t channel_ID, int layer, float coef_ang,  std::ofstream &out)
{
	// PRODUZIONE DI SEGNALI DA PIXEL ADIACENTI A QUELLO COLPITO: SE LA TRACCIA HA COEFFICIENTE ANGOLARE
	// MOLTO GRANDE (IN MODULO) GENERO SEGNALE ANCHE DAL PIXEL A DX (O SX SE COEFF. E' NEGATIVO)
	if(coef_ang>0.50)
	{
		uint32_t smearing_dx = channel_ID+1;
		word parola(layer, smearing_dx);
		uint32_t dx = parola.encoder();
		out.write((char* ) &dx, sizeof(dx));
	}
	else if(coef_ang < -0.50)
	{
		uint32_t smearing_sx = channel_ID-1;
		word parola_sx(layer, smearing_sx);
		uint32_t sx = parola_sx.encoder();
		out.write((char* ) &sx, sizeof(sx));	
	}
}

void hit(int piano, float coef_ang, float offset, float dist_interazione, float dist_piani, float dimensioni_tracker, int num_pixel, float &x, uint32_t &channel_ID, uint32_t &layer, std::ofstream &out)
{
	float dist = dist_interazione + (piano-1)*dist_piani;
	x = coef_ang*dist + offset;

	if(x < (dimensioni_tracker/2) and x > (-dimensioni_tracker/2))
	{	
		// STAMPO SU FILE IL CANALE DEL PIXEL CORRISPONDENTE ALLA COORDINATA x
		std::cout << "x = " << x << "\tdist = " << dist << std::endl;
		channel_ID = (x+dimensioni_tracker/2)*(num_pixel/dimensioni_tracker);
		// 4 BIT LIBERI: USABILI PER AGGIUNGERE CANALI AL TRACKER O PER IDENTIFICAZIONE 
		// DEL TIPO DI TRACKER (PIXEL, STRAWTUBE, MICROSTRIP ECC.)
		word parola(piano, channel_ID);
		uint32_t new_w = parola.encoder();
		out.write((char*) &new_w, sizeof(new_w));
		// AGGIUNGO SMEARING: INSERISCO SEGNALI SUI CANALI channel_ID+-1
		smearing(channel_ID, piano, coef_ang, out);	
	}
}	

void background(int media_bkg, int piano, float dimensioni_tracker, uint32_t &channel_ID_bkg, int num_pixel, std::ofstream &out)
{
	TRandom3* rnd_num_bkg = new TRandom3();
	rnd_num_bkg->SetSeed();
	int num_bkg = rnd_num_bkg->Gaus(media_bkg, 0);
	if(num_bkg<0){num_bkg=0; std::cout << "Nessun evento fondo" << std::endl;}
	for(int j=0; j<num_bkg; j++)
	{
		uint32_t channel_ID_bkg;
		TRandom3* rnd_bkg = new TRandom3();
		rnd_bkg->SetSeed();
		float coord_bkg = rnd_bkg->Uniform((-dimensioni_tracker/2), dimensioni_tracker/2);
		channel_ID_bkg = (coord_bkg+dimensioni_tracker/2)*(num_pixel/dimensioni_tracker);
		word bkg(piano, channel_ID_bkg);
		uint32_t w_bkg = bkg.encoder();
		out.write((char* ) &w_bkg, sizeof(w_bkg));
	}	
}


int main(int argc, char* argv[])
{	
	// LETTURA FILE DI CONFIGURAZIONE 
	std::string nomi_param;
	std::vector<float> parametri;
	std::ifstream infile("config.txt");
	float par_value;
	while(infile >> nomi_param >> par_value)
	{
		std::cout << nomi_param << "\t" << par_value << std::endl;
		parametri.push_back(par_value);
	}
	infile.close();
	// INSERISCO PARAMETRI DEL RIVELATORE E DI GENERAZIONE SEGNALI
	int num_piani= int(parametri[0]);
	float dist_interazione = parametri[1];
	float dist_piani = parametri[2];
	float dimensioni_tracker = parametri[3];
	float dim_pixel = parametri[4];
	int num_pixel = int((dimensioni_tracker*1000)/dim_pixel);
	std::cout << "numero di pixel: " << num_pixel << std::endl;
	std::cout << "dimensioni pixel: " << (dimensioni_tracker/num_pixel)*1000 << "micron" << std::endl;
	int num_eventi = int(parametri[5]);
	int media_sig = int(parametri[6]);
	int media_bkg = int(parametri[7]);
	int crossing_time = int(parametri[8]);
	int on_off_background = int(parametri[9]);
	float max_offset = parametri[10];
	int count=0;

	// CREAZIONE DEL FILE SU CUI SCRIVERE I DATA
	std::ofstream out("file.dat", std::ios::binary);
	// CREAZIONE FILE SU CUI SALVARE I COEFFICIENTI DELLA RETTA DA RICOSTURIRE
	std::ofstream data_retta("data_retta.dat");
	// QUI INIZIA LA PRODUZIONE DI DATA: RICHIAMO LE FUNZIONI PER LA GENERAZIONE DEI SEGNALI SUI PIANI
	std::cout << "***\n INIZIO GENERAZIONE \n***" << std::endl;
	uint32_t start=0xBABACACA;
	out.write( (char*) &start, sizeof(start));
	uint32_t header_size = 4;
	out.write((char*) &header_size, sizeof(header_size));
	uint32_t data_format=0x000000AA;
	out.write( (char*) &data_format, sizeof(data_format));
	uint32_t checkword=0xDEADCAFE;
	uint32_t source_ID = 0x0FFAACC0;
	out.write((char*) &source_ID, sizeof(source_ID));

	int event_ID=1;
	while(event_ID<=num_eventi)
	{
		std::cout << "*** \tEvento numero " << event_ID << "\t *** " << std::endl;
		out.write( (char*) &checkword, sizeof(checkword));
		out.write( (char*) &event_ID, sizeof(event_ID));
	// PER OGNI EVENTO GENERO UN NUMERO DI SEGNALI CHE SEGUE UNA DISTRIBUIZONE GAUSSIANA
		int num_sig;
		eventi(media_sig, num_sig);
		std::cout << "Numero di eventi di segnale: " << num_sig << std::endl;
		// GENERO num_sig EVENTI E LI SCRIVO SU FILE
		for(int i=1; i<=num_sig; i++)
		{
			float offset;
			float coef_ang;
			// GENERO PARAMETRI m e q PER CIASCUNA DELLE num_sig TRACCE E LI SALVO SU FILE "data_retta.txt"
			retta(max_offset, dimensioni_tracker, dist_interazione, offset, coef_ang, num_piani, dist_piani);
			std::cout << "segnale numero " << i << ":" << "\t" << "(m, q) = " << "(" << coef_ang << "," << offset << ")" << std::endl;
			data_retta << "Event_ID: " << event_ID << "\tm: " << coef_ang << "\t q: " << offset << std::endl;
			// PER CIASCUN PIANO DATA UNA TRACCA CALCOLO COORDINATE HIT -> PIXEL_ID, LAYER_ID E STAMPO
			for(int piano=1; piano<=num_piani; piano++)
			{
				float x;
				uint32_t channel_ID;
				uint32_t layer;
				// PER CIASCUNA TRACCIA CALCOLO IL PUNTO DI INTERAZIONE CON L'i-esimo PIANO TRACKER
				hit(piano, coef_ang, offset, dist_interazione, dist_piani, dimensioni_tracker, num_pixel, x, channel_ID, layer, out);
				// SULL' i-ESIMO PIANO DEL TRACKER AGGIUNGO SEGNALE DI FONDO GENERATO IN MODO UNIFORME
				uint32_t channel_ID_bkg;
				if(on_off_background==1){background(media_bkg, piano, dimensioni_tracker, channel_ID_bkg, num_pixel, out);}
			}
			
		}
	uint32_t fine_evento = 0x00FF00FF;
	out.write((char*) &fine_evento, sizeof(fine_evento));
	event_ID++;	
	}



	// TERMINATA TRASMISSIONE DEI DATI, INSERISCO CHECKWORD DI CHIUSURA FILE	
	uint32_t fine=0x0A11DEAD;
	out.write( (char*) &fine, sizeof(fine));
	out.close();
	data_retta.close();
	return 0;
}



