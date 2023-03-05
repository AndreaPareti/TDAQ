#include<iostream>
#include<cstdio>
#include<stdlib.h>
#include<iomanip>
#include<stdint.h>
#include<inttypes.h>
#include<string.h>
#include<vector>
#include<fstream>
#include"header.h"
const unsigned int num_bins_a = 1024;	//dump a 1445
const unsigned int num_bins_b = 1024;	//dump a 1445
const unsigned int threads_side = 31; 	// 32x32=1024 masssimo di thread in una block
#include<cuda.h>
#include<cuda_runtime.h>

#define gpuErrchk(ans) {gpuAssert( (ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
	if(code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
	}
}


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
	//if(event_ID==0x00ff00ff){return;}
	event_ID = int(event_ID);
	std::cout << "\nEvento_numero: " << std::dec << event_ID << std::endl;
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

__global__ void hough_gpu(matrix d_m, float* d_vec_bin, float* d_vec_x, float* d_vec_y, int size_xy)
{

	int index = blockIdx.x*blockDim.x + threadIdx.x;
	for(int k=0; k<size_xy; k++)
	{
		float elm = (d_vec_x[k]-d_vec_bin[index])/d_vec_y[k];
		d_m.e[k*d_m.w + index] = elm;
	}
} 


__global__ void hough_gpu_new(matrix d_m, float* d_vec_bin, float* d_vec_x, float* d_vec_y, int size_xy)
{
	int block_id = blockIdx.x;
	int index = threadIdx.x;
	float elm = (d_vec_x[block_id]-d_vec_bin[index])/d_vec_y[block_id];
	d_m.e[block_id*d_m.w + index] = elm;

} 

void hough_cpu(matrix &m, float* vettore_di_bin, float *vec_x, float *vec_y)
{
	for(int k=0; k<m.h; k++) // i indicatore della riga
		for(int j=0; j<m.w; j++) // j indicarore colonna
		{
			if(vec_x[k]==0){;}
			else m.e[k*m.w + j] = (vec_x[k]-vettore_di_bin[j])/vec_y[k];
			//std::cout << m.e[k*m.w + j] << " boh " ;
			//std:: cout << (vec_x[k]-vettore_di_bin[j])/vec_y[k] << " ";
		}
}



__global__ void ricerca_max_gpu(unsigned int *d_histo, float* d_vec_bin, float* d_vec_coeff, int event_ID, int threshold)
{
	int bin_a=threadIdx.x+blockIdx.x*blockDim.x;
	__shared__ int array_max[num_bins_a];
	int location=1;
	__shared__ int loc[num_bins_a];
	// cerco massimo lungo colonna bin offset -> Salvo valore massimo in array_max e bin coef. ang. in loc[]
	for(int c=0; c<num_bins_b; c++)
	{
		if(d_histo[bin_a*num_bins_b+c] > d_histo[bin_a*num_bins_b+location])
		{
			location=c;
		}
		loc[bin_a]=location;
		array_max[bin_a] = d_histo[bin_a*num_bins_b+location];
	}
	__syncthreads();

	if(bin_a>1 and bin_a<num_bins_a-1 and array_max[bin_a]>threshold and array_max[bin_a-1]<array_max[bin_a] and array_max[bin_a]>array_max[bin_a+1]){

                printf("Evento %d \t Massimo %d \tCoeffAng %.4f \toffset %.4f\n", event_ID, array_max[bin_a], d_vec_coeff[bin_a], d_vec_bin[loc[bin_a]]);
                }

}
__global__ void ricerca_max_gpu_locale2(unsigned int *d_histo, float *d_vec_bin, float* d_vec_coeff, int event_ID, int threshold)
{
	int bin_b=threadIdx.x + blockIdx.x*blockDim.x;
	int bin_a=threadIdx.y + blockIdx.y*blockDim.y;
	unsigned int channelID = bin_a*num_bins_b + bin_b;

	unsigned int dx = bin_a*num_bins_b+(bin_b+1);
	unsigned int sx = bin_a*num_bins_b+(bin_b-1);
	unsigned int up = (bin_a+1)*num_bins_b+bin_b;
	unsigned int down = (bin_a-1)*num_bins_b+bin_b;
	
	unsigned int dx2 = bin_a*num_bins_b+(bin_b+2);
	unsigned int sx2 = bin_a*num_bins_b+(bin_b-2);
	unsigned int up2 = (bin_a+2)*num_bins_b+bin_b;
	unsigned int down2 = (bin_a-2)*num_bins_b+bin_b;
	
	unsigned int ud = (bin_a+1)*num_bins_b+(bin_b+1);
	unsigned int us = (bin_a+1)*num_bins_b+(bin_b-1);
	unsigned int dd = (bin_a-1)*num_bins_b+(bin_b+1);
	unsigned int ds = (bin_a-1)*num_bins_b+(bin_b-1);
	// CASO GENERALE, VENGONO ESCLUSE LE THREAD AI BORDI DELLA BLOCK
	if(threadIdx.x>1 and threadIdx.x<threads_side-2 and threadIdx.y>1 and threadIdx.y<threads_side-2)
	{
		if(d_histo[channelID]>threshold and d_histo[channelID]>=d_histo[up] and d_histo[channelID]>=d_histo[down] and d_histo[channelID]>=d_histo[dx] and d_histo[channelID]>=d_histo[sx] and d_histo[channelID]>d_histo[up2] and d_histo[channelID]>d_histo[down2] and d_histo[channelID]>d_histo[dx2] and d_histo[channelID]>d_histo[sx2])
			{
				if(d_histo[channelID]>d_histo[ud] and d_histo[channelID]>d_histo[us] and d_histo[channelID]>d_histo[dd] and d_histo[channelID]>d_histo[ds])
					printf("Evento %d \tMassimo %d \tCoeffAng %.3f \toffset %.3f\n", event_ID, d_histo[channelID], d_vec_coeff[bin_a], d_vec_bin[bin_b]);
			}
	}
	
	// CONSIDERO LE THREAD LUNGO I BORDI DELLA BLOCK
	else if(threadIdx.x == 0 and threadIdx.y>1 and threadIdx.y<threads_side-2)
	{
		if(d_histo[channelID]>threshold and d_histo[channelID]>d_histo[up] and d_histo[channelID]>d_histo[up2] and d_histo[channelID]>d_histo[down] and d_histo[channelID]>d_histo[down2] and d_histo[channelID]>d_histo[dx] and d_histo[channelID]>d_histo[dx2])
				printf("Evento %d \tMassimo %d \tCoeffAng %.3f \toffset %.3f\n", event_ID, d_histo[channelID], d_vec_coeff[bin_b], d_vec_bin[bin_a]);
	}

	else if(threadIdx.x == threads_side and threadIdx.y>1 and threadIdx.y<threads_side-2)
	{
		if(d_histo[channelID]>threshold and d_histo[channelID]>d_histo[up] and d_histo[channelID]>d_histo[up2] and d_histo[channelID]>d_histo[down] and d_histo[channelID]>d_histo[down2] and d_histo[channelID]>d_histo[sx] and d_histo[channelID]>d_histo[sx2])
				printf("Evento %d \tMassimo %d \tCoeffAng %.3f \toffset %.3f\n", event_ID, d_histo[channelID], d_vec_coeff[bin_b], d_vec_bin[bin_a]);			
	}


	else if(threadIdx.y == 0 and threadIdx.x>1 and threadIdx.x<threads_side-2)
	{
		if(d_histo[channelID]>threshold and d_histo[channelID]>d_histo[up] and d_histo[channelID]>d_histo[up2] and d_histo[channelID]>d_histo[sx] and d_histo[channelID]>d_histo[sx2] and d_histo[channelID]>d_histo[dx] and d_histo[channelID]>d_histo[dx2])
				printf("Evento %d \tMassimo %d \tCoeffAng %.3f \toffset %.3f\n", event_ID, d_histo[channelID], d_vec_coeff[bin_b], d_vec_bin[bin_a]);
	}

	else if(threadIdx.y ==threads_side and threadIdx.x>1 and threadIdx.x<threads_side)
	{
		if(d_histo[channelID]>threshold and d_histo[channelID]>d_histo[down] and d_histo[channelID]>d_histo[down2] and d_histo[channelID]>d_histo[dx] and d_histo[channelID]>d_histo[dx2] and d_histo[channelID]>d_histo[sx] and d_histo[sx2])
				printf("Evento %d \tMassimo %d \tCoeffAng %.3f \toffset %.3f\n", event_ID, d_histo[channelID], d_vec_coeff[bin_b], d_vec_bin[bin_a]);
	}
	
	else if(threadIdx.x==1 and threadIdx.y>1 and threadIdx.y<threads_side-2)
	{
		if(d_histo[channelID] > threshold and d_histo[channelID]>=d_histo[sx] and d_histo[channelID]>=d_histo[up] and d_histo[channelID]>=d_histo[down] and d_histo[channelID]>=d_histo[down] and d_histo[channelID]>d_histo[up2] and d_histo[channelID]>d_histo[down2] and d_histo[channelID]>d_histo[dx2] and d_histo[channelID]>d_histo[ud] and d_histo[channelID]>d_histo[us] and d_histo[channelID]>d_histo[dd] and d_histo[channelID]>d_histo[ds])
				printf("Evento %d \tMassimo %d \tCoeffAng %.3f \toffset %.3f\n", event_ID, d_histo[channelID], d_vec_coeff[bin_b], d_vec_bin[bin_a]);
	}
	
	else if(threadIdx.x == threads_side-1 and threadIdx.y>1 and threadIdx.y<threads_side-2)
		if(d_histo[channelID]>threshold and d_histo[channelID]>=d_histo[up] and d_histo[channelID]>=d_histo[down] and d_histo[channelID]>=d_histo[sx]  and d_histo[channelID]>=d_histo[dx] and d_histo[channelID]>d_histo[up2] and d_histo[channelID]>d_histo[down2] and d_histo[channelID]>d_histo[sx2] and d_histo[channelID]>d_histo[ud] and d_histo[channelID]>d_histo[us] and d_histo[channelID]>d_histo[dd] and d_histo[channelID]>d_histo[ds])
				printf("Evento %d \tMassimo %d \tCoeffAng %.3f \toffset %.3f\n", event_ID, d_histo[channelID], d_vec_coeff[bin_b], d_vec_bin[bin_a]);			
	
	else if(threadIdx.y==1 and threadIdx.x>1 and threadIdx.y<threads_side-2)
		if(d_histo[channelID]>threshold and d_histo[channelID]>=d_histo[up] and d_histo[channelID]>=d_histo[dx] and d_histo[channelID]>=d_histo[sx] and d_histo[channelID]>=d_histo[down] and d_histo[channelID]>d_histo[up2] and d_histo[channelID]>d_histo[dx2] and d_histo[channelID]>d_histo[sx2] and d_histo[channelID]>d_histo[ud] and d_histo[channelID]>d_histo[us] and d_histo[channelID]>d_histo[dd] and d_histo[channelID]>d_histo[ds])
				printf("Evento %d \tMassimo %d \tCoeffAng %.3f \toffset %.3f\n", event_ID, d_histo[channelID], d_vec_coeff[bin_b], d_vec_bin[bin_a]);
				
	else if(threadIdx.y==threads_side-1 and threadIdx.x>1 and threadIdx.y<threads_side-2)
		if(d_histo[channelID]>threshold and d_histo[channelID]>=d_histo[up] and d_histo[channelID]>=d_histo[dx] and d_histo[channelID]>=d_histo[sx] and d_histo[channelID]>=d_histo[down] and d_histo[channelID]>d_histo[dx2] and d_histo[channelID]>d_histo[sx2] and d_histo[channelID]>d_histo[down2] and 
d_histo[channelID]>d_histo[ud] and d_histo[channelID]>d_histo[us] and d_histo[channelID]>d_histo[dd] and d_histo[channelID]>d_histo[ds])	
				printf("Evento %d \tMassimo %d \tCoeffAng %.3f \toffset %.3f\n", event_ID, d_histo[channelID], d_vec_coeff[bin_b], d_vec_bin[bin_a]);			
	
	
	else
		if(d_histo[channelID]>threshold)
			printf("Ciaone mi hai mancato\n");
	
}	



__global__ void histo_fill_gpu(matrix d_m, int num_bins_a, int size_xy, unsigned int *histo, float *d_vec_bin, float *d_vec_coeff)
{
	int index = threadIdx.x + blockIdx.x*blockDim.x;
	for(int i=0; i<size_xy; i++)
	{
		for(int j=0; j<num_bins_a; j++)
		{
			if(d_m.e[i*d_m.w+index] < d_vec_coeff[j])
			{
				atomicAdd(&histo[j*num_bins_a + index], 1);
				break;
			}
		}
	
	}
}

__global__ void prova_fill(matrix d_m, unsigned int *histo, float *d_vec_bin, float *d_vec_coeff)
{
	// Per ogni hit (totale size_xy eventi) uso un block di 1024 thread
	// per ciascun bin offset (max 1024 bin) uso una thread
	// ciascuna thread compara un canale di d_m e riempie frequenza istogramma
	int hit_id = blockIdx.x;
	int bin_offset=threadIdx.x;
	for(int bin_coeff=0; bin_coeff<num_bins_a; bin_coeff++)
	{
		if(d_m.e[hit_id*d_m.w + bin_offset] < d_vec_coeff[bin_coeff])
		{
			atomicAdd(&histo[bin_coeff*num_bins_b + bin_offset], 1);
			break;
		}

	}
}


int main()
{
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
	int num_pixel = int((dimensioni_tracker*1000)/dim_pixel);
	int threshold = int(parametri[11]);

	float max_coeff_ang = 3;

	// INIZIO LETTURA DI FILE DATA
	// PARTO LEGGENDO FULL HEADER
	std::ifstream data("file.dat", std::ios::binary);
	full_header(data);

	// RIEMPIO UN VETTORE I CUI ELEMENTI SONO I BIN DELL'ISTOGRAMMA IN CUI VIENE SUDDIVISO L'OFFSET
	float *vettore_di_bin, *vettore_coeff_ang;
	vettore_di_bin = ( float*) malloc(num_bins_b*sizeof(float) );
	vettore_coeff_ang = (float*) malloc(num_bins_a*sizeof(float) );
	for(int i=0; i<num_bins_b; i++)
	{
		vettore_di_bin[i] = -max_offset + i*(2*max_offset/num_bins_b);
		//std::cout << vettore_di_bin[i] << " " ;
	}
	std::cout << std::endl;
	for(int j=0; j<num_bins_a; j++)
	{
		vettore_coeff_ang[j] = -max_coeff_ang + j*(2*max_coeff_ang/num_bins_a);
		//std::cout << vettore_coeff_ang[j] << " ";
	}
	std::cout << std::endl;
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
			// LEGGO IL PAYLOAD E RIEMPIO UN VETTORE EVENT(EVENT_ID, LAYER, CHANNEL_ID)
			fragment(data, dimensioni_tracker, num_pixel, vec_eventi);
			const int size_xy = vec_eventi.size();
			uint32_t event_ID = vec_eventi[1].event_ID();

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

			// ### RIEMPIMENTO MATRICE USANDO CPU ###
			//hough_cpu(m, vettore_di_bin, vec_x, vec_y);

			// ### RIEMPIMENTO MATRICE USANDO GPU ### 
			// alloco spazio in memoria per matrice d_m
			matrix d_m;
			d_m.w = m.w;
			d_m.h = m.h;
			size_t size_matrix = d_m.w * d_m.h * sizeof(float);
			cudaError_t err = cudaMalloc(&d_m.e, size_matrix);
			cudaError_t err1 = cudaMemcpy(d_m.e, m.e, size_matrix, cudaMemcpyHostToDevice);
			
			cudaEvent_t start_hough, stop_hough, start_histo, stop_histo, start_max, stop_max;
			cudaEventCreate(&start_hough);
			cudaEventCreate(&stop_hough);
			cudaEventCreate(&start_histo);
			cudaEventCreate(&stop_histo);
			cudaEventCreate(&start_max);
			cudaEventCreate(&stop_max);

			float *d_vec_x, *d_vec_y, *d_vec_bin; 

			// Copy Coordinate Vectors from Host do Device
			size_t fragment_size = size_xy*sizeof(float);
			size_t bins_size = num_bins_b*sizeof(float);
			cudaMalloc( (float **)&d_vec_x, fragment_size);
			cudaMalloc( (float **)&d_vec_y, fragment_size);
			cudaMalloc( (float **)&d_vec_bin, bins_size);
			cudaMemcpy(d_vec_x, vec_x, fragment_size, cudaMemcpyHostToDevice);
			cudaMemcpy(d_vec_y, vec_y, fragment_size, cudaMemcpyHostToDevice); 
			cudaMemcpy(d_vec_bin, vettore_di_bin, bins_size, cudaMemcpyHostToDevice);
			//partenza misura tempo 
			cudaEventRecord(start_hough, 0);
			// launch hough algorith
			//hough_gpu<<<num_bins_b/64,64>>>(d_m, d_vec_bin, d_vec_x, d_vec_y, size_xy);
			//hough_gpu<<<1, num_bins_b>>>(d_m, d_vec_bin, d_vec_x, d_vec_y, size_xy);
			hough_gpu_new<<<size_xy, num_bins_b>>>(d_m, d_vec_bin, d_vec_x, d_vec_y, size_xy);
			cudaError_t KernelError=cudaGetLastError();
			cudaDeviceSynchronize();
			float t_hough;
			cudaEventRecord(stop_hough);
			cudaEventSynchronize(stop_hough);

			cudaEventElapsedTime(&t_hough, start_hough, stop_hough);
			// Copy matrix back to m
			cudaMemcpy(m.e, d_m.e, size_matrix, cudaMemcpyDeviceToHost);


			cudaFree(d_vec_x);
			cudaFree(d_vec_y);
			// PREPRARO ARRAY BIDIMENSIONALE IL CUI CONTENUTO SONO LE FREQUENZE DELL'ISTOGRAMMA
			
			//histo_fill_cpu(istogramma,  array_a, size_xy, m, vettore_di_bin);
			
			// Allocazione spazio memoria vettore di bin coefficiete angolare
			float *d_vec_coeff;
			gpuErrchk(cudaMalloc((float**)&d_vec_coeff, num_bins_a*sizeof(float)));
			gpuErrchk(cudaMemcpy(d_vec_coeff, vettore_coeff_ang,num_bins_a*sizeof(float), cudaMemcpyHostToDevice));
			// Allocazione spazio di memoria per istogramma
			unsigned int *d_histo;
			unsigned int *histo;
			histo = (unsigned int*) malloc(num_bins_a*num_bins_b*sizeof(unsigned int) ); 
			gpuErrchk(cudaMalloc((unsigned int**)&d_histo, num_bins_a*num_bins_b*sizeof(unsigned int)));
			cudaEventRecord(start_histo);

			const int num_hits=size_xy;
			dim3 dimGrid(num_hits,1,1);
			dim3 dimBlock(num_bins_b,1, 1);
			prova_fill<<<dimGrid,dimBlock>>>(d_m, d_histo, d_vec_bin, d_vec_coeff);
			//  Riempimento istogramma  usando GPU
			//histo_fill_gpu<<<num_bins_b/128, 128>>>(d_m, num_bins_a, size_xy, d_histo, d_vec_bin, d_vec_coeff);
			cudaDeviceSynchronize();
			cudaEventRecord(stop_histo);
			float t_histo;
			cudaEventSynchronize(stop_histo);
			cudaEventElapsedTime(&t_histo, start_histo, stop_histo);

			gpuErrchk(cudaMemcpy(histo, d_histo, num_bins_a*num_bins_b*sizeof(unsigned int), cudaMemcpyDeviceToHost));
			cudaFree(d_m.e);
			
			cudaEventRecord(start_max);			
			//ricerca_max_gpu<<<num_bins_a/(num_bins_a/2),(num_bins_a/2) >>>(d_histo, d_vec_bin, d_vec_coeff, event_ID, threshold);

			int lato_a_grid = num_bins_a/threads_side;
			int lato_b_grid = num_bins_b/threads_side;
			dim3 dimBlock_histo(threads_side, threads_side, 1);
			dim3 dimGrid_histo(lato_a_grid, lato_b_grid,1);
			ricerca_max_gpu_locale2<<<dimGrid_histo, dimBlock_histo>>>(d_histo, d_vec_bin, d_vec_coeff, event_ID, threshold);

			cudaDeviceSynchronize();

			cudaFree(d_vec_bin);
			cudaFree(d_histo);
			cudaFree(d_vec_coeff);
			//cudaFree(event_ID);
			//cudaFree(threshold); 
			
			cudaEventRecord(stop_max);
			float t_max;
			cudaEventSynchronize(stop_max);
			cudaEventElapsedTime(&t_max, start_max, stop_max);
			printf("Time Hough: %.6f ms\n", t_hough);
			printf("Time Histo: %.6f ms\n", t_histo);
			printf("Time Max: %.6f ms\n", t_max);


			std::cout << "fine Evento" << event_ID << std::endl;
			std::cout << "Size XY = " << size_xy << std::endl;
		}



	}	


	std::cout << "finito!" << std::endl;
	return 0;
} 
