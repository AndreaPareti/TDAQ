#include<stdint.h>
#include<iostream>
#include<vector>
#include<string>
#include<stdio.h>

class event{
public:
	event(uint32_t word, int event_ID);
	//event(int event_ID, uint32_t channel_ID, uint32_t layer);
	~event(){};
	uint32_t channel_ID() {return m_word & 0x000FFFFF;}
	uint32_t layer() {return m_word>>24 & 0x000000FF;}	
	int event_ID() {return m_eventid;}
	uint32_t encoder();
	uint32_t decoder();
	friend std::ostream &operator<<(std::ostream &stream, event &evento);

private:
	uint32_t m_word;
	int m_eventid;
	//uint32_t m_channel_id;
	//uint32_t m_layer;
};


class word{
public:
	word(int layer, uint32_t channel_ID);
	~word(){};
	uint32_t encoder();
	uint32_t decoder();
	friend std::ostream &operator<<(std::ostream &stream, word &word);


private:
	uint32_t m_channel_id;
	int m_layer;

};

class histo{
public:
	histo(int n_channel_a, int n_channel_b, float min_a, float min_b, float max_a, float max_b);
	void dump(); 		// stampa contenuto canale
	int entries();
	bool fill(float val);
	float chValue(int channel);


private:
	std::vector<int> m_channel_a;
	std::vector<int> m_channel_b;
	//std::string m_name;
	float m_inf_a;
	float m_inf_b;
	float m_sup_a;
	float m_sup_b;
	float m_binsize_a;
	float m_binsize_b;

};

struct matrix
{
	unsigned long int w;		//Width
	unsigned long int h;		//Height
	float* e;	//element
	void fill();
};

