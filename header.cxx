#include"header.h"
#include<iostream>
#include<stdlib.h>
#include<stdint.h>
#include<inttypes.h>

#include<algorithm>
#include<thread>
#include<fstream>
#include<typeinfo>





event::event(uint32_t word, int event_ID)
{
	m_word = word;
	m_eventid = event_ID;
}

std::ostream &operator<<(std::ostream &stream, event &evento)
{
	stream << "[" << evento.event_ID() << "," << evento.layer() << "," << evento.channel_ID() << "]";
	return stream;
}

word::word(int layer, uint32_t channel_ID)
{
	m_channel_id = channel_ID;
	m_layer = layer;
	uint32_t free_space;
}
	
uint32_t word::encoder()
{
	uint32_t free_space = 0x00000000;
	uint32_t encoded_word = uint32_t(m_layer);
	encoded_word <<= 4;
	encoded_word |= free_space;
	encoded_word <<= 20;
	encoded_word |= m_channel_id;
	return	encoded_word;
}
 
std::ostream &operator<<(std::ostream &stream, word &parola)
{
	stream << parola.encoder();
	return stream;
}

histo::histo(int num_bin_a, int num_bin_b, float min_a, float min_b, float max_a, float max_b)
{
	m_channel_a.assign(num_bin_a, 0);
	m_channel_b.assign(num_bin_b, 0);

	m_inf_a = min_a;
	m_inf_b = min_b;
	m_sup_a = max_a;
	m_sup_b = max_b;

	m_binsize_a = (m_sup_a - m_inf_a)/num_bin_a;
	m_binsize_b = (m_sup_b - m_inf_b)/num_bin_b;
}

void histo::dump()
{
	for(int i=0; i<m_channel_a.size(); i++)
	{
		float ch_inf_a = m_inf_a + m_binsize_a*i;
		float ch_sup_a = ch_inf_a + m_binsize_a;

		for(int j=0; j<m_channel_b.size(); j++)
		{
			float ch_inf_b = m_inf_b + m_binsize_b*j;
			float ch_sup_b = ch_inf_b + m_binsize_b;
			std::cout << "[" << i << j << "]\t" << ch_inf_b << "-" << ch_sup_b << "\t" << m_channel_b[j] << m_channel_a[i] << std::endl;
		}
	}
}




