import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

data=pd.read_csv('tempi_read_data_singlethread.txt', sep='\t', usecols=[0,1,2,3,4], header=0, names=['t_hough', 't_fill', 't_max_locale', 't_max_alt', 'num_hits'])

print(data)

plt.title("Tempo riempimento histogramma")
input1=data['t_fill']
input2=data['num_hits']
plt.plot(input2, input1)
plt.show()
