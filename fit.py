import pandas as pd
import numpy as np
import sys

data_mc=pd.read_csv('data_retta.dat', sep='\t', header=None, names=['event_ID', 'm', 'q'])
data_reco=pd.read_csv('reco_retta.dat', sep='\t', header=None, names=['reco_eventID', 'reco_max', 'reco_m', 'reco_q'])


num_events_mc=data_mc['event_ID'].max()
num_events_reco=data_reco['reco_eventID'].max()
print(num_events_mc)
print(num_events_reco)
'''
if(num_events_mc != num_events_reco):
	print("Il numero di eventi ricostruiti non corrisponde al numero di eventi generati! ")
	sys.exit()
'''

missed=0
counter_reco=[]
tot_eventi=len(data_mc)
for i in range(1,num_events_mc+1):								# seleziono gli eventi in base all'eventID
	count=0

	print("Lettura evento %d: " %i)
	df_mc=pd.DataFrame(data_mc[data_mc['event_ID']==i])					
	df_reco=pd.DataFrame(data_reco[data_reco['reco_eventID']==i])
	df_mc=df_mc.sort_values(by='m')												# li ordino per comodita'
	df_reco=df_reco.sort_values(by='reco_m')
	check=[]
	mc_m=df_mc['m']														
	mc_q=df_mc['q']
	reco_m=df_reco['reco_m']
	reco_q=df_reco['reco_q']										# loop sulle tracce generate
	for track in df_mc.itertuples(): 							# track e' un oggetto di classe tuple
		track=pd.DataFrame(data=track).T							# conversione di track in pandas dataframe
		track.columns=['index', 'event_ID', 'm', 'q']
		m = track['m'].iat[0]
		q = track['q'].iat[0]										
		if( reco_m.between(m-0.1, m+0.1).any() and reco_q.between(q-0.8, q+0.8).any() ):
			check.append(True)										# se esiste una traccia ricostruita con m e q 
		else:																# abbastanza vicini a quelli generati 
			check.append(False)										# ritorno True, altrimenti False
	npcheck=np.array(check)
	print(npcheck)
	if( np.all(npcheck) ):
		print( "Tutte le tracce sono state ricostruite correttamente")
	else:
		print( "Ci sono ", np.size(npcheck) - np.count_nonzero(npcheck), " tracce non ricostruite")
	missed+=( np.size(npcheck) - np.count_nonzero(npcheck) )




	df_reco=df_reco.sort_values('reco_max', ascending=False).drop_duplicates('reco_m').sort_index()
	counter_reco.append(df_reco)
	#print("Numero di eventi generati: ", len(df_mc), "\tNumero di eventi ricostruiti: ", len(df_reco))
	
	df_mc=df_mc.sort_values(by='m')	
	df_reco=df_reco.sort_values(by='reco_m')
	#print(df_mc)
	#print(df_reco)
	

print("Totale eventi generati: ", tot_eventi)
print("numero eventi non ricostruiti: ", missed)
print("Numero totale di eventi ricostruiti: ", len(data_reco))
	

