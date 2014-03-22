import numpy as np

def make_name(SN_Array):

	min_phase,min_red,min_dm15,min_mb,min_bminusv= 10000000
	max_phase,max_red,max_dm15,max_mb,max_bminusv= -10000000
	phase=[]
	red=[]
	dm15=[]
	mb=[]
	bminusv=[]

	for SN in SN_Array:
		if ((SN.phase<min_phase) & (np.isnan(SN.phase)==False)):
			min_phase=SN.phase
		if ((SN.phase>max_phase) & (np.isnan(SN.phase)==False)):
			max_phase=SN.phase
		if ((SN.redshift<min_red) & (np.isnan(SN.redshift)==False)):
			min_red=SN.redshift
		if ((SN.redshift>max_red) & (np.isnan(SN.redshift)==False)):
			max_red=SN.redshift
		if ((SN.dm15<min_dm15) & (np.isnan(SN.dm15)==False)):
			min_dm15=SN.dm15
		if ((SN.dm15>max_dm15) & (np.isnan(SN.dm15)==False)):
			max_dm15=SN.dm15
		if ((SN.m_b<min_mb) & (np.isnan(SN.m_b)==False)):
			min_mb=SN.m_b
		if ((SN.m_b>max_mb) & (np.isnan(SN.m_b)==False)):
			max_mb=SN.m_b
		if ((SN.B_minus_v<min_bminusv) & (np.isnan(SN.B_minus_v)==False)):
			min_bminusv=SN.B_minus_v
		if ((SN.B_minus_v>max_bminusv) & (np.isnan(SN.B_minus_v)==False)):
			max_bminusv=SN.B_minus_v
		
		if(np.isnan(SN.phase)==False):
			phase.append(SN.phase)
		if(np.isnan(SN.redshift)==False):
			red.append(SN.redshift)
		if(np.isnan(SN.dm15)==False):
			dm15.append(SN.dm15)
		if(np.isnan(SN.m_b)==False):
			mb.append(SN.m_b)
		if(np.isnan(SN.B_minus_v)==False):
			bminusv.append(SN.B_minus_v)

	avg_phase=np.average(phase)
	avg_red=np.average(red)
	avg_dm15=np.average(dm15)
	avg_mb=np.average(mb)
	avg_bminusv=np.average(bminusv)

	f_name1='composite,'+min_phase+'.'+max_phase+'.'+min_red+'.'+max_red+'.'+min_dm15+'.'+max_dm15+'.'+min_mb+'.'+max_mb+'.'+min_bminusv+'.'+max_bminusv
	f_name2='--'+avg_phase+'.'+avg_red+'.'+avg_dm15+'.'+avg_mb+'.'+avg_bminusv+'--'+len(SN_Array)+'SN'
	f_name=f_name1+f_name2
	return f_name


