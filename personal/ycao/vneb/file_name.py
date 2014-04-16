import numpy as np

def make_name(SN_Array):

	min_phase = 10000000
	min_red = 10000000
	min_dm15 = 10000000
	min_mb = 10000000
	min_bminusv = 10000000
	max_phase = -10000000
	max_red = -10000000
	max_dm15 = -10000000
	max_mb = -10000000
	max_bminusv = -10000000
	phase=[]
	red=[]
	dm15=[]
	mb=[]
	bminusv=[]

	for SN in SN_Array:
		if ((SN.phase<min_phase) & (SN.phase !=  None)):
			min_phase=SN.phase
		if ((SN.phase>max_phase) & (SN.phase != None)):
			max_phase=SN.phase
		if ((SN.redshift<min_red) & (SN.redshift != None)):
			min_red=SN.redshift
		if ((SN.redshift>max_red) & (SN.redshift != None)):
			max_red=SN.redshift
		if ((SN.dm15<min_dm15) & (SN.dm15 != None)):
			min_dm15=SN.dm15
		if ((SN.dm15>max_dm15) & (SN.dm15 != None)):
			max_dm15=SN.dm15
		if ((SN.m_b<min_mb) & (SN.m_b != None)):
			min_mb=SN.m_b
		if ((SN.m_b>max_mb) & (SN.m_b != None)):
			max_mb=SN.m_b
		if ((SN.B_minus_v<min_bminusv) & (SN.B_minus_v != None)):
			min_bminusv=SN.B_minus_v
		if ((SN.B_minus_v>max_bminusv) & (SN.B_minus_v != None)):
			max_bminusv=SN.B_minus_v
		
		if(SN.phase != None):
			phase.append(SN.phase)
		if(SN.redshift != None):
			red.append(SN.redshift)
		if(SN.dm15 != None):
			dm15.append(SN.dm15)
		if(SN.m_b != None):
			mb.append(SN.m_b)
		if(SN.B_minus_v != None):
			bminusv.append(SN.B_minus_v)

	avg_phase=np.average(phase)
	avg_red=np.average(red)
	avg_dm15=np.average(dm15)
	avg_mb=np.average(mb)
	avg_bminusv=np.average(bminusv)

	f_name1='composite,'+str(min_phase)+'.'+str(max_phase)+'.'+str(min_red)+'.'+str(max_red)+'.'+str(min_dm15)+'.'+str(max_dm15)+'.'+str(min_mb)+'.'+str(max_mb)+'.'+str(min_bminusv)+'.'+str(max_bminusv)
	f_name2='--'+str(avg_phase)+'.'+str(avg_red)+'.'+str(avg_dm15)+'.'+str(avg_mb)+'.'+str(avg_bminusv)+'--'+str(len(SN_Array))+'SNe'
	f_name=f_name1+f_name2
	return f_name


