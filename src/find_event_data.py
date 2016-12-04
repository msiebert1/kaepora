import composite as comp

def find_data(event):
	event = "'" + event + "'"
	query = "SELECT * FROM Supernovae where SN = " + event
	print query
	SN_Array = comp.grab(query)
	return SN_Array

