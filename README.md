astr596
=======
Spectral template code using an SQL database for Type Ia Supernova observations.

Version specific dependencies:
msgpack-python version 0.4.6
msgpack-numpy version 0.3.5 or 0.3.6

Example Query:
python query_db.py nb "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where phase between 36 and 44 and velocity between -98 and -12" "SELECT * from Supernovae inner join Photometry ON Supernovae.SN = Photometry.SN where phase between 36 and 44 and velocity >= -12"

	First argument "b" additionally estimates errors via bootstrap resampling for template spectra. "nb" generates template spectra without errors (much faster).