import build_spectral_table
import build_event_table
import add_ryan_hst_data
import add_swift_uv_data
import db_maintenance

if __name__ == "__main__":
	# Builds the main spectral table. Homogenizes data by deredshifting,
	# correcting for MW extinction, generating variance spectra, and interpolating.
	# Spectrum specific metadata are also stored in this table.
	build_spectral_table.main()

	# adds hst spectral data separately
	add_ryan_hst_data.main()

	# adds swift spectral data separately
	add_swift_uv_data.main()

	# Builds the event table with available metadata. Light curves
	# from the OSC are also included
	build_event_table.main()

	# Since the scripts above were written by many different contributors, they are
	# difficult to debug. Addtionally, rebuilding the spectral table can take several hours.
	# Instead this script fixes various problems that they introduce. Most importantly, 
	# phases are adjusts to a common t_max for each event. 
	db_maintenance.main()