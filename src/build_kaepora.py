import build_spectral_table
import build_event_table
import add_ryan_hst_data
import add_swift_uv_data
import db_maintenance

if __name__ == "__main__":
	build_spectral_table.main()
	add_ryan_hst_data.main()
	add_swift_uv_data.main()
	build_event_table.main()
	db_maintenance.main()