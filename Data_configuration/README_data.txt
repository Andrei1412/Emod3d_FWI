1. Create station lists in cartesian only ("STATION_dh_4km.txt") and in both cartesian, lon/lat ("STATION_utm_dh_4km.txt") tp put in StaInfo:

	run create_broadband_station_list_10BB.py

2. Create an event list:

 	a)run select_events_with_CMT_solution_GeoNet_TEMPLATE.py #to have a list of events with given space/ time and Mw conditions
	b)run constrain_selected_event_geometrically_save_srf_list.py # to constrain the number of events included in the inversion (with most geometrical balance)
	c)run create_utm_source_list_vs_SRF_folder.py # to generate the source list file "SOURCE_SRF.txt" to put in Model/Models and srf-file folder 
	----(SRF_13s in this case, then put this folder in Kernels/iter/iter1 for inversion run in Maui)
	----srf-file found in Andrei_Sources_20191209 database of 422 events OR generated from mahuika!

3. Download broadband data for events listed in "SOURCE_SRF.txt"

	run get_geonet_broadband_data.py # to download data and store in Vel_ob_200s_HH_13s, then put this folder in Kernels/ for inversion run in Maui
