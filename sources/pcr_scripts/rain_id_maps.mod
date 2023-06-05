#! --degrees --clone mask.map --lddin --matrixtable
#-------------------------------------------
# prepare id.map based on rainfall source

binding
	#INPUT MAPS
	id_in = id_raw.map;
	id_tbl = rain_id.tbl;

	#OUTPUT MAPS
	id_out = id.map;
	
initial

report id_out = lookupscalar(id_tbl, 1, id_in);
