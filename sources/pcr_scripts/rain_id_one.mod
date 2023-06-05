#! --degrees --clone mask.map --lddin --matrixtable
#-------------------------------------------
# prepare id.map based on rainfall source

binding
	#INPUT MAPS
	catch = catchment.map;

	#OUTPUT MAPS
	id_out = id.map;
	
initial

report id_out = nominal(if(catch, 1));
