#! --degrees --clone mask.map --lddin --matrixtable
#-------------------------------------------
# prepare input for tcrp model by Takken 2001
# from the chan and dem a channel ldd must be made
# the channel map must be made boolean with 0 outside
# create roughness maps with lookup table

binding
	#INPUT MAPS
	dem = dem.map;
	chan = channels.map;
	fields = fields.map;
	ro = roughness.tbl;

	#OUTPUT MAPS
	chanldd = chanldd.map;
	channum = channum.map;
	chanout = chan.map;
	roughnes = roughnes.map;

initial
report  channum = cover(chan, 0);	
report	chanout = boolean(if(channum > 0, 1, 0));
		chandem = if(chanout, dem);
		chandem = lddcreatedem(chandem, 1e31, 1e31, 1e31, 1e31);
report	chanldd = lddcreate(chandem, 1e31, 1e31, 1e31, 1e31);

report	roughnes = lookupscalar(ro, 1, fields);
