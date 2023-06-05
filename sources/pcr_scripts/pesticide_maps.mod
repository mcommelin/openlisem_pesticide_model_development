#! --clone mask.map --lddin --matrixtable
#------------------------------
# Produce dataset of fixed pesticide variables for a lisem event

binding

#INPUT
	landuse = lu.map;
	pesttbl = pest.tbl;
	
#OUTPUT
	# maps for pesticide data
	pestwat = pcmixwat.map;
	pestsed = pcmixsoil.map;
	mixdep = pestmixdep.map;
	pestsoildep = pestsoildep1.map;
	pcsoil1 = pcsoil1.map;
	
initial

	# pesticide maps
report pestwat = lookupscalar(pesttbl, 2, landuse);
report pestsed = lookupscalar(pesttbl, 1, landuse);
report mixdep = lookupscalar(pesttbl, 3, landuse);
report pestsoildep = lookupscalar(pesttbl, 4, landuse);
report pcsoil1 = lookupscalar(pesttbl, 5, landuse);