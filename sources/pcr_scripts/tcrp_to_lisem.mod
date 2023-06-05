#! --clone mask.map --lddin --matrixtable
#------------------------------
# derive catchment boundaries from tcrp result
# produce lisem map database based on that.

binding

	#INPUT
	lddres = lddres.map;
	outpoint = outpoint.map;
	sloperes = sloperes.map;
	fields = fields.map;
	dem = dem.map;
	thalweg = thalweg.map;
	roughness = roughnes.map;
	roughtbl = roughness2020.tbl;
	
	
	#OUTPUT
	catch = catchment.map;
	outp = outlet.map;
	demout = dem_catch.map;
	fieldsout = landuse.map;
	slopeout = slope.map;
	lddout = ldd_catch.map;
	lddtop = lddtopo_catch.map;

initial

		outpointB = boolean(if(cover(outpoint, 2) == 0, 1, 0));
report 	outp = outpointB;
report	catch = catchment(lddres, outp);
report	demout = if(catch, dem);
report	fieldsout = if(catch, fields);
report  slopeout = if(catch, sloperes);
	    lddtmp = if(catch, lddres);
report	lddout = if(outpointB, 5, lddtmp);
        lddtopo = lddcreate(dem, 1e31, 1e31, 1e31, 1e31);
		catchtopo = catchment(lddtopo, outpointB);
		lddtmp = if(catchtopo, lddtopo);
report  lddtop = if(outpointB, 5, lddtmp);

