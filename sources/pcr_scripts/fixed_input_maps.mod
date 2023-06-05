#! --clone mask.map --lddin --matrixtable
#------------------------------
# Produce dataset of fixed variables for a lisem event

# 1 make homogeneous maps: thetas1, thetas2, psi1, psi2, ksat2, thetai2
# 2 make landuse specific maps: lai, smax, ch, cover, rr
# 3 create buffer volume in thalweg by adjusting rr
# 4 adjust variables for compacted wheeltracks
# 5 smooth channel slope

binding

#INPUT
	lu_in = landuse.map;
	catch = catchment.map;
	croptbl = fixed.tbl;
	homtbl = homogeneous.tbl;
	catchldd = ldd_catch.map;
	thalweg = thalweg.map;
	roughness = roughnes.map;
	compacted = compacted.map;	
	
	channum = channum.map;
	slope = slope.map;
	outpoint = outlet.map;
	tilldir = tilldir.map;
	
#OUTPUT
	# helper maps
	one = one.map;			# map with catchment = 1
	zero = zero.map; 		# map with catchment = 0
	rr_temp = rr_tmp.map;

	# homogeneous maps
	thetas1 = thetas1.map;
	thetas2 = thetas2.map;
	psi2 = psi2.map;
	ksat2 = ksat2.map;
	thetai2 = thetai2.map;
	soildep1 = soildep1.map;
	soildep2 = soildep2.map;
	
	#landuse specific maps
	landuse = lu.map;
	lai = lai.map;
	ch = ch.map;
	cover = per.map;
	smax = smax.map;
	rr = rr.map;	
	
	slope_chan = slopechan.map;
	stream = stream.map;
	
	ksat = ksat1.map;
	thetai = thetai1.map;
	n = n.map;
	ksatcomp = ksatcomp.map;
	thetascomp = porecomp.map;
	
	# maps for dynamic or sediment transport
	aggr = aggrstab.map;
	
initial

	# helper maps
report zero = scalar(lu_in) * 0;
report one = zero + 1;
	area = nominal(one);

	# adjust channel lu class
	landuse = if(lu_in == 99, 3, lu_in);
report	stream = streamorder(catchldd);
	landuse = if(channum == 1 and stream > 4, 99, landuse);
report landuse = if(lu_in == 8, 8, landuse);

	# homogeneous maps
	thetas1 = lookupscalar(homtbl, 1, area);
report thetas2 = lookupscalar(homtbl, 2, area);
report ksat2 = lookupscalar(homtbl, 3, area);
report thetai2 = lookupscalar(homtbl, 4, area);
report soildep1 = lookupscalar(homtbl, 5, area);
report soildep2 = lookupscalar(homtbl, 6, area);
report psi2 = lookupscalar(homtbl, 8, area);

	#landuse specific maps
report	lai = lookupscalar(croptbl, 1, landuse);
report	smax = lookupscalar(croptbl, 2, landuse);
report	ch = lookupscalar(croptbl, 3, landuse);	
	cover = lookupscalar(croptbl, 4, landuse);
report	rr_temp = lookupscalar(croptbl, 5, landuse);
	aggr = lookupscalar(croptbl, 6, landuse);
	# no splash erosion on asphalt roads
report aggr = if(landuse == 8, 500, aggr);

	# adjust variables for compacted wheeltracks
	# rr, thetas, cover

	comp_frac = if(compacted > 1, 1, compacted);
    rr = (comp_frac * 0.5 * rr_temp) + ((1 - comp_frac) * rr_temp);
report thetas1 = (comp_frac * 0.9 * thetas1) + ((1-comp_frac) * thetas1);
	cover = (1 - comp_frac) * cover;
	
	# give roads a full cover to prevent splash erosion
report cover = if(landuse == 8, 1, cover);

# buffer volume in thalweg (adjust rr)
# select cells which have the following 3 chcracteristics:
# 	1: selected by tcrp as thalweg cell
#	2: have an upstream area of more than 1000 m2
#	3: are not within 3 cell of a tillage direction border 
# the last criteria is because tcrp also select a lot of cells which are not
# a thalweg in reality, these are all within 3 cells distance to a border between
# tillage directions.

	dx = celllength();
	accu = accuflux(catchldd, dx * dx);
	upsign = if(accu > 1000, 1, 0);
	thal = cover(thalweg, 0);
	thalaccu = if(thal and upsign, 1, 0);
	thalaccu = cover(thalaccu, 0);

	bordercell = windowdiversity(ordinal(tilldir), 10);
	thalselect = if(bordercell == 1 and thalaccu == 1, 1, 0);

# also add buffer volume for the upstream cells (more realistic representation of
# total buffer volume by potato ridges)
	thalupstr = downstream(catchldd, thalselect);
	rradj = if(thalupstr, roughness *0.4, 0);
	rradj = if(thalselect, roughness * 0.8, rradj);
 	thalwegrr = if(catch, rradj);
report	rr = rr + thalwegrr;

	# smooth channel slope (based on height difference from top to flume)
	# elevation flume = 79.5, elevation start ditch = 86.1
	# length ditch = 300 m --> mean slope = 0.022
	channel_slope = lookupscalar(homtbl, 7, area);
report slope_chan = if(stream > 4 and channum == 1, channel_slope, slope);
