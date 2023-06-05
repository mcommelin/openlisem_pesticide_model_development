#! --clone mask.map --lddin --matrixtable
#------------------------------
# Produce dataset of dynamic variables for a lisem event

# 1 make landuse specifc maps: n, ksat1, thetai1
# 2 adjust variables for compacted wheeltracks

binding

#INPUT
	landuse = landuse.map;
	catch = catchment.map;
	inftbl = dynamic.tbl;		# dynamic table not used yet
	wheeltbl = wheels.tbl;		# ksat for wheeltracks
	catchldd = ldd_catch.map;
	compacted = compacted.map;	
	rr = rr_tmp.map;
	
	# helper maps
	one = one.map;			# map with catchment = 1
	zero = zero.map; 		# map with catchment = 0
	
#OUTPUT
	
	#landuse specific maps
	n = n.map;
	ksat1 = ksat1.map;
	psi1 = psi1.map;
	thetai1 = thetai1.map;
	coh = coh.map;
	d50 = d50.map;
	comp_coh = cohcomp.map;
	
initial

	#landuse specific maps
	ksat1 = lookupscalar(inftbl, 1, landuse);
	n = lookupscalar(inftbl, 2, landuse);
report	thetai1 = lookupscalar(inftbl, 3, landuse);	
	psi1 = lookupscalar(inftbl, 4, landuse);


	# adjust variables for compacted wheeltracks

	comp_frac = if(compacted > 1, 1, compacted);
	ksat_wheel = lookupscalar(wheeltbl, 1, one);
	n_wheel = rr/100 * lookupscalar(wheeltbl, 4, one);
	psi_wheel = lookupscalar(wheeltbl, 2, one);
report	ksat1 = comp_frac * ksat_wheel + (1-comp_frac) * ksat1;
report  n = comp_frac * n_wheel + (1-comp_frac)*n;
report	psi1 = comp_frac * psi_wheel + (1-comp_frac)*psi1;

	# erosion related variables
	coh = lookupscalar(inftbl, 5, landuse);
	coh_wheel = lookupscalar(wheeltbl, 3, one);
	coh_wheel = coh_wheel / (4.7 - 30 * thetai1 + 60 * thetai1**2);
	# increase cohesion for arable land with flow concentration
	accu_coh = if(coh < 0, 100, coh);
	accu = accuflux(catchldd, 1/abs(accu_coh));
	accu = accu/areamaximum(accu, nominal(one));
	accu = if(accu > 0.1, 1, 0);
report	comp_coh = accu;
report	coh = if(accu eq 1 and landuse ne 99 and landuse ne 8, 7.85, comp_frac * coh_wheel + (1-comp_frac) * coh);
report  d50 = lookupscalar(inftbl, 6, landuse);

