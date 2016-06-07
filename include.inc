
	integer*4 max_pnts			!max # of E' points.
	parameter  (max_pnts = 150)

	integer*4 max_sets			!max # of data sets.
	parameter (max_sets = 8)

	real*8 param1(40), param2(40), current_slope, current_offset
	real*8 innt, innp, fermip(197), epsn(197), aa, bb, cc, dd
	real*8 param3(40), param4(40), manual_e1(3000), manual_e2(3000)
	real*8 manual_deg(3000), manual_anum(3000)

	common /fits/ param1, param2, current_slope, current_offset,
     &   innp, innt, fermip, epsn, aa, bb, cc, dd, param3, param4,
     &  manual_e1, manual_e2, manual_deg, manual_anum
