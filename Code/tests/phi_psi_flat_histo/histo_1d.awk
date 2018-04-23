#!/usr/bin/awk -f

BEGIN {
	#min, max, nbin and col can be passed by the user
	min=100000000
	max=-100000000
	tmin=100000000
	tmax=-100000000
	nbin=100
	col=1
	n=0;
}

$1!~/#/ {
	x[n]=$col
	if(max<min)
	{
		if($col<tmin) tmin=$col
		if($col>tmax) tmax=$col
	}
	n++;
}
END{
	if(max<min){
		max=tmax
		min=tmin
	}
	max=max+0.001
	delta=(max-min)/(nbin)
	area=0
	for(i=0;i<n;i++)
	{
		#values below min will have ix<0
		#values above max will have ix>nbin
		ix=int((x[i]-min)/delta)
		hist[ix]++;
	#	print min,max,delta,x[i],ix,hist[ix]
		area++
	}
	for(ix=0;ix<nbin;ix++)
	{
		printf("%f %f\n",min+delta*(ix),hist[ix]/area)
	}
}

