pwf - a simple Gauss-Seidal / Newton-Raphon load flow program
========================================

## Usage

pwf [-agrcnlyj] [-lim] [+tol] [-ooutfile] [infile...]

-a	adjusted load flow
-g	G-S method
-r	rectangular N-R method
-c	write to CDF file
-n	print the nodal information
-l	print the line flows
-y	print the YBus
-j	print the Jacobian
-lim	set the iteration limit to lim(INT)
+tol	set the precision to tol(DOUBLE)
-ofilename
	set the output file to filename
infile...
	CDF input files

## Description

The pwf calculates load flows from infiles, and write the results to outfile.
If no options given, pwf reads from stdin and write to stdout, using the 
default settings: limit 150, tolerance .001, and polar form N-R method.
