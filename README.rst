Peak calling by CWT and subsequent filtering
========

An experimental method of CLIP-seq peak calling.
Peak-shaped bits of signal are highlighted by scipy's CWT.
Various control methods are then applied to create different peak lists.

Uses one gtf library file, bed files of reads for negative controls (a minimum of a negative IP and an RNA-seq dataset), CLIP bed files and CLIP bedgraph files.
Paths to files are defined in a config file.

testbed/, testwigs/, lib/mtdna.gtf and configfile are example input files.

Example: ::

	$ python main.py -c configfile -x

The -x ensures that new peaks are called instead of loading existing data.

Outputs tables of calls by different methods to data/experiment_name/cwt_calls/.

