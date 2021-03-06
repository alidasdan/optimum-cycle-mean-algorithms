optimum-cycle-mean-algorithms
==============================

This package contains the software that implements many optimum
cycle mean algorithms as detailed in [Da04]. I will refer to the
software as CYCLE_MEAN for ease of reference.

## OPTIMUM CYCLE MEAN PROBLEM

Consider a cyclic graph where every edge has a single number
associated with it, called its weight. The weight of a cycle is equal
to the total weight of the edges along the cycle. The length of a
cycle is equal to the number of edges along the cycle. The mean of a
cycle is the total weight divided by the total length. The mean is
like the average weight of the cycle; in other words, it is like
saying what would be the edge weight of this cycle if each of its
edges had the same weight.

Refer to my optimum cycle ratio algorithms package to learn about a
generalization of cycle means to cycle ratios.

Finding the shortest (in weight) cycle is tractable but finding the
longest (in weight) cycle is NP-hard. It is interesting that finding
the cycle whose mean is optimum, i.e., maximum or minimum, is
tractable. The algorithms in this package solve this problem.

This problem is fundamental to analyzing the performance of discrete
event systems. This is another way of saying if you need to find the
optimum performance of a system, say, the optimum speed an electronic
circuit can run at or the optimum capacity a railway network can
carry, you will need the algorithms implemented in this package.

## MORE INTRODUCTION

I originally developed CYCLE_MEAN during my PhD study (which was
PhinisheD in 1999) but re-implemented a couple of times to get the
current efficient versions. The principles guiding my development
effort were simplicity and efficiency. I originally wrote CYCLE_MEAN
in C++ on Solaris operating system. I later ported it to Linux. I
expect that it can also run on other operating systems with minor
modifications, if any.

CYCLE_MEAN is available on an "as is" basis. I do not say or imply
that it will be useful for whatever you want to do with it. It may
also contain bugs, and I assume no responsibility for any potential
problems associated with its use. You can use CYCLE_MEAN free of
charge in academic research and teaching. For any commercial use,
contact Ali Dasdan at ali_dasdan@yahoo.com. See the COPYRIGHT section
below.

## HOW TO BUILD

Under the 'src' directory, type 'make' (or 'gmake') to build all
executables. An executable named ALGO can be generated for any file
named 'ad_alg_ALGO.cc'. The corresponding make target is ALGO. The
executables are all have .x extension. You can also build the same
target by using any prefix of the name. For example, you can build the
executable 'yto.x' (of the Young-Tarjan-Orlin's algorithm) by typing
'make yto' or 'make yt' or 'make y'.

With no targets following the make command, the following executables
will be generated:
- 'burns.x'  (Excluded due to its slowness)
- 'howard.x'
- 'ko.x'
- 'lawler.x'
- 'szymanski.x'
- 'tarjan.x'
- 'valiter.x'
- 'yto.x'

## HOW TO RUN

Under the 'src' directory, type the name of one of the executables in
your command line to get the usage information. For example, typing
'yto.x' prints out the following information:

```
> yto.x
Usage: yto.x
   [input_file]     file to read input graph (MUST BE 1ST ARG)
   [-m/ode 0/1/2]   read or generate -- see ad_util.cc for details
   [-v/ersion 0/1]  min or max version
   [-n nruns]       number of runs to perform
   [-o offset]      subtract offset from every edge weight
   [-p/aram n m]    num nodes and edges for graph generation
   [-d/ist 0/1/2]   distribution to use -- see ad_util.c for details
   [-w/eight w1 w2] min and max weight bounds
   [-s seed]        random number generator seed
   [-f dump_file]   file to dump output
Below are what is known at this point.
	mode= 0
	input file= 
	version= min
	offset= 0
	num runs= 1
	n= 0
	m= 0
	dist= uniform
	[w1:w2]= [ 1 : 300 ]
	seed= -1
	dump file= 
```

The simplest non-trivial usage is the executable name followed by the
input file name that contains the input graph. For example, 'yto.x
sample.d' prints out the following output:

```
> yto.x sample.d
time to read input graph=       0.00
time to find components=       0.00
run_no= 0
final min_lambda=       40.0 time=       0.00
```

This output shows that the minimum cycle mean of the graph described
in 'sample.d' is 40.0. To get the maximum cycle mean of the graph,
run the same command followed by '-v 0', which should produce
50.0. Note that the minimum version is the default. Also note that the
output also shows how many seconds each main step of the program took.

These flags should be self explanatory but as the usage information
shows, you can do a couple of powerful manipulations with these
flags. For example, you can regenerate arc weights using a number of
supported distributions or you can do multiple runs for runtime
measurement purposes. 

For more information on the input flags, see the code and Makefile.

## HOW TO TEST

Under the 'src' directory, type 'make test' to test each executable on
the sample.d file. The result will be a 'pass' or a 'fail'.

To see the results of all executables on all .d files, see the files
'all-min-runs.txt' and 'all-max-runs.txt' under
'github/alidasdan/graph-benchmarks'.

## HOW TO CLEAN

Under the 'src' directory, type 'make clean'.

## INPUT FILE FORMAT

The input file format is the DIMACS format. See the file sample.d for
a sample input file, which is also explained below.

```
> cat sample.d
p sample-253926760 4 7
a 1 2 40 9
a 2 1 60 17
a 2 3 50 8
a 3 1 30 24
a 4 3 60 22
a 2 4 70 14
a 4 1 30 20
```

Here the 'p' line (the 'problem' line) states that our graph, named
'sample', has 4 nodes or vertices and 7 arcs or directed edges. The arc
weights are generated from a (usually uniform) random number generator
initialized by a seed of '253926760'.

The 'a' lines (the 'arc' lines) following the 'p' line list each arc
as from a source node to a target node with two integer weights: an
arc weight followed by a transit time. For example, the first 'a' line
indicates an arc from node '1' to node '2' with a weight '40' and a
transit time '9'. Note that the node ids start from 1 instead of 0.

For cycle mean algorithms, the transit time is ignored or can be
thought of as equal to '1'. Any lines marked with a 'c' is a comment.

The file sample.pdf (generated from sample.dot using the dot tool in
the graphviz package) shows a picture of the graph in sample.d.

You can use the graphs under 'github/alidasdan/graph-benchmarks' as
input graphs to CYCLE_MEAN algorithms. You can also use the scripts
under 'github/alidasdan/graph-benchmarks/scripts' to change the arc
weights, e.g., regenerate using different random number generator
seeds or using different distributions. 

## REFERENCE

Please cite this reference if you use my programs in your research
work. A preprint is accessible from the 'doc' directory.

```
@article{Da04,
 author = {Ali Dasdan},
 title = {Experimental analysis of the fastest optimum cycle ratio and mean algorithms},
 journal = {ACM Transactions on Design Automation of Electronic Systems},
 volume = {9},
 number = {4},
 year = {2004},
 issn = {1084-4309},
 pages = {385--418},
 doi = {http://doi.acm.org/10.1145/1027084.1027085},
 publisher = {ACM Press},
 address = {New York, NY, USA},
 }
```

## COPYRIGHT

COPYRIGHT C 1999 - Ali Dasdan

Permission to use for non-commercial purposes is granted provided that
proper acknowledgments are given. For a commercial licence, contact
Ali Dasdan at ali_dasdan@yahoo.com.

This software is provided on an "as is" basis, without warranties or
conditions of any kind, either express or implied including, without
limitation, any warranties or conditions of title, non-infringement,
merchantability or fitness for a particular purpose.

## END OF FILE
