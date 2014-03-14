  GARD (graded autocatalysis replication domain) model for the Origin of Life
-------------------------------------------------------------------------------
Cite this as: [Markovitch O. and Lancet D. (2012) Artificial Life 18:3; http://www.mitpressjournals.org/doi/abs/10.1162/artl_a_00064].
Contact: Omer Markovitch, omermar@gmail.com, http://sites.google.com/site/omermar .
Additional sites: www.weizmann.ac.il/molgen/members/lancet.html, http://en.wikipedia.org/wiki/Gard_model .

GARD is a general kinetic model for homeostatic-growth and fission of compositional-assemblies, with specific application towards lipids.
Additional information can be found in the aforementioned publication and its supporting information, and by contacting me.

This short text will give a quick guide on how to run GARD on MATLAB:
A1) Unzip the files into a single directory.
A2) Add this path into MATLAB, either by File-->Set Path, or by "addpath(path,'<put the full path here>');"

B1) Create a parameter structure with "p=tgs_parameters_v10();", the default values are identical to those in the publication (except for "p.GEN=5000"). Note that the beta matrix is not generated with this command.
B2) Set a random seed, that will also determine the beta, with: "p.seed=[1 1 1]". The first number is relevant to beta.
B3) Generate a beta network: "p.Beta=tgs_newbeta_v10(p);".
B4) Now you can run a GARD simulation, including the clustering for composomes/compotypes: "o=tgs_agard_v10(p, 1);".
B5) Plotting a 'similarity carpet': "c=tgs_carpet_v10(o.trace); title('GARD'); xlabel('Generation'); ylabel('Generation');".
B6) Some description of the output data structure:
	B6i) "o.trace" is a NG by GEN matrix holding the composition of the assembly at each generation. Each composition is a NG long vector of size Nmax (Nmax=p.splitsize*p.NG).
	B6ii) 'o.tags' tells for each generation to which compotype it belongs to, or zero for drift. Try: "hist(o.tags, [0:1:size(o.comps,2)]); xlabel('Compotype index'); ylabel('Frequency'); title('GARD');".
	B6iii) 'o.comps' is a matrix holding all the compotypes compositions, un-normalized. For example, if only a single compotype found then then this is a NG by 1 matrix.
	B6iv) 'c' is a NG by NG matrix holding all the H values.

C1) Perform a 'biased GARD' simulation, to study the selection response of the most frequent GARD compotype: "ob=biased_gard_v10(p, 1.1, 1);".
C2) Some of description of the output data structure:
	C2i) For "ob.regular" see (B6).
	C2ii) "ob.h" is the H between the target in the regular and in the biased GARD.
	C2iii) "ob.target" is the composition vector of the target from the regular simulation.
	C2iv) "ob.targetfreq" and "ob.targetfreqbiased" are how many generations the target occured in the regular and the biased simulations.
	C2v) "ob.it" and "ob.itbiased" are the indices of the target compotype in the regular and the biased simulations.

D1) Perform a 'population GARD' simulation, for a population of size 'p.gen' and for a duration of 2000 splits: "opop=population_gard_nmin_v10(p, 2000);".
D2) Some of description of the data structure:
	D2i) "opop.trace" is a NG by p.gen matrix holding the population after the 2000 splits.
D3) Perform a 'population selection GARD': "opopb=population_gard_nmin_v10(p, 2000, trgt);". Where 'trgt' is a NG by 1 vector holding the composition of the selection-target.

E) "pmc=pmc_beta_v10(o.Beta);" is the mutual catalysis power of this beta network.

F1) "[timearray, ctnorm]=correlate_carpet_v10(o.trace, [], 1);" calculates and plot the similarity autocorrelation for this GARD simulation.
F2) "exp1=omer_curvefit_v10(timearray, ctnorm, display);" fits the autocorrelation function with a first order exponent. Note that in the paper first the tail is smoothed before fitting.

That is it to begin with. Feel free to read the paper and contact me for further information and/or help. Omer.

Acknowledgement: The Lipid World and GARD related research was invented and is headed by Prof. Doron Lancet from the Department of Molecular Genetics at Weizmann Institute of Science. Throughout the years the following contributed and continued the work: Daniel Segre, Barak Shenhav, Aia Oz, Hamutal Arbel, Aron Inger and Omer Markovitch.


