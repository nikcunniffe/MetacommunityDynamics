# MetacommunityDynamics
Supporting code for Metacommunity dynamics and the spurious detection of species associations in co-occurrence analyses by Vincent Calcagno, Nik J. Cunniffe and Frédéric M. Hamelin. 

The repository contains five files
  - generalisedMetapopulations.R: the library to simulate patch occupancy matrices according to our generalised metapopulation model, and to do various calculations
  - metapopulationDynamics.Rmd: a R Markdown which creates all figures in the paper, to show how the library can be used
  - metapopulationDynamics.html: the knitted output of the R Markdown file
  - refs.bib: the bibtex file referred to in the R Markdown file and required to make it knit
  - metapopulationDynamicsTest.Rmd: a R Markdown file which is *identical* to the other markdown, except it has testRun <- TRUE on L38. This has settings for the various stochastic algorithms which allow the code to run much faster (but to not obtain statistically valid results)
  - metapopulationDynamicsTest.html: the knitted output of the test markdown file
  
