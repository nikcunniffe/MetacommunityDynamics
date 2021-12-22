# MetacommunityDynamics
Supporting code for "Metacommunity dynamics and the detection of species associations in co-occurrence analyses: why patch disturbance matters" by Vincent Calcagno, Nik J. Cunniffe and Frédéric M. Hamelin. 

Aside from this README.md, the repository contains four files
  - generalisedMetapopulations.R: the library to simulate patch occupancy matrices according to our generalised metapopulation model, and to do various calculations;
  - metapopulationDynamics.Rmd: a R Markdown which creates all figures in the paper, to show how the library can be used;
  - metapopulationDynamics.html: the knitted output of the R Markdown;
  - refs.bib: the bibtex file referred to in the R Markdown file and required to make it knit.
  
Note that on L38 of the .Rmd file there is a boolean switch that can be used to toggle between a short testing run and a full length run (the latter takes ~20 hours or so on a reasonably powerful desktop).
