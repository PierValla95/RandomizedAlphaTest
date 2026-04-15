Each folder allows the user to reproduce one of the tables of the Monte Carlo analysis in the main body of Massacci et al. (2026). 

Within each folder, scripts are named according to the test that they run, using the same acronyms as in the paper. Script "_smallT" are for T<= 500 while results for T=1000 and T=2000 can be replicated with the scripts "_largeT".

Before running the scripts, data must be generated and saved using the scripts "SimulateGaussianData.R", "SimulateTData.R", and "SimulateGARCHData.R". For all these script, setting Type="size" simulates data under the null, while data under the alternative are simulated by setting Type="power". Data for the large sample sizes T=1000 and T=2000 are generated in the analysis scripts so as to avoid saving extremely large files (roughly 12Gb, using .Rdata files). 

Results using the output of the above functions may differ across different R versions. Experiments with the output of the above scripts showed qualitatively similar results. Our data were generated with the December 2024 version of R and are available upon request for perfect replication.


