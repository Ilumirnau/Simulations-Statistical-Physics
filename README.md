<strong> RW-1D </strong> inlcudes a code that will allow you to illustrate a one-dimensional random walk (RW) algorithm. The final results are two pdf documents:

-One for the plot of all the positions along the M steps of the N trajectories together with the probability distribution of the final positions for these trajectories (the file name will be ('..._total_fig.pdf').

-The other one will contain the data of the Gaussian distribution of the final positions of the simulated walks and the expected theoretical values (the file name will be ('..._distr_params.pdf').

-This algorithm allows the study of 1-dimensional diffusion and includes the effect of drift in the RW, as a result of a biased step probability.


Ising includes the file ising.py and it is a Montecarlo simulation implementing the Metropolis algorithm of a 20x20 matrix of up&down interacting spins. The result will be:

-Plots of the system as the temperature increases to visualize its evolution

-A gif of all the plots together and a plot of the magnetization evolution

-A text file with energies and magnetization will also be created independently so that data for further calculations of parameters like critical exponents.
