1D-RW is a code that will allow you to illustrate a one-dimensional random walk algorithm. The final result are two pdf documents:

-One for the plot of all the positions along the M steps of the N trajectories together with the probability distribution of the final positions for these trajectories (the file name will be ('..._total_fig.pdf').

-The other one will contain the data of the Gaussian distribution of the numerical simulation and the expected values from the theoretical point of view (the file name will be ('..._distr_params.pdf').


ising.py is a Montecarlo simulation via the Metropolis algorithm of a 20x20 matrix of up&down spins. The result will be:

-Plots of the system as the temperature increases to visualize its evolution

-A gif of all the plots together and a plot of the magnetization evolution

-A text file with energies and magnetization will also be created independently so that data for further calculation of parameters like critical exponents can be calculated afterwards
