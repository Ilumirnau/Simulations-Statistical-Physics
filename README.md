[Hard-sphere-Gas](/Hard-sphere-Gas/) is a folder that contains the simulation of a hard-sphere gas model trapped in a cubic container of hard walls.
<ul>
<li> Dynamics of the system. </li>
<li>  Change of the initial conditions for specific temperatures of the gas. </li>
<li>  Thermal expansion coefficient.</li>
<li>  Pressure-volume relation in free expansion (mass-less walls).</li>

</ul>

[Ising](/Ising/) includes the file ising.py and it is a Montecarlo simulation implementing the Metropolis algorithm of a 20x20 matrix of up&down interacting spins. The result will be:
<ul>
<li>Plots of the system as the temperature increases to visualize its evolution.</li>

<li>A gif of all the plots together and a plot of the magnetization evolution.</li>

<li>A text file with energy and magnetization data will also be created independently so that data for further calculations of parameters like critical exponents or fluctuations.</li>
</ul>

With the parameters in the code, there should be a phase transition when the factor (temperature parameter) is close to 2.3, some images of the system at different temperatures around this value are included as an example.

[RW-1D](/RW-1D/) includes a code that will allow you to illustrate a one-dimensional random walk (RW) algorithm. The final results are two pdf documents:
<ul>
<li>One for the plot of all the positions along the M steps of the N trajectories together with the probability distribution of the final positions for these trajectories (the file name will be ('..._total_fig.pdf').</li>

<li>The other one will contain the data of the Gaussian distribution of the final positions of the simulated walks and the expected theoretical values (the file name will be ('..._distr_params.pdf').</li>

<li>This algorithm allows the study of 1-dimensional diffusion and includes the effect of drift in the RW, as a result of a biased step probability.</li>
</ul>
