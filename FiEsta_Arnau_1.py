import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity
from scipy.optimize import curve_fit
import random
import time
tic = time.perf_counter()

# Function to write time in a fancy way
def fancy_t(time):
    if time/3600 > 1:
        hores = int(time/3600)
        minuts = int((time/3600 - hores)*3600/60)
        segons = int((time/60 - minuts)*60*100)/100
        return str(hores)+' hores '+str(minuts)+' min '+str(segons)+' s'
    elif time/60>1:
        minuts = int(time/60)
        segons = int((time/60 - minuts)*60*100)/100
        return str(minuts)+' min '+str(segons)+' s'
    else:
        return str(int(time*100)/100)+' s'

# Function to calculate the Gaussian with constants a, b, and c
def gaussian(x, a, b, c):
    return a*np.exp(-np.power(x - b, 2)/(2*c)) 

# Function to round number "valor" so that it has "num_xs" decimals
def xs(valor, num_xs):
    num = valor*10**(num_xs)
    if num < 0:
        return -int(num)/10**(num_xs)
    else:
        return int(num)/10**(num_xs)
    
# Function to simulate and plot the 1D random walk and the probability distribution of trajectories in the final position
def rw_1d(l, p, N, M): #p (probability of walking +l) N (trajectories) M(steps)
    camins = np.zeros((N,M))
    mitjana = [0]
    for i in range(1,M):
        step_mean = mitjana[i-1]
        for j in range(N):
            r = random.random()
            if r <= p: 
                camins[j][i] = camins[j][i-1] + l
                step_mean += 1/N
            else: 
                camins[j][i] = camins[j][i-1] - l
                step_mean -= 1/N
        mitjana.append(step_mean)
    dada_h = [camins[i][-1] for i in range(N)]
    x = np.linspace(min(dada_h)*1.05,max(dada_h)*1.05, int(M*1.05))[:, np.newaxis]
    x2 = np.linspace(min(dada_h)*1.05,max(dada_h)*1.05, int(M*1.05))
    fig, axs = plt.subplots(1, 2, sharey = True, gridspec_kw={'width_ratios': [3, 1]}, figsize=(10,5))
    fig.subplots_adjust(wspace=.01)
    for j in range(N):
        axs[0].plot(range(M), camins[j][:], color = 'grey', linewidth = .7, alpha = 0.7, drawstyle = 'steps-pre')
    axs[0].plot(range(M), mitjana, color = 'r', linewidth = .8, drawstyle = 'steps-pre', label = 'Mean position per step')
    axs[0].set_ylabel('Position', fontsize=8)
    axs[0].set_xlabel('Steps', fontsize=8)
    axs[0].set_title('Positions of the 1D-RW', fontsize=10)
    for tick in axs[0].xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    hist, bin_edges = np.histogram(dada_h, bins=int((max(dada_h)-min(dada_h))*.5), 
                                   range = (min(dada_h),max(dada_h)))
    axs[1].hist(dada_h, bins=int((max(dada_h)-min(dada_h))*.5), range = (min(dada_h),max(dada_h)), 
                density=True, color = 'grey', orientation = 'horizontal')
    dada_h = np.array(dada_h)
    re_dada = dada_h.reshape(-1, 1)
    kd_dada_h = KernelDensity(kernel = 'gaussian', bandwidth = 5).fit(re_dada)
    kd_vals_dada_h = np.exp(kd_dada_h.score_samples(x))
    axs[1].plot(kd_vals_dada_h, x, color='orange', label= 'Probability Distribution')
    axs[1].legend(loc = 'lower right', prop={'size': 6})
    pars, cov = curve_fit(f=gaussian, xdata=x2, ydata=kd_vals_dada_h, 
                          p0=[1/np.sqrt(2*np.pi * l**2 * (M - (2*p-1)**2)), M*(2*p - 1)*l, l**2 * M - (2*p-1)**2], 
                          bounds=(-np.inf, np.inf))
    stdevs = np.sqrt(np.diag(cov))
    print('val_num')
    for i in range(len(pars)):
        print('n.',i,')',xs(pars[i],5),'+/-',xs(stdevs[i],7))
    print('val_teo \n', 't. 0) ',xs(1/np.sqrt(2*np.pi * l**2 * (M - (2*p-1)**2)),4), 
          '\n t. 1) ', xs(M*(2*p - 1)*l,2), 
          '\n t. 2) ',xs(l**2 * M - (2*p-1)**2, 2))
    axs[1].set_xlabel('Amount of trajectories', fontsize=8)
    axs[1].set_title('Distribution of final positions', fontsize=10)
    for tick in axs[1].xaxis.get_major_ticks():
        tick.label.set_fontsize(8)
    #fig.suptitle('N='+str(N)+'   M='+str(M)+'   p='+str(p))
    plt.savefig('rw_1d_'+str(N)+'_'+str(M)+'_'+str(p)+'_total_fig.pdf')
    plt.figure()
    plt.axis('off')
    plt.text(-.1, .9, r'Distribution: $\rho (x)=a \cdot e^\frac{(x-b)^2}{2c}$')
    plt.text(-.1,.7, 'Values from the numerical simulation: ')
    plt.text(-.1,.6, r'$a=$'+str(xs(pars[0],6))+r'$\pm$'+str(xs(stdevs[0],6))+'       '+r'$b=$'+str(xs(pars[1],3))+r'$\pm$'+str(xs(stdevs[1],3))+'       '+r'$c=$'+str(xs(pars[2],2))+r'$\pm$'+str(xs(stdevs[2],2)))
    plt.text(-.1,.4, 'Values from the theoretical approach: ')
    plt.text(-.1,.3, r'$b=\langle x \rangle = (2p-1)Ml=$'+str(xs(M*(2*p - 1)*l,2))+'       '+r'$c=\sigma^2 = l^2[M-(2p-1)^2]=$'+str(xs(l**2 * M - (2*p-1)**2, 2)))
    plt.text(-.1,.2, r'$a=\frac{1}{\sqrt{2\pi c}}=$'+str(xs(1/np.sqrt(2*np.pi * l**2 * (M - (2*p-1)**2)),3)))
    plt.savefig('rw_1d_'+str(N)+'_'+str(M)+'_'+str(p)+'_distr_params.pdf')
l,p,N,M = 1, 0.6, 3000, 300 #data to plot an example
rw_1d(l,p,N,M)
toc = time.perf_counter()
print('Temps final:', fancy_t(toc-tic))
