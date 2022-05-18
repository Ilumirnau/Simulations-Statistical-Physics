# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 23:28:42 2021

@author: arnau
"""

import numpy as np
import random
import matplotlib.pyplot as plt

#per tractar dades amb xifres significatives
def xs(valor, num_xs):
    num = valor*10**(num_xs)
    return int(num)/10**(num_xs) 

def magnetitzacio(M, system):
    m = np.sum(system)
    return m/M**2/spin

def energia(M, spin, system):
    spin = abs(spin)
    e = 0
    for c in range(M):
        for d in range(M):
            left=c-1 #c-1
            right=c+1 #c+1
            top=d-1 #d-1
            bottom=d+1 #d+1
            e += -J/spin**2 * (system[d, left%M]+system[d, right%M]+system[top%M, c]+system[bottom%M, c])*system[d,c]
    return e/M**2

#guardar dades que costen de computar en fitxer txt per posterior tractament
def guardar_dades(nom_fitxer, x, y):
  import numpy as np
  a_file = open(str(nom_fitxer)+".txt", "w")
  dades = [[x], [y]]
  for row in dades:
    np.savetxt(a_file, row)
  a_file.close()
  
#calcular la mitjana dels x=items últims valors d'una llista
def last_mean(vector, items):
    l = len(vector)
    v = vector[int(l - items):l]
    return v

#contar caselles per energies
def mb_check(apunts, system, row, col):
    left = col-1
    right = col+1
    top = row-1
    bottom = row+1
    c = col
    d = row
    ee = -J/spin**2 * (system[d, left%M]+system[d, right%M]+system[top%M, c]+system[bottom%M, c])*system[d,c]
    try: apunts[ee] +=1
    except : apunts[ee] = 1
    return apunts

def metropolis(M, N, factor, spin, system):
    spin = abs(spin)
    k = J
    #epsilon= kt * factor
    #initial system definition
    if system.all() == np.zeros((M,M)).all():
        for i in range(M):
            for j in range(M):
                random_num = random.randint(0, 2*spin) #index llista
                spins = np.linspace(-spin, spin, int(2*spin + 1)) #llista de possibles valors spin
                spin_i = spins[random_num] #agafar index random de la llista de spins
                system[i,j] = spin_i #establir spin a la posició ij
      
    #magnetització inicial
    magnetitzacio1=magnetitzacio(M, system)
    #energia inicial
    energia1 = energia(M, spin, system)
    #print(b)
    p=0 #probabilitat
    left=0 #c-1
    right=0 #c+1
    top=0 #d-1
    bottom=0 #d+1
  
    #evolucio magnetitzacio i energia
    magnetitzacio_t = [magnetitzacio1]
    energia_t = [energia1]
    for i in range (N-1):
        e_f=0 #energia final
        e_i=0 #energia inicial
        c=random.randint(0,M-1)
        d=random.randint(0,M-1)
        left=c-1 #c-1
        right=c+1 #c+1
        top=d-1 #d-1
        bottom=d+1 #d+1
        #initial energy
        #used when considering deltaE the energy difference of the systems
        """micro_sys_i = np.empty((5, 5))
        for j in range(top-1, bottom+2):
            for k in range(left-1, right+2):
                micro_sys_i[j-(top-1)][k-(left-1)] = system[j%M][k%M]
        e_i = energia(5, spin, micro_sys_i)"""
        #deltaE the energy difference of the specific position
        deltaE = 2 * J/spin**2 * system[d, c] * (system[d, left%M]+system[d, right%M]+system[top%M, c]+system[bottom%M, c])
        #spin change
        system[d,c] *= -1
        #final energy
        """micro_sys_f = np.empty((5, 5))
        for j in range(top-1, bottom+2):
            for k in range(left-1, right+2):
                micro_sys_f[j-(top-1)][k-(left-1)] = system[j%M][k%M]
        e_f = energia(5, spin, micro_sys_f)"""
        if deltaE>0 and random.uniform(0,1)>=np.exp(-deltaE/(k*factor)):
            system[d,c] *= -1
                
        magnetitzacio_i = magnetitzacio(M, system)
        magnetitzacio_t.append(magnetitzacio_i)
        
        energia_i = energia(M, spin, system)
        energia_t.append(energia_i)
    m_final = np.mean(magnetitzacio_t[int(0.6*N):N])
    u_final = np.mean(energia_t[int(0.6*N):N])
    return magnetitzacio_t, energia_t, system, m_final, u_final

def m_plot(x, y, factor, N):
    import matplotlib.pyplot as plt
    plt.style.use('dark_background')
    plt.plot(x, y, linewidth=0.3)
    plt.ylabel(r'M/$s_i$')
    plt.xlabel('Steps (N)')
    plt.title(str(factor))
    plt.grid(alpha=0.5)
    plt.xlim(0, N)
    plt.savefig('magnetitzacio_'+str(xs(factor, 4))+'.jpg', dpi=1000)
    plt.show()

def u_plot(x, y, factor, N):
    import matplotlib.pyplot as plt
    plt.style.use('dark_background')
    plt.plot(x, y, linewidth=0.3)
    plt.ylabel(r'U/$\epsilon$')
    plt.xlabel('Steps (N)')
    plt.title(str(factor))
    plt.grid(alpha=0.5)
    plt.xlim(0, N+1)
    plt.savefig('energia_'+str(xs(factor, 4))+'.jpg', dpi=1000)
    plt.show()
  
def system_plot(system, spin, factor, N, M):
    import matplotlib.pyplot as plt
    xx = np.linspace(1/2, M+1/2, M+1)
    yy = np.linspace(1/2, M+1/2, M+1)
    X, Y = np.meshgrid(xx, yy)
    plt.pcolormesh(X, Y, system, cmap = 'Spectral', vmin = -spin, vmax = spin, edgecolors = 'white', linewidth = 0.3)
    plt.title('spin='+str(spin)+', factor='+str(xs(factor, 4))+', steps='+str(N))
    xlabels = []
    for i in range(1, M+1):
        if i%2!=0:
            xlabels.append(str( ))
        else:
            xlabels.append(str(int(i)))
    plt.xticks([i for i in range(1, M+1)], xlabels)
    ylabels = []
    for i in range(1, M+1):
        if i%2!=0:
            ylabels.append(str( ))
        else:
            ylabels.append(str(int(i)))
    plt.yticks([i for i in range(1, M+1)], ylabels)
    #plt.grid(alpha=0.5, markevery = 5)
    cbar = plt.colorbar()
    spins = np.linspace(-spin, spin, int(2*spin + 1))
    cbar.set_ticks(spins)
    cbar.set_ticklabels([str(i) for i in spins])
    cbar.set_label('Spin')
    plt.savefig('system_'+str(xs(factor, 4))+'.jpg', dpi=1000)
    plt.show()
    
import numpy as np
import matplotlib.pyplot as plt
M=20 #mida matriu
N=10 * M**3
spin = 1/2
J = 1

def evo_m_and_e():
    system = spin * np.zeros((M,M))
    factor = 20
    y_m, y_e, system, m, u = metropolis(M, N, factor, spin, system)
    plt.style.use('bmh')
    plt.plot(range(N), y_m, linewidth = 1, color= 'teal')
    plt.xlabel('N intents')
    plt.ylabel(r'$m (T=20)$')
    plt.ylim((-1, 1))
    plt.savefig('mtempst=20.png', dpi = 700)
    plt.show()
    plt.plot(range(N), y_e, linewidth = 1, color= 'teal')
    plt.ylim((-1, 1))
    plt.xlabel('N intents')
    plt.ylabel(r'$u (T = 20)$')
    plt.savefig('etempst=20.png', dpi = 700)
    plt.show()
evo_m_and_e()

def figures_and_video_complete():
    import numpy as np
    import matplotlib.pyplot as plt
    M=20 #mida matriu
    N=[int(5*M**3), int(M**3)] #numero d'intents
    factor = 1E08 #ratio E/kT  #T<<1 => factor>>1  # factor \propto beta
    factors1 = np.linspace(0.25, 1.39, int(30))
    factors2 = np.linspace(1.4, 3, int(200))
    factors3 = np.linspace(3.01, 3.2, int(20))
    factors = [i for i in factors1] + [i for i in factors2] + [i for i in factors3]
    spin = 1/2
    J = 1
    
    system = spin * np.ones((M,M))
    
    magnetizations = []
    energies = []
    for i in factors:
        if i == factors[0]:
            n = N[0]
        else:
            n = N[1]
        x = range(n)
        y_m, y_e, system, m, u = metropolis(M, n, i, spin, system)
        plt.style.use('bmh')
        plt.plot(range(n), y_m)
        plt.savefig('mtemps.jpg')
        #y_m=[i/spin for i in y_m]
        magnetizations.append(m)
        energies.append(u)
        #m_plot(x, y_m, i, N)
        #u_plot(x, y_e, i, N)
        system_plot(system, spin, i, n, M)
    
    guardar_dades('magnetitzacio_vs_T', factors, magnetizations)
    plt.plot(factors, magnetizations)
    plt.ylabel('m')
    #plt.gca().invert_xaxis()
    plt.xlabel(r'$T$')
    plt.savefig('m_evolution.jpg', dpi=400)
    plt.show()
    
    guardar_dades('energies_vs_T', factors, energies)
    plt.plot(factors, energies)
    plt.ylabel('u')
    #plt.gca().invert_xaxis()
    plt.xlabel(r'$T$')
    plt.savefig('u_evolution.jpg', dpi=400)
    plt.show()
    
    
    import imageio
    import os
    images=[]
    #print(os.listdir())
    for i in factors:
        for filename in os.listdir():
            if filename.startswith("system_"+str(xs(i, 4))+".jpg"):
                images.append(imageio.imread(filename))
            else: continue
    imageio.mimsave('system evolution video.gif', images, fps=5)
    
figures_and_video_complete()

def histogram_4():
    import numpy as np
    M=20 #mida matriu
    N=M**3 #numero d'intents
    factor = 3
    spin = 1/2
    J = 1
    
    system = np.zeros((M,M))
    
    magnetizations = []
    energies = []
    y_m, y_e, system, m, u = metropolis(M, N, factor, spin, system)
    system_plot(system, spin, factor, N, M)
    apunts = {}
    for i in range(M):
        for j in range(M):
            apunts = mb_check(apunts, system, i, j)
    bins = [i for i in list(apunts.keys())]
    counts = [i/M**2 for i in list(apunts.values())]
    
    # Function to calculate the Gaussian with constants a, b, and c
    def gaussian(E, T):
        k=J
        return np.exp(-E**2/(k*T))
    E = np.linspace(min(bins)*1.5, max(bins)*1.5, 100)
    #plt.plot(E, [gaussian(i, factor) for i in E], 
                    #color = 'green', linewidth = .7, label = 'Theroetical Distribution')
    plt.style.use('bmh')
    plt.bar(bins, counts, zorder=2, color= 'teal')
    plt.grid(alpha=0.7, zorder=1)
    plt.xlabel('Energies')
    plt.ylabel('Freqüència')
    plt.savefig('histo3.jpg', dpi=900)
    plt.show()

histogram_4()