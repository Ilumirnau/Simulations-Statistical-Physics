from vpython import *
import numpy as np
import matplotlib.pyplot as plt
#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood


#initial params
Natoms = 400  # change this to have more or fewer atoms

# Typical values
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container
mass = 4E-3/6E23 # helium mass
Ratom = 0.03 # wildly exaggerated size of helium atom, valor inicial 0.03
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature
dt = 1E-5
r = 0.005

#mostra t bonic
def fancy_t(time):
    if time/60>1:
        minuts = int(time/60)
        segons = int((time/60 - minuts)*60*100)/100
        return str(minuts)+' min '+str(segons)+' s'
    else:
        return str(int(time*100)/100)+' s'
        

deltav = 100 # binning for v histogram

def checkCollisions(apos):
    hitlist = []
    r2 = 2*Ratom
    r2 *= r2
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i) :
            aj = apos[j]
            dr = ai - aj
            if mag2(dr) < r2: hitlist.append([i,j])
    return hitlist


###########################################################################################################################
def ViP(mass, k, T, Natoms, Ratom, L, it_final, canvis, T_lim):
    Lx = L
    p = []
    apos = []
    from scipy.stats import maxwell
    p_sp = maxwell.rvs(scale=(k * T / mass)**(1/2),size=Natoms)
    #valors inicials
    for i in range(Natoms): #crea posicions random
        x = L*random()-L/2
        y = L*random()-L/2
        z = L*random()-L/2
        apos.append(vec(x,y,z))
        theta = pi*random() #crea vectors de moment random a partir del moment mitjà
        phi = 2*pi*random()
        px = p_sp[i]*sin(theta)*cos(phi)
        py = p_sp[i]*sin(theta)*sin(phi)
        pz = p_sp[i]*cos(theta)
        p.append(vector(px,py,pz))
    
    #temperature check
    velocities1 = [i**2 for i in p_sp]
    v_mean1 = np.mean(velocities1)
    T_meani = v_mean1 * mass / (3*k)
    T_means1 = [T_meani]
    
    p = [i*mass for i in p]
    #loop
    
    resultats=[]
    canvi = 0
    while canvi<= canvis:
        canvi += 1
        voltes = 0
        Lx+=L/20
        P_loop = []
        #temperature check
        while voltes<it_final:
            voltes+=1
            if abs(T_meani-T)>T_lim: break
            for i in range(Natoms): apos[i] = apos[i] + (p[i]/mass)*dt
            # Check for collisions
            hitlist = checkCollisions(apos)
        
            # If any collisions took place, update momenta of the two atoms
            for ij in hitlist:
                i = ij[0]
                j = ij[1]
                ptot = p[i]+p[j]
                posi = apos[i]
                posj = apos[j]
                vi = p[i]/mass
                vj = p[j]/mass
                vrel = vj-vi
                a = vrel.mag2
                if a == 0: continue;  # exactly same velocities
                rrel = posi-posj
                if rrel.mag > Ratom: continue # one atom went all the way through another
            
                # theta is the angle between vrel and rrel:
                dx = dot(rrel, vrel.hat)       # rrel.mag*cos(theta)
                dy = cross(rrel, vrel.hat).mag # rrel.mag*sin(theta)
                # alpha is the angle of the triangle composed of rrel, path of atom j, and a line
                #   from the center of atom i to the center of atom j where atome j hits atom i:
                alpha = asin(dy/(2*Ratom)) 
                d = (2*Ratom)*cos(alpha)-dx # distance traveled into the atom from first contact
                deltat = d/vrel.mag         # time spent moving from first contact to position inside atom
                
                posi = posi-vi*deltat # back up to contact configuration
                posj = posj-vj*deltat
                mtot = 2*mass
                pcmi = p[i]-ptot*mass/mtot # transform momenta to cm frame
                pcmj = p[j]-ptot*mass/mtot
                rrel = norm(rrel)
                pcmi = pcmi-2*pcmi.dot(rrel)*rrel # bounce in cm frame
                pcmj = pcmj-2*pcmj.dot(rrel)*rrel
                p[i] = pcmi+ptot*mass/mtot # transform momenta back to lab frame
                p[j] = pcmj+ptot*mass/mtot
                apos[i] = posi+(p[i]/mass)*deltat # move forward deltat in time
                apos[j] = posj+(p[j]/mass)*deltat
            velocities1 = [i.mag2/mass**2 for i in p]
            v_mean1 = np.mean(velocities1)
            T_mean1 = v_mean1 * mass / (3*k)
            T_means1.append(T_mean1)
            #wallhitcheck
            P_partial = []
            for i in range(Natoms):
                loc = apos[i]
                if loc.x < -L/2: 
                    p[i].x =  abs(p[i].x) 
                    if voltes > 75: P_partial.append(2* abs(p[i].x) /dt)
                if loc.x > Lx/2:
                    p[i].x =  -abs(p[i].x)
                    if voltes > 75: P_partial.append(2* abs(p[i].x) /dt)
                
                if abs(loc.y) > L/2:
                    if loc.y < 0: 
                        p[i].y = abs(p[i].y)
                        if voltes > 75: P_partial.append(2* abs(p[i].y) /dt)
                    else: 
                        p[i].y =  -abs(p[i].y)
                        if voltes > 75: P_partial.append(2* abs(p[i].y) /dt)
                
                if abs(loc.z) > L/2:
                    if loc.z < 0: 
                        p[i].z =  abs(p[i].z)
                        if voltes > 75: P_partial.append(2* abs(p[i].z) /dt)
                    else: 
                        p[i].z =  -abs(p[i].z)
                        if voltes > 75: P_partial.append(2* abs(p[i].z) /dt)
            if voltes > 75: P_loop.append(sum(P_partial) / ((4*L**2)+ (2*L*Lx)))
        #print(P_loop)
        if abs(T_meani-T)>T_lim: 
            canvi -= 1
            Lx-=L/20
            break
        print('Canvi', canvi, ' | Temperatura', T_meani, 'K')
        if voltes == it_final: P_error = (abs(Natoms*k*T_meani/(L**2*(L+Lx)/2) - np.mean(P_loop))/(Natoms*k*T_meani/(L**2*(L+Lx)/2))*100) 
        resultats.append((np.mean(P_loop), L*L*(L+Lx)/2, T_meani, P_error, np.std(P_loop)))
              
    if abs(T_meani-T)>T_lim: return [(0, 0, T_meani)]
    else: 
        plt.xlabel('Iteracions')
        plt.ylabel('Temperatura (K)')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-4,4), useOffset = False, useMathText=True)
        plt.tick_params(axis='x', direction='in')
        plt.tick_params(axis='y', direction='in')
        plt.plot(np.linspace(0, it_final*canvis, len(T_means1)), T_means1, color='coral')
        plt.savefig('T_Q4.pdf')
        plt.close()
        return resultats

class Regr:
  def xs(self, valor, num_xs):
    return int(valor*10**(num_xs))/10**(num_xs)
class Regr(Regr):
  def pend(self, a, b):
    import numpy as np
    from sklearn.linear_model import LinearRegression
    x = np.array(a).reshape((-1, 1))
    y = np.array(b)
    model = LinearRegression().fit(x, y)
    r_sq = model.score(x, y)
    return self.xs(float(model.coef_),25)
  def i_pend(self, a, b):
    import numpy as np
    from sklearn.linear_model import LinearRegression
    x = np.array(a).reshape((-1, 1))
    y = np.array(b)
    S_xx=np.std(x)**2*(len(x)-1)
    S_yy=np.std(y)**2*(len(y)-1)
    model = LinearRegression().fit(x, y)
    S_x_y=np.sqrt((S_yy-(S_xx*float(model.coef_)**2))/(len(x)-2))
    inc_pen = S_x_y/np.sqrt(S_xx)
    return self.xs(inc_pen, 25)
  def cf_reg(self, a, b):
    import numpy as np
    from sklearn.linear_model import LinearRegression
    x = np.array(a).reshape((-1, 1))
    y = np.array(b)
    model = LinearRegression().fit(x, y)
    r_sq = model.score(x, y)
    return self.xs(r_sq, 5)
  def graf(self, a, b, e_x, e_y, tit_a, tit_b):
    import matplotlib.pyplot as plt
    from sklearn.linear_model import LinearRegression
    import numpy as np
    x = np.array(a).reshape((-1, 1))
    y = np.array(b)
    model = LinearRegression().fit(x, y)
    y_pred = model.predict(x)
    fig = plt.figure(figsize=(6.5,4))
    ax = fig.add_subplot(111)
    ax.set(xlabel=tit_a, ylabel=tit_b)
    plt.tick_params(axis='x', direction='in')
    plt.tick_params(axis='y', direction='in')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(-1,2), useOffset = False, useMathText=True) #or style plain
    for tick in ax.yaxis.get_major_ticks():
      tick.label.set_fontsize(8)  #mida lletres eix y
    ax.plot(x,y_pred, color='coral', lw=2, linestyle=(0, (5, 4)), zorder=1)
    ax.scatter(x,y, color= 'teal', lw=0.3, zorder=3)
    for i in range(len(x)):
      ax.errorbar(x[i], y[i], xerr=e_x[i], yerr=e_y[i],capsize=3 ,color='k', alpha =0.8,zorder=2)
    plt.savefig('regrQ4.pdf')
    plt.close()
it_final=225
canvis=20
T_lim = 0.2
pendent_lim = 6.0
funcio = ViP(mass, k, T, Natoms, Ratom, L, it_final, canvis, T_lim)
while abs(funcio[0][2]-T)>T_lim or (abs(Regr().pend([1/funcio[i][1] for i in range(len(funcio))], [funcio[i][0] for i in range(len(funcio))])-np.mean([Natoms*k*funcio[i][2] for i in range(len(funcio))]))/np.mean([Natoms*k*funcio[i][2] for i in range(len(funcio))])*100)>pendent_lim:
    funcio = ViP(mass, k, T, Natoms, Ratom, L, it_final, canvis, T_lim)
#print(funcio2)
errx = [0.0 for i in range(len(funcio))]
erry = [0.0 for i in range(len(funcio))]
erryG = [funcio[i][4] for i in range(len(funcio))]
mitjana_errors = np.mean(erryG)
print("La mitjana d'errors en la pressió és de", mitjana_errors, 'Pa')
print('Pendent:',Regr().pend([1/funcio[i][1] for i in range(len(funcio))], [funcio[i][0] for i in range(len(funcio))]))
print('Incertesa pendent:', Regr().i_pend([1/funcio[i][1] for i in range(len(funcio))], [funcio[i][0] for i in range(len(funcio))]))
print('Pendent teòric =', np.mean([Natoms*k*funcio[i][2] for i in range(len(funcio))]))
print('% Error pendent =', abs(Regr().pend([1/funcio[i][1] for i in range(len(funcio))], [funcio[i][0] for i in range(len(funcio))])-np.mean([Natoms*k*funcio[i][2] for i in range(len(funcio))]))/np.mean([Natoms*k*funcio[i][2] for i in range(len(funcio))])*100)
print('Coeficient de regressió:',Regr().cf_reg([1/funcio[i][1] for i in range(len(funcio))], [funcio[i][0] for i in range(len(funcio))]))
Regr().graf([1/funcio[i][1] for i in range(len(funcio))], [funcio[i][0] for i in range(len(funcio))], errx, erry, '1/V', 'P')

#taula dades
from tabulate import tabulate
#import latextable
titol = ['1/V', 'Pressió (Pa)']
rows = [titol]
for i in range(len(funcio)):
    afegit = [1/funcio[i][1], funcio[i][0]]
    rows.append(afegit)
#print('Tabulate Table:')
#print(tabulate(rows, headers='firstrow'))
print('\nTabulate Latex:')
print(tabulate(rows, headers='firstrow', tablefmt='latex'))
