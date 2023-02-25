from vpython import *
import numpy as np
import matplotlib.pyplot as plt


#win = 500
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

d = L/2+Ratom
r = 0.005


"""
deltav = 100 # binning for v histogram

def barx(v):
    return int(v/deltav) # index into bars array


#print(accum)
def interchange(v1, v2):  # remove from v1 bar, add to v2 bar
    barx1 = barx(v1)
    barx2 = barx(v2)
    if barx1 == barx2:  return
    if barx1 >= len(histo) or barx2 >= len(histo): return
    histo[barx1] -= 1
    histo[barx2] += 1
 """ 
def checkCollisions():
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


def PiT(mass, k, T, Natoms, Ratom, it_final, T_lim):
    global apos
    global L
    L = 1
    #box
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
        theta = pi*random() #crea vectors de velocitat random a partir del moment mitjà
        phi = 2*pi*random()
        px = p_sp[i]*sin(theta)*cos(phi)
        py = p_sp[i]*sin(theta)*sin(phi)
        pz = p_sp[i]*cos(theta)
        p.append(vector(px,py,pz))
    
    
    #temperature check
    velocities1 = [i**2 for i in p_sp]
    v_mean1 = np.mean(velocities1)
    T_mean1 = v_mean1 * mass / k
    T_means1 = []
    T_means1.append(T_mean1/3)
    
    
    p = [i*mass for i in p]
    P_loop = []
    #loop
    voltes = -1
    P_error = 0
    if abs(T_means1[0]-T) < T_lim: print(f"Fent PiT per T={T}")
    while voltes <= it_final:
        voltes+=1
        #print(voltes)
        if abs(T_means1[0]-T)>T_lim: 
            P_error=0
            break
        for i in range(Natoms): apos[i] = apos[i] + (p[i]/mass)*dt
        # Check for collisions
        hitlist = checkCollisions()
    
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
        #wallhitcheck
        F_partial = []
        
        for i in range(Natoms):
            loc = apos[i]
            if abs(loc.x) > L/2:
                if loc.x < 0: 
                    p[i].x =  abs(p[i].x) 
                    F_partial.append(2* abs(p[i].x) /dt)
                else: 
                    p[i].x =  -abs(p[i].x)
                    F_partial.append(2* abs(p[i].x) /dt)
            
            if abs(loc.y) > L/2:
                if loc.y < 0: 
                    p[i].y = abs(p[i].y)
                    F_partial.append(2* abs(p[i].y) /dt)
                else: 
                    p[i].y =  -abs(p[i].y)
                    F_partial.append(2* abs(p[i].y) /dt)
            if abs(loc.z) > L/2:
                if loc.z < 0: 
                    p[i].z =  abs(p[i].z)
                    F_partial.append(2* abs(p[i].z) /dt)
                else: 
                    p[i].z =  -abs(p[i].z)
                    F_partial.append(2* abs(p[i].z) /dt)
        P_loop.append(sum(F_partial) / (6*L**2))
        P_error = 0
        if voltes == 100: 
            P_error = (abs(Natoms*k*np.mean(T_means1)/L**3 - np.mean(P_loop))/(Natoms*k*np.mean(T_means1)/L**3)*100)
            print('100P_error=', P_error)
        #if P_error > P_lim : break
        
        if voltes % 20 == 0:
            velocities1 = [i.mag2/mass**2 for i in p]
            v_mean1 = np.mean(velocities1)
            T_mean1 = v_mean1 * mass / (1.4E-23)
            T_means1.append(T_mean1/3)
    if abs(T_means1[0]-T) > T_lim: 
        return (0, 0, T_means1[0], P_error)
    P_error = (abs(Natoms*k*np.mean(T_means1)/L**3 - np.mean(P_loop))/(Natoms*k*np.mean(T_means1)/L**3)*100)
    print('P_error_final=', P_error)       
    plt.xlabel('Iteracions')
    plt.ylabel(u"Pressió (Pa)")
    plt.tick_params(axis='x', direction='in')
    plt.tick_params(axis='y', direction='in')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useOffset = False, useMathText=True)
    plt.plot(np.linspace(0, voltes, len(P_loop)), P_loop)
    plt.savefig(f'Pressio T {T} K.pdf')
    plt.close()
    print(f"PiT fet per T={T}")
    return (np.mean(T_means1), np.mean(P_loop), T_means1[0], P_error, np.std(P_loop))
###########################################################################################################################
def ViP(mass, k, T, Natoms, Ratom, L, it_final, T_lim):
    global apos
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
    T_meansi = v_mean1 * mass / (3*k)
    
    p = [i*mass for i in p]
    P_loop=[]
    #loop
    voltes = -1
    P_error = 0
    if abs(T_meansi-T)<T_lim: print(f"Fent lnViP per V={L**3}")
    while voltes<=it_final:
        voltes+=1
        if abs(T_meansi-T)>T_lim: 
            break
        for i in range(Natoms): apos[i] = apos[i] + (p[i]/mass)*dt
        # Check for collisions
        hitlist = checkCollisions()
    
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
        #wallhitcheck
        P_partial = []
        for i in range(Natoms):
            loc = apos[i]
            if abs(loc.x) > L/2:
                if loc.x < 0: 
                    p[i].x =  abs(p[i].x) 
                    P_partial.append(2* abs(p[i].x) /dt)
                else: 
                    p[i].x =  -abs(p[i].x)
                    P_partial.append(2* abs(p[i].x) /dt)
            
            if abs(loc.y) > L/2:
                if loc.y < 0: 
                    p[i].y = abs(p[i].y)
                    P_partial.append(2* abs(p[i].y) /dt)
                else: 
                    p[i].y =  -abs(p[i].y)
                    P_partial.append(2* abs(p[i].y) /dt)
            
            if abs(loc.z) > L/2:
                if loc.z < 0: 
                    p[i].z =  abs(p[i].z)
                    P_partial.append(2* abs(p[i].z) /dt)
                else: 
                    p[i].z =  -abs(p[i].z)
                    P_partial.append(2* abs(p[i].z) /dt)
        P_loop.append(sum(P_partial) / (6*L**2))
        if voltes == 100: P_error = (abs(Natoms*k*np.mean(T_meansi)/L**3 - np.mean(P_loop))/(Natoms*k*np.mean(T_meansi)/L**3)*100)
    if voltes == it_final: P_error = (abs(Natoms*k*np.mean(T_meansi)/L**3 - np.mean(P_loop))/(Natoms*k*np.mean(T_meansi)/L**3)*100)
    if abs(T_meansi-T) > T_lim: return (0, 0, T_meansi)
    else: 
        print(f"lnViP fet per V={L**3}")
        return (np.mean(P_loop), -np.log(L**3), T_meansi, P_error, np.std(P_loop))
#falta canviar codi per començar de nou a la que se sap que T no val

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
  def graf(self, a, b, e_x, e_y, tit_a, tit_b, nom):
    import matplotlib.pyplot as plt
    from sklearn.linear_model import LinearRegression
    import numpy as np
    x = np.array(a).reshape((-1, 1))
    y = np.array(b)
    model = LinearRegression().fit(x, y)
    y_pred = model.predict(x)
    fig = plt.figure(figsize=(6.5,4))
    ax = fig.add_subplot(111)
    ax.set(xlabel=f"{tit_a}", ylabel=f"{tit_b}")
    plt.tick_params(axis='x', direction='in')
    plt.tick_params(axis='y', direction='in')
    if y[0]/1E-16 < 10:
        plt.ticklabel_format(axis='y', style='sci', scilimits=(-18,-18), useOffset = False, useMathText=True)
    if x[0]/1E-16 < 10:
        plt.ticklabel_format(axis='x', style='sci', scilimits=(-18,-18), useOffset = False, useMathText=True)
    for tick in ax.yaxis.get_major_ticks():
      tick.label.set_fontsize(8)  #mida lletres eix y
    ax.plot(x,y_pred, color='coral', lw=2, linestyle=(0, (5, 4)), zorder=1)
    ax.scatter(x,y, color= 'teal', lw=0.3, zorder=3)
    for i in range(len(x)):
      ax.errorbar(x[i], y[i], xerr=e_x[i], yerr=e_y[i],capsize=3 ,color='k', alpha =0.8,zorder=2)
    plt.savefig(f"simulacio{nom}Q3.pdf")
    plt.close()

it_final = 225
T_lim = 0.3
alphes = []
inv_T = []
for j in range(1,6):
    T=100*j
    inv_T.append(1/T)

    beta_ = []
    errxbeta = []
    errybeta = []
    errybetaG = []
    for i in range(-6,6,1):
        funcio = PiT(mass, k, T+i, Natoms, Ratom, it_final, T_lim)
        while abs(funcio[2]-(T+i))>T_lim: 
            funcio = PiT(mass, k, T+i, Natoms, Ratom, it_final, T_lim)
        beta_.append((funcio[0],funcio[1],funcio[2]))
        #errxbeta.append(funcio[2]-(T+i))
        errybetaG.append(funcio[4])
        errxbeta.append(0)
        errybeta.append(0)
    betax = [beta_[i][0] for i in range(len(beta_))]
    betay = [beta_[i][1] for i in range(len(beta_))]
    pendent_beta = Regr().pend(betax, betay)
    inc_pendent_beta = Regr().i_pend(betax, betay)
    coef_regr_beta = Regr().cf_reg(betax, betay)
    k_T = []
    errxk_T = []
    errxk_TG = []
    erryk_T = []
    for i in range(-6, 6, 1):
        funcio2 = ViP(mass, k, T, Natoms, Ratom, L+(i/50), it_final, T_lim)
        while abs(funcio2[2]-T)>T_lim:
            funcio2 = ViP(mass, k, T, Natoms, Ratom, L+(i/100), it_final, T_lim)
        k_T.append((funcio2[0], funcio2[1]))
        errxk_TG.append(funcio2[4])
        erryk_T.append(0)
        errxk_T.append(0)
    k_Tx = [k_T[i][0] for i in range(len(k_T))]
    k_Ty = [k_T[i][1] for i in range(len(k_T))]
    pendent_k_T = Regr().pend(k_Tx, k_Ty)
    inc_pendent_k_T = Regr().i_pend(k_Tx, k_Ty)
    coef_regr_k_T = Regr().cf_reg(k_Tx, k_Ty)
    coef_dilat = pendent_beta * pendent_k_T
    error_coef_dilat = abs((pendent_beta * pendent_k_T)-(1/T))/(1/T) * 100
    print('Pendent beta:', pendent_beta)
    print('Incertesa pendent beta:', inc_pendent_beta)
    #print('Ordenada a l\'origen:',Regr().ordor(beta_T, beta_P))
    print(u'Coeficient de regressió beta:',coef_regr_beta)
    Regr().graf(betax, betay, errxbeta, errybeta, 'Temperatura (K)', u'Pressió (Pa)', j)
    
    print('Pendent k_T:',pendent_k_T)
    print('Incertesa pendent k_T:', inc_pendent_k_T)
    print(u'Coeficient de regressió k_T:', coef_regr_k_T)
    Regr().graf(k_Tx, k_Ty, errxk_T, erryk_T, u'Pressió (Pa)', '-ln V', 7*j)
    
    print(u'Coef. dilatació =', coef_dilat)
    alphes.append((coef_dilat, error_coef_dilat))
    print('%err_alpha =', error_coef_dilat)
    #plt.savefig(figure_name)
    #print(Natoms*k/L**3)
    print(np.mean(errxk_TG))
    print(np.mean(errybetaG))
alphesy = [alphes[i][0] for i in range(len(alphes))]
errors = [0 for i in range(len(alphes))]
print('Pendent alpha vs 1/T:', Regr().pend(inv_T, alphesy))
print('inc. pendent alphes:', Regr().pend(inv_T, alphesy))
Regr().graf(inv_T, alphesy, errors, errors, '1/T', r'$\alpha$', 10000000)
#taula dades
from tabulate import tabulate
#import latextable
titol = ['1/T', 'alpha']
rows = [titol]
for i in range(len(alphesy)):
    afegit = [inv_T[i], alphesy[i]]
    rows.append(afegit)
#print('Tabulate Table:')
#print(tabulate(rows, headers='firstrow'))
print('\nTabulate Latex:')
print(tabulate(rows, headers='firstrow', tablefmt='latex'))
