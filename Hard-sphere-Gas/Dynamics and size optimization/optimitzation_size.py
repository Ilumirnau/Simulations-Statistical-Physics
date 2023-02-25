from vpython import *
import numpy as np
import matplotlib.pyplot as plt

win = 500
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

deltav = 100 # binning for v histogram

def find_min(llista):
    minim = 100
    index = 'not found'
    for i in range(len(llista)):
        if llista[i] < minim:
            minim = llista[i]
            index = i
    return (index, minim)

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

def ViP(mass, k, T, Natoms, Ratom, L, it_final):
    global apos
    global histo
    global win
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
    
    lis = []
    dic = {}
    indexs = []
    for w in range(Natoms):
        indexs.append(barx(p_sp[w]))
    indexs.sort()
    
    for index in indexs:
        if index not in lis:
            lis.append(index)
            dic[str(index)] = 1
        else:
            dic[str(index)] = dic.get(str(index), 0) + 1
    #print(dic)
    
    nhisto = int(4500/deltav)
    histo = []

    for key in dic: histo.append(dic[key])
    for i in range(nhisto):
        if str(i) not in dic:
            histo.insert(i, 0.0)
    #print(histo, len(histo))
    #print(sum(histo), Natoms)
    gg = graph( width=.9*win, height=0.8*win, xmax=3000, align='left',
        xtitle='speed, m/s', ytitle='Number of atoms', ymax=Natoms*deltav/1000)
    
    theory = gcurve( color=color.cyan )
    dv = 10
    
    vel_distr = np.linspace(0,3001+dv,len(range((3001+dv))))
    params = maxwell.fit(p_sp, floc=0)
    probdist = Natoms *deltav* maxwell.pdf(vel_distr, *params)
    for v in range(0,3001+dv,dv):
        theory.plot(v, probdist[v])
        
    
    accum = []
    for i in range(int(3000/deltav)): accum.append([deltav*(i+.5),0.0])
    vdist = gvbars(color=color.red, delta=deltav )
    #temperature check
    velocities1 = [i**2 for i in p_sp]
    v_mean1 = np.mean(velocities1)
    T_meansi = v_mean1 * mass / (3*k)
    
    p = [i*mass for i in p]
    P_loop=[]
    #loop
    voltes = -1
    snap = 0
    while voltes<=it_final:
        voltes+=1
        for i in range(len(accum)): accum[i][1] = (snap*accum[i][1] + histo[i])/(snap+1)
        if snap % 5 == 0:
            vdist.data = accum
        snap += 1
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
            interchange(vi.mag, p[i].mag/mass)
            interchange(vj.mag, p[j].mag/mass)
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
    print(f"infoP fet per Ratom={Ratom}")
    return (np.mean(P_loop), T_meansi, P_error, np.std(P_loop))
#falta canviar codi per començar de nou a la que se sap que T no val
it_final = 250
minims_finals = []
for j in range(1,10):
    errors_totals = []
    radis = []
    for i in range(1,10,1):
        print('Fent Ratom =', Ratom + i*.005, f'| Volta {j}')
        info_P = ViP(mass, k, T, Natoms, Ratom + i*.005, L, it_final)
        errors_totals.append(info_P[2])
        radis.append(Ratom + i*.005)
    minims_finals.append(find_min(errors_totals))
    plt.xlabel('Radis')
    plt.ylabel('% Error Pressió')
    plt.plot(radis, errors_totals)
    plt.savefig(f'P_errors_volta_{j}.pdf')
print(minims_finals)
contador_indexs = {i[0]:minims_finals.count(i[0]) for i in minims_finals}
print(contador_indexs)
print(Ratom + 3*.005, Ratom + 4*.005)
