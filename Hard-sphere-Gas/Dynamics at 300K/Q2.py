from vpython import *
#GlowScript 3.0 VPython

# Hard-sphere gas.

# Bruce Sherwood


win = 500

Natoms = 400  # change this to have more or fewer atoms

# Typical values
L = 1 # container is a cube L on a side
gray = color.gray(0.7) # color of edges of container
mass = 4E-3/6E23 # helium mass
Ratom = 0.04 # wildly exaggerated size of helium atom, valor inicial 0.03
k = 1.4E-23 # Boltzmann constant
T = 300 # around room temperature
dt = 1E-5

animation = canvas( width=win, height=win, align='left')
animation.range = L


animation.title = 'A "hard-sphere" gas'
s = """  Theoretical and averaged speed distributions (meters/sec).
  Initially all atoms have the same speed, but collisions
  change the speeds of the colliding atoms. One of the atoms is
  marked and leaves a trail so you can follow its path.
  
"""
animation.caption = s

#mostra t bonic
def fancy_t(time):
    if time/60>1:
        minuts = int(time/60)
        segons = int((time/60 - minuts)*60*100)/100
        return str(minuts)+' min '+str(segons)+' s'
    else:
        return str(int(time*100)/100)+' s'
        
#interruptor
running = True
def Run(b):
        global running
        running = not running
        if running: b.text = "Stop"
        else: b.text = "Stopped"
    
button(text="Stop", pos=animation.caption_anchor, bind=Run)

#crono i foto
import time
tic = time.perf_counter()
timecounter=0

def Captura(b):
    toc = time.perf_counter()
    timeout = toc - tic
    animation.capture('scrshot_'+fancy_t(timeout))
button(text="Capture", pos=animation.caption_anchor, bind=Captura)

def Tempo(b):
    global timecounter
    toc = time.perf_counter()
    timeout = toc - tic
    output.text = 'Current time: '+fancy_t(timeout)
    #animation.capture(f"scrshot_t{timeout:0.4f}s")

button(text="Time", pos=animation.caption_anchor, bind=Tempo)
output = wtext(pos=animation.caption_anchor, text=f"Counting...")

d = L/2+Ratom
r = 0.005
#box
boxbottom = curve(color=gray, radius=r)
boxbottom.append([vector(-d,-d,-d), vector(-d,-d,d), vector(d,-d,d), vector(d,-d,-d), vector(-d,-d,-d)])
boxtop = curve(color=gray, radius=r)
boxtop.append([vector(-d,d,-d), vector(-d,d,d), vector(d,d,d), vector(d,d,-d), vector(-d,d,-d)])
vert1 = curve(color=gray, radius=r)
vert2 = curve(color=gray, radius=r)
vert3 = curve(color=gray, radius=r)
vert4 = curve(color=gray, radius=r)
vert1.append([vector(-d,-d,-d), vector(-d,d,-d)])
vert2.append([vector(-d,-d,d), vector(-d,d,d)])
vert3.append([vector(d,-d,d), vector(d,d,d)])
vert4.append([vector(d,-d,-d), vector(d,d,-d)])

Atoms = []
p = []
apos = []
from scipy.stats import maxwell
p_sp = maxwell.rvs(scale=(k * T / mass)**(1/2),size=Natoms)
#valors inicials
for i in range(Natoms): #crea posicions random
    x = L*random()-L/2
    y = L*random()-L/2
    z = L*random()-L/2
    if i == 0:
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.cyan,
                            make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else: Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z))
    theta = pi*random() #crea vectors de moment random a partir del moment mitjà
    phi = 2*pi*random()
    px = p_sp[i]*sin(theta)*cos(phi)
    py = p_sp[i]*sin(theta)*sin(phi)
    pz = p_sp[i]*cos(theta)
    p.append(vector(px,py,pz))
#print(p[100], p[200], p[300])

deltav = 100 # binning for v histogram

def barx(v):
    return int(v/deltav) # index into bars array

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
#for i in range(nhisto): histo.append(dic) #canviar aquesta part si es canvia les distribucions inicials
#histo[barx(pavg/mass)] = Natoms
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

import numpy as np
vel_distr = np.linspace(0,3001+dv,len(range((3001+dv))))
params = maxwell.fit(p_sp, floc=0)
probdist = Natoms *deltav* maxwell.pdf(vel_distr, *params)
for v in range(0,3001+dv,dv):
    theory.plot(v, probdist[v])
    

accum = []
for i in range(int(3000/deltav)): accum.append([deltav*(i+.5),0.0])
vdist = gvbars(color=color.red, delta=deltav )
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

snap = 0 # number of histogram snapshots to average

#temperature check
velocities1 = [i**2 for i in p_sp]
v_mean1 = np.mean(velocities1)
T_mean1 = v_mean1 * mass / k
T_means1 = []
T_means1.append(T_mean1/3)
print(T_means1)
velocitiesx = []
velocitiesy = []
velocitiesz = []
for i in p:
    velocityx = i.x
    velocitiesx.append(velocityx**2)
    velocityy = i.y
    velocitiesy.append(velocityy**2)
    velocityz = i.z
    velocitiesz.append(velocityz**2) 
v_meanx = np.mean(velocitiesx)
T_meanx = v_meanx * mass / k
T_meansx = []
T_meansx.append(T_meanx)

v_meany = np.mean(velocitiesy)
T_meany = v_meany * mass / k
T_meansy = []
T_meansy.append(T_meany)

v_meanz = np.mean(velocitiesz)
T_meanz = v_meanz * mass / k
T_meansz = []
T_meansz.append(T_meanz)


times = 0
p = [i*mass for i in p]
#loop


while True:
    if running:
        rate(300)
        times+=1
        # Accumulate and average histogram snapshots
        for i in range(len(accum)): accum[i][1] = (snap*accum[i][1] + histo[i])/(snap+1)
        if snap % 5 == 0:
            vdist.data = accum
        snap += 1
        
        # Update all positions
        for i in range(Natoms): Atoms[i].pos = apos[i] = apos[i] + (p[i]/mass)*dt
        
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
        
        for i in range(Natoms):
            loc = apos[i]
            if abs(loc.x) > L/2:
                if loc.x < 0: p[i].x =  abs(p[i].x)
                else: p[i].x =  -abs(p[i].x)
            
            if abs(loc.y) > L/2:
                if loc.y < 0: p[i].y = abs(p[i].y)
                else: p[i].y =  -abs(p[i].y)
            
            if abs(loc.z) > L/2:
                if loc.z < 0: p[i].z =  abs(p[i].z)
                else: p[i].z =  -abs(p[i].z)
        velocities1 = [i**2 for i in p_sp]
        v_mean1 = np.mean(velocities1)
        T_mean1 = v_mean1 * mass / (1.4E-23)
        T_means1.append(T_mean1/3)
        
        velocitiesx = []
        velocitiesy = []
        velocitiesz = []
        for i in p:
            velocityx = i.x/mass
            velocitiesx.append(velocityx**2) 
            velocityy = i.y/mass
            velocitiesy.append(velocityy**2)
            velocityz = i.z/mass
            velocitiesz.append(velocityz**2)
        v_meanx = np.mean(velocitiesx)
        T_meanx = v_meanx * mass / k
        T_meansx.append(T_meanx)

        v_meany = np.mean(velocitiesy)
        T_meany = v_meany * mass / k
        T_meansy.append(T_meany)
             
        v_meanz = np.mean(velocitiesz)
        T_meanz = v_meanz * mass / k
        T_meansz.append(T_meanz)
        
    else:
        break
toc = time.perf_counter()
timeout = toc - tic
output.text = "Total time: "+fancy_t(timeout)
animation.capture(f"scrshot_final_"+fancy_t(timeout))
import matplotlib.pyplot as plt
x1 = np.linspace(0, timeout, len(T_means1))  
xx = np.linspace(0, timeout, len(T_meansx))
xy = np.linspace(0, timeout, len(T_meansy))  
xz = np.linspace(0, timeout, len(T_meansz))          
y1 = T_means1
yx = T_meansx
yy = T_meansy
yz = T_meansz
plt.title('Temps total de '+fancy_t(timeout)+', '+str(times)+' iteracions')
plt.xlabel(r'$t$ (s)')
plt.ylabel(r'$T$ (K)')
plt.plot(x1, y1, label=r'$T(\langle v^2 \rangle/3)$')
plt.plot(xx, yx, label=r'$T(\langle v_x^2 \rangle)$')
plt.plot(xy, yy, label=r'$T(\langle v_y^2 \rangle)$')
plt.plot(xz, yz, label=r'$T(\langle v_z^2 \rangle)$')
plt.legend()
figure_name='T_expected_'+fancy_t(timeout)+'.pdf'
plt.savefig(figure_name)
