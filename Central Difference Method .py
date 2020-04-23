#!/usr/bin/env python
# coding: utf-8

# In[1]:


#%matplotlib notebook
from IPython.display import Image
import matplotlib.pyplot as plt
import numpy as np
#import scipy.linalg as linalg


# # Central Difference Method

# In[2]:


# properties
m = 0.2533 #k-s^2/in.
k = 10 #k/in
Tn = 1 #s
wn = 2*np.pi/Tn
#c = 0.1592
zeta = 0.0
c = zeta*2*m*wn
c


# In[3]:


# loading function
def pDef(t):
    if t < 0.600001:
        p = (10/0.6*t)
    else:
        p = 0.0
    return p


# In[4]:


# numerical properties
dt = 0.1 #s
T = 3 #s end time
nsteps = int(round(T/dt))
t = np.linspace(0.0,T,nsteps+1)


# In[5]:


# constants
kh = m/dt**2 + c/(2*dt)
a = m/dt**2 - c/(2*dt)
b = k-2*m/dt**2


# In[6]:


p = np.zeros((nsteps+1))
for j in range(0,nsteps):
    p[j] = pDef(t[j])

p


# In[7]:


# initialize
u0 = 0
v0 = 0
a0 = (p[0]-c*v0-k*u0)/m
un = u0-dt*v0+dt**2/2*a0

# store
u = np.zeros((nsteps+1))
u[0] = u0
ph = p[0]-a*un-b*u[0]
u[1] = ph/kh
u


# In[8]:


# loop over time steps
for j in range(0,nsteps):
    ph = p[j]-a*u[j-1]-b*u[j]
    u[j+1] = ph/kh
u


# In[9]:


# plot
plt.plot(t, p, 'b-o', label=r'interpolated loading')
plt.legend(loc='upper right')
plt.ylabel('$p$ kips')
plt.xlabel('$t$ [s]')
plt.show()


# In[10]:


# define more parameters
td = 0.6

# define exact solution
te = np.linspace(0.0,T,100)
ue = np.zeros((100))
for j in range(1,100):
    if te[j]<td:
        ue[j] = te[j]/td-(Tn/(2*np.pi*td))*np.sin(2*np.pi*te[j]/Tn)
    else:
        ue[j] = np.cos(2*np.pi/Tn*(te[j]-td))+ (Tn/(2*np.pi*td))*np.sin(2*np.pi*(te[j]-td)/Tn)-(Tn/(2*np.pi*td))*np.sin(2*np.pi*te[j]/Tn)

ue


# In[11]:


# plot
plt.plot(te, ue, 'r', label=r'exact loading')
plt.plot(t, u, 'b-o', label=r'approximated loading')
plt.legend(loc='upper right')
plt.ylabel('$u$ in.')
plt.xlabel('$t$ [s]')
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




