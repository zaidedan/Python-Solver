"""
Solve (u_xx + u_yy)*nu_u = u_t on a uniform grid
and (v_xx + v_yy)*nu_v = v_t on a uniform grid

with exact solution 
u, v = exp(-nu*pi^2*(a^2+b^2)t)*sin(a*pi*x)*sin(b*pi*y)

for initial condition
u, v = sin(a*pi*x)*sin(b*pi*y)

for integer choices a and b
this allows for a lot of flexibility 
in how things want to be defined

lets pick 
u_a = 1, u_b = 1, nu_u = 1
v_a = 4, v_b = 2, nu_v = 0.5
tf = 1, the final time

These are independent of eachother, but are sufficient to check the
unsteady solver is working

The discretization is the same as in poisson2D, with an additional
du/dt = R(U), handled by odeint.

the main difference is that the time of a block (.t) is used in the 
boundary blocks, which are defined by a single source returning a dictionary

boundary block state = S(boundary block)

This is very much a work in progress and not actively used, but it could be...

The code does a rough accuracy estimation,
and outputs the accuracy for both u and v. The accuracy should be around 2 if 
the code is performing correctly.

"""
import math as math
import sys
import src.blocks as b
import src.flux as f
import src.problem as p
import src.source as s
import numpy as np

""" 
lets define a uniform square mesh on [-1, 1] x [1, 1]
and create boundary blocks as we go,
initializing based on the exact solution, and naming the block by its coordinates
"""
u_a = 1
u_b = 1
nu_u = 1
v_a = 1
v_b = 1
nu_v = 0.5
tf = 1

""" Sources defined here """
def time(B,P):
  u = math.exp(-nu_u*math.pi*math.pi*(u_a*u_a+u_b*u_b)*B.t)*math.sin(u_a*math.pi*B.p['x'])*math.sin(u_b*math.pi*B.p['y'])
  v = math.exp(-nu_v*math.pi*math.pi*(v_a*v_a+v_b*v_b)*B.t)*math.sin(v_a*math.pi*B.p['x'])*math.sin(v_b*math.pi*B.p['y'])
  return {'u':u,'v':v}

""" Fluxes defined here """
def difference(B,N,P):
  return dict((s,(N[s]-B[s])/P['d']) for s in B.state)

def diffusion2D(N):
  d = 2./float(N) # spacing, delta 
  # initialize with exact solution at t = 0
  B = []
  for i in range(-1,N+1):
    for j in range(-1,N+1):
      x = (i*d+d/2)
      y = (j*d+d/2)
      # block has a simple name
      name = '('+str(i)+','+str(j)+')'
      B.append(b.Block(name,
        # initialize to exact solution
        {'u':math.sin(u_a*math.pi*x)*math.sin(u_b*math.pi*y),
         'v':math.sin(v_a*math.pi*x)*math.sin(v_b*math.pi*y)}, 
         None, # no parameter functions
         {'x':x,'y':y})) 

  # Flux geometry 
  P = {'d':d*d}
  n = N+2 # add two for the boundaries
  for i in range(1,n-1):
    for j in range(1,n-1):
      # Add fluxes, looping around the edges
      for k in [(i-1)*n+j, i*n+j-1,(i+1)*n+j, i*n+j+1]:
        B[i*n+j].addFlux(f.Flux(B[k],difference,P))

  interiorBlocks = [B[i*n+j] for i in range(1,n-1) for j in range(1,n-1)]

  # sort out the boundary blocks
  boundaryBlocks = []
  bcRange = [j         for j in range(1,n-1)] + \
            [(n-1)*n+j for j in range(1,n-1)] + \
            [i*n       for i in range(0,n)] + \
            [i*n+n-1   for i in range(0,n)] 

  for k in bcRange:
    B[k].addSource(s.Source(time,None,B[k].name))
    boundaryBlocks.append(B[k])

  # solve the problem on the interior blocks,
  # pass in the boundary blocks so they can be updated
  # at each needed time.
  # solve the unstead problem at the following timesteps
  P = p.Problem(interiorBlocks,boundaryBlocks)
  P.solveUnst(np.linspace(0,tf,10))
  
  # calculate the error for accuracy checking
  Eu = 0
  Ev = 0
  for bb in interiorBlocks:
    (x,y) = (bb.p['x'],bb.p['y'])
    ue = math.exp(-nu_u*math.pi*math.pi*(u_a*u_a+u_b*u_b)*tf)*math.sin(u_a*math.pi*x)*math.sin(u_b*math.pi*y)
    ve = math.exp(-nu_v*math.pi*math.pi*(v_a*v_a+v_b*v_b)*tf)*math.sin(v_a*math.pi*x)*math.sin(v_b*math.pi*y)
    Eu += (bb['u']-ue)**2
    Ev += (bb['v']-ve)**2

  return (math.sqrt(Eu/(n-2)/(n-2)),math.sqrt(Ev)/(n-2)/(n-2))

def test():
  n = 3
  Error = [diffusion2D(n),diffusion2D(n*2)]
  Rate = [(math.log(Error[1][0])-math.log(Error[0][0]))/(math.log(2./(2*n))-math.log(2./(n))),
  (math.log(Error[1][1])-math.log(Error[0][1]))/(math.log(2./(2*n))-math.log(2./(n)))]
  return Rate
