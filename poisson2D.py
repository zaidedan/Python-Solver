"""
Solve u_xx + u_yy = f(x,y) on a uniform grid
  and v_xx + v_yy = f(x,y) on a uniform grid

with known solution

u = exp(xy)         f['u'] = (x^2+y^2)*exp(xy)
v = exp((x^2+y^2))  f['v'] = (4*(x^2+y^2+1))*exp((x^2+y^2)) 

This is discretized using fluxes, for example, denoting u^{i,j} as the cell
centered at (x^{i,j},y^{i,j})

u^{i+1/2,j}_x = u^{i+1,j}-u^{i,j})/(Delta x)

for u^{i+1/2,j} the midpoint on an edge of the grid

We can form fluxes across the four edges of a cell to
sum their contributions to u^{i,j}, since
u^{i,j}_xx = u^{i+1/2,j}_x-u^{i-1/2,j}_x)/(Delta x)
And put it all together to get

u^{i+1/2,j}_x-u^{i-1/2,j}_x)/(Delta x)+ u^{i,j+1/2}_y-u^{i,j-1/2}_y)/(Delta y) - (x^{i,j}^2+y^{i,j}^2)*exp(x^{i,j}y^{i,j})

Our flux functions are then a simple difference, 
where P is a dictionary of parameters, in this case, with P['d'] = Delta x (or Delta y),
the grid spacing. Both 'u' and 'v' have similar forms, and so the same function is used

def difference(B,N,P):
  return dict((s,(N[s]-B[s])/P['d']) for s in B.state)


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


""" Fluxes defined here """
def difference(B,N,P):
  return dict((s,(N[s]-B[s])/P['d']) for s in B.state)

def poisson2D(N):
  """ 
  lets define a uniform square mesh on [-1, 1] x [1, 1]
  and create boundary blocks as we go,
  initializing based on the exact solution, and naming the block by its coordinates
  """
  d = 2./float(N) # spacing, Delta x = Delta y
  # For N blocks in each direction, N+2 blocks since we have a constant block on each end surrounding the domain
  # a ghost cell approach, there are no parameters in this case, or parameter functions needed
  # the center of the cell is at (i*d+d/2,j*d+d/2)
  B = []
  for i in range(-1,N+1):
    for j in range(-1,N+1):
      x = (i*d+d/2)
      y = (j*d+d/2)
      # block has a simple name
      name = '('+str(i)+','+str(j)+')'
      B.append(b.Block(name,
        {'u':math.exp(x*y),
         'v':math.exp(x**2+y**2)}, # initialize to exact solution
         None, # no parameter functions
         {'x':x,'y':y})) # parameter of coordinates
  # interior sources, -f['u'] and -f['v']
  for block in B:
    (x,y) = (block.p['x'],block.p['y'])
    block.addSource(s.Source(s.constant, # standard function defined in source.py
      {'u':-(x*x+y*y)*math.exp(x*y),
      'v':-4.0*(x*x+y*y+1.0)*math.exp(x*x+y*y)}, # parameters
      'constant')) # name

  # Flux geometry 
  P = {'d':d*d} # divide by delta x^2 in the end
  n = N+2 # add two for the boundaries
  for i in range(1,N+1):
    for j in range(1,N+1):
      # Add fluxes, figuring out which neighbors to connect to
      for k in [(i-1)*n+j, i*n+j-1,(i+1)*n+j, i*n+j+1]:
        B[i*n+j].addFlux(f.Flux(B[k],difference,P))

  interiorBlocks = [B[i*n+j] for i in range(1,N+1) for j in range(1,N+1)]

  # solve the problem on the interior blocks
  P = p.Problem(interiorBlocks)
  P.solve()
  # compute the L-2 error against the exact solution for both variables
  Eu = math.sqrt(sum([(math.exp(block.p['x']*block.p['y'])-block['u'])**2 for block in interiorBlocks])/(n-2)/(n-2))
  Ev = math.sqrt(sum([(math.exp(block.p['x']**2+block.p['y']**2)-block['v'])**2 for block in interiorBlocks])/(n-2)/(n-2))
  return (Eu,Ev)

def test():
  n = 3
  Error = [poisson2D(n),poisson2D(n*2)]
  # do a quick check of convergence rate of error, should be > 2
  Rate = [(math.log(Error[1][0])-math.log(Error[0][0]))/(math.log(2./(2*n))-math.log(2./(n))),
  (math.log(Error[1][1])-math.log(Error[0][1]))/(math.log(2./(2*n))-math.log(2./(n)))]
  return Rate
