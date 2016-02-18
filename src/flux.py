""" 
flux.py contains a list of flux functions 

This is meant to serve as a wrapper
for all flux functions

Each flux computes the heat transfer rate between two blocks
Each flux has:
  The block it belongs to (.B)
  The block its connected to (.N)
  The type (.f), this is a pointer to the function
  The geometry between the left and right state (.P) (or parameters)

Each flux function is f(B,N,P), a function of a block,
its neighbor, and its parameters. For example, consider
a finite difference scheme, du/dx = (u_i-u_{i-1})/dx
we could write the flux as

def fluxFunction(B,N,P):
  f = (B['u']-N['u'])*(P['dx'])
  return {'u':f}

u1 = Block('u1',{u:'1'})
u2 = Block('u2',{u:'1'})

u1_flux = Flux(u2,fluxFunction,{'dx':1},'u_flux')

u1.addFlux(u1_flux)

"""

class Flux(object):
  """
  Flux object. Each flux is effectively a boundary condition

  __init__:   Flux Constructor

  input(s):   (N) Neighboring block
              (f) Flux function name
              (P) Parameters
              (name) identifying name

  output(s):  None
  """
  def __init__(self,N,f,P=None,name=''):
    self.B = None # this will be set when its added to the block
    self.N = N
    self.F = f 
    self.P = P
    self.name = name

  """
  The following is a wrapper for flux function choices defined with
  Fluxes should be define as such
  input(s):   (B) Block
              (N) Neighboring block
              (P) Parameters (optional) 
  output(s):  dict with entries for the flux for each state contributed to

  Fluxes can have blocks with different physical states
  provided the flux function passes contributions to some of the states

  """

  def flux(self):
    return self.F(self.B,self.N,self.P)

  """
  __repr__    overloading for print command

  input(s):   None
  output(s):  name
  """
  def __repr__(self):
    return name+" "+self.F.__name__