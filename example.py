"""

This file contains both implementations of the example system
in examples.tex

C_f(u_1 - u_0) + C_su_1 = 0 
C_f(u_2 - u_1) + C_su_2 = 0
C_f(v_1 - u_2) - g(v_1)v_1 = 0

for unknowns u_1, u_2, and v_1
and boundary block u_0

u_0 = C_s = C_f = 1, g(U) = U
are chosen so the solution is unique

Three implementations are presented, the first two are 
described in the aformentioned document.

"""

import sys
import src.blocks as b
import src.flux as f
import src.problem as p
import src.source as s
import numpy as np

"""
Fluxes are here, the u flux returns a value
corresponding to the 'u' variable, and
the parameter dictionary with 'C_f' in it
is used to demonstrate it

The v flux is similar
"""
def F_u(B,N,P):
  return {'u':P['C_f']*(B['u']-N['u'])}

def F_v(B,N,P):
  return {'v':(B['v']-N['u'])}

"""
Sources are here, again, dictionaries are
returned. The g(U) function is defined as
well, returning simply the 'v' value in
the block, and is called from the source
using the parameter dictionary, but in this
case, the entry takes in an argument, the block
"""
def S_u(B,P):
  return {'u':P['C_s']*B['u']}

def g(B):
  return B['v']

def S_v(B,P):
  return {'v':B['v']*P['g'](B)}


"""
This is the first implementation, using a block
to form the boundary
"""
def implementationOne():
  # constants
  u0 = 1.0
  C_f = 1.0
  C_s = 1.0

  # create the blocks, with no parameters or parameter functions
  b0 = b.Block('0', # name
              {'u':u0}) # boundary state, this is fixed
  b1 = b.Block('1', # name
              {'u':0,'v':0}) # states, and their initial guesses
  b2 = b.Block('2', # name
              {'u':0}) # states, and their initial guesses

  # create the fluxes
  # flux from 0 to 1
  F_P = {'C_f':C_f} # dictionary of parameters for the flux
  # flux is from b0, function is F_u
  f_u01 = f.Flux(b0,F_u,F_P,'f_u01')
  # flux is from b1, function is F_u
  f_u12 = f.Flux(b1,F_u,F_P,'f_u12')
  # flux is from b2, function is F_v
  # parameter argument is optional, if not used
  # we don't technically need a name, they really exist
  # for debugging, and if we want to do something that
  # requires us to remove the flux, here it is defaulted to ''
  f_v21 = f.Flux(b2,F_v)
  # add these fluxes to the right blocks
  b1.addFlux(f_u01)
  b1.addFlux(f_v21)
  b2.addFlux(f_u12)

  # create the sources
  S_uP = {'C_s':C_s}
  S_vP = {'g':g}
  s_u1 = s.Source(S_u,S_uP,'s_u1')
  s_u2 = s.Source(S_u,S_uP,'s_u2')
  # again name argument is optional
  s_v1 = s.Source(S_v,S_vP)
  b1.addSource(s_u1)
  b2.addSource(s_u2)

  # create and solve the problem
  problem = p.Problem([b1,b2])
  problem.solve()

  # the blocks now contain their own values,
  # which are the correct values
  print "Running first implementation of example code"
  print b1
  print b2

"""
This is the second implementation, using a source
to form the boundary
"""

""" the boundary source, with u_0 as a parameter """
def S_0(B,P):
  return {'u':P['C_f']*(B['u']-P['u_0'])}

def implementationTwo():
  # constants
  u0 = 1.0
  C_f = 1.0
  C_s = 1.0

  # create the blocks, with no parameters or parameter functions
  b1 = b.Block('1', # name
              {'u':0,'v':0}) # states, and their initial guesses
  b2 = b.Block('2', # name
              {'u':0}) # states, and their initial guesses

  # create the fluxes
  # flux from 0 to 1
  F_P = {'C_f':C_f} # dictionary of parameters for the flux
  # flux is from b1, function is F_u
  f_u12 = f.Flux(b1,F_u,F_P,'f_u12')
  # flux is from b2, function is F_v
  # name argument is optional
  # parameter argument is optional too, if not used
  f_v21 = f.Flux(b2,F_v)
  # add these fluxes to the right blocks
  b1.addFlux(f_v21)
  b2.addFlux(f_u12)

  # create the sources
  # these are parameters for the new boundary source
  S_0P = {'C_f':C_f, 'u_0':u0}
  S_uP = {'C_s':C_s}
  S_vP = {'g':g}
  s_u0 = s.Source(S_0,S_0P,'boundary')
  s_u1 = s.Source(S_u,S_uP,'s_u1')
  s_u2 = s.Source(S_u,S_uP,'s_u2')
  # again name argument is optional
  s_v1 = s.Source(S_v,S_vP)
  # add the boundary source in with the other ones
  b1.addSource(s_u0)
  b1.addSource(s_u1)
  b2.addSource(s_u2)

  # create and solve the problem
  problem = p.Problem([b1,b2])
  problem.solve()

  # the blocks now contain their own values,
  # which are the correct values
  print "Running second implementation of example code"
  print b1
  print b2

"""
This is the third implementation, using a block
to form the boundary but using Block parameters and 
parameter functions
"""

"""
The u flux returns a value
corresponding to the 'u' variable, and
assumes the block has a parameter dictionary 
containing C_f
"""
def F_u3(B,N,P):
  return {'u':B.p['C_f']*(B['u']-N['u'])}

"""
The u source now has C_s in B.p, again,
it must exist

the v source has a parameter function, as
if the block had a material with some function
It must also exist
"""
def S_u3(B,P):
  return {'u':B.p['C_s']*B['u']}

def S_v3(B,P):
  return {'v':B['v']*B.P['g'](B)}

def implementationThree():
  # constants
  u0 = 1.0
  C_f = 1.0
  C_s = 1.0

  # create the blocks, with no parameters or parameter functions
  b0 = b.Block('0', # name
              {'u':u0}) # boundary state, this is fixed
  b1 = b.Block('1', # name
              {'u':0,'v':0}, # states, and their initial guesses
              {'g':g}, # parameter functions
              {'C_f':C_f,'C_s':C_s}) # parameter constants
  b2 = b.Block('2', # name
              {'u':0}, # states, and their initial guesses
              None, # no parameter functions are needed in this block
              {'C_f':C_f,'C_s':C_s}) # parameter constants

  # create the fluxes
  # flux from 0 to 1
  # flux is from b0, function is F_u
  # no parameter dictionaries needed
  f_u01 = f.Flux(b0,F_u3)
  # flux is from b1, function is F_u
  f_u12 = f.Flux(b1,F_u3)
  # flux is from b2, function is F_v
  # name argument is optional
  # parameter argument is optional too, if not used
  f_v21 = f.Flux(b2,F_v)
  # add these fluxes to the right blocks
  b1.addFlux(f_u01)
  b1.addFlux(f_v21)
  b2.addFlux(f_u12)

  # create the sources
  # no parameters or names needed,
  # they come from the plots
  s_u1 = s.Source(S_u3)
  s_u2 = s.Source(S_u3)
  s_v1 = s.Source(S_v3)
  b1.addSource(s_u1)
  b2.addSource(s_u2)

  # create and solve the problem
  problem = p.Problem([b1,b2])
  problem.solve()

  # the blocks now contain their own values,
  # which are the correct values
  print "Running third implementation of example code"
  print b1
  print b2

if __name__ == "__main__":
  implementationOne()
  implementationTwo()
  implementationThree()
