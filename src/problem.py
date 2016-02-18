"""
problem.py contains the Problem Class

Each Problem has:
  (.b) a list of the blocks to solve for in the problem
  (.bc) a list of the blocks used as boundary blocks
  (.mapping) the mapping between the local and global systems
  (._band) a pair (a,b) corresponding to the known bandedness of 
    the problem, if the matrix structure of the global problem is known
    This depends a lot on the order of the list of blocks
    if specified, this can speed up fsolve by orders of magnitude

This class solves R(U) = F(U,U_N) + S(U) = 0
by assembling the global system and solving it with fsolve

"""

"""
fsolve is used for solving steady state
ode is used for solving transient
"""
from scipy.optimize import fsolve
from scipy.integrate import odeint
from collections import OrderedDict
from sys import exit
from numpy import array
try:
  from numdifftools import nd
  JACOBIAN = True
except ImportError:
  JACOBIAN = False

class Problem(object):
  """ 
  Problem Class

  __init__:   Problem Constructor

  input(s):   (blocks) relevant blocks
              (boundaries) blocks on the boundaries
                these are needed for unsteady problems
  output(s):  None
  """
  def __init__(self,blocks,boundaries = [],**parameters):
    self.b = blocks
    # uniqueness in block names is required
    namelist = set([b.name for b in blocks])
    if(len(namelist) != len(blocks)):
      exit("multiple blocks have the same name")

    self.bc = boundaries
    self.mapping = [(i, k) for i, b in enumerate(blocks) \
      for k in b.state.keys()]

    # Bandedness (default to None)
    self._band = None
  """
  __repr__    overloading for print command

  input(s):   None
  output(s):  other information
  """
  def __repr__(self):
    return "problem with "+str(len(b))+" blocks"

  """
  __getitem__, __setitem__: get/set blocks by integer key in list

  input(s):   access blocks by index
  output(s):  as shown
  """
  def __getitem__(self, key):
    return self.b[key]
  def __setitem__(self, key, value):
    self.b[key] = value

  """
  setBand:    sets the bandedness of fsolve (see fsolve help)

  input(s):   (int,int), pair of number of sub and super diagonals
  output(s):  none
  """
  def setBand(self,band):
    self._band = band

  """
  update:     Updates the blocks by unwrapping the new solution

  input(s):   (solution) global array of floats corresponding to mapping
  output(s):  None
  """     
  def update(self,solution):
    for ix, (i,k) in enumerate(self.mapping):
      self.b[i][k] = solution[ix]
    for bc in self.bc:
      for s in bc.state:
        bc[s] = sum([S.source()[s] for S in bc.S])

  def updateUnst(self,t):
    for b in self.b + self.bc:  
      b.t = t

  def getSolutionVec(self):
    solution = [None]*len(self.mapping)
    for ix, (i,k) in enumerate(self.mapping):
      solution[ix] = self.b[i][k]
    return solution

  """
  r:          Global residual function r(solution)
  rVec:       Same as r, but in numpy array format (used elsewhere)

  input(s):    (solution) global array of floats corresponding to mapping
  output(s):  R(solution) global array of floats corresponding to residual

  updates solution first, then computes
  should be passed into another function
  """
  def r(self,solution):
    self.update(solution)
    return [self.b[i].R()[v] for i,v in self.mapping]

  def rVec(self,solution):
    return array(self.r(solution))

  """
  rUnst:      Unsteady version of above

  input(s):   (solution) global array of floats corresponding to mapping
              (t) time to evaluate solution at
  output(s):  R(solution) global list of floats corresponding to residual

  updates solution first, then computes
  should be passed into another function
  """
  def rUnst(self,solution,t):
    self.updateUnst(t)
    self.update(solution)
    return [self.b[i].R()[v]/self.b[i].T(self.b[i])[v] for i,v in self.mapping]

  """
  solve:      wrapper for chosen (non)linear solver

  input(s):   None
  output(s):  None

  unwraps blocks, passes into solver, finishes by updating blocks one last time
  """
  def solve(self,t=0):
    solution = fsolve(self.r,self.getSolutionVec(),band=self._band)

  """
  jacobian:   computes numerical jacobian at a given solution

  input(s):   None
  output(s):  Jacobian Matrix

  this is useful for looking at matrix structure, among other things
  """
  def jacobian(self):
    if JACOBIAN is False:
      exit("cannot compute numerical jacobian without numdifftools package")
    Jfun = nd.Jacobian(self.rVec)
    return Jfun(self.getSolutionVec())

  """
  solveUnst:  solve the transient problem

  input(s):   (ti) initial time
              (tf) final time
              (n)  number of timesteps
  output(s):  Solution at each timestep

  unwraps blocks, passes into solver, finishes by updating blocks one last time
  """
  def solveUnst(self,t):
    solution = [None]*len(self.mapping)
    # This has the unsteady part
    # Solver, just live and let live  
    for ix, (i,k) in enumerate(self.mapping):
      solution[ix] = self.b[i][k]
    soln = odeint(self.rUnst, solution, t,hmax=(t[-1]-t[0])/len(t), \
      rtol = 1e-4, atol = 1e-4)

    # final update
    self.updateUnst(t[-1])
    self.update(soln[-1,:])

    # Lets collect all the steps
    fullSolution = dict([(b.name + '_'+s,[]) for b in self.b for s in b.state])
    for j in range(0,len(t)):
      for ix, (i,k) in enumerate(self.mapping):
        fullSolution[self.b[i].name+'_'+k].append(soln[j,ix])
    fullSolution['t'] = t
    return fullSolution
