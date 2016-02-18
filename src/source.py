"""
source.py contains the Source class

Each Source has:
  (.S) source function, which returns a dictionary
  (.P) The parameters associated with it
  (.name) a string associated with the source

Each source function is f(B,P), a function of a block
and its parameters. For example, consider a block with 
two variables, u and v with the source for v, 
S(v) = 5*(u-v), one possible implementation is

def sourceFunction(B,P):
  s = (B['u']-B['v'])*P['Constant']
  return {'u':s}

u = Block('u1',{'u':'1'})

u_source = Source(sourceFunction,{'Constant',5},'u_source')

u.addSource(u_source)

"""

""" commonly used sources """
def constant(B,P):
  return dict((key,P[key]) for key in P.keys())

def linear(B,P):
  return dict((key,B[key]*P[key]) for key in P.keys())

class Source(object):
  """ 
  Source Object. Similar to a flux, but only needing one state

  __init__:   Source Constructor

  input(s):   (s) string corresponding to function name
              (parameters) optional dictionary with arguments for the 
              source functions
  output(s):  None
  """
  def __init__(self,s,P=None,name=''):
    self.B = None
    self.S = s
    self.P = P
    self.name = name

  """
  The following is a wrapper for flux function choices defined with
  Sources should be define as such
  input(s):   (B) Block
              (P) Parameters (optional) 
  output(s):  dict with entries for the flux for each state contributed to

  """
  def source(self):
    return self.S(self.B,self.P)

  """
  __repr__    overloading for print command

  input(s):   None
  output(s):  name
  """
  def __repr__(self):
    return name