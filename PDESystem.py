class PDESystem():
  def __init__(self, modele, **kwargs):
    if(modele=="Advection"):
      self.a = kwargs['a']
      self.f = lambda U : self.a*U
      self.J = lambda U : np.array([[a]])