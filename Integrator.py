class Integrator():
  def __init__(self,schema,  CFL=1.0):
    if(schema=="LaxFriedrich"):
      print("Vous avez choisi d'intégrer via le schéma de Lax Friedrich avec CFL="+str(CFL))
      self.step = lambda myModel : LaxFriedrich1D(myModel, CFL)
      self.nb_pas_max = 500
      self.n = 0
      
  def integrate(self, myModel, Tf):
    while(myModel.T<Tf and self.n<self.nb_pas_max):
      self.n+=1
      self.step(myModel)
    if(self.n==self.nb_pas_max):
      print("Attention, nombre maximal de pas atteint")

# Ne fonctionne qu'avec un pas d'espace fixe
# Conditions de bord périodiques

def LaxFriedrich1D(myModel, CFL):
  U_ancien = myModel.U
  f_U = myModel.f(myModel.U)
  speed = myModel.maxSoundSpeed()
  U_neuf = np.zeros(myModel.U.shape)

  # Pas de temps et d'espace
  delta_x = myModel.dx
  delta_T = CFL*delta_x/speed

  # Flux numériques en i-1/2
  flux_numeriques_gauche = np.zeros(myModel.U.shape)
  flux_numeriques_gauche[1:,:] = 0.5*(f_U[:-1,:]+f_U[1:,:])-0.5*delta_x/delta_T*(U_ancien[1:,:]-U_ancien[:-1,:])
  flux_numeriques_gauche[0,:] = 0.5*(f_U[-1,:]+f_U[0,:])-0.5*delta_x/delta_T*(U_ancien[0,:]-U_ancien[-1,:])

  # Calcul de la solution à l'instant T+delta_T
  U_neuf[0:-1,:]=U_ancien[0:-1,:]+0.5*delta_T/delta_x*(flux_numeriques_gauche[:-1,:]-flux_numeriques_gauche[1:,:])
  U_neuf[-1,:]=U_ancien[-1,:]+0.5*delta_T/delta_x*(flux_numeriques_gauche[-1,:]-flux_numeriques_gauche[0,:])

  # On update le modèle
  myModel.T+=delta_T
  myModel.U=U_neuf 