import numpy as np
from Model import *

class Integrator():
  def __init__(self,schema,  CFL=1.0):
    if(schema=="LaxFriedrich"):
      print("Vous avez choisi d'intégrer via le schéma de Lax Friedrich avec CFL="+str(CFL))
      self.step = lambda myModel : Lax_wendroff(myModel, CFL)
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

def Lax_wendroff(myModel, CFL):
  U_ancien = myModel.U
  f_U = myModel.f(U_ancien)
  f_D = np.zeros(f_U.shape)
  f_G = np.zeros(f_U.shape)
  f_D[:-1]=f_U[1:]
  f_D[-1]=f_U[0]
  f_G[1:]=f_U[:-1]
  f_G[0]=f_U[-1]

  U_mid = np.zeros(myModel.U.shape)
  U_mid[1:]=0.5*(U_ancien[:-1]+U_ancien[1:])
  U_mid[0]=0.5*(U_ancien[-1]+U_ancien[0])

  A_G = np.zeros((myModel.n_mailles,myModel.n,myModel.n))
  A_G = myModel.J(U_mid)
  A_D = np.zeros((myModel.n_mailles,myModel.n,myModel.n))
  A_D[:-1] = A_G[1:]
  A_D[-1] = A_G[0]

  speed = myModel.maxSoundSpeed()
  dx = myModel.dx
  dt = CFL*dx/speed

  # On update le modèle
  myModel.T+=dt
  for i in range(myModel.n_mailles):
    myModel.U[i]+=0.5*dt/dx*(f_G[i]-f_D[i])+0.5*dt**2/dx**2*(A_D[i].dot(f_D[i]-f_U[i])-A_G[i].dot(f_U[i]-f_G[i]))