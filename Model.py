import numpy as np
import matplotlib.pyplot as plt

class Model():
  def __init__(self, longueur, n_mailles, dimension, myPDESystem):
    self.longueur = longueur
    self.n_mailles = n_mailles # plutôt pair svp
    self.n=dimension
    self.f=myPDESystem.f
    self.J=myPDESystem.J
    self.x = self.longueur*np.array(range(self.n_mailles+1))/self.n_mailles
    self.dx = self.longueur/self.n_mailles
    self.x_ctr = self.x[:-1]+self.dx/2
    self.U=np.zeros((self.n_mailles,self.n))
    self.T=0

  def maxSoundSpeed(self):
    maxi = -1
    for i in range(self.n_mailles):
      maxi = max(maxi, np.max(np.abs(np.linalg.eigvals(self.J(self.U[i])[0]))))
    self.max_speed = maxi
    return maxi
  
  def plot(self):
    if self.n>1:
        y = np.sin(self.x_ctr ** 2)    
        fig, axs = plt.subplots(self.n)
        fig.suptitle('Etat du système au temps '+str(self.T))
        for i in range(self.n):
            axs[i].plot(self.x_ctr, self.U[:,i])
        plt.show()
            
    else:
      y = np.sin(self.x_ctr ** 2)    
      fig, axs = plt.subplots()
      fig.suptitle('Etat du système au temps '+str(self.T))
      axs.plot(self.x_ctr, self.U[:,0])
      plt.show()
