import numpy as np
import math

class PDESystem():
    def __init__(self, modele, **kwargs):
        if(modele=="Advection"):
            self.a = kwargs['a']
            self.f = lambda U : self.a*U
            #self.J = lambda U : np.array([[self.a]])
            self.J = lambda U : self.a*np.ones((max(U.shape[0],1),1,1))
            self.n_dim = 1
            self.modele = "Advection"
            self.iKnowExactSol = False
    
    def initialize_RP(self, myModel, UL, UR, x0):
        self.UL = UL
        self.UR = UR
        self.x0 = x0
        self.IC="Riemann"
        self.L = myModel.longueur

        if(self.modele=="Advection"):
            self.iKnowExactSol = True

        if self.n_dim==1:
            print(UL)
            for i in range(myModel.n_mailles):
                if(myModel.x_ctr[i]<x0):
                    myModel.U[i]=UL
                else:
                    myModel.U[i]=UR
        else:
            for i in range(myModel.n_mailles):
                if(myModel.x_ctr[i]<x0):
                    myModel.U[i,:]=UL
                else:
                    myModel.U[i,:]=UR

    def exact_sol(self, myModel, t,shape):
        if self.iKnowExactSol:
            if (self.modele=="Advection" and self.IC=="Riemann"):
                L = self.L
                n_parc = math.floor(self.a*t/L)
                _x0 = self.x0+self.a*t - n_parc*L
                print(_x0)
                res = np.zeros(shape)
                if self.n_dim==1:
                    for i in range(myModel.n_mailles):
                        if((_x0>self.x0 and myModel.x_ctr[i]<_x0 and myModel.x_ctr[i]>_x0-self.x0) or (_x0<self.x0 and (myModel.x_ctr[i]<_x0 or myModel.x_ctr[i]>L-self.x0+_x0))):
                            res[i]=self.UL
                        else:
                            res[i]=self.UR
                else:
                    for i in range(myModel.n_mailles):
                        if((_x0>self.x0 and myModel.x_ctr[i]<_x0 and myModel.x_ctr[i]>_x0-self.x0) or (_x0<self.x0 and (myModel.x_ctr[i]<_x0 or myModel.x_ctr[i]>L-self.x0+_x0))):
                            res[i,:]=self.UL
                        else:
                            res[i,:]=self.UR
                return res



    