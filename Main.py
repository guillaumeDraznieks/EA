from PDESystem import *
from Model import *
from Integrator import *

myPDESystem = PDESystem("Advection", a=10)
myModel = Model(1,100,1,myPDESystem)
myPDESystem.initialize_RP(myModel,1,0,0.2)
myIntegrator = Integrator("LaxFriedrich", CFL=0.5)

myIntegrator.integrate(myModel,0.03)

myModel.plot()