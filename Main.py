from PDESystem import *
from Model import *
from Integrator import *

a=10
myPDESystem = PDESystem("Advection", a=10)
myModel = Model(1,100,2,myPDESystem)

for i in range(99):
  myModel.U[i,0]=np.sin(2*np.pi*myModel.x[i])
  if(i<50):
    myModel.U[i,1]=1

myModel.plot()

myIntegrator = Integrator("LaxFriedrich", CFL=1.5)
myIntegrator.integrate(myModel,0.03)

myModel.plot()

print("bonjour")