import pickle
from numpy import zeros,ones,array
from numpy.linalg import lstsq,norm
from scipy.optimize import minimize


f = open(u'/Users/vsahakian/anza/models/pckl/regr_0.0_3.3_4.5_6.5_resid_2663.62423294.pckl')
p=pickle.load(f)
f.close()

# run traditional least squares
mL2,res,rank,s=lstsq(p.G,p.d)


# Objective function (the thing to be minimzed, i.e. ||GM-d||
def objective(m,G=p.G,d=p.d):
    from numpy.linalg import norm
    residuals=G.dot(m)-d
    L2norm=norm(residuals)
    return L2norm

#Check we get the same residuals with objective function than what came out of lstsq
obj_res=objective(mL2,p.G,p.d)
print 'lstsq residual is %.2f, objective residuals is %.2f' % (res**0.5,obj_res)




#Perform unconstrained inversion with scipy.minimize

#Inital guess of the solution
m0=ones(mL2.shape)
#Run the fucker
optimization_result=minimize(objective,m0, method='SLSQP')
#What are the coefficients?
mmin=optimization_result.x






#Now add constraints to inversion, still use scipy.minimize
param_constraints= ({'type': 'ineq','fun' : lambda m: array(-m[2])},
                    {'type': 'ineq','fun' : lambda m: array(-m[7])},
                    {'type': 'ineq','fun' : lambda m: array(-m[12])})
                    
#Bounds on the parameters, basically all are free except a3,a8 and a13
bnds=((-100,100),(-100,100),(-100,100),(0.9,1.1),(-100,100),
      (-100,100),(-100,100),(-100,100),(0.9,1.1),(-100,100),
      (-100,100),(-100,100),(-100,100),(0.9,1.1),(-100,100))


#Run the sonbitch
optimization_result=minimize(objective,m0, method='SLSQP',constraints=param_constraints,bounds=bnds)
mmin_const=optimization_result.x





#Spit out each different result
print '\n'
print 'lstsq result is:'
print mL2

print '\n'
print 'minimize uncosntrained is'
print mmin

print '\n'
print 'minimize with constraints is'
print mmin_const