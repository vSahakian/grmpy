from numpy import load
import inversion as inv
from scipy.linalg import lstsq

G1500=load('/Users/vsahakian/Desktop/tmpG1500.npy')
G81=load('/Users/vsahakian/Desktop/tmpG81.npy')

d1500=load('/Users/vsahakian/Desktop/tmpd1500.npy')
d81=load('/Users/vsahakian/Desktop/tmpd8.1.npy')

m1500,residual1500,rank,singular_vals=lstsq(G1500,d1500)
m81,residual81,rank,singular_vals=lstsq(G81,d81) 
