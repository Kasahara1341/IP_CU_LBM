import numpy as np

def integrated4(a,index,t):
    b = max(a) + t
    a0 = a[(0+index)%4]
    a1 = a[(1+index)%4]
    a2 = a[(2+index)%4]
    a3 = a[(3+index)%4]
    c = 3*b**4 - 4*(a1+a2+a3)*b**3 \
        + 6*(a1*a2 + a2*a3 + a3*a1)*b**2 - 12*a1*a2*a3*b
    c /= 12*(a0 - a1)*(a0 - a2)*(a0 - a3)
    return c

def integrated2(a0,a1,b):
    c = (b/2-a1)*b/(a0-a1)
    return c
a = [10,7.5,2.5,0]
b = 2.5

coeff = [0,0,0,0]
# coeff[0] = integrated2(a[0],a[1],b) - integrated2(a[0],a[1],2)
# coeff[1] = integrated2(a[1],a[0],b) - integrated2(a[1],a[0],2)
for i in range(4):
    coeff[i] = integrated4(a,i,b) - integrated4(a,i,a[0])
# coeff[1] = integrated4(a[1],a[2],a[3],a[0],b) - integrated4(a[1],a[2],a[3],a[0],4)
# coeff[0] = integrated4(a[0],a[1],a[2],a[3],b) - integrated4(a[0],a[1],a[2],a[3],4)
# coeff[3] = integrated4(a[3],a[0],a[1],a[2],b) - integrated4(a[3],a[0],a[1],a[2],4)
# coeff[2] = integrated4(a[2],a[3],a[0],a[1],b) - integrated4(a[2],a[3],a[0],a[1],4)
for i in range(4):
    # coeff[i] *= 24.0
    pass
print(coeff)