from scipy.optimize import fsolve
import math

def equations(p,*kargs):
    al, be = p
    
    a,b,x,y = kargs
    
    return (-be * x/(b*b) + al*y/(a*a)  - al* be * (1/(a*a) - 1/(b*b)) , 
            al *al/(a*a) + be*be/(b*b) - 1 )


def equation1(p):
    x, y = p
    return x+y**2-4

def equation2(p):
    x, y = p
    return math.exp(x) + x*y - 3

def equationsBis(p,*kargs):
    t, de = p
    
    a,b,x,y = kargs
    
    return (x-(b+de)*math.cos(t) , 
            y-(a+de)*math.sin(t) )


a  = 6.0
b  = 6.0
x  = 8.0
y  = 8.0

args = (a,b,x,y)
xx,yy =  fsolve(equationsBis, (2., -5.),args=args)
print equationsBis((xx, yy),*args)