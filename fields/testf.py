#
#
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
#
#import the modulus
import fields as tf

order = 1

x  = np.linspace(-10,10,num=100)
dx = x[1] - x[0]
x2 = x*x

res = tf.deriv1D(x2, order, dx)
res2 = tf.deriv1D(x2, 2, dx)

res2_2 = tf.deriv1D(2*x,order,dx)

plt.plot(x,x2,label = 'x2')
plt.plot(x,2*x,label = '2x')
plt.plot(x,res,label='deriv')
plt.plot(x,res2,label='dd')
plt.plot(x,res2_2,label='dd_bis')
plt.grid()
plt.legend(loc='best')
plt.show()