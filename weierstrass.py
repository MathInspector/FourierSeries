"""
The Weierstrass approximation is one of the most important theorems in analysis.  It asserts more than completeness, it asserts the possibility of uniform convergence.

The functions in this file can be used for creating a visualization of both the 1-D Weierstrass and 2-D approximation formulas.  
"""
import numpy as np
from scipy import integrate

def J_n(u, n):
	"""
	Used for computing the denominator of P_n
	"""
	return (1 - u**2)**n

def P_n(x,f,n):
	"""
	computes the partial sums of the weierstrass approximation of an arbitrary function on a single variable f(x)

	Since this approximation only works for functions on a domain [0,1], there are some helper functions in for shrinking functions before plugging them into J_n.

	P_n = 1/(2J_n) integral_0^1 {(1 - (u-x))^n du}

	Parameters
	----------
	x : array
	    the points to evaluate the partial sums at

	f: callable
		the function to approximate

	n: int
		the number of terms of the partial sum to compute

	Returns
	-------
	y: array
		the value of the partial sum
	
	Examples
	--------

	>>>  u = np.linspace(0,1)
	>>>  P_n(u, shrink(u, deg4,  6, -2), 11)
	[ 1.08333333  0.56826812  0.0860291  -0.35055781 -0.73094778 -1.04720151
	 ...
	  0.56826812  1.08333333]

	Animation
	>>> app.animate("P_n", start=1, stop=100, delay=0.1, step=1, argname="n")
	"""	
	scale = x[len(x)-1] - x[0]
	translate = x[0]
	u = np.linspace(0,1,num=len(x))
	
	def J_numerator(u, f, n, x):
		return f(x) * (1 - (u - x)**2)**n 

	J2 = 1 / integrate.quad(J_n, -1, 1, args=n)[0]
	values = []
	for i in x:
		val = J2 * integrate.quad(J_numerator, 0, 1, args=(f,n,i))[0]
		values.append(val)
	return np.array(values)

def P2_n(f,n,num=50):
	"""
	Computes a 3d surface whose height is given by the 2d Weierstrass approximation

	This function is currently not working properly, and has a lot of visual artifacts.

	Please fix!
	"""
	J2 = 1 / integrate.quad(J_n, -1, 1, args=n)[0]**2
	r = np.linspace(0, 1, num=num)
	theta = np.linspace(0, 2*np.pi - 0.001, num=num)
	
	def phi_n(r,theta,f,n):
		J2 = 1 / integrate.quad(J_n, -1, 1, args=n)[0]
		
		def J2_numerator(u, f, n, r, theta):
			return f(theta) * ((1 - (u - r * np.cos(theta))**2) * (1 - (u - r * np.sin(theta))**2))**n 
		
		return r*np.cos(theta), integrate.quad(J2_numerator, 0, 1, args=(f,n,r,theta))[0], r*np.sin(theta)
	
	values = []
	for i in r:
		temp = []
		for j in theta:
			temp.append(phi_n(i,j,f,n))
		values.append(temp)

	values2 = []
	for j in theta:
		temp = []
		for i in r:
			temp.append(phi_n(i,j,f,n))
		values2.append(temp)
	return values, values2

#####################
### TEST FUNCTION ###
#####################
def deg4(x, scale=1, translation=0):
	"""
	4th degree polynomial with specific coefficients which looks like the letter W
	looks best when plotted on the range -2 < x < 4
	"""
	poly_coeff = [-3,12,-2,-4,1]
	return scale * np.polynomial.polynomial.polyval(x + translation, poly_coeff)

def wavy_line(x):
	"""
	Just a random function for testing animations
	"""
	return np.sin(x)**3 + (0.3*x)**2 + 0.5*x

#########################
### UTILITY FUNCTIONS ###
#########################
def shrink(u, f, scale, translate):
	"""
	scale is how much to multiply the closed interval [0,1] by
	translate is how much to shift the new interval left

	for example, if scale is 6 and translate is -2
	first [0,1] --> [0,6] by multiplying by 6
	then  [0,6] --> [-2,4] by subracting by 2
	"""
	return lambda u: f(u * scale + translate) / scale

def shrunk(u,f,scale,translate):
	"""
	This is useful showing a target for the weierstrass approximation of a shrunk down function
	It's the same as shrink, except instad of returning a lambda function, it returns the values
	"""
	return shrink(u,f,scale,translate)(u)