"""
The functions in this file are useful for computing the fourier expansion of an arbitrary function, and for animating the convergence of the partial sums.

There are many possible performance optimizations to make to these implementations.  Please submit a PR if you improve the code.
"""
import numpy as np
from scipy import integrate

def fourier_coeff(fn, n, lower_bound=-np.pi, upper_bound=np.pi, return_complex=True):
	"""
	computes the fourier coefficients of an arbitrary function

	Parameters
	----------
	fn : callable
	    a real valued function of a real variable

	n: int
		the number of coefficients to return

	lower_bound: float
		lower bound of integration for computing the coefficients

	upper_bound: float
		upper bound of integration for computing the coefficients

	return_complex(optional): bool
		whether or not to return the complex coefficients (default is True)

	Returns
	-------
	c_n: list
		when return_complex is True, returns a list of complex numbers
	
	a0, a_n, b_n: tuple
		when return_complex is False it returns a tuple which contains three items (a0, a_n, a_b) where a_n and a_b are lists of length n


	Examples
	--------

	>>> fourier_coeff(lambda x: x**2, 3)
	[6.5797362673929065, (-3.9999999999999996+0j), (1.0000000000000004+0j)]
	"""
	a_n, b_n = [], []

	a0 = integrate.quad(lambda u: 1/np.pi * fn(u), lower_bound, upper_bound)[0]

	for i in range(1,n):
		a = integrate.quad(lambda u: 1/np.pi * fn(u)*np.cos(i*u), lower_bound, upper_bound)[0]
		b = integrate.quad(lambda u: 1/np.pi * fn(u)*np.sin(i*u), lower_bound, upper_bound)[0]
		a_n.append(a)
		b_n.append(b)


	if return_complex:
		c_n = [a0]
		for i in range(0,n-1):
			c_n.append(a_n[i] + 1j*b_n[i])
		return c_n
	return a0, a_n, b_n

def partial_fourier_series(x, fn, N):
	"""
	computes the partial fourier series an arbitrary function

	f(x) = 1/2 a0 + a1 cos(x) + b1 sin(x) + a2 cos(2x) + b2 sin(2x) + ...

	Parameters
	----------
	x : array
	    the points to evaluate the function at

	fn: callable
		an arbitrary real function of a real variable to approximate

	N: int
		the number of terms of the partial sum to compute

	Returns
	-------
	y: array
		the values of the partial sum

	Examples
	--------

	>>> X = np.linspace(-np.pi, np.pi, num=100)
	>>> partial_fourier_series(X, wavy_line, 3)
	[8.28986813 8.27376965 8.22563609 8.14595086 8.03551232 7.89542286
	 ...
	 8.14595086 8.22563609 8.27376965 8.28986813]

	Animation
	>>> app.animate("partial_fourier_series", start=3, stop=25, delay=0.1, step=1, argname="N")
	"""	
	xmin = x[0]
	xmax = x[len(x) - 1]
	a0, a, b = fourier_coeff(fn, N, xmin, xmax, False)
	result = np.full(len(x), 0.5*a0)
	# result = np.zeros(len(x))
	for n in range(0,N-1):
		result += a[n]*np.cos((n+1)*x) + b[n]*np.sin((n+1)*x)
	return result