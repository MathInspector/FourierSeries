"""
The discontinuous h-function lies at the heart of the proof that Fourier series converge for a wide class of functions.  The h-function is discontinuous at zero, but it has one of the simplest possible fourier expansions possible.  This makes it an excellent introduction to computational fourier series.
"""

import numpy as np

def h(x):
	"""
	computes the discontinuous h-function

	`h(x) = 0.5 * (pi - x)`
	`h(0) = 0`
	`h(x + 2pi) = h(x)`

	Parameters
	----------
	x : array
	    the points to evaluate the h-function at

	Returns
	-------
	y: array
		the value of h(x)
	
	Examples
	--------

	>>> X = np.linspace(-np.pi, np.pi, num=100)
	>>> h(X)
	[-0.         -0.03173326 -0.06346652 -0.09519978 -0.12693304 -0.1586663
	 ...
	  0.09519978  0.06346652  0.03173326  0.        ]
	"""
	result = []
	for i in x:
		if i > 0:
			result.append(0.5*(np.pi - i))
		elif i ==0:
			result.append(0)
		else:
			result.append(-0.5*(np.pi + i))
	return np.array(result)


def h_partial_sum(x, N):
	"""
	computes the partial sums of the fourier series for the discontinuous h-function

	h(x,N) = sin(x) + 1/2 sin(2x) + ... + 1/N sin(Nx)

	Parameters
	----------
	x : array
	    the points to evaluate the partial sums at

	N: int
		the number of terms of the partial sum to compute

	Returns
	-------
	y: array
		the value of the partial sum
	
	Examples
	--------

	>>> X = np.linspace(-np.pi, np.pi, num=100)
	>>> h_partial_sum(X, 11)
	[ 0.         -0.00229248 -0.01716756 -0.05187561 -0.10524304 -0.16812812
	 ...
	  0.05187561  0.01716756  0.00229248  0.        ]

	Animation
	>>> app.animate("h_partial_sum", start=1, stop=100, delay=0.1, step=1, argname="N")
	"""
	result = np.zeros(len(x))
	for n in range(1,N):
		result += 1/n * np.sin((n*x))
	return np.array(result)

def g(x):
	"""
	computes the modified version of the discontinuous h-function

	`g(x) = h(x)(1 - cos(x))`

	Parameters
	----------
	x : array
	    the points to evaluate the g-function at

	Returns
	-------
	y: array
		the value of g(x)
	
	Examples
	--------

	>>> X = np.linspace(-np.pi, np.pi, num=100)
	>>> g(X)
	[-0.00000000e+00 -6.34026289e-02 -1.26422436e-01 -1.88679171e-01
	 ...
	  1.88679171e-01  1.26422436e-01  6.34026289e-02  0.00000000e+00]
	  """

	result = []
	for i in x:
		# print (i)
		if i > 0:
			result.append(0.5*(np.pi - i)*(1 - np.cos(i)))
		elif i ==0:
			result.append(0)
		else:
			result.append(-0.5*(np.pi + i)*(1 - np.cos(i)))
	return np.array(result)

def g_partial_sum(x, N):
	"""
	computes the partial sums of the fourier series for the g-function

	g(x,N) = (1 - 1/2) sin(x) + (1/2 - 0.5(1/2 + 1/3)) sin(2x) + ... + (1/N - 0.5(1/(N-1) + 1/(N+1))) sin(Nx)

	Parameters
	----------
	x : array
	    the points to evaluate the partial sums at

	N: int
		the number of terms of the partial sum to compute

	Returns
	-------
	y: array
		the value of the partial sum
	
	Examples
	--------

	>>> X = np.linspace(-np.pi, np.pi, num=100)
	>>> g_partial_sum(X, 11)
	[-1.23021338e-16 -6.36700872e-02 -1.26842046e-01 -1.89069618e-01
	 ...
	  1.89069618e-01  1.26842046e-01  6.36700872e-02  1.23021338e-16]
	
	Animation
	>>> app.animate("g_partial_sum", start=1, stop=100, delay=0.2, step=1, argname="N")
	"""	
	result = np.zeros(len(x))
	for n in range(1,N):
		if n == 1:
			result += (1 - 0.5*(1/2)) * np.sin((n*x))
		else:
			result += (1/n - 0.5*(1/(n-1) + 1/(n+1))) * np.sin((n*x))
	return np.array(result)
