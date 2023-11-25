from mpmath import *
from numpy import array
import time

def F(arg : array):
	return 9*arg[0][0]**2 + arg[1][0]**2

def g1(arg : array):
	return (arg[0, 0] - 3)**2 + (arg[1, 0] - 2)**2 - 1

def g2(arg : array):
	return arg[0, 0] + arg[1, 0] - 5

def H_extern(arg : array):
	return max(0, g1(arg))**2 + max(0, g2(arg))**2

def vert_descent(start_point : array, F, step : mpf = 0.5, milling_coef : mpf = 0.5, epsilon : mpf = 1e-16, dim_num : int = 2):
	start_time = time.time()
	count_of_iterat = 0

	old_min_value_F = F(start_point)
	old_cycle_start_point = start_point.copy()

	while True:
		changing_vector = False
		for iterat in range(dim_num):
			new_point = start_point.copy()
			count_of_iterat += 1

			new_point[iterat, 0] += step
			new_value_F = F(new_point)
			if new_value_F < old_min_value_F: 
				old_min_value_F = new_value_F
				start_point = new_point.copy()
				changing_vector = True
				continue

			new_point[iterat, 0] -= 2*step
			new_value_F = F(new_point)
			if new_value_F < old_min_value_F: 
				old_min_value_F = new_value_F
				start_point = new_point.copy()
				changing_vector = True
				continue

		if not changing_vector: step *= milling_coef
		if norm([start_point[i, 0] - old_cycle_start_point[i, 0] for i in range(dim_num)]) < epsilon and step < epsilon * milling_coef: break
		old_cycle_start_point = start_point.copy()

	return start_point, count_of_iterat, round(time.time() - start_time, 3)

def external_penalty_method(start_point : array, epsilon : mpf = 1e-16):

	start_time = time.time()
	count_of_iterat = 0

	while H_extern(start_point) > epsilon:
		phi = lambda x: F(x) + 10 ** count_of_iterat * H_extern(x)
		start_point = vert_descent(start_point, phi, epsilon=epsilon)[0]
		count_of_iterat += 1

	return start_point, count_of_iterat, round(time.time() - start_time, 3)

'''
def internal_penalty_method(start_point : array, epsilon : mpf = 1e-16):

	start_time = time.time()
	count_of_iterat = 0

	while H_extern(start_point) > epsilon:
		phi = lambda x: F(x) + 10 ** count_of_iterat * H_extern(x)
		start_point = vert_descent(start_point, phi, epsilon=epsilon)[0]
		count_of_iterat += 1

	return start_point, count_of_iterat, round(time.time() - start_time, 3)
'''

start_point = array([[1.9], [1.1]])
print(external_penalty_method(start_point, epsilon=1e-8))