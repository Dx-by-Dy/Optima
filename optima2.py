import matplotlib.pyplot as plt
from mpmath import *
import time

mp.dps = 50
def F(arg : mpf):
	global count_of_func_call
	count_of_func_call += 1
	#return -sin(arg)
	#return arg**4 + arg
	return exp(arg**2)*sin(arg)**2/cos(arg)

def F_prime(arg : mpf):
	global count_of_func_prime_call
	count_of_func_prime_call += 1
	#return -cos(arg)
	#return 4*arg**3 + 1
	return diff(lambda arg: exp(arg**2)*sin(arg)**2/cos(arg), arg)

def F_prime_2(arg : mpf):
	global count_of_func_prime_2_call
	count_of_func_prime_2_call += 1
	#return sin(arg)
	#return 12*arg**2
	return diff(lambda arg: exp(arg**2)*sin(arg)**2/cos(arg), arg, 2)

def golden_ratio(start_a : mpf, start_b : mpf, epsilon : mpf, type_of_search : str = "min"):
	count_of_iterations = 1

	left_phi = (3 - pow(5, 0.5)) / 2
	right_phi = (pow(5, 0.5) - 1) / 2
	left_point = start_a + left_phi * (start_b - start_a)
	right_point = start_a + right_phi * (start_b - start_a)

	left_value = F(left_point)
	right_value = F(right_point)

	while start_b - start_a > 2 * epsilon:
		count_of_iterations += 1

		if (type_of_search == "min" and left_value < right_value) or (type_of_search == "max" and left_value > right_value):
			start_b = right_point
			right_point = left_point
			right_value = left_value
			left_point = start_a + left_phi * (start_b - start_a)
			left_value = F(left_point)
		else:	
			start_a = left_point
			left_point = right_point
			left_value = right_value
			right_point = start_a + right_phi * (start_b - start_a)
			right_value = F(right_point)

	return start_a, start_b, count_of_iterations

def tangent_method(start_a : mpf, start_b : mpf, epsilon : mpf, type_of_search : str = "min"):
	start_time = time.time()

	start_a, start_b, count_of_iterations_golden_ratio = golden_ratio(start_a, start_b, 1e-1, type_of_search)
	count_of_iterations_tangent_method = 1

	left_value = F(start_a)
	right_value = F(start_b)
	left_value_prime = F_prime(start_a)
	right_value_prime = F_prime(start_b)

	while start_b - start_a > 2 * epsilon:
		count_of_iterations_tangent_method += 1

		new_point = (right_value - left_value + left_value_prime * start_a - right_value_prime * start_b) / (left_value_prime - right_value_prime)
		new_value = F(new_point)
		new_value_prime = F_prime(new_point)

		if (type_of_search == "min" and new_value_prime > 0) or (type_of_search == "max" and new_value_prime < 0):
			start_b = new_point
			right_value = new_value
			right_value_prime = new_value_prime
		else:	
			start_a = new_point
			left_value = new_value
			left_value_prime = new_value_prime

	return ((start_a + start_b) / 2, round(time.time() - start_time, 3), count_of_iterations_golden_ratio, count_of_iterations_tangent_method)

def Newton_Rav_method(start_a : mpf, start_b : mpf, epsilon : mpf, type_of_search : str = "min", type_of_method : str = "prime_2"):

	start_time = time.time()

	start_a, start_b, count_of_iterations_golden_ratio = golden_ratio(start_a, start_b, 1e-1, type_of_search)
	count_of_iterations_Newton_Rav_method = 1

	if type_of_method == "prime_2":
		value_prime = F_prime(start_a)
		value_prime_2 = F_prime_2(start_a)

		while abs(value_prime) > 2 * epsilon:
			count_of_iterations_Newton_Rav_method += 1

			start_a = start_a - value_prime / value_prime_2
			value_prime = F_prime(start_a)
			value_prime_2 = F_prime_2(start_a)

	elif type_of_method == "prime":
		first_value_prime = F_prime(start_a)
		second_value_prime = F_prime(start_b)

		while start_b - start_a > 2 * epsilon:
			count_of_iterations_Newton_Rav_method += 1

			new_point = start_a - (first_value_prime * (start_b - start_a)) / (second_value_prime - first_value_prime)
			start_b = start_a
			second_value_prime = first_value_prime
			start_a = new_point
			first_value_prime = F_prime(start_a)

	return (start_a, round((time.time() - start_time), 3), count_of_iterations_golden_ratio, count_of_iterations_Newton_Rav_method)


def grafics() -> None:
	global count_of_func_call, count_of_func_prime_call, count_of_func_prime_2_call

	plt.ion()

	ax, ax_text = [], []
	data_of_time, data_of_iter_gold, data_of_iter_tang, data_of_func_call, data_of_func_prime_call = [], [], [], [], []
	for i in range(2, 11):
		ax += [i]

		count_of_func_call = 0
		count_of_func_prime_call = 0
		epsilon = pow(10, -i)

		result, time_of_work, count_of_iterations_golden_ratio, count_of_iterations_tangent_method = tangent_method(-pi/2, pi/2, epsilon)
		ax_text += ["$ { 10}^{" + str(-i) + "} (" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter_gold += [count_of_iterations_golden_ratio]
		data_of_iter_tang += [count_of_iterations_tangent_method]
		data_of_func_call += [count_of_func_call]
		data_of_func_prime_call += [count_of_func_prime_call]

	axis[0].clear()
	axis[0].plot(ax, data_of_iter_gold, "g", alpha = 0.7, label = "Количество итераций метода золотого сечения")
	axis[0].plot(ax, data_of_iter_tang, "r", label = "Количество итераций метода касательных")
	axis[0].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции")
	axis[0].plot(ax, data_of_func_prime_call, "b", label = "Количество вызовов производной функции")
	axis[0].set_title("Метод касательных")
	axis[0].set_xlabel("Степень погрешности")
	axis[0].set_xticks(ax, ax_text)
	axis[0].grid()
	axis[0].legend()

	ax, ax_text = [], []
	data_of_time, data_of_iter_gold, data_of_iter_New_Rav, data_of_func_prime_call, data_of_func_prime_2_call = [], [], [], [], []
	data_of_time_modif, data_of_iter_New_Rav_modif, data_of_func_prime_call_modif = [], [], []
	for i in range(2, 11):
		ax += [i]

		count_of_func_prime_call = 0
		count_of_func_prime_2_call = 0
		epsilon = pow(10, -i)

		result, time_of_work, count_of_iterations_golden_ratio, count_of_iterations_Newton_Rav_method = Newton_Rav_method(-pi/2, pi/2, epsilon)
		data_of_time += [time_of_work]
		data_of_iter_gold += [count_of_iterations_golden_ratio]
		data_of_iter_New_Rav += [count_of_iterations_Newton_Rav_method]

		count_of_func_prime_call = 0
		result, time_of_work, count_of_iterations_golden_ratio, count_of_iterations_Newton_Rav_method = Newton_Rav_method(0, pi, epsilon, type_of_method = "prime")
		data_of_time_modif += [time_of_work]
		data_of_iter_New_Rav_modif += [count_of_iterations_Newton_Rav_method]
		data_of_func_prime_call_modif += [count_of_func_prime_call]

		ax_text += ["$ { 10}^{" + str(-i) + "} (" + str(data_of_time[-1]) +"с)(" + str(data_of_time_modif[-1]) +"с)$"]

	axis[1].clear()
	axis[1].plot(ax, data_of_iter_gold, "g", alpha = 0.7, label = "Количество итераций метода золотого сечения")
	axis[1].plot(ax, data_of_iter_New_Rav, "k", label = "Количество итераций метода Ньютона-Рассела, количество вызовов первой и второй производной функции")
	axis[1].plot(ax, data_of_iter_New_Rav_modif, "b", alpha = 0.7, label = "Количество итераций модифицированного метода Ньютона-Рассела")
	axis[1].plot(ax, data_of_func_prime_call_modif, "--r", label = "Количество вызовов первой производной функции (мод. метод)")
	axis[1].set_title("Метод Ньютона-Рассела")
	axis[1].set_xlabel("Степень погрешности")
	axis[1].set_xticks(ax, ax_text)
	axis[1].grid()
	axis[1].legend()

	plt.draw()
	plt.ioff()

count_of_func_call, count_of_func_prime_call, count_of_func_prime_2_call = 0, 0, 0
'''
count_of_iterations = 0
epsilon = pow(10, -20)
start_time = time.time()

points = golden_ratio(0, pi, 1e-1)
result = tangent_method(points[0], points[1], epsilon)
time_of_work = round(time.time() - start_time, 3)

print(f"Время работы: {time_of_work}")
print(f"Количество итераций: {count_of_iterations}")
print(f"Степень требуемой точности: {round(log(epsilon, 10), 1)}")
print(f"Степень истинной погрешности: {round(log(abs(pi/2 - result), 10), 1)}")
'''

fg, axis = plt.subplots(nrows= 2, ncols= 1, figsize=(15, 8))
fg.subplots_adjust(bottom=0.1, left=0.06, right=0.98, top=0.95, hspace=0.3)

grafics()

plt.show()