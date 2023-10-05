import matplotlib.pyplot as plt
from mpmath import *
import time

mp.dps = 200
def F(arg : list[mpf]):
	global count_of_func_call
	count_of_func_call += 1
	return 9*arg[0]**2 + arg[1]**2

def gradF(arg : list[mpf]):
	global count_of_gradF_call
	count_of_gradF_call += 1
	return [18*arg[0], 2*arg[1]]

def norm(arg : list[mpf]):
	return pow(sum([pow(coord, 2) for coord in arg]), 0.5)

def copy_point(arg : list[mpf]):
	new_point = []
	for item in arg:
		new_point += [item]
	return new_point

def sub_lists(list1 : list[mpf], list2 : list[mpf]):
	new_list = []
	for id_item in range(len(list1)):
		new_list += [list1[id_item] - list2[id_item]]
	return new_list

def mult_list_on_num(arg : list[mpf], num : mpf):
	return [coord*num for coord in arg]

def vert_descent(start_point : list[mpf], step : mpf = 1, milling_coef : mpf = 0.5, epsilon : mpf = 1e-16, dim_num : int = 2):
	start_time = time.time()
	count_of_iterat = 0

	old_min_value_F = F(start_point)
	old_cycle_start_point = copy_point(start_point)

	while True:
		changing_vector = False
		for iterat in range(dim_num):
			new_point = copy_point(start_point)
			count_of_iterat += 1

			new_point[iterat] += step
			new_value_F = F(new_point)
			if new_value_F < old_min_value_F: 
				old_min_value_F = new_value_F
				start_point = copy_point(new_point)
				changing_vector = True
				continue

			new_point[iterat] -= 2*step
			new_value_F = F(new_point)
			if new_value_F < old_min_value_F: 
				old_min_value_F = new_value_F
				start_point = copy_point(new_point)
				changing_vector = True
				continue

		if not changing_vector: step *= milling_coef
		if norm([start_point[i] - old_cycle_start_point[i] for i in range(dim_num)]) < epsilon and step < epsilon * milling_coef: break
		old_cycle_start_point = copy_point(start_point)

	return start_point, count_of_iterat, round(time.time() - start_time, 3)

def gradient_descent(start_point : list[mpf], step : mpf = 1, milling_coef : mpf = 0.5, epsilon : mpf = 1e-16, delta : mpf = 0.5, milling_steps = False, \
						harmopnic_series = False):
	start_time = time.time()
	count_of_iterat = 0

	default_step = step
	gradF_value = gradF(start_point)
	harmopnic_series_coef = mpf(2)

	while norm(gradF_value) > epsilon:
		count_of_iterat += 1

		new_point = sub_lists(start_point, mult_list_on_num(gradF_value, step))
		if F(new_point) - F(start_point) < -step*delta*norm(gradF_value)**2:
			if milling_steps: step = default_step
			start_point = copy_point(new_point)
			gradF_value = gradF(start_point)
		else:
			if harmopnic_series: 
				step = 1/harmopnic_series_coef
				harmopnic_series_coef += 1
			else: step *= milling_coef

	return start_point, count_of_iterat, round(time.time() - start_time, 3)

def grafics() -> None:
	global count_of_func_call, count_of_gradF_call

	plt.ion()

	ax, ax_text = [], []
	data_of_time, data_of_iter, data_of_func_call = [], [], []
	for i in range(100, 160, 10):
		ax += [i]

		count_of_func_call, count_of_gradF_call = 0, 0
		epsilon = pow(10, -i)

		result, count_of_iterations, time_of_work = vert_descent(start_point = [100, 100], epsilon = epsilon)
		ax_text += ["$ { 10}^{" + str(-i) + "} (" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter += [count_of_iterations]
		data_of_func_call += [count_of_func_call]

	axis[0, 0].clear()
	axis[0, 0].plot(ax, data_of_iter, "g", label = "Количество итераций")
	axis[0, 0].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции")
	axis[0, 0].set_title("Метод покоординатного спуска")
	axis[0, 0].set_xlabel("Степень погрешности")
	axis[0, 0].set_xticks(ax, ax_text)
	axis[0, 0].grid()
	axis[0, 0].legend()

	ax, ax_text = [], []
	data_of_time, data_of_iter, data_of_func_call, data_of_gradF_call = [], [], [], []
	for i in range(100, 160, 10):
		ax += [i]

		count_of_func_call, count_of_gradF_call = 0, 0
		epsilon = pow(10, -i)

		result, count_of_iterations, time_of_work = gradient_descent(start_point = [100, 100], epsilon = epsilon, milling_steps = True)
		ax_text += ["$ { 10}^{" + str(-i) + "} (" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter += [count_of_iterations]
		data_of_func_call += [count_of_func_call]
		data_of_gradF_call += [count_of_gradF_call]

	axis[0, 1].clear()
	axis[0, 1].plot(ax, data_of_iter, "g", label = "Количество итераций")
	axis[0, 1].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции")
	axis[0, 1].plot(ax, data_of_gradF_call, "r", label = "Количество вызовов градиента функции")
	axis[0, 1].set_title("Метод градиентного спуска с дроблением шага")
	axis[0, 1].set_xlabel("Степень погрешности")
	axis[0, 1].set_xticks(ax, ax_text)
	axis[0, 1].grid()
	axis[0, 1].legend()

	ax, ax_text = [], []
	data_of_time, data_of_iter, data_of_func_call, data_of_gradF_call = [], [], [], []
	for i in range(100, 160, 10):
		ax += [i]

		count_of_func_call, count_of_gradF_call = 0, 0
		epsilon = pow(10, -i)

		result, count_of_iterations, time_of_work = gradient_descent(start_point = [100, 100], epsilon = epsilon, milling_coef = 0.99)
		ax_text += ["$ { 10}^{" + str(-i) + "} (" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter += [count_of_iterations]
		data_of_func_call += [count_of_func_call]
		data_of_gradF_call += [count_of_gradF_call]

	axis[1, 0].clear()
	axis[1, 0].plot(ax, data_of_iter, "g", label = "Количество итераций")
	axis[1, 0].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции")
	axis[1, 0].plot(ax, data_of_gradF_call, "r", label = "Количество вызовов градиента функции")
	axis[1, 0].set_title("Метод градиентного спуска с постоянным шагом")
	axis[1, 0].set_xlabel("Степень погрешности")
	axis[1, 0].set_xticks(ax, ax_text)
	axis[1, 0].grid()
	axis[1, 0].legend()

	ax, ax_text = [], []
	data_of_time, data_of_iter, data_of_func_call, data_of_gradF_call = [], [], [], []
	for i in range(100, 160, 10):
		ax += [i]

		count_of_func_call, count_of_gradF_call = 0, 0
		epsilon = pow(10, -i)

		result, count_of_iterations, time_of_work = gradient_descent(start_point = [100, 100], epsilon = epsilon, milling_coef = 0.99, harmopnic_series = True)
		ax_text += ["$ { 10}^{" + str(-i) + "} (" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter += [count_of_iterations]
		data_of_func_call += [count_of_func_call]
		data_of_gradF_call += [count_of_gradF_call]

	axis[1, 1].clear()
	axis[1, 1].plot(ax, data_of_iter, "g", label = "Количество итераций")
	axis[1, 1].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции")
	axis[1, 1].plot(ax, data_of_gradF_call, "r", label = "Количество вызовов градиента функции")
	axis[1, 1].set_title("Метод градиентного спуска с постоянным шагом и расходящимся рядом")
	axis[1, 1].set_xlabel("Степень погрешности")
	axis[1, 1].set_xticks(ax, ax_text)
	axis[1, 1].grid()
	axis[1, 1].legend()

	plt.draw()
	plt.ioff()

count_of_func_call, count_of_gradF_call = 0, 0

fg, axis = plt.subplots(nrows= 2, ncols= 2, figsize=(15, 8))
fg.subplots_adjust(bottom=0.1, left=0.06, right=0.98, top=0.95, hspace=0.3, wspace = 0.15)

grafics()

plt.show()

'''
print(vert_descent(start_point = [1, 1], step = 10, milling_coef = 0.3, epsilon = 1e-150))
print(gradient_descent(start_point = [1, 1], step = 10, milling_coef = 0.3, epsilon = 1e-150, delta = 0.5, milling_steps = True))
print(gradient_descent(start_point = [1, 1], step = 10, milling_coef = 0.3, epsilon = 1e-150, delta = 0.5, milling_steps = False))
print(gradient_descent(start_point = [1, 1], step = 10, milling_coef = 0.3, epsilon = 1e-150, delta = 0.5, milling_steps = False, harmopnic_series = True))
'''