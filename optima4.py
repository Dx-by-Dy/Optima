import matplotlib.pyplot as plt
from mpmath import *
import time

mp.dps = 50
def F(arg : list[mpf]):
	global count_of_func_call
	count_of_func_call += 1
	return 9*arg[0]**2 + arg[1]**2
	#return (1 - arg[0])**2 + 100*(arg[1] - arg[0]**2)**2

def gradF(arg : list[mpf]):
	global count_of_gradF_call
	count_of_gradF_call += 1
	return [18*arg[0], 2*arg[1]]
	#return [-2 + 2*arg[0] - 400*arg[0]*arg[1] + 400*arg[0]**3, 200*arg[1] - 200*arg[0]**2]

def invGessianF_mul_gradF(arg : list[mpf]):
	global count_of_invGessianF_call
	count_of_invGessianF_call += 1
	return [arg[0], arg[1]]

def norm(arg : list[mpf]):
	return pow(sum([pow(coord, 2) for coord in arg]), 0.5)

def sub_lists(list1 : list[mpf], list2 : list[mpf]):
	new_list = []
	for id_item in range(len(list1)):
		new_list += [list1[id_item] - list2[id_item]]
	return new_list

def mult_list_on_num(arg : list[mpf], num : mpf):
	return [coord*num for coord in arg]

def golden_ratio(start_a : mpf, start_b : mpf, epsilon : mpf, F, type_of_search : str = "min"):
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

	return (start_b - start_a) / 2, count_of_iterations

def busted_gradient_descent(start_point : list[mpf], epsilon : mpf = 1e-16, deg_method = 10):

	start_time = time.time()
	count_of_iterat = 0

	gradF_value = gradF(start_point)
	old_point = start_point

	while norm(gradF_value) > epsilon:
		count_of_iterat += 1

		for iterat in range(deg_method):
			count_of_iterat += 1

			D1F = lambda x: F(sub_lists(start_point, mult_list_on_num(gradF_value, x)))

			step, _ = golden_ratio(0, 100, 1e-1, D1F)
			start_point = sub_lists(start_point, mult_list_on_num(gradF_value, step))
			gradF_value = gradF(start_point)

		if deg_method > 1:
			D1F = lambda x: F(sub_lists(start_point, mult_list_on_num(sub_lists(start_point, old_point), x)))

			step, _ = golden_ratio(0, 100, 1e-1, D1F)
			start_point = sub_lists(start_point, mult_list_on_num(sub_lists(start_point, old_point), step))
			gradF_value = gradF(start_point)

			old_point = start_point

	return start_point, count_of_iterat, round(time.time() - start_time, 3)

def gully_method(start_point : list[mpf], epsilon : mpf = 1e-16, delta : mpf = 1e-18, deg_method = 10):

	start_time = time.time()
	count_of_iterat = 0

	dim = len(start_point)
	second_point = sub_lists(start_point, [delta for i in range(dim)])

	gradF_value_first_point = gradF(start_point)
	gradF_value_second_point = gradF(second_point)

	while norm(gradF_value_first_point) > epsilon:
		count_of_iterat += 1

		second_point = sub_lists(start_point, [delta for i in range(dim)])
		gradF_value_second_point = gradF(second_point)

		for iterat in range(deg_method):
			count_of_iterat += 1

			D1F = lambda x: F(sub_lists(start_point, mult_list_on_num(gradF_value_first_point, x)))

			step, _ = golden_ratio(0, 100, 1e-1, D1F)
			start_point = sub_lists(start_point, mult_list_on_num(gradF_value_first_point, step))
			gradF_value_first_point = gradF(start_point)

			D1F = lambda x: F(sub_lists(second_point, mult_list_on_num(gradF_value_second_point, x)))

			step, _ = golden_ratio(0, 100, 1e-1, D1F)
			second_point = sub_lists(second_point, mult_list_on_num(gradF_value_second_point, step))
			gradF_value_second_point = gradF(second_point)

		D1F = lambda x: F(sub_lists(start_point, mult_list_on_num(sub_lists(start_point, second_point), x)))

		step, _ = golden_ratio(0, 100, 1e-1, D1F)
		start_point = sub_lists(start_point, mult_list_on_num(sub_lists(start_point, second_point), step))
		gradF_value_first_point = gradF(start_point)

	return start_point, count_of_iterat, round(time.time() - start_time, 3)

def modified_Newton_method_ND(start_point : list[mpf], epsilon : mpf = 1e-16):

	start_time = time.time()
	count_of_iterat = 0

	old_point = start_point
	start_point = sub_lists(start_point, invGessianF_mul_gradF(start_point))

	while norm(sub_lists(start_point, old_point)) > epsilon:
		count_of_iterat += 1

		D1F = lambda x: F(sub_lists(start_point, mult_list_on_num(invGessianF_mul_gradF(start_point), x)))
		step, _ = golden_ratio(0, 100, 1e-1, D1F)

		old_point = start_point
		start_point = sub_lists(start_point, mult_list_on_num(invGessianF_mul_gradF(start_point), step))

	return start_point, count_of_iterat, round(time.time() - start_time, 3)

def grafics() -> None:
	global count_of_func_call, count_of_gradF_call, count_of_invGessianF_call

	plt.ion()

	ax, ax_text = [], []
	data_of_time, data_of_iter, data_of_func_call, data_of_gradF_call = [], [], [], []
	for i in range(10, 16):
		ax += [i]

		count_of_func_call, count_of_gradF_call = 0, 0
		epsilon = pow(10, -i)

		result, count_of_iterations, time_of_work = busted_gradient_descent(start_point = [mpf(1), mpf(1)], epsilon = epsilon)
		ax_text += ["$ { 10}^{" + str(-i) + "} (" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter += [count_of_iterations]
		data_of_func_call += [count_of_func_call / 10]
		data_of_gradF_call += [count_of_gradF_call]

	axis[0, 0].clear()
	axis[0, 0].plot(ax, data_of_iter, "g", label = "Количество итераций")
	axis[0, 0].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции (/10)")
	axis[0, 0].plot(ax, data_of_gradF_call, "r", label = "Количество вызовов градиента функции")
	axis[0, 0].set_title("Метод ускоренного градиентного спуска 10 порядка")
	axis[0, 0].set_xlabel("Степень погрешности")
	axis[0, 0].set_xticks(ax, ax_text)
	axis[0, 0].grid()
	axis[0, 0].legend()


	ax, ax_text = [], []
	data_of_time, data_of_iter, data_of_func_call, data_of_gradF_call = [], [], [], []
	for i in range(10, 16):
		ax += [i]

		count_of_func_call, count_of_gradF_call = 0, 0
		epsilon = pow(10, -i)

		result, count_of_iterations, time_of_work = gully_method(start_point = [mpf(1), mpf(1)], epsilon = epsilon)
		ax_text += ["$ { 10}^{" + str(-i) + "} (" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter += [count_of_iterations]
		data_of_func_call += [count_of_func_call / 10]
		data_of_gradF_call += [count_of_gradF_call]

	axis[0, 1].clear()
	axis[0, 1].plot(ax, data_of_iter, "g", label = "Количество итераций")
	axis[0, 1].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции (/10)")
	axis[0, 1].plot(ax, data_of_gradF_call, "r", label = "Количество вызовов градиента функции")
	axis[0, 1].set_title("Овражный метод")
	axis[0, 1].set_xlabel("Степень погрешности")
	axis[0, 1].set_xticks(ax, ax_text)
	axis[0, 1].grid()
	axis[0, 1].legend()


	ax, ax_text = [], []
	data_of_time, data_of_iter, data_of_func_call, data_of_gradF_call = [], [], [], []
	for i in range(1, 6):
		ax += [i]

		count_of_func_call, count_of_gradF_call = 0, 0

		result, count_of_iterations, time_of_work = busted_gradient_descent(start_point = [mpf(1), mpf(1)], deg_method = i)
		ax_text += ["$" + str(i) + "(" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter += [count_of_iterations]
		data_of_func_call += [count_of_func_call / 10]
		data_of_gradF_call += [count_of_gradF_call]

	axis[1, 0].clear()
	axis[1, 0].plot(ax, data_of_iter, "g", label = "Количество итераций")
	axis[1, 0].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции (/10)")
	axis[1, 0].plot(ax, data_of_gradF_call, "r", label = "Количество вызовов градиента функции")
	axis[1, 0].set_title("Метод ускоренного градиентного спуска 10 порядка")
	axis[1, 0].set_xlabel("Степень метода")
	axis[1, 0].set_xticks(ax, ax_text)
	axis[1, 0].grid()
	axis[1, 0].legend()


	ax, ax_text = [], []
	data_of_time, data_of_iter, data_of_func_call, data_of_gradF_call = [], [], [], []
	for i in range(1, 6):
		ax += [i]

		count_of_func_call, count_of_gradF_call = 0, 0

		result, count_of_iterations, time_of_work = gully_method(start_point = [mpf(1), mpf(1)], deg_method = i)
		ax_text += ["$" + str(i) + "(" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter += [count_of_iterations]
		data_of_func_call += [count_of_func_call / 10]
		data_of_gradF_call += [count_of_gradF_call]

	axis[1, 1].clear()
	axis[1, 1].plot(ax, data_of_iter, "g", label = "Количество итераций")
	axis[1, 1].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции (/10)")
	axis[1, 1].plot(ax, data_of_gradF_call, "r", label = "Количество вызовов градиента функции")
	axis[1, 1].set_title("Овражный метод")
	axis[1, 1].set_xlabel("Степень метода")
	axis[1, 1].set_xticks(ax, ax_text)
	axis[1, 1].grid()
	axis[1, 1].legend()


	plt.draw()
	plt.ioff()


count_of_func_call, count_of_gradF_call, count_of_invGessianF_call = 0, 0, 0

fg, axis = plt.subplots(nrows= 2, ncols= 2, figsize=(15, 8))
fg.subplots_adjust(bottom=0.1, left=0.06, right=0.98, top=0.95, hspace=0.3, wspace = 0.15)

grafics()

plt.show()

'''
print(busted_gradient_descent(start_point = [mpf(1), mpf(1)], deg_method = 10))
print(gully_method(start_point = [mpf(1), mpf(1)], deg_method = 10))
print(modified_Newton_method_ND(start_point = [mpf(1), mpf(1)]))
'''