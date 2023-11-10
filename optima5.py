import matplotlib.pyplot as plt
from mpmath import *
from numpy import array, eye, matmul, transpose
import time


def F(arg : array):
	global count_of_func_call
	count_of_func_call += 1
	return 9*arg[0][0]**4 + arg[1][0]**2
	#return (1 - arg[0])**2 + 100*(arg[1] - arg[0]**2)**2


def gradF(arg : array):
	global count_of_gradF_call
	count_of_gradF_call += 1
	return array([[36*arg[0][0]**3], [2*arg[1][0]]])
	#return [-2 + 2*arg[0] - 400*arg[0]*arg[1] + 400*arg[0]**3, 200*arg[1] - 200*arg[0]**2]


def norm(arg : array):
	return pow(sum([pow(coord[0], 2) for coord in arg]), 0.5)


def golden_ratio(start_a : mpf, start_b : mpf, epsilon : mpf, F, type_of_search : str = "min"):
	count_of_iterations = 1

	left_phi = (3 - pow(mpf(5), 0.5)) / 2
	right_phi = (pow(mpf(5), 0.5) - 1) / 2
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

	return (start_b + start_a) / 2, count_of_iterations


def quasi_Newton_method_ND(start_point : array, epsilon : mpf = 1e-16):

	start_time = time.time()
	count_of_iterat = 0

	old_point_2 = start_point
	grad_point_2 = gradF(old_point_2)

	D1F = lambda x: F(old_point_2 - x * grad_point_2)
	step, _ = golden_ratio(0, 1, 1e-12, D1F)
	old_point_1 = old_point_2 - step * grad_point_2
	grad_point_1 = gradF(old_point_1)

	dim = len(start_point)
	H = eye(dim)

	while norm(grad_point_1) > epsilon:
		count_of_iterat += 1

		if count_of_iterat % (10 * dim) != 0:
			delta = old_point_1 - old_point_2
			gamma = grad_point_1 - grad_point_2
			tran_delta_minus_H_gamma = transpose(delta - matmul(H, gamma))
			corr_matrix = matmul(transpose(tran_delta_minus_H_gamma), tran_delta_minus_H_gamma) / matmul(tran_delta_minus_H_gamma, gamma)
			H = H + corr_matrix

		else:
			H = eye(dim)

		H_mul_grad = matmul(H, grad_point_1)
		D1F = lambda x: F(old_point_1 - x * H_mul_grad)
		step, _ = golden_ratio(0, 1, 1e-12, D1F)

		old_point_2 = old_point_1
		old_point_1 = old_point_1 - step * H_mul_grad
		grad_point_2 = grad_point_1
		grad_point_1 = gradF(old_point_1)

	return old_point_1, count_of_iterat, round(time.time() - start_time, 3)

def Fletcher_Reeves_method(start_point : array, epsilon : mpf = 1e-16):

	start_time = time.time()
	count_of_iterat = 0

	dim = len(start_point)

	old_point = start_point
	old_grad = gradF(old_point)

	vec_d = -old_grad

	D1F = lambda x: F(old_point + x * vec_d)
	step, _ = golden_ratio(0, 1, 1e-12, D1F)
	new_point = old_point + step * vec_d
	new_grad = gradF(new_point)

	while norm(new_grad) > epsilon:
		count_of_iterat += 1

		if count_of_iterat % (dim * 10) == 0:
			vec_d = -new_grad
		else:
			vec_d = -new_grad + (norm(new_grad) / norm(old_grad))**2 * vec_d

		D1F = lambda x: F(new_point + x * vec_d)
		step, _ = golden_ratio(0, 10000, 1e-12, D1F)

		old_point = new_point
		new_point = new_point + step * vec_d
		old_grad = new_grad
		new_grad = gradF(new_point)

	return new_point, count_of_iterat, round(time.time() - start_time, 3)

def grafics() -> None:
	global count_of_func_call, count_of_gradF_call

	plt.ion()

	ax, ax_text = [], []
	data_of_time, data_of_iter, data_of_func_call, data_of_gradF_call = [], [], [], []
	for i in range(10, 17):
		ax += [i]

		count_of_func_call, count_of_gradF_call = 0, 0
		epsilon = pow(10, -i)

		result, count_of_iterations, time_of_work = quasi_Newton_method_ND(start_point = array([[mpf(1)], [mpf(1)]]), epsilon = epsilon)
		ax_text += ["$ { 10}^{" + str(-i) + "} (" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter += [count_of_iterations]
		data_of_func_call += [count_of_func_call / 10]
		data_of_gradF_call += [count_of_gradF_call]

	axis[0].clear()
	axis[0].plot(ax, data_of_iter, "g", label = "Количество итераций")
	axis[0].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции (/10)")
	axis[0].plot(ax, data_of_gradF_call, "r", label = "Количество вызовов градиента функции")
	axis[0].set_title("Квази-Нютоновский метод")
	axis[0].set_xlabel("Степень погрешности")
	axis[0].set_xticks(ax, ax_text)
	axis[0].grid()
	axis[0].legend()

	
	ax, ax_text = [], []
	data_of_time, data_of_iter, data_of_func_call, data_of_gradF_call = [], [], [], []
	for i in range(10, 17):
		ax += [i]

		count_of_func_call, count_of_gradF_call = 0, 0
		epsilon = pow(10, -i)

		result, count_of_iterations, time_of_work = Fletcher_Reeves_method(start_point = array([[mpf(1)], [mpf(1)]]), epsilon = epsilon)
		ax_text += ["$ { 10}^{" + str(-i) + "} (" + str(time_of_work) +"с.)$"]
		data_of_time += [time_of_work]
		data_of_iter += [count_of_iterations]
		data_of_func_call += [count_of_func_call / 10]
		data_of_gradF_call += [count_of_gradF_call]

	axis[1].clear()
	axis[1].plot(ax, data_of_iter, "g", label = "Количество итераций")
	axis[1].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции (/10)")
	axis[1].plot(ax, data_of_gradF_call, "r", label = "Количество вызовов градиента функции")
	axis[1].set_title("Метод Флетчара-Ривза")
	axis[1].set_xlabel("Степень погрешности")
	axis[1].set_xticks(ax, ax_text)
	axis[1].grid()
	axis[1].legend()


	plt.draw()
	plt.ioff()

mp.dps = 30
start_point = array([[5], [8]])

count_of_func_call, count_of_gradF_call = 0, 0

fg, axis = plt.subplots(nrows= 2, ncols= 1, figsize=(15, 8))
fg.subplots_adjust(bottom=0.1, left=0.06, right=0.98, top=0.95, hspace=0.3, wspace = 0.15)

grafics()

plt.show()

'''
print(quasi_Newton_method_ND(start_point))
print(Fletcher_Reeves_method(start_point))
'''