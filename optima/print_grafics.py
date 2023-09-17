import matplotlib.pyplot as plt
from mpmath import *
from os import startfile

def deg_mpf(arg : mpf) -> nstr:
	return "{ 10}^{" + str(int(floor(log(arg, 10)))) +"}с."

def grafics() -> None:

	file = open("passiv_search_data.txt")
	ax = []
	data_of_time = []
	data_of_iterat = []
	data_of_func_call = []
	idx = 0

	for line in file:
		if idx == 5: 
			idx = 0
			continue
		if idx == 0: data_of_time += [mpf(line)]
		elif idx == 1: data_of_iterat += [int(line)]
		elif idx == 2: data_of_func_call += [int(line)]
		elif idx == 4: ax += [-int(log(float(line), 10))]
		idx += 1

	ax_text = ["$ { 10}^{" + str(-ax[i]) + "} (" + deg_mpf(data_of_time[i])+")$" for i in range(len(ax))]

	plt.ion()
	axis[0, 0].clear()
	axis[0, 0].plot(ax, data_of_iterat, "r", label = "Количество итераций")
	axis[0, 0].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции")
	axis[0, 0].set_title("Метод пассивного поиска")
	axis[0, 0].set_xlabel("Степень погрешности")
	axis[0, 0].set_xticks(ax, ax_text)
	axis[0, 0].grid()
	axis[0, 0].legend()

	file = open("dihotomia_search_data.txt")
	ax = []
	data_of_time = []
	data_of_iterat = []
	data_of_func_call = []
	idx = 0
	
	for line in file:
		if idx == 5: 
			idx = 0
			continue
		if idx == 0: data_of_time += [mpf(line)]
		elif idx == 1: data_of_iterat += [int(line)]
		elif idx == 2: data_of_func_call += [int(line)]
		elif idx == 4: ax += [-int(log(float(line), 10))]
		idx += 1

	ax_text = ["$ { 10}^{" + str(-ax[i]) + "} (" + deg_mpf(data_of_time[i])+")$" for i in range(len(ax))]

	axis[1, 0].clear()
	axis[1, 0].plot(ax, data_of_iterat, "r", label = "Количество итераций")
	axis[1, 0].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции")
	axis[1, 0].set_title("Метод половинного поиска")
	axis[1, 0].set_xlabel("Степень погрешности")
	axis[1, 0].set_xticks(ax, ax_text)
	axis[1, 0].set_yticks(data_of_func_call, data_of_func_call)
	axis[1, 0].grid()
	axis[1, 0].legend()

	file = open("golden_ratio_search_data.txt")
	ax = []
	data_of_time = []
	data_of_iterat = []
	data_of_func_call = []
	idx = 0
	
	for line in file:
		if idx == 5: 
			idx = 0
			continue
		if idx == 0: data_of_time += [mpf(line)]
		elif idx == 1: data_of_iterat += [int(line)]
		elif idx == 2: data_of_func_call += [int(line)]
		elif idx == 4: ax += [-int(log(float(line), 10))]
		idx += 1

	ax_text = ["$ { 10}^{" + str(-ax[i]) + "} (" + deg_mpf(data_of_time[i])+")$" for i in range(len(ax))]

	axis[0, 1].clear()
	axis[0, 1].plot(ax, data_of_iterat, "r", label = "Количество итераций")
	axis[0, 1].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции")
	axis[0, 1].set_title("Метод золотого сечения")
	axis[0, 1].set_xlabel("Степень погрешности")
	axis[0, 1].set_xticks(ax, ax_text)
	axis[0, 1].set_yticks(data_of_func_call, data_of_func_call)
	axis[0, 1].grid()
	axis[0, 1].legend()

	file = open("fibonacci_search_data.txt")
	ax = []
	data_of_time = []
	data_of_iterat = []
	data_of_func_call = []
	idx = 0
	
	for line in file:
		if idx == 5: 
			idx = 0
			continue
		if idx == 0: data_of_time += [mpf(line)]
		elif idx == 1: data_of_iterat += [int(line)]
		elif idx == 2: data_of_func_call += [int(line)]
		elif idx == 4: ax += [-int(log(float(line), 10))]
		idx += 1

	ax_text = ["$ { 10}^{" + str(-ax[i]) + "} (" + deg_mpf(data_of_time[i])+")$" for i in range(len(ax))]

	axis[1, 1].clear()
	axis[1, 1].plot(ax, data_of_iterat, "r", label = "Количество итераций")
	axis[1, 1].plot(ax, data_of_func_call, "k", label = "Количество вызовов функции")
	axis[1, 1].set_title("Метод чисел Фибоначчи")
	axis[1, 1].set_xlabel("Степень погрешности")
	axis[1, 1].set_xticks(ax, ax_text)
	axis[1, 1].set_yticks(data_of_func_call, data_of_func_call)
	axis[1, 1].grid()
	axis[1, 1].legend()

	plt.draw()
	plt.ioff()

fg, axis = plt.subplots(nrows= 2, ncols= 2, figsize=(15, 8))
fg.subplots_adjust(bottom=0.1, left=0.06, right=0.98, top=0.95, hspace=0.3)

grafics()

plt.show()