#include <iostream>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <chrono>
#define PI 3.14159265358979
#define NANOSEC chrono::duration_cast<chrono::nanoseconds>(chrono::system_clock::now().time_since_epoch())

using namespace std;

int count_F_call;
int count_of_iterations;
auto start_time = NANOSEC;
int deg_of_err_rate = 6;
string type_F = "-sin(x)";
long double epsilon = pow(10, -deg_of_err_rate);

long double msin(long double arg) {
    arg -= (int)(arg / (2 * PI)) * 2 * PI;
    long double value_term;
    long double result = 0;

    for (int term = 0; term < 14; ++term) {
        value_term = 1;
        for (int mult = 0; mult < 2 * term + 1; ++mult) {
            value_term *= arg / ((double)mult + 1);
        }
        if (term % 2 == 0) result += value_term;
        else result -= value_term;
    }

    return result;
}

long double F(long double arg) {
    count_F_call += 1;
    if (type_F == "-sin(x)") return -msin(arg);
    else if (type_F == "1/x") return 1 / arg;
}

void set_type_F(string type) {
    type_F = type;
}

void set_deg_err_rate(int deg) {
    deg_of_err_rate = deg;
    epsilon = pow(10, -deg);
}

void print_method_data(string type_of_method, long double start_a, long double start_b, long double extr_point) {
    auto end_time = NANOSEC;
    if (type_of_method == "passiv_search") cout <<            "------- Метод пассивного поиска ------" << endl;
    else if (type_of_method == "dihotomia_search") cout <<    "------ Метод половинного поиска ------" << endl;
    else if (type_of_method == "golden_ratio_search") cout << "------- Метод золотого сечения -------" << endl;
    else if (type_of_method == "fibonacci_search") cout <<    "------- Метод чисел Фибоначчи --------" << endl;

    cout.precision(7);
    cout << "Функция:                     " << type_F << " (" << start_a << "; " << start_b << ")" << endl;
    cout << "Время работы программы:      " << (end_time.count() - start_time.count()) / pow(10, 9) << " сек." << endl;
    cout.precision(deg_of_err_rate + 1);
    cout << "Количество итераций:         " << count_of_iterations << endl;
    cout << "Количество вызовов функции:  " << count_F_call << endl;
    cout << "Точка экстремума функции:    " << extr_point << endl;
    cout.precision(deg_of_err_rate);
    cout << "Погрешность:                 " << epsilon << endl << endl;
}

long double passiv_search(long double start_a, long double start_b, long double epsilon, string type_of_search = "min") {
    count_F_call = 0;
    count_of_iterations = 1;
    start_time = NANOSEC;

    long double old_F_value = F(start_a);
    long double new_F_value;
    start_a += epsilon;

    while (start_a < start_b){
        count_of_iterations += 1;
        new_F_value = F(start_a);

        if (type_of_search == "min") {
            if (new_F_value < old_F_value) old_F_value = new_F_value;
            else return start_a - epsilon;
        }
        else if (type_of_search == "max") {
            if (new_F_value > old_F_value) old_F_value = new_F_value;
            else return start_a - epsilon;
        }
        start_a += epsilon;
    }

    return start_b;
}

long double dihotomia_search(long double start_a, long double start_b, long double epsilon, string type_of_search = "min") {
    count_F_call = 0;
    count_of_iterations = 0;
    start_time = NANOSEC;

    long double delta = epsilon / 8;
    long double left_point;
    long double right_point;
    long double left_value;
    long double right_value;

    while (start_b - start_a > 2 * epsilon) {
        count_of_iterations += 1;
        left_point = (start_b + start_a) / 2 - delta;
        right_point = (start_b + start_a) / 2 + delta;
        left_value = F(left_point);
        right_value = F(right_point);

        if (type_of_search == "min") {
            if (left_value < right_value) start_b = right_point;
            else start_a = left_point;
        }
        else if (type_of_search == "max") {
            if (left_value < right_value) start_a = left_point;
            else start_b = right_point;
        }
    }

    return (start_b + start_a) / 2;
}

long double golden_ratio_search(long double start_a, long double start_b, long double epsilon, bool mirror_point = 0, string type_of_search = "min") {
    count_F_call = 0;
    count_of_iterations = 1;
    start_time = NANOSEC;
    
    long double left_phi = (3 - pow(5, 0.5)) / 2;
    long double right_phi = (pow(5, 0.5) - 1) / 2;
    long double left_point = start_a + left_phi * (start_b - start_a);
    long double right_point = start_a + right_phi * (start_b - start_a);
    long double left_value = F(left_point);
    long double right_value = F(right_point);

    while (start_b - start_a > 2 * epsilon) {
        count_of_iterations += 1;

        if (type_of_search == "min") {
            if (left_value < right_value) {
                start_b = right_point;
                right_point = left_point;
                right_value = left_value;
                if (mirror_point) left_point = start_a + start_b - right_point;
                else left_point = start_a + left_phi * (start_b - start_a);
                left_value = F(left_point);
            }
            else {
                start_a = left_point;
                left_point = right_point;
                left_value = right_value;
                if (mirror_point) right_point = start_b - (left_point - start_a);
                else right_point = start_a + right_phi * (start_b - start_a);
                right_value = F(right_point);
            }
        }
        else if (type_of_search == "max") {
            if (left_value < right_value) {
                start_a = left_point;
                left_point = right_point;
                left_value = right_value;
                if (mirror_point) right_point = start_b - (left_point - start_a);
                else right_point = start_a + right_phi * (start_b - start_a);
                right_value = F(right_point);
            }
            else {
                start_b = right_point;
                right_point = left_point;
                right_value = left_value;
                if (mirror_point) left_point = start_a + start_b - right_point;
                else left_point = start_a + left_phi * (start_b - start_a);
                left_value = F(left_point);
            }
        }
    }

    return (start_b + start_a) / 2;
}

vector<long int> get_fibonacci_numbers(long double barrier) {
    vector<long int> numbers{ 1, 1 };
    while (*(numbers.end() - 1) < barrier) {
        numbers.push_back(*(numbers.end() - 1) + *(numbers.end() - 2));
    }
    numbers.push_back(*(numbers.end() - 1) + *(numbers.end() - 2));
    return numbers;
}

long double fibonacci_search(long double start_a, long double start_b, long double epsilon, string type_of_search = "min") {
    count_F_call = 0;
    count_of_iterations = 1;
    start_time = NANOSEC;
    
    vector<long int> fib_num = get_fibonacci_numbers(11 * (start_b - start_a) / (20 * epsilon));
    int count_iter = fib_num.size() - 2;

    long double left_point = start_a + (double)fib_num[count_iter - 1] / fib_num[count_iter + 1] * (start_b - start_a);
    long double right_point = start_a + (double)fib_num[count_iter] / fib_num[count_iter + 1] * (start_b - start_a);
    long double left_value = F(left_point);
    long double right_value = F(right_point);

    while (start_b - start_a > 2 * epsilon) {
        count_of_iterations += 1;

        if (type_of_search == "min") {
            if (left_value < right_value) {
                start_b = right_point;
                right_point = left_point;
                right_value = left_value;
                left_point = start_a + (double)fib_num[count_iter - count_of_iterations] / fib_num[count_iter - count_of_iterations + 2] * (start_b - start_a);
                if (count_iter == count_of_iterations) left_point -= (start_b - start_a) / (10 * fib_num[count_iter + 1]);
                left_value = F(left_point);
            }
            else {
                start_a = left_point;
                left_point = right_point;
                left_value = right_value;
                right_point = start_a + (double)fib_num[count_iter - count_of_iterations + 1] / fib_num[count_iter - count_of_iterations + 2] * (start_b - start_a);
                if (count_iter == count_of_iterations) right_point += (start_b - start_a) / (10 * fib_num[count_iter + 1]);
                right_value = F(right_point);
            }
        }
        else if (type_of_search == "max") {
            if (left_value < right_value) {
                start_a = left_point;
                left_point = right_point;
                left_value = right_value;
                right_point = start_a + (double)fib_num[count_iter - count_of_iterations + 1] / fib_num[count_iter - count_of_iterations + 2] * (start_b - start_a);
                if (count_iter == count_of_iterations) right_point += (start_b - start_a) / (10 * fib_num[count_iter + 1]);
                right_value = F(right_point);
            }
            else {
                start_b = right_point;
                right_point = left_point;
                right_value = left_value;
                left_point = start_a + (double)fib_num[count_iter - count_of_iterations] / fib_num[count_iter - count_of_iterations + 2] * (start_b - start_a);
                if (count_iter == count_of_iterations) left_point -= (start_b - start_a) / (10 * fib_num[count_iter + 1]);
                left_value = F(left_point);
            }
        }
    }

    return (start_b + start_a) / 2;
}

void set_info_in_file() {
    set_type_F("-sin(x)");
    long double start_a = 0;
    long double start_b = PI;

    ofstream out;
    out.open("C:\\Users\\user\\source\\repos\\optima\\passiv_search_data.txt");
    out.setf(std::ios::fixed);
    out.precision(7);

    for (int i = 1; i <= 6; i++) {
        set_deg_err_rate(i);
        long double extr_point = passiv_search(start_a, start_b, epsilon);
        out << (NANOSEC.count() - start_time.count()) / pow(10, 9) << endl;
        out << count_of_iterations << endl;
        out << count_F_call << endl;
        out << extr_point << endl;
        out << epsilon << endl;
        out << '-' << endl;
    }

    out.close();

    out.open("C:\\Users\\user\\source\\repos\\optima\\dihotomia_search_data.txt");
    out.setf(std::ios::fixed);
    out.precision(8);

    for (int i = 1; i <= 7; i++) {
        set_deg_err_rate(i);
        long double extr_point = dihotomia_search(start_a, start_b, epsilon);
        out << (NANOSEC.count() - start_time.count()) / pow(10, 9) << endl;
        out << count_of_iterations << endl;
        out << count_F_call << endl;
        out << extr_point << endl;
        out << epsilon << endl;
        out << '-' << endl;
    }

    out.close();

    out.open("C:\\Users\\user\\source\\repos\\optima\\golden_ratio_search_data.txt");
    out.setf(std::ios::fixed);
    out.precision(8);

    for (int i = 1; i <= 7; i++) {
        set_deg_err_rate(i);
        long double extr_point = golden_ratio_search(start_a, start_b, epsilon);
        out << (NANOSEC.count() - start_time.count()) / pow(10, 9) << endl;
        out << count_of_iterations << endl;
        out << count_F_call << endl;
        out << extr_point << endl;
        out << epsilon << endl;
        out << '-' << endl;
    }

    out.close();

    out.open("C:\\Users\\user\\source\\repos\\optima\\fibonacci_search_data.txt");
    out.setf(std::ios::fixed);
    out.precision(8);

    for (int i = 1; i <= 7; i++) {
        set_deg_err_rate(i);
        long double extr_point = fibonacci_search(start_a, start_b, epsilon);
        out << (NANOSEC.count() - start_time.count()) / pow(10, 9) << endl;
        out << count_of_iterations << endl;
        out << count_F_call << endl;
        out << extr_point << endl;
        out << epsilon << endl;
        out << '-' << endl;
    }

    out.close();
}

int main()
{
    setlocale(LC_CTYPE, "Russian");
    cout.setf(std::ios::fixed);

    long double start_a = 0;
    long double start_b = PI;

    //set_deg_err_rate(3);
    //print_method_data("passiv_search", start_a, start_b, passiv_search(start_a, start_b, epsilon));
    //set_deg_err_rate(7);
    //print_method_data("dihotomia_search", start_a, start_b, dihotomia_search(start_a, start_b, epsilon));
    //print_method_data("golden_ratio_search", start_a, start_b, golden_ratio_search(start_a, start_b, epsilon));
    //print_method_data("fibonacci_search", start_a, start_b, fibonacci_search(start_a, start_b, epsilon));

    set_type_F("1/x");
    start_a = 1;
    start_b = 1000;
    //print_method_data("fibonacci_search", start_a, start_b, fibonacci_search(start_a, start_b, epsilon));

    cout << "-------- Без отражения точки -------" << endl;
    print_method_data("golden_ratio_search", start_a, start_b, golden_ratio_search(start_a, start_b, epsilon));
    cout << "-------- C отражением точки --------" << endl;
    print_method_data("golden_ratio_search", start_a, start_b, golden_ratio_search(start_a, start_b, epsilon, 1));

    //set_info_in_file();
    system("C:\\Users\\user\\source\\repos\\optima\\print_grafics.py");
    
    return 0;
}
