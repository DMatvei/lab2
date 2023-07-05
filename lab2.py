import numpy as np
import sympy as sp
from tabulate import tabulate

import formuls as f 

# start data
a = 0.5
b = 1.0
h = (b-a)/10
arr_x = [round(a + i*h, 2) for i in range(11)]
x_1star = 0.77
x_2star = 0.52
x_3star = 0.97
x_4star = 0.73

# function
fun_x = lambda x: x**2 - np.sin(x)
t = lambda x: (x - f.choice_x0(x, arr_x)) / h



arr_y = f.create_table(fun_x, arr_x)
table_dif_y = f.dif_table(arr_y)

# cheking sins
x = sp.Symbol('x')
fun = x**2 - sp.sin(x)
dif_n, crit_point_min, crit_point_max = f.calc_dif(fun, x, len(table_dif_y[0]) + 1)

fac = 1 / f.factorial(len(table_dif_y[0]) + 1)
min_R = crit_point_min * fac
max_R = crit_point_max * fac

# table output

table_data = {}
table_data['x'] = arr_x
for i in range(len(table_dif_y)):
    table_data[f'f{i}'] = table_dif_y[i]
table = tabulate(table_data, headers="keys", tablefmt="fancy_grid")
print(table)


#  output
print('L(x**) = ' , f.choice_Gaus_or_Newtons(x_2star, h, arr_x, table_dif_y))
print('f(x**) = ', fun_x(x_2star))
print(f'R(x**) = ', "{:0.20f}".format(abs(fun_x(x_2star) - f.choice_Gaus_or_Newtons(x_2star, h, arr_x, table_dif_y))))

print(f'min R = {"{:0.20f}".format(abs(max_R * f.calc_omega(x_2star, arr_x)))}')
print(f'max R = {"{:0.20f}".format(abs(min_R * f.calc_omega(x_2star, arr_x)))}')
print( 'R_min < R < R_max is ', abs(max_R * f.calc_omega(x_2star, arr_x)) < abs(fun_x(x_2star) - f.choice_Gaus_or_Newtons(x_2star, h, arr_x, table_dif_y)) < abs(min_R * f.calc_omega(x_2star, arr_x)))
print()


print('L(x***) = ' ,f.choice_Gaus_or_Newtons(x_3star, h, arr_x, table_dif_y))
print('f(x***) = ', fun_x(x_3star))
print(f'R(x***) = ', "{:0.20f}".format(abs(fun_x(x_3star) - f.choice_Gaus_or_Newtons(x_3star, h, arr_x, table_dif_y))))

print(f'min R = {"{:0.20f}".format(abs(max_R * f.calc_omega(x_3star, arr_x)))}')
print(f'max R = {"{:0.20f}".format(abs(min_R * f.calc_omega(x_3star, arr_x)))}')
print( 'R_min < R < R_max is ', abs(max_R * f.calc_omega(x_3star, arr_x)) < abs(fun_x(x_3star) - f.choice_Gaus_or_Newtons(x_3star, h, arr_x, table_dif_y)) < abs(min_R * f.calc_omega(x_3star, arr_x)))
print()


print('L(x****) = ' ,f.choice_Gaus_or_Newtons(x_4star, h, arr_x, table_dif_y))
print('f(x****) = ', fun_x(x_4star))
print(f'R(x****) = ', "{:0.20f}".format(abs(fun_x(x_4star) - f.choice_Gaus_or_Newtons(x_4star, h, arr_x, table_dif_y))))

print(f'min R = {"{:0.20f}".format(abs(max_R * f.calc_omega(x_4star, arr_x)))}')
print(f'max R = {"{:0.20f}".format(abs(min_R * f.calc_omega(x_4star, arr_x)))}')
print( 'R_min < R < R_max is ', abs(max_R * f.calc_omega(x_4star, arr_x)) < abs(fun_x(x_4star) - f.choice_Gaus_or_Newtons(x_4star, h, arr_x, table_dif_y)) < abs(min_R * f.calc_omega(x_4star, arr_x)))


# print(*arr_y, sep='\n')
# print('dif table')
# print(*table_dif_y, sep='\n \n')