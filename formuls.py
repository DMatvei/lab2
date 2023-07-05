import sympy as sp
import numpy as np


def choiceX0(x, arr_X):
    if abs(x - arr_X[0]) < abs(x - arr_X[-1] ):
        x0 = arr_X[0]
    else: 
        x0 = arr_X[-1]
    return x0

def choice_x0(x, arr_X):
    i = 0
    while (x > arr_X[i]):
        i += 1
        
    if i != 0:
        if (x - arr_X[i - 1]) < (arr_X[i] - x):
            x0 = arr_X[i - 1]
            i = i - 1
        else: 
            x0 = arr_X[i]
    else:
        x0 = arr_X[i]

    return x0, i


def choice_Gaus_or_Newtons(x, h, arr_x,  difTable):
    ''' принимает x, вычисляет t и по модулю выбирает, какую формулулу применять 

        параметры difTable нужнa для передачи в дальнейшую функцию для вычисления
    ''' 
    t = (x - choiceX0(x, arr_x)) / h
    
    if round(abs(t), 2) < 1:
        
        return choice_Newtons_formulas(t, difTable)
    else:
        t = (x - choice_x0(x, arr_x)[0]) / h
        print('Gauss')
        return first_second_gf(t, choice_x0(x, arr_x)[-1], difTable) 


def choice_Newtons_formulas(t, difTable):
    '''используя значение t выбирает первую или вторую формулу Ньютона 

    При 0 < t < 1 выбирает 1ую формулу Ньютона, при -1 < t < 0 выбирает 2ую формулу Ньютона
    '''
    if 0 < round(t, 1) < 1:
        print('First Newton')
        return first_nf(t, difTable)
    elif -1 < round(t, 1) < 0:
       print('Second Newton ', t)
       return  second_nf(t, difTable)
    else:
        print('t не в диапазоне (-1, 1) \ 0')
    return



def first_nf(t, difTable):
    '''1ая формула Ньютона'''
    result = 0
    w = 1
    
    for i in range(len(difTable[0])):
        result += difTable[i][0] * w
        w *= (t-i) / (i + 1)
        
    return result

def second_nf(t, difTable):
    
    result = 0
    w = 1
    
    for i in range(len(difTable[0])):
        
        result += difTable[i][-1] * w
        w *= (t + i) / (i + 1)
        
    return result



def first_second_gf(t, index, difTable):   
    result = difTable[0][index]

    count = 1
    # если условие выполнятся, то это вторая функция гауса
    if -0.5 <= round(t, 2) < 0: count = 2
    factorial = 1
    for i in range(1, len(difTable)):
        factorial /= i 
        result += difTable[i][index - (count // 2)] * factorial * sup_fun_for_gf(t, count - 1, i)
        count+= 1
        
    return result




def sup_fun_for_gf(t, count, m):
    coef = 1
    q = t + (count // 2)
    for i in range(m):
        coef *= q - i
    return coef

def dif_table(arr_y):

    dif_y_list = []
    dif_y_list.append(arr_y)

    for i in range(len(arr_y) - 1):
        current_arr = dif_y_list[i]
        add_arr = [(current_arr[j+1] - current_arr[j]) for j in range(len(current_arr) - 1)]
        dif_y_list.append(add_arr)

    return dif_y_list


def create_table(fun, arr_x):

    arr_y = [fun(arr_x[i]) for i in range(len(arr_x))]
    return arr_y
    

def calc_dif(fun, x, n): 
    diff_fun = fun

    for i in range(n):
        diff_fun = sp.diff(diff_fun, x)

    crit_point_max = sp.maximum(diff_fun, x, sp.Interval(0.5, 1.0).closure)
    crit_point_min = sp.minimum(diff_fun, x, sp.Interval(0.5, 1.0).closure)
    
    return diff_fun, crit_point_min, crit_point_max

def calc_omega(x, arr_x):
    
    w = 1
    for i in arr_x:
        w *= x - i

    return w

def factorial(n):
    result = 1
    for i in range(1, n+1):
        result *= i

    return result

def main():
    pass


if __name__ == "__main__": main()







