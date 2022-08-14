import math

def factorial_f(num):
    if num == 0:
        return 1
    else:
        return num*factorial_f(num-1)
    

factorial_f(5)
math.factorial(5)==factorial_f(5)