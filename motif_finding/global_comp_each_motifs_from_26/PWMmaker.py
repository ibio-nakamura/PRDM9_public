#Reading PWM in test file, make it to np matrix.
import numpy as np


def PWMmaker(file_name):
    PWM = []
    with open(file_name) as f:
        lines = f.read()
        splited_list = lines.split()

    for i in range(0, len(splited_list)-3, 4):
        row = [float(splited_list[i]), float(splited_list[i+1]), float(splited_list[i+2]), float(splited_list[i+3])] 
        PWM.append(row)

    return np.array(PWM)
