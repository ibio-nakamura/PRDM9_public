from collections import OrderedDict
import sys


with open('multi_PWM.txt') as f:
    lines = f.readlines()


dir = './dev_PWM/'
current_filter = []
for i, line in enumerate(lines):
    if line[0]=='>':
        if i!=0:
            with open(dir+current_key+'.txt', 'w') as f:
                string = ''
                for i, row in enumerate(current_filter):
                    if i == 13:
                        string += row.rstrip('\n')
                    else:
                        string += row
                f.write(string)
            current_key = line.lstrip('>').rstrip('\n')
            current_filter = []
            continue
        else:
            current_key = line.lstrip('>').rstrip('\n')
            continue
    else:
        current_filter.append(line)

