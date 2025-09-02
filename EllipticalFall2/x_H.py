import matplotlib.pyplot as plt
import csv
import numpy as np

def make_csv(csv_file):
    f = csv.reader(open(csv_file+".csv", 'r'), delimiter=',', doublequote=True, lineterminator='\r\n',
                   quotechar='"', skipinitialspace=True)
    C_D = []
    C_time = []
    for row in f:
        if row[0] != 'x_H':
            C_D.append(float(row[1]))
            C_time.append(float(row[0]))
    plt.figure()
    plt.plot(C_time, C_D, marker='o', linestyle='-', color='b')
    plt.ylabel(r'$y/H$')
    plt.xlabel(r'$x/H$')
    # plt.title('Drag Coefficient C_D over Time')
    plt.grid()
    # plt.xlim(200,300)
    # plt.ylim(0,2.3)
    plt.savefig(csv_file+'.png')

make_csv("x_H")
make_csv("x_theta")