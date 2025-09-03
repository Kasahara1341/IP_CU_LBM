import matplotlib.pyplot as plt
import csv
import numpy as np

C_D       =[]
C_time    =[]
C_D_ref   =[]
C_time_ref=[]

csv_file = 'x_H.csv'
f = csv.reader(open(csv_file, 'r'), delimiter=',', doublequote=True, lineterminator='\r\n',
               quotechar='"', skipinitialspace=True)
for row in f:
    if row[0] != 'x_H':
        C_D.append(float(row[1]))
        C_time.append(float(row[0]))
csv_file = 'Suzuki_1.1y.csv'
f = csv.reader(open(csv_file, 'r'), delimiter=',', doublequote=True, lineterminator='\r\n',
               quotechar='"', skipinitialspace=True)
for row in f:
    if row[0] != 'x_H':
        C_D_ref.append(float(row[1]))
        C_time_ref.append(float(row[0]))
plt.plot(C_time, C_D, marker='o', linestyle='-', color='b', label="my LBM")
plt.plot(C_time_ref, C_D_ref, color='gray', label="Suzuki & Yoshino(2018)")

plt.ylabel(r'$y/H$')
plt.xlabel(r'$x/H$')
plt.title('Drag Coefficient C_D over Time')
plt.legend()
plt.grid()
# plt.xlim(200,300)
# plt.ylim(0,2.3)
plt.savefig('C_D_plot.png')