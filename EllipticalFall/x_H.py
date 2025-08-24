import matplotlib.pyplot as plt
import csv
import numpy as np

csv_file = 'x_H.csv'
f = csv.reader(open(csv_file, 'r'), delimiter=',', doublequote=True, lineterminator='\r\n',
               quotechar='"', skipinitialspace=True)
C_D = []
C_time = []
for row in f:
    if row[0] != 'x_H':
        C_D.append(float(row[1]))
        C_time.append(float(row[0])-2.5)
plt.plot(C_time, C_D, marker='o', linestyle='-', color='b')
plt.ylabel(r'$y/H$')
plt.xlabel(r'$x/H$')
plt.title('Drag Coefficient C_D over Time')
plt.grid()
# plt.xlim(200,300)
# plt.ylim(0,2.3)
plt.savefig('C_D_plot.png')
plt.show()