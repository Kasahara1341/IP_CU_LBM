import matplotlib.pyplot as plt
import csv
import numpy as np

csv_file = 'C_D.csv'
f = csv.reader(open(csv_file, 'r'), delimiter=',', doublequote=True, lineterminator='\r\n',
               quotechar='"', skipinitialspace=True)
C_D = []
C_time = []
for row in f:
    if row[0] != 'C_D':
        C_D.append(float(row[0]))
        C_time.append(float(row[1]))
plt.plot(C_time, C_D, marker='o', linestyle='-', color='b')
plt.xlabel('Time (s)')
plt.ylabel('C_D')
plt.title('Drag Coefficient C_D over Time')
plt.grid()
# plt.xlim(300,350)
# plt.ylim(1.5,1.55)
# plt.xlim(150,400)
# plt.ylim(1.2,1.6)
plt.savefig('C_D_plot.png')