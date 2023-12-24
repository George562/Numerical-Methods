import matplotlib.pyplot as plt
import os

fig = plt.figure()
ax = fig.add_subplot()
ax.set(xlim=[0, 105], ylim=[0, 3.5])

dim = [i for i in range(5, 55)]
ratio = []
cur = [0] * 10

for i in dim:
    print(" i = ", i)
    for j in range(len(cur)):
        os.system(f"gen.exe {i} {i} 500000 > LUm.txt")
        os.system(f"gen.exe 500 {i} 500000 > LUb.txt")
        os.system("run.exe < input.txt > output.txt")
        with open("output.txt", "r") as file:
            file.readline()
            cur[j] = (float(file.readline()[23:-9]) / float(file.readline()[23:-9])) ** -1
    ratio.append(sum(cur) / len(cur))

ax.plot(dim, ratio, marker='o', color='green', alpha=0.3)
plt.show()