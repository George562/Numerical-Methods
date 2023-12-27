import matplotlib.pyplot as plt
import os

fig = plt.figure()
ax = fig.add_subplot()
ax.set(xlim=[0, 305], ylim=[0, 4.5])

dim = [i for i in range(15, 305, 15)]
ratio = []
cur = [0] * 5

for i in dim:
    print("i = ", i)
    for j in range(len(cur)):
        os.system(f"iter_gen.exe {i} 50000 > Iterm.txt")
        os.system(f"gen.exe 50 {i} 50000 > Iterb.txt")
        os.system("run.exe < input.txt > output.txt")
        with open("output.txt", "r") as file:
            file.readline()
            file.readline()
            cur[j] = (float(file.readline()[42:-9]) / float(file.readline()[23:-9])) ** -1
    ratio.append(sum(cur) / len(cur))

ax.plot(dim, ratio, marker='o', color='green', alpha=0.3)
plt.show()