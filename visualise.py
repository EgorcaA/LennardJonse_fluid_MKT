import matplotlib.pyplot as plt

arr = []
with open("scalars.txt") as f:
    for line in f:
        #print(float(x) for x in coloumn.split())
        arr.append(line.split())
en = []
pot = []
kin = []
for i in range(50):
    en.append(float(arr[i][1]))
    pot.append(float(arr[i][2]))
    kin.append(float(arr[i][3]))



brr = []
with open("sost.txt") as f:
    for line in f:
        #print(float(x) for x in coloumn.split())
        brr.append(line.split())
p = []
d = []
c = []
for i in range(len(brr)):
    p.append(float(brr[i][0]))
    d.append(float(brr[i][1]))
    c.append(float(brr[i][2]))

'''
err = []
with open("errs.txt") as f:
    for line in f:
        #print(float(x) for x in coloumn.split())
        err.append(float(line))
'''


ddx = []
ddv = []



fig, (ax, bx, cx) = plt.subplots(1, 3)  # Create a figure containing a single axes

ax.plot(range(50), en,'o-', label='full energy' )  # Plot some data on the axes.
ax.plot(range(50), pot,'x-', label = 'potential evergy')
ax.plot(range(50), kin,'.-', label =' kinetic energy')
ax.set_xlabel('timesteps')  # Add an x-label to the axes.
ax.set_ylabel('energy')  # Add a y-label to the axes.
ax.set_title("Energies")
ax.legend()  # Add a legend.



bx.plot(p, label='p' )  # Plot some data on the axes.

bx.set_xlabel('timesteps')  # Add an x-label to the axes.
bx.set_ylabel('')  # Add a y-label to the axes.
bx.set_title("pressure")
bx.legend()  # Add a legend.


cx.plot(d, label = '<KE/TE>')
cx.plot(c, label = '<pk/p>')
cx.set_title("div en")
cx.legend()  # Add a legend.

plt.show()
'''
    print(sum(en)/len(en), "en")
    print(sum(p[-200:])/(200), "p")
    print(sum(d[-200:])/(200), "<KE/TE>")
    print(sum(c[-200:])/(200), "pk/p")
   
'''
