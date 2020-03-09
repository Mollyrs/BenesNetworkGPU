from random import seed
from random import random
a=[]
n=0
seed(512)
print(int(random()*256), int(random()*256), int(random()*256))
while(n<256):
    if (int(random()*256) not in a):
        a.append(int(random()*256))
        n += 1
f = open("demofile3.dat", "w")
for i in a:
    print(format(i, 'x'))
    f.write(format(i, 'x'))
f.close()
#     x= rand()
#     if x not in a :

