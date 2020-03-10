from random import seed
from random import random
a=[]
n=0
seed(1024)
# print(int(random()*256), int(random()*256), int(random()*256))
while(n<128):
    if (int(abs(random()*128)) not in a):
        a.append(int(abs(random()*128)))
        n += 1
n=0
f = open("demofile3.dat", "w")
for i in a:
    print(i)
    n += 1
    if (n%10==0):
        f.write("\n")
    f.write(chr(i))
    
f.close()
#     x= rand()
#     if x not in a :

