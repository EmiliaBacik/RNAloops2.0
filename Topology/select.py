import sys

FS = ";"

a = {}  
b = {}  
c = {}  
d = {} 

for line in sys.stdin:
    line = line.rstrip("\n")
    fields = line.split(FS)

   
    f2 = fields[1]
    f4 = fields[3]
    f6 = fields[5]

   
    try:
        f6_val = float(f6)
    except ValueError:
        continue

    if f4 not in a or a[f4] > f6_val:
        a[f4] = f6_val
        b[f4] = line
        d[f4] = f2
    elif a[f4] == f6_val and f2 != d[f4]:
        b[f4] += "\n" + line

    c[f4] = c.get(f4, 0) + 1

for key in a:
    for row in b[key].split("\n"):
        print(f"{row};{c[key]}")
