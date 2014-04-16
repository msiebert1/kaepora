import glob

b = glob.glob('*')

for i in range(len(b)):
    if (b[i][-4:] == ".flm"):
        print b[i]
