import numpy as np
import sqlite3 as sq3

with sq3.connect('SNe.db') as conn:
    conn.text_factory = str
    cur = conn.cursor()
    cur.execute('SELECT DISTINCT SN from Supernovae')
    names =[r[0] for r in cur.fetchall()]

with open('snlist.dat') as f:
    lines = f.readlines()
    sndict = {}
    for line in lines:
        if not line.startswith('#'):
            data = line.split()
            if data:
                name = data[0]
                indices= [i for i,x in enumerate(data) if x.startswith(('19', '20'))]
                
                if indices:
                    index = indices[len(indices) - 1]
                    data = data[index:]
                else:
                    data = None

                if data and len(data) > 1:
                    discoverer = ' '.join(data[1:])
                else:
                    discoverer = None
                sndict[name] = discoverer

"""
discdict = {}
for name in names:
    try:
        disc = sndict[name]
    except:
        print 'No discovery data found for ' + name + '.'
        disc = None
    discdict[name] = disc
"""
#words to ask about: panstarrs, la sagra
#changed so that keywords correspond to untargeted
keywords = ['catalina', 'Drake et al', 'master', 'ptf', 'quest', 'rotse', 'palomar', 'sdss', 'sloan', 'la sagra', 'panstarrs', 'pan-starrs']
targeted = []
untargeted = []
for name in sndict:
    disc = sndict[name]
    if disc:
        if any(word in disc.lower() for word in keywords):
            untargeted.append(name)
            print name, 'was discovered independently.'
        else:
            targeted.append(name)
            print name, ' is from a targeted survey.'

for i in untargeted:
    print sndict[i]


print len(targeted)
print len(untargeted)

