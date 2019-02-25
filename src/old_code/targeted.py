import numpy as np
import sqlite3 as sq3

def old_get_names():
    with sq3.connect('SNe.db') as conn:
        conn.text_factory = str
        cur = conn.cursor()
        cur.execute('SELECT DISTINCT SN from Supernovae')
        names =[r[0] for r in cur.fetchall()]
    return names

def disc_dict():
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
        return sndict

def targ_lists():
    #names = old_get_names()
    sndict = disc_dict()
    #words to ask about: panstarrs, la sagra, EROS, deep lens survey team, Neat/woods-vassey et al
    #changed so that keywords correspond to untargeted
    keywords = ['catalina', 'Drake et al', 'master', 'ptf', 'quest', 'rotse', 'palomar', 'la sagra', 'panstarrs', 'pan-starrs', 'sloan', 'sdss', 'EROS']
    targwords =['loss', 'supernova', 'lotoss' ]
    targeted = []
    untargeted = []
    unknown = []
    for name in sndict:
        disc = sndict[name]
        if disc:
            if any(word in disc.lower() for word in keywords):
                untargeted.append(name)
            elif any(word in disc.lower() for word in targwords):
                targeted.append(name)
            else:
                unknown.append(name)
    """
    for item in targeted:
        if not any(word in sndict[item].lower() for word in targwords):
            print item, sndict[item]

    for item in untargeted:
        print item
    """
    return targeted, untargeted, unknown
