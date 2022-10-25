import os
import re
from chemistry.atom import Atom
from chemistry.molecule import Molecule

def parse_gjf(filename):
    if os.path.splitext(filename)[1] not in ('.gjf', '.com'):
        raise Exception('%s is not a gjf file.' % filename)
    
    file = open(filename)
    newline = lambda: file.readline().strip()
    
    commands = []
    line = newline()
    while line.startswith('%'):
        commands.append(line)
        line = newline()
    
    if not line.startswith('#'):
        raise Exception('Invalid gjf file.')
    
    route = ''
    while line:
        if route:
            route = route + ' ' + line
        else:
            route = line
        line = newline()
    
    title = []
    line = newline()
    while line:
        title.append(line)
        line = newline()
    
    charge, multiplicity = [int(_) for _ in newline().split()]
    
    atoms = []
    layers = {'H':0, 'M':1, 'L':2}
    line = newline()

    poscoords = min([i for i, seg in enumerate(line.split()) 
        if '.' in seg])
    while line:
        l = line.split()
        symbol = l[0]
        coords = [float(_) for _ in l[poscoords:poscoords+3]]
        try:
            layer = layers[l[poscoords+3]]
        except IndexError:
            layer = 0
        atoms.append(Atom(symbol, coords, layer))
        line = newline()
    
    if 'connectivity' in route:
        line = newline()
        # if not line:
            # line = newline()
        length = len(atoms)
        connectivity = [[0 for i in range(length)] for j in range(length)]
        for i in range(length):
            if not re.match(r'\d+( \d+ \d+\.\d+)*', line) or \
                int(line.split()[0]) != i + 1:
                raise Exception('Invalid connectivity.')
            l = line.split()
            for atomid, bondorder in zip(l[1::2], l[2::2]):
                connectivity[int(l[0]) - 1][int(atomid) - 1] = float(bondorder)
                connectivity[int(atomid) - 1][int(l[0]) - 1] = float(bondorder)
            line = newline()
    else:
        connectivity = []
    
    rest = file.read()
    
    file.close()
    
    # return commands, route, title, mol, rest
    return {'molecule': Molecule(atoms, charge, multiplicity, connectivity), 
        'command': '\n'.join(commands), 'route': route, 
        'title': '\n'.join(title), 'rest': rest}

if __name__ == '__main__':
    import sys
    d = parse_gjf(sys.argv[1])
    m = d['molecule']
    print(m.charge, m.multiplicity)
    # for atom in m.atoms:
        # print(atom.coords)
    # for line in m.connectivity:
        # print(' '.join([str(_) for _ in line]))
    print(d['command'])
    print(d['route'])
    print(d['title'])
    print(len(d['molecule'].atoms))
    print(d['rest'])

