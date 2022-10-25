import re
import numpy as np
from AO import AtomicOrbital
from collections import OrderedDict

def fromstring(string):
    pat = re.compile(r'( +-?0\.[0-9]+)(-[1-9][0-9]{2})')
    string = pat.sub(r'\1E\2', string)
    return np.fromstring(string, sep=' ')
    # return np.fromiter(string.split(), dtype='float')

def tri2square(tri, dim = None):
    assert tri.ndim == 1
    if not dim:
        dim = int((tri.shape[0]*2) ** 0.5)
    try:
        assert tri.shape == (dim*(dim+1)/2,)
    except AssertionError:
        print(tri.shape, dim)
        assert tri.shape == (dim*(dim+1)/2,)
    square = np.zeros((dim, dim))
    id = np.tril_indices(dim)
    square[id] = tri
    square[id[::-1]] = tri
    return square

def parse_mat(text, dim_bs):
    data = fromstring(text)
    # data.shape = (-1, dim_bs)
    print(data.shape)

def parse_47(filename, silent=False, h11flag=False):
    f = open(filename)
    text = f.read()
    f.close()
    
    d = dict()
    
    sections = [s.strip() for s in text.split('$END')]
    assert sections[-1] == ''
    for section in sections[:-1]:
        assert section.startswith('$')
        try:
            module, data_string = section.split(None, 1)
        except ValueError:
            module = section.split(None, 1)[0]
            data_string = ''
        d[module] = data_string
    # print(d.keys())
    # GENNBO, NBO, COORD, BASIS, CONTRACT
    # OVERLAP, DENSITY, FOCK, LCAOMO, DIPOLE
    
    data = dict()
    
    # INFO
        # Read $NBO
        # Read $CONTRACT
    
    # $GENNBO
    keywords = [kwd.strip().upper() for kwd in d['$GENNBO'].split()]
    assert 'BODM' in keywords
    openshell = 'OPEN' in keywords
    if openshell:
        dim_spin = 2
    else:
        dim_spin = 1
    
    # $COORD
        # Read NO. ecp electrons
    lines = [line.strip() for line in d['$COORD'].splitlines()]
    title = lines.pop(0)
    ddfff = lambda l: (int(l[0]), int(l[1]), 
        float(l[2]), float(l[3]), float(l[4]))
    atoms = [ddfff(line.split()) for line in lines]
    N_atoms = len(atoms)
    data['atoms'] = [atom[0:2] for atom in atoms]
    data['coords'] = [atom[2:5] for atom in atoms]
    
    # $BASIS
        # Read centers & labels
    c, l = re.match(r'CENTER =((?:\s+\d+)+)\s+LABEL =((?:\s+\d+)+)', 
        d['$BASIS']).groups()
    centers = [int(_) for _ in c.strip().split()]
    labels = [int(_) for _ in l.strip().split()]
    dim_bs = len(centers)
    assert len(labels) == dim_bs
    data['basis'] = [AtomicOrbital(center, label) 
        for center, label in zip(centers, labels)]
    if not silent:
        print('Basic info loaded.')
    
    # MAT
        # Read $OVERLAP
        # Read $DENSITY
        # Read $FOCK
        # Read $LCAOMO
        # Read $DIPOLE
        # Determine dimension
    
    data['overlap'] = tri2square(fromstring(d['$OVERLAP']), dim_bs)
    if not silent:
        print('Overlap matrix loaded.')
    
    data['trafomo'] = fromstring(d['$LCAOMO']).reshape(dim_spin, -1, dim_bs)
    data['dim'] = data['trafomo'].shape
    dim = data['dim'][1]
    if not silent:
        print('LCAOMO matrix loaded.')
    
    # data['population'] = np.array([tri2square(u, dim) for u in fromstring(d['$DENSITY']).reshape(dim_spin, -1)])
    data['population'] = np.array([tri2square(u, dim_bs) for u in fromstring(d['$DENSITY']).reshape(dim_spin, -1)])
    # assert data['population'].shape == (dim_spin, dim_bs, dim)
    assert data['population'].shape == (dim_spin, dim_bs, dim_bs)
    if not silent:
        print('Population matrix loaded.')
    
    # data['density'] = np.array([data['overlap'].dot(pop).dot(data['overlap']) for pop in data['population']])
    data['density'] = data['overlap'] @ data['population'] @ data['overlap']
    assert data['density'].shape == (dim_spin, dim_bs, dim_bs)
    if not silent:
        print('Density matrix calculated.')
    
    if '$FOCK' in d:
        # data['fock'] = tri2square(fromstring(d['$FOCK']), dim_bs)
        data['fock'] = np.array([tri2square(u, dim) for u in fromstring(d['$FOCK']).reshape(dim_spin, -1)])
        assert data['fock'].shape == (dim_spin, dim_bs, dim_bs)
        if not silent:
            print('Fock matrix loaded.')
    if '$DIPOLE' in d:
        data['dipole'] = np.array([tri2square(u, dim_bs) for u in 
            fromstring(d['$DIPOLE']).reshape(3, -1)])
        assert data['dipole'].shape == (3, dim_bs, dim_bs)
        if not silent:
            print('Dipole matrix loaded.')
    
    if '$NUCLEAR' in d:
        data['nuclear'] = tri2square(fromstring(d['$NUCLEAR']), dim_bs)
        if not silent:
            print('Nuclear matrix loaded.')
    
    if '$KINETIC' in d:
        data['kinetic'] = tri2square(fromstring(d['$KINETIC']), dim_bs)
        if not silent:
            print('Kinetic matrix loaded.')
    
    if '$C10' in d:
        data['C10'] = fromstring(d['$C10']).reshape(3, -1, dim_bs)
        # data['C10'] = fromstring(d['$C10']).reshape(3, dim_bs, -1).transpose(0,2,1)
        if not silent:
            print('C10 matrix loaded.')
    
    if '$H01' in d:
        data['H01'] = np.array([tri2square(u, dim_bs) for u in 
            fromstring(d['$H01']).reshape(N_atoms*3, -1)])\
            .reshape(N_atoms, 3, dim_bs, dim_bs)
        if not silent:
            print('H01 matrix loaded.')
    
    if '$H11' in d and h11flag == True:
        data['H11'] = np.array([tri2square(u, dim_bs) for u in 
            fromstring(d['$H11']).reshape(N_atoms*3*3, -1)])\
            .reshape(N_atoms, 3, 3, dim_bs, dim_bs)
        if not silent:
            print('H11 matrix loaded.')
    
    for kw in d.keys():
        if kw in ('$GENNBO', '$NBO', '$COORD', '$BASIS', '$CONTRACT', 
            '$OVERLAP', '$DENSITY', '$FOCK', '$LCAOMO', '$DIPOLE', 
            '$NUCLEAR', '$KINETIC', 
            '$C10', '$H01', '$H11'):
            pass
        else:
            print(kw, )
            parse_mat(d[kw], dim_bs)
    
    return data

def collect_basis(basis):
    atoms = OrderedDict()
    for oid, AO in enumerate(basis):
        atoms.setdefault(AO.center, OrderedDict())
        atoms[AO.center].setdefault(AO.l, [])
        atoms[AO.center][AO.l].append(oid)
    return atoms

def main():
    import sys
    data = parse_47(sys.argv[1])
    np.set_printoptions(suppress=True)
    for k, v in data.items():
        try:
            v.shape
            print(k, v.shape, v.dtype)
        except:
            print(k, v)
    print(collect_basis(data['basis']))

if __name__ == '__main__':
    main()

