import numpy as np
from collections import Counter
from atom import Atom

VERSION = 'V1.2'

class Molecule(object):
    '''Class Molecule defines a collection of atoms.'''
    
    def __init__(self, atoms = [], charge = 0, multiplicity = 1, 
        connectivity = [], gridinfo = (None, None, None), title = ''):
        '''input: atoms (Atom)
        output: Molecule'''
        
        self.setAtoms(atoms)
        self.setConnectivity(connectivity)

        self.setCharge(charge)
        self.setMultiplicity(multiplicity)
        self.setTitle(title)

        self.setGridinfo(gridinfo)
    
    def setAtoms(self, atoms):
        '''Set atoms for a molecule afterwards.'''
        
        self.atoms = []
        self.connectivity = []
        for at in atoms:
            self.addAtom(at)
    
    def addAtom(self, at):
        assert type(at) is Atom
        self.atoms.append(at)
        for conn_i in self.connectivity:
            conn_i.append(0)
        self.connectivity.append([0 for _ in self.atoms])
    
    def removeAtoms(self, symbol):
        self.atoms = [at for at in self.atoms 
            if at.symbol != symbol.capitalize()]
    
    def setCharge(self, charge):
        self.charge = int(charge)
    
    def setMultiplicity(self, multiplicity):
        self.multiplicity = int(multiplicity)
    
    def setConnectivity(self, connectivity):
        if not connectivity:
            self.connectivity = [[0 for i, ai in enumerate(self.atoms)] 
                for j, aj in enumerate(self.atoms)]
        else:
            assert len(connectivity) == len(self.atoms)
            self.connectivity = connectivity
    
    def connect(self, i, j):
        self.connectivity[i][j] = 1
        self.connectivity[j][i] = 1
    
    def setGridinfo(self, gridinfo):
        self.gridinit, self.gridsize, self.griddata = gridinfo
        if self.gridinit is None:
            self.gridinit = np.array((0, 0, 0))
        if self.gridsize is None:
            self.gridsize = np.array(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
    
    def setTitle(self, title):
        self.title = title

    def __repr__(self):
        return ''.join(['%s%d' % (k, v) for k, v in Counter([at.symbol for at in self.atoms]).items()])


