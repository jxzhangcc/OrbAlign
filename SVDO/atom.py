from constants import ATOMLIST

class Atom(object):
    '''Class Atom collects the atomic number of an atom,
    as well as its coordinates for general use.'''
    
    def __init__(self, an_as, coordinates, layer = 0):
        '''input: an (int) or symbol (str), coordinates (float * 3)'''
        
        try:
            self.an = int(an_as)
            self.symbol = ATOMLIST.an2symbol(self.an)
        except ValueError:
            self.symbol = str(an_as)
            self.an = ATOMLIST.symbol2an(self.symbol)
        
        assert len(coordinates) == 3
        self.coords = tuple(coordinates)
        self.layer = layer
    
    def setLayer(self, layer):
        self.layer = layer
    
    def output_form(self, layer=False):
        return ' %-15s' % self.symbol + \
            '%14.8f%14.8f%14.8f' % tuple(self.coords) + \
            (' %s' % {0:'H', 1:'M', 2:'L'}[self.layer] if layer else '')
    
    def __eq__(self, atom):
        return (self.an == atom.an) and (self.coords == atom.coords)

    def __str__(self):
        return '%-2s' % self.symbol + ' %5.2f %5.2f %5.2f' % self.coords
                                              
    def __repr__(self):                       
        return '%-2s' % self.symbol + ' %5.2f %5.2f %5.2f' % self.coords

