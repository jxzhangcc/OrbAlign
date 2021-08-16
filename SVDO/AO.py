class AtomicOrbital(object):
    def __init__(self, center, label=None, n=0, l=0, m=0, pure=True, val=None):
        self.center = center
        if label:
            self.n, self.l, self.m, self.pure = self.convert_label(label)
        else:
            self.n, self.l, self.m, self.pure = n, l, m, pure
        self.val = val

    @staticmethod
    def convert_label(label):
        if label is None:
            return (0, 0, 0, None)
        label = int(label)
        return (0, label // 100, label % 50, label // 100 <= 1 or (label % 100) // 50 >= 1)
    
    def set_n(self, n):
        self.n = n
    
    def label(self):
        return AOLABEL.nlmlabel(self.n, self.l, self.m)
    
    def __str__(self):
        return str(self.center) + '-' + AOLABEL.nlmlabel(self.n, self.l, self.m)
    
    def shift(self, delta):
        self.center += delta
        return self
    
    def __eq__(self, ao2):
        return (self.center == ao2.center and self.n == ao2.n 
            and self.l == ao2.l and self.m == ao2.m and self.pure == ao2.pure)

    def __repr__(self):
        return self.__str__()

class AOLABEL():
    Llabels = ['s', 'p', 'd', 'f', 'g', 'i']
    Mlabels = [[''], 
            ['x', 'y', 'z'], 
            ['xy', 'xz', 'yz', 'x2y2', 'z2'], 
            ['0', 'c1', 's1', 'c2', 's2', 'c3', 's3']]
    Vlabels = ['Cor', 'Val', 'Ryd']

    @classmethod
    def nlmlabel(self, n, l, m):
        return (str(n) if n else '') + AOLABEL.Llabels[l] + str(m)

    @classmethod
    def labellist(self, c=None, n=None, l=None, m=None, val=None):
        res = []
        if c is not None: res.append('%d' % c)
        if n is not None: res.append('%d' % n)
        if l is not None:
            res.append('%s' % AOLABEL.Llabels[l])
            if m is not None: res.append('%s' % AOLABEL.Mlabels[l][m-1])
        if val is not None: res.append('%s' % AOLABEL.Vlabels[val])
        return res


