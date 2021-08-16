import os
import sys
import numpy as np
from functools import reduce
from numpy.linalg import eigh, svd, inv
from scipy.linalg import sqrtm
from indexparser import parse_index
from parse47 import parse_47
from genfchk import quicksave
from autoSVDO import parse_id

np.set_printoptions(suppress=True, precision=3)

def Lowdin_ortho(S):
    return sqrtm(inv(S))

def gen_SVDO(mfn, filename0, fnid, Stol=1e-6, Ntol=3.0e-5):
    print('Raw filenames:', mfn, filename0, fnid, sep='\n')
    fragmentation, charges, muls, labels, filenames = parse_id(fnid)
    print('Ref filenames:', filenames, sep='\n')
    f0 = parse_47(filename0, silent=True)
    ufs = dict([(filename, parse_47(filename, silent=True)) for filename in set(filenames)])
    fs = [ufs[filename] for filename in filenames]
    print('All files loaded in.')
    print()

    oidss = []
    for frag, f, m in zip(fragmentation, fs, muls):
        atids = parse_index(frag)
        if not np.allclose((np.array(f['coords']) - [f0['coords'][_-1] for _ in atids]).std(axis=0), 0, atol=2e-6):
            print(frag)
        assert np.allclose((np.array(f['coords']) - [f0['coords'][_-1] for _ in atids]).std(axis=0), 0, atol=2e-6)
        oids = np.array(reduce(lambda x,y: x+y, [[k for k, ao in enumerate(f0['basis'])
            if ao.center == at] for at in atids]))
        assert np.allclose(f0['overlap'][oids, oids[:,None]], f['overlap'], atol=Stol)
        oidss.append(oids)
        if m < 0:
            f['trafomo'] = f['trafomo'][::-1]
            f['density'] = f['density'][::-1]
            f['population'] = f['population'][::-1]
            f['fock'] = f['fock'][::-1]

    oidst = np.hstack(oidss)
    dim_mo0, dim_bs0 = f0['dim'][1:3]
    dim_mos, dim_bss = np.transpose([f['dim'][1:3] for f in fs])
    dim_mot, dim_bst = dim_mos.sum(), dim_bss.sum()

    maxspin = max(f0['dim'][0], max([f['dim'][0] for f in fs]))
    if maxspin == 1:
        pass
    elif maxspin == 2:
        print('Replicating matrices.')
        for weight, kw in zip((1., 2., 2., 1.), ('trafomo', 'population', 'density', 'fock')):
            if kw in f0 and f0[kw].shape[0] == 1:
                f0[kw] = np.stack((f0[kw][0], f0[kw][0])) / weight
            for f in fs:
                if kw in f and f[kw].shape[0] == 1:
                    f[kw] = np.stack((f[kw][0], f[kw][0])) / weight
    else:
        raise UserWarning('Wrong dimension of spin!!!')

    nmoao_to_save = np.zeros((maxspin, dim_mo0, dim_bs0))
    evalmo_to_save = np.zeros((maxspin, dim_mo0))
    nfoao_to_save = np.zeros((maxspin, dim_mo0, dim_bs0))
    evalfo_to_save = np.zeros((maxspin, dim_mo0))
    for spin in range(maxspin):
        if maxspin == 2:
            if spin == 0:
                print('Alpha space')
            else:
                print('Beta space')
        Nepfs = [f['trafomo'][spin].dot(f['density'][spin]).dot(f['trafomo'][spin].T).trace() * maxspin / 2 for f in fs]
        Neps = np.round(Nepfs).astype(int)
        assert np.allclose(Neps, Nepfs)
        Nepf0 = f0['trafomo'][spin].dot(f0['density'][spin]).dot(f0['trafomo'][spin].T).trace() * maxspin / 2
        Nep0 = int(np.round(Nepf0))
        assert np.isclose(Nep0, Nepf0)
        print('#electrons in each fragment:', Neps, '#electrons in whole molecule:', Nep0, sep='\n')
        
        S0 = f0['overlap']
        dmao0 = f0['density'][spin]
        assert np.allclose(S0, S0.T)
        faoao = np.zeros((dim_bst, dim_bs0))
        faoao[np.arange(dim_bst),oidst] = 1
        xoao = faoao
        Sfo = xoao.dot(S0).dot(xoao.T)
        
        moends = np.cumsum(Neps)
        mostarts = np.hstack(([0], moends[:-1]))
        aoends = np.cumsum(dim_bss)
        aostarts = np.hstack(([0], aoends[:-1]))
        ofoao = np.zeros((moends[-1], aoends[-1]))
        fiter = zip(fs, Neps, mostarts, moends, aostarts, aoends)
        for f, Nep, mostart, moend, aostart, aoend in fiter:
            ofoao[mostart:moend,aostart:aoend] = f['trafomo'][spin][:Nep]
        xoao = ofoao.dot(xoao)
        Sfo = xoao.dot(S0).dot(xoao.T)
    
        orthofo = Lowdin_ortho(Sfo)
        xoao = orthofo.dot(xoao)
        Sfo = xoao.dot(S0).dot(xoao.T)
    
        moao = f0['trafomo'][spin]
        Smofo = moao[:Nep0].dot(S0).dot(xoao.T)
        U, D2, V = svd(Smofo, full_matrices=True)
        print('Overlap:', D2, sep='\n')

        K = D2.shape[0]
        Norbs_removed = (D2 < 1e-4).sum()
        print('%d orbitals have singular values smaller than 1e-4 and are moved to the other group.' % Norbs_removed)
        K -= Norbs_removed
        D2 = D2[:K]
        print()
        
        print('Real orbitals:')
        nmoao = U.T.dot(moao[:Nep0])
        print('#orbitals to diagonalize by Fock matrix:', nmoao[K:].shape)
        if 'fock' in f0:
            matao = f0['fock'][spin]
        else:
            matao = f0['density'][spin]
        matnmo = nmoao[K:].dot(matao).dot(nmoao[K:].T)
        nmoao[K:] = eigh(matnmo)[1].T.dot(nmoao[K:])
        matnmo = nmoao[K:].dot(matao).dot(nmoao[K:].T)
        print('Diagonalized Fock matrix:', np.diag(matnmo), sep='\n')
        nmoao_to_save[spin][:Nep0] = nmoao
        evalmo_to_save[spin][:Nep0] = np.hstack((D2, [-1]*(Nep0-K)))
        print()
        
        print('Reference orbitals:')
        nfoao = V.dot(xoao)
        print('#orbitals to diagonalize by Fock matrix:', nfoao[K:].shape)
        if 'fock' in f0:
            matao = f0['fock'][spin]
        else:
            matao = np.zeros_like(dmao0)
            for f, oids in zip(fs, oidss):
                matao[oids,oids[:,None]] += f['population'][spin]
            matao = S0.dot(matao).dot(S0)
        matnfo = nfoao[K:].dot(matao).dot(nfoao[K:].T)
        nfoao[K:] = eigh(matnfo)[1].T.dot(nfoao[K:])
        matnfo = nfoao[K:].dot(matao).dot(nfoao[K:].T)
        print('Diagonalized Fock matrix', np.diag(matnfo), sep='\n')
        nfoao_to_save[spin][:Neps.sum()] = nfoao
        evalfo_to_save[spin][:Neps.sum()] = np.hstack((D2, [-1]*(Neps.sum()-K)))
        print()

    quicksave(mfn, nmoao_to_save, evalmo_to_save, suffix='_svdmo')
    quicksave(mfn, nfoao_to_save, evalfo_to_save, suffix='_svdfo')

def main():
    mfn, fn0, fnid = sys.argv[1:4]
    gen_SVDO(mfn, fn0, fnid)
    print()

if __name__ == '__main__':
    main()


