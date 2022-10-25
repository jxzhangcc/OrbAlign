import sys
import numpy as np
from functools import reduce
from numpy.linalg import eigh, svd, inv
from scipy.linalg import sqrtm
from indexparser import parse_index
from parse47 import parse_47
from genfchk import quicksave
from autofrag import parse_id

Lowdin_ortho = lambda S: sqrtm(inv(S))
array_str = lambda arr: np.array2string(arr, suppress_small=True, precision=6, max_line_width=100, 
    formatter={'str_kind': lambda s: '%9s' % s.center(7)})

class Molecule():
    def __init__(self, fnfchk=None, name=''):
        self.mfn = fnfchk
        self.name = name

    def load47(self, fn47, d47=None, silent=True, Dtol=1e-6):
        self.fn47 = fn47
        if d47 is None:
            d47 = parse_47(fn47, silent=silent)
        assert isinstance(d47, dict)

        self.atoms = np.array(d47['atoms'])
        self.coords = np.array(d47['coords'])
        self.basis = d47['basis']
        self.sao = d47['overlap']

        self.mo_ao = d47['trafomo']
        (self.nspin, self.nmo, self.nbasis) = self.mo_ao.shape
        assert self.nspin in (1, 2)
        self.dmao = d47['density']
        assert self.dmao.shape == (self.nspin, self.nbasis, self.nbasis)
        self.fmao = d47.get('fock', None)
        assert self.fmao is None or self.fmao.shape == (self.nspin, self.nbasis, self.nbasis)

        nelec = (self.mo_ao @ self.dmao @ self.mo_ao.transpose(0,2,1)).trace(axis1=1, axis2=2)
        nocc = nelec * self.nspin / 2
        self.nocc = np.round(nocc).astype(int)
        assert np.allclose(nocc, self.nocc)
        # assert HF/DFT calculations
        for spin in range(self.nspin):
            CS = self.mo_ao[spin,:self.nocc[spin]] @ self.sao
            D = CS.T @ CS
            try:
                assert np.allclose(self.dmao[spin] / 2 * self.nspin, D, atol=Dtol)
            except AssertionError:
                print(abs(self.dmao[spin] / 2 * self.nspin - D).max())

        return self

    def separate_spin(self):
        if self.mo_ao.shape[0] == 2:
            print('%s is open-shell.' % self.name)
            return
        print('%s is closed-shell. Matrices replicated.' % self.name)
        self.mo_ao = np.array((self.mo_ao[0], self.mo_ao[0]))
        self.dmao = np.array((self.dmao[0], self.dmao[0])) / 2.
        self.nocc = (self.nocc[0], self.nocc[0])
        if self.fmao is not None:
            self.fmao = np.array((self.fmao[0], self.fmao[0]))

class Fragment(Molecule):
    def __init__(self, parent, atids, oids, name=''):
        super().__init__(name=name)
        self.parent = parent
        self.atids = np.array(atids, dtype=int)
        self.oids = np.array(oids)

        self.atoms = self.parent.atoms[self.atids-1]
        self.coords = self.parent.coords[self.atids-1]
        self.basis = [_ for _ in self.parent.basis if _.center in atids]
        self.sao = self.parent.sao[oids[:,None], oids]

    def load47(self, fn47, d47=None, silent=True, Stol=1e-6, ctol=2e-6):
        atoms0 = self.atoms
        coords0 = self.coords
        basis0 = self.basis
        sao0 = self.sao

        super().load47(fn47, d47, silent)

        try:
            assert (atoms0 == self.atoms).all()
            shift = coords0.mean(axis=0) - self.coords.mean(axis=0)
            assert np.allclose(self.coords, coords0 - shift, atol=ctol)
        except AssertionError as e:
            print('Fragment atoms do not match with the complex. Please check.')
            raise e
        try:
            assert len(self.basis) == len(basis0)
            for b1, b2 in zip(basis0, self.basis):
                assert b1.center == self.atids[b2.center-1]
                assert b1.label() == b2.label()
            assert np.allclose(self.sao, sao0, atol=Stol)
        except AssertionError as e:
            print('Fragment basis does not match with the complex. Please check.')
            raise e

    def flipspin(self):
        self.mo_ao = self.mo_ao[::-1]
        self.dmao = self.dmao[::-1]
        if self.fmao is not None:
            self.fmao = self.fmao[::-1]

class FAMO(Molecule):
    def __init__(self, fnfchk, fn47, fnid):
        super().__init__(fnfchk, 'complex')
        print('Load 47 file for complex:', fn47)
        self.load47(fn47)

        # Build fragments
        print('Load fragmentation:', fnid)
        fragmentation, charges, muls, labels, filenames = parse_id(fnid)
        print('Load 47 files for fragments:', filenames)
        ufs = {}
        for filename in set(filenames):
            try:
                ufs[filename] = parse_47(filename, silent=True)
            except:
                print(filename)
                ufs[filename] = parse_47(filename, silent=True)
        print('All files loaded in.')
        print()
        self.frags = []
        for atids_str, ch, mul, label, fn in zip(fragmentation, charges, muls, labels, filenames):
            atids = parse_index(atids_str)
            oids = np.array(reduce(lambda x,y: x+y, [[k for k, ao in enumerate(self.basis)
                if ao.center == at] for at in atids]))
            frag = Fragment(self, atids, oids, label)
            try:
                frag.load47(fn, ufs[fn])
                if mul[0] < 0:
                    frag.flipspin()

                nuclear_charge = frag.atoms[:,1].sum()
                assert ch[0] == nuclear_charge - frag.nocc.sum() * 2 / frag.nspin
                assert mul[0] == (1 if frag.nocc.shape == (1,) else frag.nocc @ (1, -1) + 1)

                if len(ch) > 1 or len(mul) > 1:
                    print('Cautious! Fragment density of %s (%s) adapted to charge/multiplicity.' % (label, fn))
                    na = (nuclear_charge - ch[-1] + (abs(mul[-1]) - 1) * np.sign(mul[-1])) / 2
                    nb = (nuclear_charge - ch[-1] - (abs(mul[-1]) - 1) * np.sign(mul[-1])) / 2
                    assert np.allclose((na, nb), np.array((na, nb)).round().astype(int))
                    na, nb = np.array((na, nb)).round().astype(int)
                    nocc = (na,) if (na == nb and frag.nspin == 1) else (na, nb)
                    if self.nspin == 2 and frag.nspin == 1:
                        frag.separate_spin()
                    if self.nspin == 1 and frag.nspin == 2:
                        self.separate_spin()
                    for spin in range(frag.nspin):
                        # # Choose orbitals by FMO order
                        # CS = frag.mo_ao[spin,:frag.nocc[spin]] @ frag.sao
                        # frag.dmao[spin] = CS.T @ CS

                        # Choose orbitals by FMO overlap with MO
                        if nocc[spin] < frag.nocc[spin]:
                            Socc = frag.mo_ao[spin,:frag.nocc[spin]] @ self.sao[frag.oids] @ \
                                self.mo_ao[spin,:self.nocc[spin]].T
                            U, D, V = np.linalg.svd(Socc)
                            frag.mo_ao[spin,:frag.nocc[spin]] = U.T @ frag.mo_ao[spin,:frag.nocc[spin]]
                            print('%d orbitals removed from fragment density with max overlap being %.2f' % 
                                (frag.nocc[spin] - nocc[spin], D[nocc[spin]:].max()))
                            frag.nocc[spin] = nocc[spin]
                            CS = frag.mo_ao[spin,:frag.nocc[spin]] @ frag.sao
                            frag.dmao[spin] = CS.T @ CS
                        elif nocc[spin] > frag.nocc[spin]:
                            raise UserWarning('Increasing fragment density is not suggested.')
                            # Svir = frag.mo_ao[spin,:frag.nocc[spin]] @ CS.T
                            # U, D, V = np.linalg.svd(Svir)
                            # print(D)
                            # frag.nocc[spin] = nocc[spin]
                        else:
                            pass

            except Exception as e:
                print(atids_str, ch, mul, label, fn)
                raise e
            self.frags.append(frag)
        print('Fragments prepared')
        print()

    def match_spin(self):
        # Replicate matrices if spin mismatch
        maxspin = max(self.nspin, max(frag.nspin for frag in self.frags))
        if maxspin == 2:
            self.separate_spin()
            for frag in self.frags:
                frag.separate_spin()
            print()

    def check_complete_fragmentation(self):
        existing_atids = np.hstack([frag.atids for frag in self.frags])
        return (np.arange(len(self.atoms)) + 1 == existing_atids).all()

    def complete_fragmentation(self):
        existing_atids = np.hstack([frag.atids for frag in self.frags])
        rest_atids = [_ for _ in range(1, len(self.atoms) + 1) if _ not in existing_atids]
        if not rest_atids:
            return
        rest_oids = np.hstack([[k for k, ao in enumerate(self.basis) if ao.center == at] for at in rest_atids])
        rest = Fragment(self, rest_atids, rest_oids, 'rest')
        self.frags.append(rest)

    def align(self):
        # original SVDO
        print('Computing Fragment-Aligned Molecular Orbital (FAMO) analysis')

        # Order AO by fragment
        aoends = np.cumsum([frag.nbasis for frag in self.frags])
        aostarts = np.hstack(([0], aoends[:-1]))
        tfnbasis = aoends[-1]
        fao_ao = np.zeros((tfnbasis, self.nbasis))
        oidss = [frag.oids for frag in self.frags]
        fao_ao[np.arange(tfnbasis), np.hstack(oidss)] = 1

        self.spcs = np.zeros((2, 4), dtype=int)
        self.famo_ao, self.famo_evals, self.famo_labels = [], [], []
        for spin in range(self.nspin):
            if self.nspin == 2:
                print('---' + 'Alpha Beta'.split()[spin] + ' space---')
                print('NO. of electrons in the complex:', self.nocc[spin])
                print('NO. of electrons in fragments:', [frag.nocc[spin] for frag in self.frags])
            else:
                print('NO. of electron pairs in the complex:', self.nocc[spin])
                print('NO. of electron pairs in fragments:', [frag.nocc[spin] for frag in self.frags])

            # Build fragment occupied MOs (FOMO)
            moends = np.cumsum([frag.nocc[spin] for frag in self.frags])
            mostarts = np.hstack(([0], moends[:-1]))
            tfnocc = moends[-1]
            fomo_fao = np.zeros((tfnocc, tfnbasis))
            for frag, mostart, moend, aostart, aoend in zip(self.frags, mostarts, moends, aostarts, aoends):
                fomo_fao[mostart:moend,aostart:aoend] = frag.mo_ao[spin][:frag.nocc[spin]]
            fomo_ao = fomo_fao @ fao_ao

            # Orthogonalize FOMOs
            sfomo = fomo_ao @ self.sao @ fomo_ao.T
            ofomo_fomo = Lowdin_ortho(sfomo)
            ofomo_ao = ofomo_fomo @ fomo_ao
            # # Without orthogonalzation of FOMOs, lead to identical results if aligned to FMO later
            # ofomo_ao = fomo_ao

            # Align
            cross_overlap = self.mo_ao[spin,:self.nocc[spin]] @ self.sao @ ofomo_ao.T
            U, D2, V = svd(cross_overlap, full_matrices=True)
            # print('Singular values:')
            # print(array_str(D2))

            ni = (D2 > (1 - 1e-2)).sum()
            nt = (D2 < 1e-4).sum()
            nd = D2.shape[0] - ni - nt
            nu = self.nocc[spin] - tfnocc
            print('Inactive    space:', ni)
            print('Deformed    space:', nd)
            print('Transferred space:', nt)
            print('Unique      space:', nu)

            self.spcs[spin] = ni, nd, nt, nu
            self.famo_ao.append(U.T @ self.mo_ao[spin,:self.nocc[spin]])
            self.famo_evals.append(np.hstack((D2, [-1] * (self.nocc[spin] - D2.shape[0]))))
            self.famo_labels.append(np.array(['I'] * ni + ['D'] * nd + ['T'] * nt + ['U'] * nu, dtype='U8'))

            self.diagonalize('TU', spin=spin)
            # self.diagonalize('I', spin=spin)
            self.align_to_FMO(spin=spin)
        print()

    def save(self, suffix='_famo'):
        famo_ao_to_save = np.zeros((self.nspin, self.nmo, self.nbasis))
        famo_evals_to_save = np.zeros((self.nspin, self.nmo))
        for spin in range(self.nspin):
            famo_ao_to_save[spin,:self.nocc[spin]] = self.famo_ao[spin]
            famo_evals_to_save[spin,:self.nocc[spin]] = self.famo_evals[spin]
        ofn = quicksave(self.mfn, famo_ao_to_save, famo_evals_to_save, suffix=suffix, overwrite=True)
        if ofn:
            print('FAMO saved to file %s' % ofn)
        print()

    def diagonalize(self, spc='TU', spin=0):
        # Diagonalize specific space by Fock matrix
        # Transferred+Unique space is diagonalized by default

        assert self.fmao is not None
        spc = spc.upper()

        # if self.nspin == 2:
        #     print('---' + 'Alpha Beta'.split()[spin] + ' space---')
        ni, nd, nt, nu = self.spcs[spin]
        if spc == 'TU':
            spcmin = ni + nd
            spcmax = self.nocc[spin]
        elif spc == 'I':
            spcmin = 0
            spcmax = ni
        else:
            raise UserWarning('Unknown space: %s' % spc)
        fmfamo = self.famo_ao[spin][spcmin:spcmax] @ self.fmao[spin] @ self.famo_ao[spin][spcmin:spcmax].T
        evals, evecs = eigh(fmfamo)
        self.famo_ao[spin][spcmin:spcmax] = evecs.T @ self.famo_ao[spin][spcmin:spcmax]
        self.famo_evals[spin][spcmin:spcmax] = evals
        self.famo_labels[spin][spcmin:spcmax] = 'TU' if nt > 0 else 'U'
        # print('%s space diagonalized by Fock matrix. MO energies:' % '+'.join(spc))
        # print(array_str(evals))
        print('%s space diagonalized by Fock matrix.' % '+'.join(spc))

    def align_to_FMO(self, spc='ID', spin=0):
        spc = spc.upper()
        assert spc == 'ID'
        ni, nd, nt, nu = self.spcs.T
        assert (nt == 0).all()
        assert (nu >= 0).all()
        assert self.fmao is not None

        aoends = np.cumsum([frag.nbasis for frag in self.frags])
        aostarts = np.hstack(([0], aoends[:-1]))
        tfnbasis = aoends[-1]
        fao_ao = np.zeros((tfnbasis, self.nbasis))
        oidss = [frag.oids for frag in self.frags]
        fao_ao[np.arange(tfnbasis), np.hstack(oidss)] = 1

        # if self.nspin == 2:
        #     print('---' + 'Alpha Beta'.split()[spin] + ' space---')
        moends = np.cumsum([frag.nocc[spin] for frag in self.frags])
        mostarts = np.hstack(([0], moends[:-1]))
        tfnocc = moends[-1]
        fomo_fao = np.zeros((tfnocc, tfnbasis))
        for frag, mostart, moend, aostart, aoend in zip(self.frags, mostarts, moends, aostarts, aoends):
            fomo_fao[mostart:moend,aostart:aoend] = frag.mo_ao[spin][:frag.nocc[spin]]
        fomo_ao = fomo_fao @ fao_ao

        for frag, mostart, moend, aostart, aoend in zip(self.frags, mostarts, moends, aostarts, aoends):
            # # Canonicalize FMO by overlap with complex MO
            # # Problem: might mix core orbitals which have close singular values (~1)
            # cross_overlap = fomo_ao[mostart:moend] @ self.sao @ self.famo_ao[spin].T
            # ...

            # Canonicalize FMO by energy
            fmfomo = fomo_ao[mostart:moend] @ self.fmao[spin] @ fomo_ao[mostart:moend].T
            evals, evecs = eigh(fmfomo)
            fomo_ao[mostart:moend] = evecs.T @ fomo_ao[mostart:moend]

        K = ni[spin] + nd[spin]
        cross_overlap = self.famo_ao[spin][:K] @ self.sao @ fomo_ao.T
        U, D2, V = svd(cross_overlap)
        self.famo_ao[spin][:K] = V.T @ U.T @ self.famo_ao[spin][:K]
        overlaps = np.diag(V.T @ U.T @ cross_overlap)
        self.famo_evals[spin][:K] = overlaps
        for kfrag, (frag, mostart, moend) in enumerate(zip(self.frags, mostarts, moends)):
            self.famo_labels[spin][mostart:moend] = 'F%d' % (kfrag + 1)
        # print('%s space aligned to FMOs. Overlaps:' % '+'.join(spc))
        # print(array_str(overlaps))
        print('%s space aligned to FMOs.' % '+'.join(spc))

    def PA(self, method='Mulliken', basis='AO'):
        try:
            self.famo_ao
            self.famo_evals
            self.famo_labels
        except AttributeError as e:
            print('FAMO analysis not properly done.')
            raise e

        assert basis in ('AO', 'FMO')
        assert method in ('Mulliken', 'Lowdin', 'NPA')
        method = 'Natural' if method == 'NPA' else method
        print('%s population analysis of FAMO in the basis of %s' % (method, basis))

        for spin in range(self.nspin):
            if self.nspin == 2:
                print('---' + 'Alpha Beta'.split()[spin] + ' space---')

            if basis == 'AO':
                C = self.famo_ao[spin]
                S = self.sao
            elif basis == 'FMO':
                C = self.famo_ao[spin]
                fmo_ao = np.eye(self.nbasis)
                for frag in self.frags:
                    try:
                        fmo_ao[frag.oids[:,None],frag.oids] = frag.mo_ao[spin]
                    except AttributeError:
                        fmo_ao[frag.oids[:,None],frag.oids] = sqrtm(inv(frag.sao))
                C = C @ inv(fmo_ao)
                S = fmo_ao @ self.sao @ fmo_ao.T
            else:
                raise UserWarning('Basis not implemented yet.')

            # P = C[:,:,None] @ C[:,None,:]

            if method == 'Mulliken':
                # pop = np.diagonal(P @ S, axis1=1, axis2=2)
                pop = C * (C @ S)
            elif method == 'Lowdin':
                oao = sqrtm(S)
                assert np.allclose(oao.imag, 0)
                oao = oao.real
                # pop = np.diagonal(oao @ P @ oao, axis1=1, axis2=2)
                pop = (C @ oao) ** 2
            elif method == 'Natural' and basis == 'AO':
                raise UserWarning('Method not implemented yet.')
            else:
                raise UserWarning('Method not implemented yet.')

            pop_frag = np.array([pop[:,frag.oids].sum(axis=1) for frag in self.frags]).T
            tpop_frag = pop_frag.sum(axis=0) * 2 / self.nspin
            nuc_frag = np.array([frag.atoms[:,1].sum() for frag in self.frags]) / self.nspin
            ch_frag = nuc_frag - tpop_frag

            print('Fragment contribution to each FAMO:')
            print('Eigenvalue'.rjust(24), ''.join('%6d %-5s' % (kfrag + 1, frag.name) 
                for kfrag, frag in enumerate(self.frags)), sep='')
            for oid in range(self.nocc[spin]):
                print('%6d%6s%12.6f' % (oid+1, self.famo_labels[spin][oid], self.famo_evals[spin][oid]), 
                    ''.join('%12.6f' % _ for _ in pop_frag[oid]), sep='')
            uspc = self.famo_labels[spin] == 'U'
            print('U_total'.center(12), ' ' * 12, ''.join('%12.6f' % _ for _ in pop_frag[uspc].sum(axis=0)), sep='')
            print('Total'.center(12), '%12.6f' % tpop_frag.sum(), ''.join('%12.6f' % _ for _ in tpop_frag), sep='')
            print('Nuclear'.center(12), '%12.6f' % nuc_frag.sum(), ''.join('%12.6f' % _ for _ in nuc_frag), sep='')
            print('Charge'.center(12), '%12.6f' % ch_frag.sum(), ''.join('%12.6f' % _ for _ in ch_frag), sep='')
            print()

def main():
    fnfchk, fn47, fnid = sys.argv[1:4]
    famo = FAMO(fnfchk, fn47, fnid)
    famo.match_spin()
    famo.align()
    famo.save()
    famo.complete_fragmentation()
    famo.PA('Mulliken', 'AO')
    # famo.PA('Mulliken', 'FMO')  # same as Mulliken + AO
    # famo.PA('Lowdin', 'AO')
    # famo.PA('Lowdin', 'FMO')
    # famo.PA('NPA', 'AO')

if __name__ == '__main__':
    main()


