import os
import sys
import re
import numpy as np
from collections import Counter
from molecule import Molecule
from atom import ATOMLIST
from parsegjf import parse_gjf
from indexparser import parse_index, rev_parse_index
from gengjf import gen_gjf

def parse_id(fnid):
    with open(fnid) as f:
        text = [l for l in f.read().strip().splitlines()]
        f.close()
    fragmentation, charges, muls, labels, links = [], [], [], [], []
    for l in text:
        if l.startswith('#'):
            continue
        ll = l.split()
        fragmentation.append(ll.pop(0))
        charges.append(int(ll.pop(0)) if ll else 0)
        muls.append(int(ll.pop(0)) if ll else 0)
        labels.append(ll.pop(0) if ll else None)
        links.append(ll.pop(0) if ll else None)
    return fragmentation, charges, muls, labels, links

def gen_ofns(fragments, charges, labels):
    tmp = []
    for frag, name in zip(fragments, labels):
        c = Counter([at.symbol for at in frag.atoms])
        assert (np.array(list(c.values())) > 0).all()
        elems = map(ATOMLIST.an2symbol, sorted(map(ATOMLIST.symbol2an, c.keys()))[::-1])
        if not name:
            name = ''.join([e + (str(c[e]) if c[e] > 1 else '') for e in elems])
        tmp.append(name)
    c = Counter(tmp)
    count = {}
    names = []
    for n in tmp:
        if c[n] > 1:
            count.setdefault(n, 0)
            count[n] += 1
            names.append('%s-%d' % (n, count[n]))
        else:
            names.append(n)
    return ['%s-%s' % (name, 'p' + str(c) if c >= 0 else str(c).replace('-', 'n')) 
        for name, c in zip(names, charges)]

def main():
    fn0, fnid = sys.argv[1], sys.argv[2]

    fragmentation, charges, muls, labels, links = parse_id(fnid)
    f0 = parse_gjf(fn0)
    mol = f0['molecule']

    nfrag = len(fragmentation)
    assert len(charges) == nfrag
    assert len(muls) == nfrag
    assert len(labels) == nfrag
    print(np.array((fragmentation, charges, muls, labels, links)).T)

    atids = np.array(sorted(parse_index(','.join(fragmentation))))
    if atids.shape != (len(mol.atoms),) or not np.allclose(atids, np.arange(1, len(mol.atoms)+1)):
        print('Warning: incomplete fragmentation detected!')

    fragments = []
    for f, c, m in zip(fragmentation, charges, muls):
        atids = parse_index(f)
        atoms = [mol.atoms[atid-1] for atid in atids]
        frag = Molecule(atoms, charge=c, multiplicity=abs(m))
        fragments.append(frag)
    names = gen_ofns(fragments, charges, labels)
    moltitle = os.path.splitext(os.path.split(fn0)[1])[0]
    for frag, name in zip(fragments, names):
        frag.title = '%s-%s' % (moltitle, name)

    for fragid, frag in enumerate(fragments):
        if links[fragid]:
            continue
        links[fragid] = frag.title.upper()+'.47'
        atids = parse_index(fragmentation[fragid])

        sections = f0['rest'].splitlines()
        rest = []
        for line in sections:
            if re.match(r'^ *([A-Z][a-z]? |\d+(-\d+)? )+0$', line):
                l = line.strip().split()
                newline = ' '
                assert l[-1] == '0'
                for item in l[:-1]:
                    if re.match(r'^[A-Z][a-z]?$', item):
                        newline += '-%s ' % item
                    elif re.match(r'^\d+(-\d+)?$', item):
                        basatid = [atids.index(_)+1 for _ in parse_index(item) if _ in atids]
                        if basatid:
                            newline += rev_parse_index(basatid).replace(',', ' ') + ' '
                    else:
                        raise UserWarning('Unrecognized basis section %s' % item)
                newline += '0'
                rest.append(newline)
            elif re.match(r'.*\$nbo.*\$end', line.lower()):
                rest.append(line.replace(moltitle, frag.title))
            else:
                rest.append(line)

        route = f0['route']
        route = route.replace('guess=read', '')
        route = route if 'nosymm' in route else route + ' nosymm'
        d = {'command': '%%chk=%s.chk' % frag.title,
            'route': route,
            'title': frag.title,
            'molecule': frag,
            'rest': '\n'.join(rest)}
        gen_gjf('%s.gjf' % frag.title, **d)

    with open(fnid+'.rev', 'w') as f:
        for fr, c, m, n, l in zip(fragmentation, charges, muls, labels, links):
            f.write('%s %d %d %s %s\n' % (fr, c, m, n, l))
        f.close()

if __name__ == '__main__':
    main()


