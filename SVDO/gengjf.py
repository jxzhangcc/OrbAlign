import os

def gen_gjf(filename, molecule, command='', route='#p', title='', rest='', 
    overwrite = False, noconnect=True, layer=False):
    if not overwrite and os.path.exists(filename):
        overwriteflag = input('Overwrite file %s? Y/N ' % filename)
        if overwriteflag.lower() != 'y':
            print(filename, 'skipped.')
            return
    
    file = open(filename, 'w')
    if command:
        file.write(command+'\n')
    if not noconnect:
        if 'geom=connectivity' in route:
            file.write(route+'\n\n')
        else:
            file.write(route+' geom=connectivity'+'\n\n')
    else:
        if 'geom=connectivity' in route:
            file.write(route.replace('geom=connectivity', '')+'\n\n')
        else:
            file.write(route+'\n\n')
    if not title:
        if molecule.title:
            title = molecule.title
        else:
            title = os.path.splitext(os.path.basename(filename))[0]
    file.write('%s\n\n' % title)
    file.write('%d %d\n' % (molecule.charge, molecule.multiplicity))
    
    for atom in molecule.atoms:
        file.write(atom.output_form(layer=layer) + '\n')
    file.write('\n')
    
    if not noconnect:
        for i, ai in enumerate(molecule.atoms):
            file.write(' %d' % (i+1))
            for j, aj in enumerate(molecule.atoms):
                if j > i and molecule.connectivity[i][j]:
                    file.write(' %d %.1f' % (j+1, molecule.connectivity[i][j]))
            file.write('\n')
        file.write('\n')
    
    if rest:
        file.write(rest+'\n')
    
    file.write('\n')
    file.close()


