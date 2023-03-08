import numpy as np
from getopt import getopt
import sys

def create_new_poscar_with_select_atomtypes(atoms_to_keep,ifile, ofile):
    lv, coord, atomtypes, atomnums = parse_poscar(ifile)[:4]
            
    new_coord=[]
    new_atomtypes=[]
    new_atomnums=[]
    for i in range(len(atomtypes)):
        if atomtypes[i] in atoms_to_keep:
            new_atomtypes.append(atomtypes[i])
            new_atomnums.append(atomnums[i])
            for j in range(sum(atomnums[:i]),sum(atomnums[:i+1])):
                new_coord.append(coord[j])
        
    write_poscar(ofile, lv, new_coord, new_atomtypes, new_atomnums, seldyn=seldyn)
    
def parse_poscar(ifile):
    with open(ifile, 'r') as file:
        lines=file.readlines()
        sf=float(lines[1])
        latticevectors=[float(lines[i].split()[j])*sf for i in range(2,5) for j in range(3)]
        latticevectors=np.array(latticevectors).reshape(3,3)
        atomtypes=lines[5].split()
        atomnums=[int(i) for i in lines[6].split()]
        if 'Direct' in lines[7] or 'Cartesian' in lines[7]:
            start=8
            mode=lines[7].split()[0]
        else:
            mode=lines[8].split()[0]
            start=9
            seldyn=[''.join(lines[i].split()[-3:]) for i in range(start,sum(atomnums)+start)]
        coord=np.array([[float(lines[i].split()[j]) for j in range(3)] for i in range(start,sum(atomnums)+start)])
        if mode!='Cartesian':
            for i in range(sum(atomnums)):
                for j in range(3):
                    while coord[i][j]>1.0 or coord[i][j]<0.0:
                        if coord[i][j]>1.0:
                            coord[i][j]-=1.0
                        elif coord[i][j]<0.0:
                            coord[i][j]+=1.0
                coord[i]=np.dot(coord[i],latticevectors)
            
    #latticevectors formatted as a 3x3 array
    #coord holds the atomic coordinates with shape ()
    try:
        return latticevectors, coord, atomtypes, atomnums, seldyn
    except NameError:
        return latticevectors, coord, atomtypes, atomnums

def write_poscar(ofile, lv, coord, atomtypes, atomnums, **args):
    with open(ofile,'w') as file:
        if 'title' in args:
            file.write(str(args['title']))
        file.write('\n1.0\n')
        for i in range(3):
            for j in range(3):
                file.write(str('{:<018f}'.format(lv[i][j])))
                if j<2:
                    file.write('  ')
            file.write('\n')
        for i in atomtypes:
            file.write('  '+str(i))
        file.write('\n')
        for i in atomnums:
            file.write('  '+str(i))
        file.write('\n')
        if 'seldyn' in args:
            file.write('Selective Dynamics\n')
        file.write('Direct\n')
        for i in range(len(coord)):
            coord[i]=np.dot(coord[i],np.linalg.inv(lv))
        for i in range(len(coord)):
            for j in range(3):
                file.write(str('{:<018f}'.format(coord[i][j])))
                if j<2:
                    file.write('  ')
            if 'seldyn' in args:
                for j in range(3):
                    file.write('  ')
                    file.write(args['seldyn'][i][j])
            file.write('\n')
    print('new POSCAR written to: '+str(ofile))
    
if __name__ == '__main__':
    ifile='./POSCAR'
    ofile='./POSCAR_new'
    short_opts='hi:o:a:'
    long_opts=['help','input=','output=','atoms=']
    try:
        opts,args=getopt(sys.argv[1:],short_opts,long_opts)
    except getopt.GetoptError:
        print('error in command line syntax')
        sys.exit(2)
    for i,j in opts:
        if i in ['-h','--help']:
            print('''
Note: for options with multiple values, seperate the values with commas
                  
these options take a value:
    -a, --atoms                   list atoms to be kept, ie Ag,Rb
help options:
    -h, --help                    display this error message
''')
            sys.exit()
        if i in ['-i','--input']:
            ifile=j
        if i in ['-o','--output']:
            ofile=j
        if i in ['-a','--atoms']:
            atoms_to_keep=[k for k in j.split(',')]
    
    create_new_poscar_with_select_atomtypes(atoms_to_keep, ifile, ofile)