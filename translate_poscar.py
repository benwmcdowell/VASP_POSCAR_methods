from numpy import array,dot
from numpy.linalg import inv
from getopt import getopt
import sys

def translate_poscar(ifile,shift,**args):
    try:
        lv, coord, atomtypes, atomnums, seldyn = parse_poscar(ifile)
    except ValueError:
        lv, coord, atomtypes, atomnums = parse_poscar(ifile)
        seldyn='none'
    
    shift=array(shift)
    if 'direct' in args:
        shift=dot(shift,lv)
    for i in range(len(coord)):
        coord[i]+=shift
    
    if seldyn=='none':
        write_poscar(ifile,lv,coord,atomtypes,atomnums)
    else:
        write_poscar(ifile,lv,coord,atomtypes,atomnums,seldyn=seldyn)
    
    
def parse_poscar(ifile):
    with open(ifile, 'r') as file:
        lines=file.readlines()
        sf=float(lines[1])
        latticevectors=[float(lines[i].split()[j])*sf for i in range(2,5) for j in range(3)]
        latticevectors=array(latticevectors).reshape(3,3)
        atomtypes=lines[5].split()
        atomnums=[int(i) for i in lines[6].split()]
        if 'Direct' in lines[7] or 'Cartesian' in lines[7]:
            start=8
            mode=lines[7].split()[0]
        else:
            mode=lines[8].split()[0]
            start=9
            seldyn=[''.join(lines[i].split()[-3:]) for i in range(start,sum(atomnums)+start)]
        coord=array([[float(lines[i].split()[j]) for j in range(3)] for i in range(start,sum(atomnums)+start)])
        if mode!='Cartesian':
            for i in range(sum(atomnums)):
                for j in range(3):
                    while coord[i][j]>1.0 or coord[i][j]<0.0:
                        if coord[i][j]>1.0:
                            coord[i][j]-=1.0
                        elif coord[i][j]<0.0:
                            coord[i][j]+=1.0
                coord[i]=dot(coord[i],latticevectors)
            
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
            coord[i]=dot(coord[i],inv(lv))
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
    
if __name__ == '__main__':
    short_opts='hd'
    long_opts=['help','direct']
    direct=False
    try:
        ifile=sys.argv[1]
        shift=sys.argv[2].split(',')
        for i in range(len(shift)):
            shift[i]=float(shift[i])
    except IndexError:
        print('missing required arguments. exiting...')
        sys.exit()
    try:
        opts,args=getopt(sys.argv[3:],short_opts,long_opts)
    except IndexError:
        print('error specifying optional arguments')
        sys.exit()
    for i,j in opts:
        if i in ['-d','--direct']:
            direct=True
        if i in ['-h','--help']:
            print('''
input options:
    -d, --direct                  specify the shift in direct coordinates. default is Cartesian
help options:
    -h, --help                    the first argument is the POSCAR file to translate
                                  the second argument is shift to apply in Cartesian coordinates: x,y,z
''')
            sys.exit()
    translate_poscar(ifile,shift,direct=direct)
