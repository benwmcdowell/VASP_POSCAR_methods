from numpy import array,dot
from numpy.linalg import inv
from getopt import getopt
import sys
 
def remove_seldyn(ifile,ofile):
    try:
        lv, coord, atomtypes, atomnums, seldyn = parse_poscar(ifile)
    except ValueError:
        lv, coord, atomtypes, atomnums = parse_poscar(ifile)
    write_poscar(ofile, lv, coord, atomtypes, atomnums)

def parse_poscar(ifile):
    with open(ifile, 'r') as file:
        lines=file.readlines()
        sf=float(lines[1])
        latticevectors=[float(lines[i].split()[j])*sf for i in range(2,5) for j in range(3)]
        latticevectors=array(latticevectors).reshape(3,3)
        atomtypes=lines[5].split()
        atomnums=[int(i) for i in lines[6].split()]
        if lines[7].split()[0] == 'Direct':
            start=8
        else:
            start=9
            seldyn=[''.join(lines[i].split()[-3:]) for i in range(start,sum(atomnums)+start)]
        coord=array([[float(lines[i].split()[j]) for j in range(3)] for i in range(start,sum(atomnums)+start)])
        for i in range(sum(atomnums)):
            coord[i]=dot(latticevectors,coord[i])
            
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
        lv=inv(lv)
        for i in range(len(coord)):
            coord[i]=dot(lv,coord[i])
            for j in range(3):
                if coord[i][j]>1.0:
                    coord[i][j]-=1.0
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
    ofile='./POSCAR_seldyn'
    short_opts='hi:o:'
    long_opts=['help','input=','output=']
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
    -i, --input                    specify an input other than ./POSCAR
    -o, --output                   speicfy an output other than ./POSCAR_seldyn
help options:
    -h, --help                    display this error message
''')
            sys.exit()
        if i in ['-i','--input']:
            ifile=j
        if i in ['-o','--output']:
            ofile=j
    try:
        remove_seldyn(ifile,ofile)
    except FileNotFoundError:
        print('error reading input. exiting...')
        sys.exit(1)