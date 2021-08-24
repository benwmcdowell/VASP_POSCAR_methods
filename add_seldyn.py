from numpy import array,dot
from numpy.linalg import inv,norm
from getopt import getopt
import sys

def add_seldyn(mode, ifile, ofile,reverse, **args):
    try:
        lv, coord, atomtypes, atomnums, seldyn = parse_poscar(ifile)
    except ValueError:
        lv, coord, atomtypes, atomnums = parse_poscar(ifile)
    if 'frozen_axes'in args:
        frozen_axes=args['frozen_axes']
        if len(frozen_axes) == 0:
            frozen_axes=[0,1,2]
    else:
        frozen_axes=[0,1,2]
    if 'direct' in args:
        direct=args['direct']
    else:
        direct=False
    if 'ranges' in args:
        ranges=args['ranges']
    frozen_atoms=[]
        
    seldyn=['' for i in range(sum(atomnums))]
    if mode == 'type':
        counter=0
        for i in atomtypes:
            if i in args['frozen_atoms']:
                for j in range(atomnums[atomtypes.index(i)]):
                    frozen_atoms.append(counter)
                    counter+=1
            else:
                counter+=atomnums[atomtypes.index(i)]
    elif mode == 'number':
        frozen_atoms=args['frozen_atoms']
        
    elif mode == 'position':
        if direct == False:
            ref=dot(lv,array(args['reference_position']))
        else:
            ref=array(args['reference_position'])
        for i in range(sum(atomnums)):
            if norm(coord[i]-ref)<args['tolerance']:
                frozen_atoms.append(i)
                
    elif mode == 'range':
        for i in range(len(coord)):
            counter=0
            for j in range(3):
                if coord[i][j]>ranges[j][0] and coord[i][j]<ranges[j][1]:
                    counter+=1
            if counter==3:
                frozen_atoms.append(i)
    else:
        print('mode not recognized. exiting...')
        sys.exit(1)
    
    #when selective dynamics is applied to atoms, they are reference starting at zero. if reading atom numbers from VESTA, subtract one from each.
    for i in range(len(coord)):
        if i in frozen_atoms:
            for j in range(3):
                if j in frozen_axes and not reverse:
                    seldyn[i]+='F'
                elif j in frozen_axes and reverse:
                    seldyn[i]+='T'
                elif reverse:
                    seldyn[i]+='F'
                else:
                    seldyn[i]+='T'
        elif reverse:
            seldyn[i]+='FFF'
        else:
            seldyn[i]='TTT'
    
    print(str(len(frozen_atoms))+' atoms selected')
    if reverse == False:
        print('atoms frozen:' + str(frozen_atoms))
    else:
        print('atoms un-frozen:' + str(frozen_atoms))
            
    write_poscar(ofile, lv, coord, atomtypes, atomnums, seldyn=seldyn)
    
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
    print('new POSCAR written to: '+str(ofile))
    
if __name__ == '__main__':
    ifile='./POSCAR'
    ofile='./POSCAR_seldyn'
    reverse=False
    direct=False
    mode='none'
    frozen_axes=[]
    short_opts='hi:m:rf:p:t:dw:o:a:'
    long_opts=['help','input=','mode=','reverse','frozen_atoms=','reference_position=','tolerance=','direct','range=','output=','axes=']
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
    -m, --mode                     optional modes are type, range, position, and number
    -p, --reference_position       specify the reference point for position mode
    -t, --tolerance                atoms within this distance of the reference point will freeze
    -f, --frozen_atoms             specify the index of atoms to be frozen or the type
    -w, --range                    for range mode, specify the range of atoms to be selected
                                   form is: x_min,x_max,y_min,y_max,z_min,z_max
                                   values must be Cartesian coordinates
    -o, --output                   speicfy an output other than ./POSCAR_seldyn
    -a, --frozen_axes              select which axes will be frozen
                                   default is all
these options take no value:
    -d, --direct                   use this to specify the reference point in Cartesian
    -r, --reverse                  specify un-frozen atoms rather than frozen
help options:
    -h, --help                    display this error message
''')
            sys.exit()
        if i in ['-i','--input']:
            ifile=j
        if i in ['-o','--output']:
            ofile=j
        if i in ['-m','--mode']:
            mode=str(j)
        if i in ['-r','--reverse']:
            reverse=True
        if i in ['-f','--frozen_atoms']:
            frozen_atoms=j.split(',')
            try:
                frozen_atoms=[int(i) for i in frozen_atoms]
            except ValueError:
                pass
        if i in ['-p','--reference_position']:
            reference_position=array([float(k) for k in j.split(',')])
        if i in ['-t','--tolerance']:
            tolerance=float(j)
        if i in ['-d','--direct']:
            direct=True
        if i in ['-w','--range']:
            ranges=j.split(',')
            ranges=[[float(ranges[i]),float(ranges[j])] for i,j in zip([0,2,4],[1,3,5])]
        if i in ['-a','-axes']:
            frozen_axes=j.split(',')
            try:
                frozen_axes=[int(i) for i in frozen_axes]
            except ValueError:
                pass
    
    if mode == 'type':
        try:
            add_seldyn(mode,ifile,ofile,reverse,frozen_atoms=frozen_atoms,frozen_axes=frozen_axes)
        except TypeError:
            print('incorrect use of arguments for mode: type')
            sys.exit(2)
    elif mode == 'number':
        try:
            add_seldyn(mode,ifile,ofile,reverse,frozen_atoms=frozen_atoms,frozen_axes=frozen_axes)
        except TypeError:
            print('incorrect use of arguments for mode: number')
            sys.exit(2)
    elif mode == 'position':
        try:
            add_seldyn(mode,ifile,ofile,reverse,reference_position=reference_position,tolerance=tolerance,direct=direct,frozen_axes=frozen_axes)
        except TypeError:
            print('incorrect use of arguments for mode: position')
            sys.exit(2)
    elif mode == 'range':
        try:
            add_seldyn(mode,ifile,ofile,reverse,ranges=ranges,frozen_axes=frozen_axes)
        except TypeError:
            print('incorrect use of arguments for mode: range')
            sys.exit(2)
    else:
        print('no mode specified. exiting...')
        sys.exit(1)
