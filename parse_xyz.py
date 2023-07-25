import numpy as np

def create_poscar_from_xyz(lv,ifile,ofile):
    coord,atomtypes,atomnums=parse_xyz(ifile)
    write_poscar(ofile,lv,coord,atomtypes,atomnums)

def parse_xyz(ifile):
    with open(ifile,'r') as f:
        lines=f.readlines()
    for i in range(len(lines)):
        lines[i]=lines[i].split()
    
    coord=np.zeros((int(lines[0][0]),3))
    atomtypes=[]
    atomnums=[]
    
    for i in range(2,2+np.shape(coord)[0]):
        if lines[i][0] in atomtypes:
            atomnums[atomtypes.index(lines[i][0])]+=1
        else:
            atomtypes.append(lines[i][0])
            atomnums.append(1)
        for j in range(3):
            coord[i-2,j]=float(lines[i][j+1])
            
    return coord,atomtypes,atomnums

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