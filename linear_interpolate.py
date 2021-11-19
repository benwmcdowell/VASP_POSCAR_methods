import os
import numpy as np
from shutil import copyfile

def interpolate_structure(template,output,adatom_type,pos,lv1,lv2,n):
    os.mkdir(output)
    os.chdir(output)
    lv,coord,atomtypes,atomnums=parse_poscar(os.path.join(template,'POSCAR'))[:4]
    lv1=np.dot(lv1,lv)
    lv2=np.dot(lv2,lv)
    pos=np.dot(pos,lv)
    seldyn=[]
    for i in range(sum(atomnums)):
        seldyn.append('FFF')
    seldyn.append('FFT')
    for i in range(n):
        for j in range(n):
            os.mkdir('_{}{}'.format(i,j))
            tempcoord=[k for k in coord]
            tempcoord.append(pos+lv1*i/(n-1)+lv2*j/(n-1))
            tempcoord=np.array(tempcoord)
            temptypes=[k for k in atomtypes]
            temptypes.append(adatom_type)
            temptypes=np.array(temptypes)
            tempnums=[k for k in atomnums]
            tempnums.append(1)
            tempnums=np.array(tempnums)
            write_poscar(os.path.join('_{}{}'.format(i,j),'POSCAR'),lv,tempcoord,temptypes,tempnums,seldyn=seldyn)
            for k in ['KPOINTS','POTCAR','INCAR']:
                copyfile(os.path.join(template,k),os.path.join('_{}{}'.format(i,j),k))
            with open(os.path.join(template,'job.sh'),'r') as file:
                lines=file.readlines()
                tempvar=lines[2].split('=')
                tempvar[1]='_{}{}'.format(i,j)
                lines[2]='='.join(tempvar)+'\n'
            with open(os.path.join('_{}{}'.format(i,j),'job.sh'),'w') as file:
                for k in lines:
                    file.write(k)
                
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