import os
import numpy as np
from shutil import copyfile

#adlayer structure can be specified either as a path to a POSCAR/CONTCAR or an atom type
#if the adlayer is read from a file, the center of mass will be placed by the 3d vector pos
def interpolate_structure(template,output,adatom,pos,lv1,lv2,n,interpolate_dim=2,adlayer_shift=np.array([1,1,1]),**args):
    try:
        os.mkdir(output)
    except FileExistsError:
        pass
    os.chdir(output)
    lv,coord,atomtypes,atomnums=parse_poscar(os.path.join(template,'POSCAR'))[:4]
    lv1=np.dot(lv1,lv)
    lv2=np.dot(lv2,lv)
    pos=np.dot(pos,lv)
    
    #if rotation is specified as an argument, the adlayer is rotated counter-clockwise by the angle specified by the argument
    #'rotation=30' will rotate the adlayer system 30 degrees counter-clockwise around the center of mass
    if 'rotation' in args:
        theta=args['rotation']/180*np.pi
        rot=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1.0]])
    else:
        rot=np.identity(3)
    
    if os.path.exists(adatom):
        adatom_lv,adatom_coord,adatom_types,adatom_nums=parse_poscar(adatom)[:4]
        com=np.zeros(3)
        for i in adatom_coord:
            com+=i/len(adatom_coord)
        adatom_coord-=com*adlayer_shift
        adatom_coord=np.dot(adatom_coord,rot)
        adatom_coord+=pos
    else:
        adatom_coord=[pos]
        adatom_types=adatom
        adatom_nums=[1]
        
    seldyn=[]
    for i in range(sum(atomnums)):
        seldyn.append('FFF')
    for i in range(len(adatom_coord)):
        seldyn.append('FFT')
    
    if interpolate_dim==2:
        for i in range(n):
            for j in range(n):
                try:
                    os.mkdir('_{}{}'.format(i,j))
                except FileExistsError:
                    pass
                
                tempcoord=[k for k in coord]
                for k in adatom_coord:
                    tempcoord.append(k+lv1*i/(n-1)+lv2*j/(n-1))
                tempcoord=np.array(tempcoord)
                
                temptypes=[k for k in atomtypes]
                for k in adatom_types:
                    temptypes.append(k)
                temptypes=np.array(temptypes)
                
                tempnums=[k for k in atomnums]
                for k in adatom_nums:
                    tempnums.append(k)
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
    if interpolate_dim==1:
        for i in range(n):
            try:
                os.mkdir('_{}'.format(i))
            except FileExistsError:
                pass
            
            tempcoord=[k for k in coord]
            for k in adatom_coord:
                tempcoord.append(k+lv1*i/(n-1))
            tempcoord=np.array(tempcoord)
            
            temptypes=[k for k in atomtypes]
            for k in adatom_types:
                temptypes.append(k)
            temptypes=np.array(temptypes)
            
            tempnums=[k for k in atomnums]
            for k in adatom_nums:
                tempnums.append(k)
            tempnums=np.array(tempnums)
            
            write_poscar(os.path.join('_{}'.format(i),'POSCAR'),lv,tempcoord,temptypes,tempnums,seldyn=seldyn)
            
            for k in ['KPOINTS','POTCAR','INCAR']:
                copyfile(os.path.join(template,k),os.path.join('_{}'.format(i),k))
            with open(os.path.join(template,'job.sh'),'r') as file:
                lines=file.readlines()
                tempvar=lines[2].split('=')
                tempvar[1]='_{}'.format(i)
                lines[2]='='.join(tempvar)+'\n'
            with open(os.path.join('_{}'.format(i),'job.sh'),'w') as file:
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