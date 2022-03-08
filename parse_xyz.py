import numpy as np

def create_poscar_from_xyz(lv):
    pass

def parse_xyz(ifile,lv,atomtypes,atomnums):
    with open(ifile,'r') as f:
        lines=f.readlines()
    for i in range(len(lines)):
        lines[i]=lines[i].split()
    
    coord=np.zeros((int(lines[0]),3))
    atomtypes=[]
    atomnums=[0]
    
    for i in range(2,2+np.shape(coord)[0]):
        if lines[i,0] in atomtypes:
            atomnums[-1]+=1
        else:
            atomtypes.append(lines[i,0])
        for j in range(3):
            coord[i,j]=float(lines[i,j+1])
            
    return coord,atomnums,atomtypes