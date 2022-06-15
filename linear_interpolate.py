import os
import numpy as np
import time
import matplotlib.pyplot as plt
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
                    
class dos_vs_pos():
    def __init__(self,filepath,npts,num_adlayer_atoms,load_doscars=True):
        start=time.time()
        self.pos=[]
        self.npts=npts
        self.fp=filepath
        for i in range(self.npts):
            os.chdir(r"Y:\Au_RbI\Ben_2\rbi_trigonal_vertical_shifts_v2\_{}".format(i))
            if load_doscars:
                dos,energies=parse_doscar(r"./DOSCAR")[:2]
                if i==0:
                    self.ldos=np.zeros((3,npts,len(energies)))
                    self.energies=np.zeros(len(energies))
                self.energies+=energies/self.npts
                for j in range(1,49):
                    for k in dos[j]:
                        self.ldos[0,i]+=k
                for j in range(-6,-3):
                    for k in dos[j]:
                        self.ldos[1,i]+=k
                for j in range(-3,0):
                    for k in dos[j]:
                        self.ldos[2,i]+=k
            lv,coord,atomtypes,atomnums=parse_poscar(r"./POSCAR")[:4]
            self.pos.append(np.average(coord[-num_adlayer_atoms:,2])-np.max(coord[:-num_adlayer_atoms,2]))
            
        self.lv=lv
        self.atomtypes=atomtypes
        self.atomnums=atomnums
        self.pos=np.array(self.pos)
        
        print('total time to read files: {} min'.format((time.time()-start)/60))
        
        self.peak_pos=[]
        self.peak_energies=[]
        
    def write_ldos(self):
        for i in range(len(self.atomtypes)):
            np.savetxt(os.path.join(self.fp,'{}_dos'.format(self.atomtypes[i])),self.ldos[i])
        np.savetxt(os.path.join(self.fp,'energies'),self.energies)
            
    def load_ldos(self):
        tempvar=[]
        for i in range(len(self.atomtypes)):
            tempvar.append(np.loadtxt(os.path.join(self.fp,'{}_dos'.format(self.atomtypes[i]))))
        self.ldos=np.zeros((len(self.atomtypes),np.shape(tempvar)[1],np.shape(tempvar)[2]))
        for i in range(len(self.atomtypes)):
            self.ldos[i]+=tempvar[i]
        self.energies=np.loadtxt(os.path.join(self.fp,'energies'))
            
    def find_ldos_peaks(self,erange,min_pos,types_to_plot):
        partial_dos=np.zeros((self.npts,len(self.energies)))
        counter=0
        for i in range(len(self.atomtypes)):
            if self.atomtypes[i] in types_to_plot:
                partial_dos+=self.ldos[i]
                counter+=1
        partial_dos/=counter
        
        self.peak_pos.append([])
        self.peak_energies.append([])
        min_pos=np.argmin(abs(self.pos-min_pos))
        for i in range(2):
            erange[i]=np.argmin(abs(self.energies-erange[i]))
        ewidth=erange[1]-erange[0]
        for i in range(min_pos,self.npts):
            self.peak_pos[-1].append(self.pos[i])
            temp_index=np.argmax(partial_dos[i,erange[0]:erange[1]])+erange[0]
            self.peak_energies[-1].append(self.energies[temp_index])
            erange[0]=int(temp_index-ewidth/2)
            erange[1]=int(temp_index+ewidth/2)
            
    def plot_ldos_peaks(self,fit=True):
        if not hasattr(self,'peak_fig'):
            self.peak_fig,self.peak_ax=plt.subplots(1,1)
            
        for i in range(len(self.peak_pos)):
            self.peak_ax.scatter(self.peak_energies[i],self.peak_pos[i])
            self.ldos_ax.scatter(self.peak_energies[i],self.peak_pos[i])
        self.peak_fig.show()
            
    def plot_dos_vs_pos(self,types_to_plot):
        partial_dos=np.zeros((self.npts,len(self.energies)))
        counter=0
        for i in range(len(self.atomtypes)):
            if self.atomtypes[i] in types_to_plot:
                partial_dos+=self.ldos[i]
                counter+=1
        partial_dos/=counter
            
        self.ldos_fig,self.ldos_ax=plt.subplots(1,1)
        self.ldos_ax.pcolormesh(np.array([self.energies for i in range(len(self.pos))]),np.array([[self.pos[i] for j in range(len(self.energies))] for i in range(len(self.pos))]),partial_dos,shading='nearest',cmap='vivid')
        self.ldos_ax.set(ylabel='substrate-adlayer seperation / $\AA$')
        self.ldos_ax.set(xlabel='energy - $E_F$ / eV')
        self.ldos_fig.show()

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
            
#reads DOSCAR
def parse_doscar(filepath):
    with open(filepath,'r') as file:
        line=file.readline().split()
        atomnum=int(line[0])
        for i in range(5):
            line=file.readline().split()
        nedos=int(line[2])
        ef=float(line[3])
        dos=[]
        energies=[]
        for i in range(atomnum+1):
            if i!=0:
                line=file.readline()
            for j in range(nedos):
                line=file.readline().split()
                if i==0:
                    energies.append(float(line[0]))
                if j==0:
                    temp_dos=[[] for k in range(len(line)-1)]
                for k in range(len(line)-1):
                    temp_dos[k].append(float(line[k+1]))
            dos.append(temp_dos)
    energies=np.array(energies)-ef

    #orbitals contains the type of orbital found in each array of the site projected dos
    num_columns=np.shape(dos[1:])[1]
    if num_columns==3:
        orbitals=['s','p','d']
    elif num_columns==6:
        orbitals=['s_up','s_down','p_up','p_down','d_up','d_down']
    elif num_columns==9:
        orbitals=['s','p_y','p_z','p_x','d_xy','d_yz','d_z2','d_xz','d_x2-y2']
    elif num_columns==18:
        orbitals=['s_up','s_down','p_y_up','p_y_down','p_z_up','p_z_down','p_x_up','p_x_down','d_xy_up','d_xy_down','d_yz_up','d_yz_down','d_z2_up','d_z2_down','d_xz_up','d_xz_down','d_x2-y2_up','d_x2-y2_down']
        
    #dos is formatted as [[total dos],[atomic_projected_dos for i in range(atomnum)]]
    #total dos has a shape of (4,nedos): [[spin up],[spin down],[integrated, spin up],[integrated spin down]]
    #atomic ldos have shapes of (6,nedos): [[i,j] for j in [spin up, spin down] for i in [s,p,d]]
    #energies has shape (1,nedos) and contains the energies that each dos should be plotted against
    return dos, energies, ef, orbitals