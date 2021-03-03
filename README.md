Assorted methods for modifying VASP structure files.

add_seldyn.py: Add selective dynamics to a VASP POSCAR file.

Different modes are: 
     - number: specify the atoms to be selected by the index at which they appear in the coordinate list (ie one less than the number in VESTA)
     - type: select atoms by chemical symbol
     - range: specify the spatial range that atoms should be selected in. written as: x_min,x_max,y_min,y_max,z_min,z_max
     - position: select atoms by their distance (tolerance) from a reference point
     
these options take a value:
    -i, --input                    specify an input other than ./POSCAR
    -m, --mode                     required, optional modes are type, range, position, and number
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
    
remove_seldyn.py: removes selective dynamics conditions
these options take a value:
    -i, --input                    specify an input other than ./POSCAR
    -o, --output                   speicfy an output other than ./POSCAR_seldyn
    
copy_seldyn.py: copies selective dynamics conditions from one POSCAR to another
When specifying the files, the selective dynamics conditions of the first argument are copied to the second.
Example:
     copy_seldyn.py copied_from_this_POSCAR to_this_POSCAR

translate_poscar.py: translates the poscar by the given shift. Atom coordinates are adjusted to reside within the unit cell.
Example:
     translate_poscar.py ./POSCAR 0.0,0.0,0.5
     
align_com.py: aligns the center of mass of the second POSCAR to that of the first
Example:
     align_com.py reference_this_POSCAR to_align_this_POSCAR


Compatible with VASP 5.4.4, Python 2.7.13 and 3.6.5
