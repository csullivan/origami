# origami
A convenient question-answer application for generating input files for the nuclear reactions code WSAW, FOLD, and DWHI. It also can take a single column based input file.


Example input file:

```
projectile                                                      # Generating WSAW input file for projectile or target system?
wsaw_projectile.inp                                             # Enter a name for the WSAW input file to be generated
6LI6LI                                                          # Enter a filename for the WSAW output file. Note that this filename is restricted to 8 characters or less
../oxbash/proj-ejec/dens.dao_li6                                # Enter path for dens projectile file from oxbash/nushell
../oxbash/proj-ejec/dens.dao_li6                                # Enter path for dens ejectile file from oxbash/nushell
target                                                          # Generating WSAW input file for projectile or target system?
wsaw_target.inp                                                 # Enter a name for the WSAW input file to be generated
12C12C                                                          # Enter a filename for the WSAW output file. Note that this filename is restricted to 8 characters or less
../oxbash/target-recoil/dens.dao_c12_sk12                       # Enter path for dens projectile file from oxbash/nushell
../oxbash/target-recoil/dens.dao_c12_sk12                       # Enter path for dens ejectile file from oxbash/nushell
../li6_origami.inp                                              # Enter a template input file to load
fold.inp                                                        # Enter a name for the FOLD input file to be generated
no                                                              # Would you like the full version of entry, enter no for the abridged version.
12C6LI                                                          # Enter a filename for the FOLD output file. Note that this filename is restricted to 8 characters or less
1.0                                                             # Change in spin projectile/ejectile: <<<<<
1.0                                                             # Change in isospin projectile/ejectile:
../oxbash/proj-ejec/a2002b020.obd                               # Enter filename/path to OXBASH obtd file for projectile/ejectile
1.0                                                             # Change in spin target/recoil: <<<<<
1.0                                                             # Change in isospin target/recoil:
../oxbash/target-recoil/a0008a220.obd                           # Enter filename/path to OXBASH obtd file for target/recoil
2                                                               # The number of (jr,jp,jt) that will be entered
011                                                             # (jr,jp,jt)
211                                                             # (jr,jp,jt)
dwhi.inp                                                        # Enter a name for the DWHI input file to be generated
1210000041000000                                                # Enter dwhi options bits (ICON(0:10) array)
X(x,y)Y                                                         # Enter reaction name string
61                                                              # Number of angles:
0                                                               # Initial angle:
0.2                                                             # Angle step size:
160                                                             # Number of partial waves for elastic
2.24                                                            # coulomb radius (multiplies TARGET mass to 1/3)
1                                                               # Incoming channel potential option (1=WS, 2=surface WS, 3=second derivative, 15 = read in potential)
-120.0                                                          # OMP volume (real)
1.27                                                            # OMP radius (real)
0.840                                                           # OMP diffuseness (real)
0.0                                                             # OMP spin orbit (real)
-34.0                                                           # OMP volume (imaginary)
1.720                                                           # OMP radius (imaginary)
0.690                                                           # OMP diffuseness (imaginary)
0.0                                                             # OMP spin orbit (imaginary)
-18.67                                                          # Enter reaction Qvalue (outgoing energy):
2.24                                                            # coulomb radius (multiplies RECOIL mass to 1/3)
1                                                               # Outgoing channel potential option (1=WS, 2=surface WS, 3=second derivative, 15 = read in potential)
-120.0                                                          # OMP volume (real)
1.27                                                            # OMP radius (real)
0.840                                                           # OMP diffuseness (real)
0.0                                                             # OMP spin orbit (real)
-34.0                                                           # OMP volume (imaginary)
1.720                                                           # OMP radius (imaginary)
0.690                                                           # OMP diffuseness (imaginary)
0.0                                                             # OMP spin orbit (imaginary)
dwhi.plot                                                       # Enter path for dwhi plot file:

```
