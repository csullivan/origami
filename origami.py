#!/usr/bin/python
# Copyright 2014-2016 Chris Sullivan
import readline, glob
import fortranformat as form
import numpy as np
import math
from Queue import *
from nuclear_qvalue import get_ame_q as qval_ame12
from sympy.physics.quantum.cg import CG
from sympy import S
import argparse


input_log = []

def get(prompt, default):
    result = raw_input("%s [%s] " % (prompt, default)) or default
    if '#' in str(result):
        if '#' in result.split()[0]:
            result = default
        else:
            result = result.split()[0].strip()
    input_log.append([result,prompt])
    return result

def save_input_log():
    print
    output = open('./inputfile_log','wb')
    for entry in input_log:
        prompt = entry[1].replace("\n"," ")
        output.write(padstr(str(entry[0]),64)+"# "+prompt+"\n")
        #print padstr(str(entry[0]),64)+"# "+prompt
    output.close()


def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]

def init_tab_complete():
    readline.set_completer_delims(' \t\n;')
    readline.parse_and_bind("tab: complete")
    readline.set_completer(complete)

def init_obtd_key_map():
    """ lookup[oxbash_key] = fold_key"""
    lookup = np.zeros(45)
    lookup[1] = 1
    lookup[2] = 3
    lookup[3] = 2
    lookup[4] = 6
    lookup[5] = 5
    lookup[6] = 4
    lookup[7] = 10
    lookup[8] = 9
    lookup[9] = 8
    lookup[10] = 7
    # more needed....
    lookup = [int(x) for x in lookup]
    return lookup

def get_zT_prefactor(Ji,Ti,Tiz,dT,dTz,Tf,Tfz):
    cg = CG(Ti,Tiz,dT,dTz,Tf,Tfz)
    return ((math.sqrt(2*dT+1))/(math.sqrt(2*Ji+1)*math.sqrt(2*Tf+1)))*float(cg.doit())

all_states_dict = {}
def get_obtd(obtd_filename,dj,state=0):
    try:
        return all_states_dict[obtd_filename][state]
    except:
        pass

    oxbash_to_fold_key_map = init_obtd_key_map()
    try:
        obtd_file = open(obtd_filename,"rb")
    except IOError:
        print "Error, file does not exist"; exit()

    obtd_list = []
    for i,line in enumerate(obtd_file):
        if '!' == line[1] or i == 0:
            continue
        line = line.split()
        obtd_list.append(line)
    marks = []
    energies = Queue()
    for i,line in enumerate(obtd_list):
        if '!' not in line[0] and len(line) == 7 and float(line[0][:-1]) == dj:
            marks.append(i)
            e = float(line[3][:-1])-float(line[4][:-1])
            energies.put(e)

    # the last mark is the end of the file
    marks.append(len(obtd_list))
    obtds = []
    all_states = []
    for n,mark in enumerate(marks[:-1]):
        fold_obtds = []
        for i in range(mark+1,marks[n+1]):
            #if len(obtd_list[i])==7:
            #    break
            if len(obtd_list[i])>3:
                obtds.append([int(obtd_list[i][0][:-1]),int(obtd_list[i][1][:-1]),float(obtd_list[i][3][:-1])])
                fold_obtds.append([oxbash_to_fold_key_map[int(obtd_list[i][0][:-1])],oxbash_to_fold_key_map[int(obtd_list[i][1][:-1])],float(obtd_list[i][3][:-1])])

        all_states.append([energies.get(),fold_obtds])

    all_states_dict[obtd_filename] = all_states
    return all_states[state]

def get_number_of_obtds_from_fold_input_file(filelist):
    obtd_count_projectile = 0
    for i,line in enumerate(filelist):
        if i>4:
            line = line.split()
            if len(line) == 5:
                obtd_count_projectile+=1
            if int(line[0])==-1:
                break
    obtd_count_target = 0
    for i,line in enumerate(filelist):
        if i>(4+obtd_count_projectile+5):
            line = line.split()
            if len(line) == 5:
                obtd_count_target+=1
            if int(line[0])==-1:
                break
    return obtd_count_projectile,obtd_count_target




def padstr(x,length=8):
    if len(x) < length:
        for i in range(0,length-len(x)):
            x += " "
    return x

class CEReactions(object):
    def __init__(self):
        self.nwsaw = 0
    def wsaw_inputfile_from_dens(self):
        system = get("Generating WSAW input file for projectile or target system?", "projectile" if self.nwsaw == 0 else "target")
        self.nwsaw += 1

        filename = get("Enter a name for the WSAW input file to be generated","wsaw_"+system+".inp")
        file = open(filename,"wb")

        self.wsaw_input_filename = filename

        ################### Line 1 ###################

        line = form.FortranRecordWriter('(A10,A10,3I5)')
        line = line.write([0.1,20.,1,150,0])
        file.write(line+'\n')

        ################### Line 2 ###################
        output_file = \
            get("Enter a filename for the WSAW output file.\nNote that this filename is restricted to 8 characters or less","B13BE13")
        if len(output_file) > 8:
            print "Output WSAW filename is > 8 characters"
            exit()
        output_file = padstr(output_file)
        line = form.FortranRecordWriter('(A8)')
        line = line.write([output_file[0:8]])
        file.write(line+'\n')

        if system == "projectile":
            self.wsaw_proj = output_file
        elif system == "target":
            self.wsaw_targ = output_file
        else:
            print "Unknown system requested for wsaw input generation"
            exit()

        ################### Line 3-EOF ###################
        projdens = get("Enter path for dens projectile file from oxbash/nushell","13b-dens.dao")
        finaldens = get("Enter path for dens ejectile file from oxbash/nushell","13be-dens.dao")
        self.parse_dens(projdens,finaldens,system,file)

        line = form.FortranRecordWriter('(I2)')
        line = line.write([-1])
        file.write(line)



    def parse_dens(self,initialdens,finaldens,system,wsawfile):
        try:
            print "######", initialdens
            initial_dens_file = open(initialdens,"rb")
            final_dens_file = open(finaldens,"rb")
        except IOError:
            print "Error, file does not exist"; exit()

        for line in initial_dens_file:
            if '* ia,iz' in line:
                line = line.split()
                ai = int(line[3])
                zi = int(line[4])
                break
        for line in final_dens_file:
            if '* ia,iz' in line:
                line = line.split()
                af = int(line[3])
                zf = int(line[4])
                break
        print ai,zi," ",af,zf
        assert(ai==af)

        if system == "projectile":
            self.a_proj = ai
            self.z_proj = zi
            self.a_ejec = af
            self.z_ejec = zf
        elif system == "target":
            self.a_target = ai
            self.z_target = zi
            self.a_recoil = af
            self.z_recoil = zf


        if zi>zf:
            initial_nucleon = 'proton'
            final_nucleon = 'neutron'
        if zi<zf:
            initial_nucleon = 'neutron'
            final_nucleon = 'proton'
        if zi==zf:
            initial_nucleon = 'proton'
            final_nucleon = 'neutron'
            # doesn't matter initial == final
        valence_nucleon_initial = []
        sp_lines = False
        for i,line in enumerate(initial_dens_file):
            if sp_lines:
                if '-----' in line:
                    sp_lines = False
                    break
                valence_nucleon_initial.append(line.split())
            if initial_nucleon+' bound state results' in line:
                sp_lines = True
        valence_nucleon_final = []
        for i,line in enumerate(final_dens_file):
            if sp_lines:
                if '-----' in line:
                    sp_lines = False
                    break
                valence_nucleon_final.append(line.split())
            if final_nucleon+' bound state results' in line:
                sp_lines = True

        self.write_sp_levels(valence_nucleon_initial,ai,zi,initial_nucleon,wsawfile)
        self.write_sp_levels(valence_nucleon_final,af,zf,final_nucleon,wsawfile)

    def write_sp_levels(self,levels_list,a,z,valence,file):
        corez = 0
        if valence == 'proton':
            corez = z - 1
        else:
            corez = z
        for i,lvl in enumerate(levels_list):
            if i == 1:
                continue # first line contains indices
            if len(lvl)==0:
                continue
            if len(lvl) != 9:
                if len(lvl)==8:
                    new_lvl = []
                    for i,col in enumerate(lvl):
                        if i == 2:
                            new_lvl.append(col[0])
                            new_lvl.append(col[1:])
                        else:
                            new_lvl.append(col)
                    lvl = new_lvl
                else:
                    print "Error in sp levels, bad level: ", lvl
                    continue
            #line = form.FortranRecordWriter('(7F10.2,I2)')
            line = form.FortranRecordWriter('(F10.3,2F10.0,3F10.2,F10.1,I2)')
            line = line.write([a-1,corez,60.,0.65,1.25,1.25,7.0])
            file.write(line+'\n')
            line = form.FortranRecordWriter('(7F10.2,I2)')
            ebind = float(lvl[4])
            if ebind == -0.2:
                ebind = 1.0
            orbital_ang = lvl[2]
            if orbital_ang == 's':
                orbital_ang = 0
            elif orbital_ang == 'p':
                orbital_ang = 1
            elif orbital_ang == 'd':
                orbital_ang = 2
            elif orbital_ang == 'f':
                orbital_ang = 3
            elif orbital_ang == 'g':
                orbital_ang = 4
            elif orbital_ang == 'h':
                orbital_ang = 5
            else:
                print "error determinging orbital angular momentum: ", orbital_ang
                exit()
            ebind = abs(ebind)

            j = lvl[3]
            j = float(j[:-2])/float(j[len(j)-1])
            line = line.write([ebind, 1, orbital_ang, float(lvl[1])-1, 1 if valence == 'proton' else 0, j, 0.5])
            file.write(line+'\n')


    def fold_inputfile_from_template(self,template_filename,manual_entry,nstate=1):
        filename = get("Enter a name for the FOLD input file to be generated","fold.inp")
        assert(template_filename != filename)

        self.fold_input_filename = filename

        filelist = []
        for fileline in open(template_filename,"rb"):
            filelist.append(fileline)
            obtd_count_proj, obtd_count_targ = get_number_of_obtds_from_fold_input_file(filelist)

        n=0 # line counter
        ################### Line 1 ###################
        line = form.FortranRecordReader('(2I5,A8)')
        output_file = line.read(filelist[n])[2]; n+=1
        ################### Line 2 ###################
        line = form.FortranRecordReader('(I5,F5.2,F10.0,F10.0,I10,I4,I4)')
        nr,ns,beam_energy_lab,a,blank,blank,blank = line.read(filelist[n]); n+=1
        ################# Projectile #################
        ################### Line 3 ###################
        line = form.FortranRecordReader('(F10.1,A1,F9.1,A1)')
        jpf,ppf,jpi,ppi = line.read(filelist[n]); n+=1
        djp = max(abs(jpf - jpi),1.0) # A guess, user can enter below
        ################### Line 4 ###################
        line = form.FortranRecordReader('(F5.1,A10,F10.1,A11)')
        tpf,tzpf,tpi,tzpi = line.read(filelist[n]); n+=1
        ################### Line 5 ###################
        line = form.FortranRecordReader('(I5,I5,F7.3)')
        ntypf_p,koptn_p,alpha_p = line.read(filelist[n]); n+=1
        ################### Line 6 ###################
        n+=obtd_count_proj+1 # extra +1 for -1 -1 terminating line
        line = form.FortranRecordReader('(A8)')
        wsaw_proj = line.read(filelist[n])[0]; n+=1
        ################### Target ###################
        ################### Line 3 ###################
        line = form.FortranRecordReader('(F10.1,A1,F9.1,A1)')
        jtf,ptf,jti,pti = line.read(filelist[n]); n+=1
        djt = max(abs(jtf - jti),1.0) # a guess
        ################### Line 4 ###################
        line = form.FortranRecordReader('(F5.1,A10,F10.1,A11)')
        ttf,tztf,tti,tzti = line.read(filelist[n]); n+=1
        ################### Line 5 ###################
        line = form.FortranRecordReader('(I5,I5,F7.3)')
        ntypf_t,koptn_t,alpha_t = line.read(filelist[n]); n+=1
        ################### Line 6 ###################
        n+=obtd_count_targ+1 # extra +1 for -1 -1 terminating line
        line = form.FortranRecordReader('(A8)')
        wsaw_targ = line.read(filelist[n])[0]; n+=1
        ################### Line 6 ###################
        line = form.FortranRecordReader('(F6.3,F9.2,F11.3,A12)')
        fnrm1,fnrm2,blank,blank = line.read(filelist[n]); n+=1
        ################### Line 7 ###################
        line = form.FortranRecordReader('(I5)')
        nform = line.read(filelist[n])[0]; n+=1


        file = open(filename,"wb")


        ################### Line 1 ###################
        output_file = get("Enter a filename for the FOLD output file.\nNote that this filename is restricted to 8 characters or less",output_file)
        output_file = padstr(output_file)
        line = form.FortranRecordWriter('(I5,I5,A8)')
        line = line.write([1,1,output_file[0:9]])
        file.write(line+'\n')

        self.fold_filename = output_file

        ################### Line 2 ###################
        if manual_entry:
            nr = int(get("Number of integration steps:\n",nr))
            ns = float(get("Step size (fm):\n",ns))
            beam_energy_lab = float(get("Bombarding energy (MeV):\n",beam_energy_lab))
            if self.nwsaw == 2:
                a = self.a_proj
            else:
                a = float(get("Projectile mass number (A):\n",a))
        line = form.FortranRecordWriter('(I5,F5.2,F10.0,F10.0,I10,I4,I4)')
        line = line.write([nr,ns,beam_energy_lab,a,1,1,1])
        file.write(line+'\n')

        self.nr = nr
        self.dr = ns
        self.beam_energy_lab = beam_energy_lab
        if self.nwsaw != 2:
            self.a_proj = a
        else:
            if self.a_proj != a:
                print "Warning projectile A in FOLD template input is different from that entered for WSAW."


        ################### Line 3 ###################
        if manual_entry:
            jpf = float(get("Spin of the ejectile:\n",jpf))
            ppf = str(get("Parity of the ejectile:\n",ppf))
            jpi = float(get("Spin of the projectile:\n",jpi))
            ppi = str(get("Parity of the projectile:\n",ppi))
        line = form.FortranRecordWriter('(F10.1,A1,F9.1,A1)')
        line = line.write([jpf,ppf,jpi,ppi])
        file.write(line+'\n')

        self.j_ejec = jpf
        self.p_ejec = ppf
        self.j_proj = jpi
        self.p_proj = ppi

        ################### Line 4 ###################
        if manual_entry:
            tpf = float(get("Isospin of the ejectile:\n",tpf))
            tzpf = str(get("Isospin projection of the ejectile:\n",tzpf.strip()))
            tpi = float(get("Isospin of the projectile:\n",tpi))
            tzpi = str(get("Isospin projection of the projectile:\n",tzpi.strip()))
        line = form.FortranRecordWriter('(F5.1,A10,F10.1,A11)')
        line = line.write([tpf,tzpf,tpi,tzpi])
        file.write(line+'\n')

        self.t_ejec = tpf
        self.tz_ejec = tzpf
        self.t_proj = tpi
        self.t_proj = tzpi

        ################### Line 5 ###################
        if manual_entry:
            ntypf_p = int(get("NTYPF (=1 static, =2 inelastic, =3 charge exchange):\n",ntypf_p))
            koptn_p = int(get("Transition amplitude type: 1 (S[T]), 2 (S[pn]), 3 (Z[T] - OXBASH OBTD), 4 (Z[pn]), 5 (Wildenthal):\n",koptn_p))
            alpha_p = float(get("Alpha (=0.000 to use WSAW radial wave functions):\n",alpha_p))
        line = form.FortranRecordWriter('(I5,I5,F7.3)')
        line = line.write([ntypf_p,koptn_p,alpha_p])
        file.write(line+'\n')
        ################### Line 6 ###################
        ################## OBTDs ###################
        djp = float(get("Change in spin projectile/ejectile:\n",djp))
        dtp = float(get("Change in isospin projectile/ejectile:\n",1.0))

        self.djp = djp
        self.dtp = dtp

        obtd_filename_p = str(get("Enter filename/path to OXBASH obtd file for projectile/ejectile","t331dp150.obd"))
        fold_p_obtds = get_obtd(obtd_filename_p,djp)
        prefactor_p = get_zT_prefactor(jpi,tpi,float(tzpi),dtp,float(tzpf)-float(tzpi),tpf,float(tzpf))
        print "Projectile prefactor = ",prefactor_p
        self.ejec_ex = fold_p_obtds[0]
        for obtd in fold_p_obtds[1]:
            line = form.FortranRecordWriter('(I5,I5,I5,F5.1,F17.6)')
            obtd[2] *= prefactor_p
            line = line.write([obtd[0],obtd[1],int(djp),0.0,obtd[2]])
            file.write(line+'\n')
        line = form.FortranRecordWriter('(I5,I5)')
        line = line.write([-1,-1])
        file.write(line+'\n')

        if self.nwsaw == 2:
            wsaw_proj = self.wsaw_proj
        else:
            if manual_entry:
                wsaw_proj = str(get("Enter filename of projectile WSAW radial wavefunction file:\n",wsaw_proj))
                wsaw_proj = padstr(wsaw_proj)

        line = form.FortranRecordWriter('(A8)')
        line = line.write([wsaw_proj])
        file.write(line+'\n')

        ################### Target ###################
        ################### Line 3 ###################
        if manual_entry:
            jtf = float(get("Spin of the recoil:\n",jtf))
            ptf = str(get("Parity of the recoil:\n",ptf))
            jti = float(get("Spin of the target:\n",jti))
            pti = str(get("Parity of the target:\n",pti))
        line = form.FortranRecordWriter('(F10.1,A1,F9.1,A1)')
        line = line.write([jtf,ptf,jti,pti])
        file.write(line+'\n')

        self.j_recoil = jtf
        self.p_recoil = ptf
        self.j_target = jti
        self.p_target = pti

        ################### Line 4 ###################
        if manual_entry:
            ttf = float(get("Isospin of the recoil:\n",ttf))
            tztf = str(get("Isospin projection of the recoil:\n",tztf.strip()))
            tti = float(get("Isospin of the target:\n",tti))
            tzti = str(get("Isospin projection of the target:\n",tzti.strip()))
        line = form.FortranRecordWriter('(F5.1,A10,F10.1,A11)')
        line = line.write([ttf,tztf,tti,tzti])
        file.write(line+'\n')

        self.t_recoil = ttf
        self.tz_recoil =tztf
        self.t_target = tti
        self.tz_target = tzti

        ################### Line 5 ###################
        if manual_entry:
            ntypf_t = int(get("NTYPF (=1 static, =2 inelastic, =3 charge exchange):\n",ntypf_t))
            koptn_t = int(get("Transition amplitude type: 1 (S[T]), 2 (S[pn]), 3 (Z[T] - OXBASH OBTD), 4 (Z[pn]), 5 (Wildenthal):\n",koptn_t))
            alpha_t = float(get("Alpha (=0.000 to use WSAW radial wave functions):\n",alpha_t))
        line = form.FortranRecordWriter('(I5,I5,F7.3)')
        line = line.write([ntypf_t,koptn_t,alpha_t])
        file.write(line+'\n')
        ################### Line 6 ###################
        ################## OBTDs ###################
        djt = float(get("Change in spin target/recoil:\n",djt))
        dtt = float(get("Change in isospin target/recoil:\n",1.0))

        self.djt = djt
        self.dtt = dtt

        obtd_filename_t = str(get("Enter filename/path to OXBASH obtd file for target/recoil","t331dp150.obd"))
        fold_t_obtds = get_obtd(obtd_filename_t,djt,state=nstate)
        prefactor_t = get_zT_prefactor(jti,tti,float(tzti),dtt,float(tztf)-float(tzti),ttf,float(tztf))
        print "Target prefactor = ",prefactor_t
        self.recoil_ex = fold_t_obtds[0]
        for obtd in fold_t_obtds[1]:
            line = form.FortranRecordWriter('(I5,I5,I5,F5.1,F17.6)')
            obtd[2] *= prefactor_t
            line = line.write([obtd[0],obtd[1],int(djt),0.0,obtd[2]])
            file.write(line+'\n')
        line = form.FortranRecordWriter('(I5,I5)')
        line = line.write([-1,-1])
        file.write(line+'\n')

        if self.nwsaw == 2:
            wsaw_targ = self.wsaw_targ
        else:
            if manual_entry:
                wsaw_targ = str(get("Enter filename of target WSAW radial wavefunction file:\n",wsaw_targ))
                wsaw_targ = padstr(wsaw_targ)

        line = form.FortranRecordWriter('(A8)')
        line = line.write([wsaw_targ])
        file.write(line+'\n')
        ################### Line 6 ###################
        if manual_entry:
            fnrm1 = float(get("Normalization of SNKE Yukawas (all ranges) for transformation from tNN to tNA:\n",fnrm1))
            fnrm2 = float(get("kA (momentum of projectile in NA frame) instead of lab -- of particular importance for very light (A<10) target nuclei:\n",fnrm2))
        line = form.FortranRecordWriter('(F6.3,F9.2,F11.3,A12)')
        line = line.write([fnrm1,fnrm2,1.000,'love_140'])
        file.write(line+'\n')
        nform = int(get("The number of (jr,jp,jt) that will be entered",nform))

        self.nform = nform

        line = form.FortranRecordWriter('(I5)')
        line = line.write([nform])
        file.write(line+'\n')

        self.formfactors = []
        for i in range(0,nform):
            jrjpjt = str(get("(jr,jp,jt)","011"))

            self.formfactors.append(jrjpjt)

            line = form.FortranRecordWriter('(I5,I5,I5,I5)')
            line = line.write([jrjpjt[0],jrjpjt[1],jrjpjt[2],-1])
            file.write(line+'\n')
            line = form.FortranRecordWriter('(F5.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)')
            line = line.write([1.0,1.0,1.0,1.0,1.0,1.0,1.0,])
            file.write(line+'\n')
            line = form.FortranRecordWriter('(F5.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)')
            line = line.write([1.0,1.0,1.0,1.0,1.0,1.0,1.0,])
            file.write(line+'\n')
        file.write('\n')

    def fold_inputfile(self,nstate=0):
        filename = get("Enter a name for the FOLD input file to be generated","fold.inp")
        if nstate != 0:
            filename = filename.split('.')
            try:
                filename = filename[0]+str(nstate)+"."+filename[1]
            except:
                filename = filename[0]+str(nstate)

        self.fold_input_filename = filename


        file = open(filename,"wb")


        ################### Line 1 ###################
        output_file = get("Enter a filename for the FOLD output file.\nNote that this filename is restricted to 8 characters or less","FORMFAC")
        output_file = padstr(output_file)
        line = form.FortranRecordWriter('(I5,I5,A8)')
        line = line.write([1,1,output_file[0:9]])
        file.write(line+'\n')

        self.fold_filename = output_file

        ################### Line 2 ###################
        nr = int(get("Number of integration steps:\n",600))
        ns = float(get("Step size (fm):\n",0.03))
        beam_energy_lab = float(get("Bombarding energy (MeV):\n",600))
        if self.nwsaw == 2:
            a = self.a_proj
        else:
            a = float(get("Projectile mass number (A):\n",6))
        line = form.FortranRecordWriter('(I5,F5.2,F10.0,F10.0,I10,I4,I4)')
        line = line.write([nr,ns,beam_energy_lab,a,1,1,1])
        file.write(line+'\n')

        self.nr = nr
        self.dr = ns
        self.beam_energy_lab = beam_energy_lab
        if self.nwsaw != 2:
            self.a_proj = a
        else:
            if self.a_proj != a:
                print "Warning projectile A in FOLD template input is different from that entered for WSAW."


        ################### Line 3 ###################
        jpf = float(get("Spin of the ejectile:\n","0.0"))
        ppf = str(get("Parity of the ejectile:\n","+"))
        jpi = float(get("Spin of the projectile:\n","1.0"))
        ppi = str(get("Parity of the projectile:\n","+"))
        line = form.FortranRecordWriter('(F10.1,A1,F9.1,A1)')
        line = line.write([jpf,ppf,jpi,ppi])
        file.write(line+'\n')

        self.j_ejec = jpf
        self.p_ejec = ppf
        self.j_proj = jpi
        self.p_proj = ppi

        ################### Line 4 ###################
        tpf = float(get("Isospin of the ejectile:\n","1.0"))
        tzpf = str(get("Isospin projection of the ejectile:\n","0.0"))
        tpi = float(get("Isospin of the projectile:\n","0.0"))
        tzpi = str(get("Isospin projection of the projectile:\n","0.0"))
        line = form.FortranRecordWriter('(F5.1,A10,F10.1,A11)')
        line = line.write([tpf,tzpf,tpi,tzpi])
        file.write(line+'\n')

        self.t_ejec = tpf
        self.tz_ejec = tzpf
        self.t_proj = tpi
        self.t_proj = tzpi

        ################### Line 5 ###################
        ntypf_p = int(get("NTYPF (=1 static, =2 inelastic, =3 charge exchange):\n","3"))
        koptn_p = int(get("Transition amplitude type: 1 (S[T]), 2 (S[pn]), 3 (Z[T] - OXBASH OBTD), 4 (Z[pn]), 5 (Wildenthal):\n",3))
        alpha_p = float(get("Alpha (=0.000 to use WSAW radial wave functions) (projectile):\n","0.0"))
        line = form.FortranRecordWriter('(I5,I5,F7.3)')
        line = line.write([ntypf_p,koptn_p,alpha_p])
        file.write(line+'\n')
        ################### Line 6 ###################
        ################## OBTDs ###################
        djp = float(get("Change in spin projectile/ejectile:\n","1.0"))
        dtp = float(get("Change in isospin projectile/ejectile:\n","1.0"))

        self.djp = djp
        self.dtp = dtp

        obtd_filename_p = str(get("Enter filename/path to OXBASH obtd file for projectile/ejectile","t331dp150.obd"))
        fold_p_obtds = get_obtd(obtd_filename_p,djp)
        prefactor_p = get_zT_prefactor(jpi,tpi,float(tzpi),dtp,float(tzpf)-float(tzpi),tpf,float(tzpf))
        print "Projectile prefactor = ",prefactor_p
        self.ejec_ex = fold_p_obtds[0]
        for obtd in fold_p_obtds[1]:
            line = form.FortranRecordWriter('(I5,I5,I5,F5.1,F17.6)')
            obtd[2] *= prefactor_p
            line = line.write([obtd[0],obtd[1],int(djp),0.0,obtd[2]])
            file.write(line+'\n')
        line = form.FortranRecordWriter('(I5,I5)')
        line = line.write([-1,-1])
        file.write(line+'\n')

        if self.nwsaw == 2:
            wsaw_proj = self.wsaw_proj
        else:
            wsaw_proj = str(get("Enter filename of projectile WSAW radial wavefunction file:\n","PROJ"))
            wsaw_proj = padstr(wsaw_proj)

        line = form.FortranRecordWriter('(A8)')
        line = line.write([wsaw_proj])
        file.write(line+'\n')

        ################### Target ###################
        ################### Line 3 ###################
        jtf = float(get("Spin of the recoil:\n","1.0"))
        ptf = str(get("Parity of the recoil:\n","+"))
        jti = float(get("Spin of the target:\n","0.0"))
        pti = str(get("Parity of the target:\n","+"))
        line = form.FortranRecordWriter('(F10.1,A1,F9.1,A1)')
        line = line.write([jtf,ptf,jti,pti])
        file.write(line+'\n')

        self.j_recoil = jtf
        self.p_recoil = ptf
        self.j_target = jti
        self.p_target = pti

        ################### Line 4 ###################
        ttf = float(get("Isospin of the recoil:\n","1.0"))
        tztf = str(get("Isospin projection of the recoil:\n","0.0"))
        tti = float(get("Isospin of the target:\n","0.0"))
        tzti = str(get("Isospin projection of the target:\n","0.0"))
        line = form.FortranRecordWriter('(F5.1,A10,F10.1,A11)')
        line = line.write([ttf,tztf,tti,tzti])
        file.write(line+'\n')

        self.t_recoil = ttf
        self.tz_recoil =tztf
        self.t_target = tti
        self.tz_target = tzti

        ################### Line 5 ###################
        ntypf_t = int(get("NTYPF (=1 static, =2 inelastic, =3 charge exchange):\n",3))
        koptn_t = int(get("Transition amplitude type: 1 (S[T]), 2 (S[pn]), 3 (Z[T] - OXBASH OBTD), 4 (Z[pn]), 5 (Wildenthal):\n",3))
        alpha_t = float(get("Alpha (=0.000 to use WSAW radial wave functions) (target):\n",0.0))
        line = form.FortranRecordWriter('(I5,I5,F7.3)')
        line = line.write([ntypf_t,koptn_t,alpha_t])
        file.write(line+'\n')
        ################### Line 6 ###################
        ################## OBTDs ###################
        djt = float(get("Change in spin target/recoil:\n",1.0))
        dtt = float(get("Change in isospin target/recoil:\n",1.0))

        self.djt = djt
        self.dtt = dtt

        obtd_filename_t = str(get("Enter filename/path to OXBASH obtd file for target/recoil","t331dp150.obd"))
        fold_t_obtds = get_obtd(obtd_filename_t,djt,state=nstate)
        prefactor_t = get_zT_prefactor(jti,tti,float(tzti),dtt,float(tztf)-float(tzti),ttf,float(tztf))
        print "Target prefactor = ",prefactor_t
        self.recoil_ex = fold_t_obtds[0]
        for obtd in fold_t_obtds[1]:
            line = form.FortranRecordWriter('(I5,I5,I5,F5.1,F17.6)')
            obtd[2] *= prefactor_t
            line = line.write([obtd[0],obtd[1],int(djt),0.0,obtd[2]])
            file.write(line+'\n')
        line = form.FortranRecordWriter('(I5,I5)')
        line = line.write([-1,-1])
        file.write(line+'\n')

        if self.nwsaw == 2:
            wsaw_targ = self.wsaw_targ
        else:
            wsaw_targ = str(get("Enter filename of target WSAW radial wavefunction file:\n","TARG"))
            wsaw_targ = padstr(wsaw_targ)

        line = form.FortranRecordWriter('(A8)')
        line = line.write([wsaw_targ])
        file.write(line+'\n')
        ################### Line 6 ###################
        fnrm1 = float(get("Normalization of SNKE Yukawas (all ranges) for transformation from tNN to tNA:\n",0.954))
        fnrm2 = float(get("kA (momentum of projectile in NA frame) instead of lab -- of particular importance for very light (A<10) target nuclei:\n",2.46))
        line = form.FortranRecordWriter('(F6.3,F9.2,F11.3,A12)')
        line = line.write([fnrm1,fnrm2,1.000,'love_140'])
        file.write(line+'\n')
        nform = int(get("The number of (jr,jp,jt) that will be entered",2))

        self.nform = nform

        line = form.FortranRecordWriter('(I5)')
        line = line.write([nform])
        file.write(line+'\n')

        self.formfactors = []
        for i in range(0,nform):
            jrjpjt = str(get("(jr,jp,jt)","011"))

            self.formfactors.append(jrjpjt)

            line = form.FortranRecordWriter('(I5,I5,I5,I5)')
            line = line.write([jrjpjt[0],jrjpjt[1],jrjpjt[2],-1])
            file.write(line+'\n')
            line = form.FortranRecordWriter('(F5.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)')
            line = line.write([1.0,1.0,1.0,1.0,1.0,1.0,1.0,])
            file.write(line+'\n')
            line = form.FortranRecordWriter('(F5.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)')
            line = line.write([1.0,1.0,1.0,1.0,1.0,1.0,1.0,])
            file.write(line+'\n')
        file.write('\n')


    def dwhi_inputfile_gen(self,nstate=0):
        filename = get("Enter a name for the DWHI input file to be generated","dwhi.inp")
        if nstate != 0:
            filename = filename.split('.')
            try:
                filename = filename[0]+str(nstate)+"."+filename[1]
            except:
                filename = filename[0]+str(nstate)

        self.fold_input_filename = filename
        file = open(filename,'wb')

        ################### Line 1 ###################
        icon = int(get("Enter dwhi options bits (ICON(0:10) array)","1210000041000000"))
        runname = get("Enter reaction name string","X(x,y)Y")
        runname = padstr(runname,52)
        line = form.FortranRecordWriter('(I16,A4,A52)')
        line = line.write([icon," ",runname])
        file.write(line+'\n')
        ################### Line 2 ###################
        line = form.FortranRecordWriter('(A8)')
        line = line.write([self.fold_filename])
        file.write(line+'\n')
        ################### Line 3 ###################
        nangles = int(get("Number of angles: ","61"))
        angle0 = int(get("Initial angle: ","0"))
        d_angle = float(get("Angle step size: ","0.2"))
        line = form.FortranRecordWriter('(3F7.4)')
        line = line.write([nangles,angle0,d_angle])
        file.write(line+'\n')
        ################### Line 4 ###################
        npartial_waves = int(get("Number of partial waves for elastic","160"))
        line = form.FortranRecordWriter('(6I3)')
        line = line.write([npartial_waves,self.nform,2*self.j_proj,2*self.j_ejec,2*self.j_target,2*self.j_recoil])
        file.write(line+'\n')
        ################### Line 5 ###################
        line = form.FortranRecordWriter('(F7.4,I5)')
        line = line.write([self.dr,self.nr])
        file.write(line+'\n')


        ################### Line 6-11 ###################
        def omp_params(channel):
            line = form.FortranRecordWriter('(5F7.2,3F7.3)')
            if channel == "Incoming":
                omp_energy = self.beam_energy_lab
                omp_ap = self.a_proj
                omp_zp = self.z_proj
                omp_at = self.a_target
                omp_zt = self.z_target
                omp_radius = float(get("Entrance channel coulomb radius (positive radius will be interpreted as r*At**1/3, negative radius will be interpreted as r*[At**1/3 + Ap**1/3])","1.00"))
                if omp_radius < 0:
                    scale = (omp_at**(1./3.) + omp_ap**(1./3.))/(omp_at**(1./3.))
                    omp_radius = abs(omp_radius)*scale
                    print "\nUsing effective coulomb radius: ", omp_radius

                self.omp_radius = omp_radius
            elif channel == "Outgoing":
                self.ejec_ex = float(get("Enter ejectile state excitation energy (MeV)",self.ejec_ex))
                self.recoil_ex = float(get("Enter recoil state excitation energy (MeV)",self.recoil_ex))
                omp_energy = float(get("Enter reaction Qvalue (outgoing energy): ",  qval_ame12(self.a_recoil,self.z_target,self.z_recoil)+qval_ame12(self.a_proj,self.z_proj,self.z_ejec)-self.recoil_ex-self.ejec_ex))
                omp_ap = self.a_ejec
                omp_zp = self.z_ejec
                omp_at = self.a_recoil
                omp_zt = self.z_recoil
                omp_radius = float(get("Exit channel coulomb radius (positive radius will be interpreted as r*Ar**1/3, negative radius will be interpreted as r*[Ar**1/3 + Ae**1/3])",str(self.omp_radius)))
                if omp_radius < 0:
                    scale = (omp_at**(1./3.) + omp_ap**(1./3.))/(omp_at**(1./3.))
                    omp_radius = abs(omp_radius)*scale
                    print "\nUsing effective coulomb radius: ", omp_radius

            line = line.write([omp_energy,omp_ap,omp_zp,omp_at,omp_zt,omp_radius,1,0])
            file.write(line+'\n')

            pot_kind = int(get(channel+" channel potential option (1=WS, 2=surface WS, 3=second derivative, 15 = read in potential)","1"))
            if pot_kind == 1:
                omp_re_v = float(get("OMP volume (real)","-50.0"))
                omp_re_r = float(get("OMP radius (real) [negative implies R_omp=r*[Ar**1/3 + Ae**1/3]","1.0"))
                if omp_re_r < 0:
                    scale = (omp_at**(1./3.) + omp_ap**(1./3.))/(omp_at**(1./3.))
                    omp_re_r = abs(omp_re_r)*scale
                    print "\nUsing effective (real) OMP radius: ",omp_re_r
                omp_re_a = float(get("OMP diffuseness (real)","0.5"))
                omp_re_ls = float(get("OMP spin orbit (real)","0.0"))
                omp_im_v = float(get("OMP volume (imaginary)","-50.0"))
                omp_im_r = float(get("OMP radius (imaginary) [negative implies R_omp=r*[Ar**1/3 + Ae**1/3]","1.0"))
                if omp_im_r < 0:
                    scale = (omp_at**(1./3.) + omp_ap**(1./3.))/(omp_at**(1./3.))
                    omp_im_r = abs(omp_im_r)*scale
                    print omp_im_r
                    print "Using effective (imaginary) OMP radius: ", omp_im_r
                omp_im_a = float(get("OMP diffuseness (imaginary)","0.5"))
                omp_im_ls = float(get("OMP spin orbit (imaginary)","0.0"))
                line = form.FortranRecordWriter('(F7.4,F7.1,3F7.3,F7.1,4F7.3)')
                line = line.write([pot_kind,omp_re_v,omp_re_r,omp_re_a,omp_re_ls,omp_im_v,omp_im_r,omp_im_a,omp_im_ls,0])
                file.write(line+'\n')
            elif pot_kind == 15:
                pot_path = get("Enter path to OMP file:","./omp.dat")
                line = form.FortranRecordWriter('(4F7.3)')
                line = line.write([15,0,0,0])
                file.write(line+'\n')
                line = form.FortranRecordWriter('(A60)')
                line = line.write([pot_path])
                file.write(line+'\n')
            else:
                print "That potential type still needs to be implemented. It would be really easy to add to this script,"
                print "and copy what is done for type 1 as a template for the other types (2)"
                exit()
            line = form.FortranRecordWriter('(F2.0)')
            line = line.write([0.])
            file.write(line+'\n')
        omp_params("Incoming")
        omp_params("Outgoing")
        for jrjpjt in self.formfactors:
            jr,jp,jt = str(jrjpjt)
            line = form.FortranRecordWriter('(3I3)')
            line = line.write([jr,2*int(jp),2*int(jt)])
            file.write(line+'\n')
            line = form.FortranRecordWriter('(4F7.4)')
            line = line.write([0,0,0,1])
            file.write(line+'\n')
        plot_file = get("Enter path for dwhi plot file:","./dwhi.plot")
        plot_file = padstr(plot_file,60)
        line = form.FortranRecordWriter('(A60)')
        line = line.write([plot_file])
        file.write(line)



if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-target_state", type=int,help="nth state of target/recoil system to use in obd file",default=0)
    args = parser.parse_args()


    init_tab_complete()


    cerxn = CEReactions()
    template_filename = get("If you would like to preload FOLD input settings, enter the path to a template fold input file:","None")
    if template_filename != "" and template_filename != "None":
        manual_entry = get("Template will be loaded, would you like the full version of entry, enter no for the (shorter) abridged version.","yes")
        if manual_entry == "yes" or manual_entry == "y":
            manual_entry = True
        else:
            manual_entry = False
        cerxn.wsaw_inputfile_from_dens()
        cerxn.wsaw_inputfile_from_dens()
        cerxn.fold_inputfile_from_template(template_filename,manual_entry,nstate=args.target_state)
    else:
        cerxn.wsaw_inputfile_from_dens()
        cerxn.wsaw_inputfile_from_dens()
        cerxn.fold_inputfile(nstate=args.target_state)
    cerxn.dwhi_inputfile_gen(nstate=args.target_state)

    save_input_log()
