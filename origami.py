import readline, glob
import fortranformat as form
import numpy as np
import math
from sympy.physics.quantum.cg import CG
from sympy import S

def get(prompt, default):
    return raw_input("%s [%s] " % (prompt, default)) or default

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


def get_obtd(obtd_filename,dj):
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
    mark = 0
    for i,line in enumerate(obtd_list):
        if len(line) == 7 and float(line[0][:-1]) == dj:
            mark = i
    obtds = []
    fold_obtds = []
    for i in range(mark+1,len(obtd_list)):
        if len(obtd_list[i])==7:
            break
        if len(obtd_list[i])>3:
            obtds.append([int(obtd_list[i][0][:-1]),int(obtd_list[i][1][:-1]),float(obtd_list[i][3][:-1])])
            fold_obtds.append([oxbash_to_fold_key_map[int(obtd_list[i][0][:-1])],oxbash_to_fold_key_map[int(obtd_list[i][1][:-1])],float(obtd_list[i][3][:-1])])
    return fold_obtds

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

def init_from_template():
    template_filename = get("Enter a template input file to load","fold.inp")
    filename = get("Enter a name for the FOLD input file to be generated","fold.inp")
    assert(template_filename != filename)
    filelist = []
    for fileline in open(template_filename,"rb"):
        filelist.append(fileline)
    obtd_count_proj, obtd_count_targ = get_number_of_obtds_from_fold_input_file(filelist)

    n=0 # line counter
    ################### Line 1 ###################
    line = form.FortranRecordReader('(I5,I5,A7)')
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
    line = form.FortranRecordReader('(A7)')
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
    line = form.FortranRecordReader('(A7)')
    wsaw_targ = line.read(filelist[n])[0]; n+=1
    ################### Line 6 ###################
    line = form.FortranRecordReader('(F6.3,F9.2,F11.3,A12)')
    fnrm1,fnrm2,blank,blank = line.read(filelist[n]); n+=1
    ################### Line 7 ###################
    line = form.FortranRecordReader('(I5)')
    nform = line.read(filelist[n])[0]; n+=1



    print "\n\n######"
    print "Template file loaded, beginning generation of new input file from template. Press enter for each line and if a change is desired, enter it."
    print "######\n\n"
    file = open(filename,"wb")    
    manual_entry = get("Would you like the full version of entry, enter no for the abridged version.","no")
    if manual_entry == "yes":
        manual_entry = True
    else:
        manual_entry = False


        
      ################### Line 1 ###################        
    output_file = get("Enter a filename for the FOLD output file.\nNote that this filename is restricted to 8 characters or less",output_file)    
    line = form.FortranRecordWriter('(I5,I5,A8)')
    line = line.write([1,1,output_file[0:9]])
    file.write(line+'\n')
    ################### Line 2 ###################
    if manual_entry:
        nr = int(get("Number of integration steps:\n",nr))
        ns = float(get("Step size (fm):\n",ns))
        beam_energy_lab = float(get("Bombarding energy (MeV):\n",beam_energy_lab))
        a = float(get("Projectile mass number (A):\n",a))
    line = form.FortranRecordWriter('(I5,F5.2,F10.0,F10.0,I10,I4,I4)')
    line = line.write([nr,ns,beam_energy_lab,a,1,1,1])
    file.write(line+'\n')
    ################### Line 3 ###################
    #if manual_entry:
    jpf = float(get("Spin of the projectile final state:\n",jpf))
    ppf = str(get("Parity of the projectile final state:\n",ppf))
    jpi = float(get("Spin of the projectile initial state:\n",jpi))
    ppi = str(get("Parity of the projectile initial state:\n",ppi))
    line = form.FortranRecordWriter('(F10.1,A1,F9.1,A1)')
    line = line.write([jpf,ppf,jpi,ppi])
    file.write(line+'\n')
    ################### Line 4 ###################
    #if manual_entry:
    tpf = float(get("Isospin of the projectile final state:\n",tpf))
    tzpf = str(get("Isospin projection of the projectile final state:\n",tzpf.strip()))
    tpi = float(get("Isospin of the projectile initial state:\n",tpi))
    tzpi = str(get("Isospin projection of the projectile initial state:\n",tzpi.strip()))
    line = form.FortranRecordWriter('(F5.1,A10,F10.1,A11)')
    line = line.write([tpf,tzpf,tpi,tzpi])
    file.write(line+'\n')
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
    djp = float(get("Change in spin of projectile: <<<<<\n",djp))
    dtp = float(get("Change in isospin of projectile:\n",1.0))
    obtd_filename_p = str(get("Enter filename/path to OXBASH obtd file","t331dp150.obd"))
    fold_p_obtds = get_obtd(obtd_filename_p,djp) 
    prefactor_p = get_zT_prefactor(jpi,tpi,float(tzpi),dtp,float(tzpf)-float(tzpi),tpf,float(tzpf))   
    print "Projectile prefactor = ",prefactor_p
    for obtd in fold_p_obtds:
        line = form.FortranRecordWriter('(I5,I5,I5,F5.1,F17.6)')
        obtd[2] *= prefactor_p
        line = line.write([obtd[0],obtd[1],int(djp),0.0,obtd[2]])
        file.write(line+'\n')
    line = form.FortranRecordWriter('(I5,I5)')
    line = line.write([-1,-1])
    file.write(line+'\n')
    if manual_entry:
        wsaw_proj = str(get("Enter filename of projectile WSAW radial wavefunction file:\n",wsaw_proj))
    line = form.FortranRecordWriter('(A7)')
    line = line.write([wsaw_proj])
    file.write(line+'\n')
    ################### Target ###################
    ################### Line 3 ###################
    #if manual_entry:
    jtf = float(get("Spin of the target final state:\n",jtf))
    ptf = str(get("Parity of the target final state:\n",ptf))
    jti = float(get("Spin of the target initial state:\n",jti))
    pti = str(get("Parity of the target initial state:\n",pti))
    line = form.FortranRecordWriter('(F10.1,A1,F9.1,A1)')
    line = line.write([jtf,ptf,jti,pti])
    file.write(line+'\n')
    ################### Line 4 ###################
    #if manual_entry:
    ttf = float(get("Isospin of the target final state:\n",ttf))
    tztf = str(get("Isospin projection of the target final state:\n",tztf.strip()))
    tti = float(get("Isospin of the target initial state:\n",tti))
    tzti = str(get("Isospin projection of the target initial state:\n",tzti.strip()))
    line = form.FortranRecordWriter('(F5.1,A10,F10.1,A11)')
    line = line.write([ttf,tztf,tti,tzti])
    file.write(line+'\n')
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
    djt = float(get("Change in spin of target: <<<<<\n",djt))
    dtt = float(get("Change in isospin of target:\n",1.0))
    obtd_filename_t = str(get("Enter filename/path to OXBASH obtd file","t331dp150.obd"))
    fold_t_obtds = get_obtd(obtd_filename_t,djt)
    prefactor_t = get_zT_prefactor(jti,tti,float(tzti),dtt,float(tztf)-float(tzti),ttf,float(tztf))   
    print "Target prefactor = ",prefactor_t
    for obtd in fold_t_obtds:
        line = form.FortranRecordWriter('(I5,I5,I5,F5.1,F17.6)')
        obtd[2] *= prefactor_t
        line = line.write([obtd[0],obtd[1],int(djt),0.0,obtd[2]])
        file.write(line+'\n')
    line = form.FortranRecordWriter('(I5,I5)')
    line = line.write([-1,-1])
    file.write(line+'\n')
    if manual_entry:
        wsaw_targ = str(get("Enter filename of target WSAW radial wavefunction file:\n",wsaw_targ))
    line = form.FortranRecordWriter('(A7)')
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
    line = form.FortranRecordWriter('(I5)')
    line = line.write([nform])
    file.write(line+'\n')
    for i in range(0,nform):
        jrjpjt = str(get("(jr,jp,jt)","011"))
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









def init_from_clean():
    print 1
        
if __name__=="__main__":
    init_tab_complete()
    bool_template = get("Would you like to load a template fold input file (yes/no)?","yes")
    if bool_template=="yes":
        init_from_template()
    else:
        init_from_clean()
        
