#!/bin/env python3

import numpy as np
import openpyxl as pyxl
import re
import math
from pathlib import Path

script_dir = Path(__file__).parent.resolve()

optics='16FEB24'
vfile=['FACET2e_value.echo','FACET2p_value.echo']

outdir='oracle_upload'
xfile='FACET2_'+optics+'.xls'

print(' ')
print('   ===============================================')
print('           FACET2 Excel File Generation')
print('   ===============================================')
print(' ')

stepnum=0;

Er=5.10998918e5     # electron rest mass (eV) ... XAL value
clight=2.99792458e8 # speed of light (m/s) ... XAL value
charge=-1           # sign of electron charge ... XAL value
pi=3.141592653
rad2deg=180/pi      # degrees per radian
T2kG=10             # kG per Tesla

# ==============================================================================
# hardwired LCLS2sc MAD/XAL stuff
# ------------------------------------------------------------------------------

# file name roots

froot = [
    'FACET2e',     #  1
    'FACET2s',    #  2
    'FACET2p',  #  3
]

seq = []
seq.append({'froot':1,'name':'CATHODE TO DL10',     'beg':'BEGINJ',     'end':'ENDINJ',     'offset':[0,0],  'prev':0,  'suml':0,  'length':0})
seq.append({'froot':1,'name':'DL10 TO L1F',         'beg':'BEGDL10',    'end':'ENDDL10',    'offset':[0,0],  'prev':1,  'suml':0,  'length':0})
seq.append({'froot':1,'name':'L1F TO BC11 CHICANE', 'beg':'BEGL1F',     'end':'BC11CBEG',   'offset':[0,-1], 'prev':2,  'suml':0,  'length':0})
seq.append({'froot':1,'name':'BC11 CHICANE',        'beg':'BC11CBEG',   'end':'BC11CEND',   'offset':[0,0],  'prev':3,  'suml':0,  'length':0})
seq.append({'froot':1,'name':'BC11 CHICANE TO L2F', 'beg':'BC11CEND',   'end':'ENDBC11_2',  'offset':[+1,0], 'prev':4,  'suml':0,  'length':0})
seq.append({'froot':1,'name':'L2F TO BC14 CHICANE', 'beg':'BEGL2F',     'end':'BEGBC14E',   'offset':[0,-1], 'prev':5,  'suml':0,  'length':0})
seq.append({'froot':1,'name':'BC14E CHICANE',       'beg':'BEGBC14E',   'end':'ENDBC14E',   'offset':[0,0],  'prev':6,  'suml':0,  'length':0})
seq.append({'froot':1,'name':'BC14 CHICANE TO L3F', 'beg':'ENDBC14E',   'end':'ENDBC14_2',  'offset':[+1,0], 'prev':7,  'suml':0,  'length':0})
seq.append({'froot':1,'name':'L3F TO SCAV',         'beg':'BEGL3F_1',   'end':'ENDL3F_1',   'offset':[0,0],  'prev':8,  'suml':0,  'length':0})
seq.append({'froot':1,'name':'SCAV TO BC20E',       'beg':'BEGL3F_2',   'end':'ENDL3F_2',   'offset':[0,0],  'prev':9,  'suml':0,  'length':0})
seq.append({'froot':1,'name':'BC20E TO FF20',       'beg':'BEGBC20',    'end':'ENDBC20',    'offset':[0,0],  'prev':10, 'suml':0,  'length':0})
seq.append({'froot':1,'name':'FF20 TO EXPT20',      'beg':'BEGFF20',    'end':'ENDFF20',    'offset':[0,0],  'prev':11, 'suml':0,  'length':0})
seq.append({'froot':1,'name':'EXPT20 TO SPECT20',   'beg':'BEGEXPT20',  'end':'ENDEXPT20',  'offset':[0,0],  'prev':12, 'suml':0,  'length':0})
seq.append({'froot':1,'name':'SPECT20 TO DUMP',     'beg':'BEGSPECT20', 'end':'ENDSPECT20', 'offset':[0,0],  'prev':13, 'suml':0,  'length':0})
% FACET2 scavenger e-
seq.append({'froot':2,'name':'SCAV TO TARGET',      'beg':'BEGSCAV',    'end':'ENDSCAV',    'offset':[0,0],  'prev':9,  'suml':0,  'length':0})
% FACET2 e+
seq.append({'froot':3,'name':'BC14P CHICANE',       'beg':'BEGBC14P',   'end':'ENDBC14P',   'offset':[0,0],  'prev':6,  'suml':0,  'length':0})

seqexcl = []

# ------------------------------------------------------------------------------
# machine areas
area = []
area.append({'name':'INJ',     'beg':'BEGINJ',     'end':'L0AFBEG',    'offset':[0,-1]})
area.append({'name':'L0F',     'beg':'L0AFBEG',    'end':'ENDINJ',     'offset':[0,0]})
area.append({'name':'DL10',    'beg':'BEGDL10',    'end':'ENDDL10',    'offset':[0,0]})
area.append({'name':'L1F',     'beg':'BEGL1F',     'end':'ENDL1F',     'offset':[0,0]})
area.append({'name':'BC11_1',  'beg':'BEGBC11_1',  'end':'ENDBC11_1',  'offset':[0,0]})
area.append({'name':'BC11_2',  'beg':'BEGBC11_2',  'end':'ENDBC11_2',  'offset':[0,0]})
area.append({'name':'L2F',     'beg':'BEGL2F',     'end':'ENDL2F',     'offset':[0,0]})
area.append({'name':'BC14_1',  'beg':'BEGBC14_1',  'end':'ENDBC14_1',  'offset':[0,0]})
area.append({'name':'BC14E',   'beg':'BEGBC14E',   'end':'ENDBC14E',   'offset':[0,0]})
area.append({'name':'BC14_2',  'beg':'BEGBC14_2',  'end':'ENDBC14_2',  'offset':[0,0]})
area.append({'name':'L3F_1',   'beg':'BEGL3F_1',   'end':'ENDL3F_2',   'offset':[0,0]})
area.append({'name':'L3F_2',   'beg':'BEGL3F_2',   'end':'ENDL3F_2',   'offset':[0,0]})
area.append({'name':'BC20',    'beg':'BEGBC20',    'end':'ENDBC20',    'offset':[0,0]})
area.append({'name':'FF20',    'beg':'BEGFF20',    'end':'ENDFF20',    'offset':[0,0]})
area.append({'name':'EXPT20',  'beg':'BEGEXPT20',  'end':'ENDEXPT20',  'offset':[0,0]})
area.append({'name':'SPECT20', 'beg':'BEGSPECT20', 'end':'ENDSPECT20', 'offset':[0,0]})
area.append({'name':'SCAV',    'beg':'BEGSCAV',    'end':'ENDSCAV',    'offset':[0,0]})
area.append({'name':'BC14P',   'beg':'BEGBC14P',   'end':'ENDBC14P',   'offset':[0,0]})

areaexcl = []

from xtffs2mat import xtffs2mat

# ------------------------------------------------------------------------------
# read the MAD output files
stepnum += 1
print('   {}) Read MAD output files ...'.format(stepnum))

K, N, L, P, A, T, E, FDN, coor, S = [], [], [], [], [], [], [], [], [], []
idf, ids, idd = [], [], []  # idf: which MAD output file an element came from
                            # ids: which XAL sequence an element belongs to
                            # idd: ordinal position in MAD output file
nf = 0

for n in range(len(seq)):
		if n in seqexcl:
			continue
    if seq[n]['froot'] != nf:
        nf = seq[n]['froot']
        fname = '{}_survey.tape'.format(froot[nf-1])
        print('Opening file {}'.format(fname))
        titl, tK, tN, tL, tP, tA, tT, tE, tFDN, tcoor, tS = xtffs2mat(fname)

    id1 = tN.index(seq[n]['beg']) + seq[n]['offset'][0]
    id2 = tN.index(seq[n]['end']) + seq[n]['offset'][1]
    id_range = slice(id1, id2 + 1)
    
    K.extend(tK[id_range])
    len_id = len(tK[id_range])
    N.extend(tN[id_range])
    L.extend(tL[id_range])
    P.extend(tP[id_range])
    A.extend(tA[id_range])
    T.extend(tT[id_range])
    E.extend(tE[id_range])
    FDN.extend(tFDN[id_range])
    coor.extend(tcoor[id_range])
    S.extend(tS[id_range])
    idf.extend([nf] * len_id)
    ids.extend([n] * len_id)
    idd.extend(list(range(id1,id2+1)))

Nelem, Nc = len(N), len(N[0])

# set "display S" to linac Z-coordinate
Sd = [x[2] for x in coor]

def intersection(x,y):
    return [v for v in x if v in y]

def strmatch(n_str,N_lst,exact=False):
    if exact:
        return [ix for ix,n_ in enumerate(N_lst) if n_.strip() == n_str.strip()]
    else:
        return [ix for ix,n_ in enumerate(N_lst) if n_.startswith(n_str)]

P2 = []

# initialize some sequence data

for s in seq:
    id1 = N.index(s['beg']) + s['offset'][0]
    id2 = N.index(s['end']) + s['offset'][1]
    s['suml'] = S[id1]
    s['length']=S[id2]-S[id1]

# assign machine areas

ida=[]
for ix,a in enumerate(area):
    id1 = N.index(a['beg']) + a['offset'][0]
    id2 = N.index(a['end']) + a['offset'][1]
    ida.extend([ix]*(id2-id1+1))

# assign area "parent" names

for a in area:
    if a['name'] in ['BC11_1','BC11_2']:
        a['parent']='BD11'
    elif a['name'] in ['BC14_1','BC14E','BC14P','BC14_2']:
        a['parent']='BC14'
    elif a['name'] in ['L3F_1','L3F_2']:
        a['parent']='L3F'
    else:
        a['parent']=a['name']

FOO MPE HERE
# special handling for rolled dump lines and A-line

def fix_dump_coords(N, P, coor):
    # Implementation of FixDumpCoords function

    # set roll angle for SXR dump line components
    id1=N.index('RODMP1S')
    id2=N.index('RODMP2S')-1
    ARODMP1S=P[id1][4]
    for i in range(id1,id2+1):
      coor[i][5]=ARODMP1S

    id1=id2+1
    id2=N.index('ENDDMPS_2')
    for i in range(id1,id2+1):
      coor[i][5]=0

    # set roll angle for HXR dump line components
    id1=N.index('RODMP1H')
    id2=N.index('RODMP2H')-1
    ARODMP1H=P[id1][4];
    for i in range(id1,id2+1):
      coor[i][5]=ARODMP1H

    id1=id2+1;
    id2=N.index('ENDDMPH_2')
    for i in range(id1,id2+1):
      coor[i][5]=0

    return coor

def fix_aline_coords(N, P, coor):
    # Implementation of FixAlineCoords function
    # set roll angle for A-line components
    id1 = N.index('ROLL2')
    id2 = N.index('ENDBSYA_2')
    AROLL2 = P[id1][4]
    for i in range(id1,id2+1):
      coor[i][5] = AROLL2

    # The following block is commented out in the original code
    '''
    id1 = N.index('BEGBSYA_1')
    id2 = N.index('ROLL2') - 1
    id_slice = slice(id1, id2 + 1)
    coor[id_slice, 5] = 0  # remove residual "creeping" roll
    '''
    return coor

def fix_sxtes_coords(N, coor):
    # Implementation of FixSXTESCoords function
    # fix BSY coordinates for selected SXTES system devices per P. Stephens
    name = [
        'MR1K3_VGC_1', 'ND1S', 'SP1K1_MONO_VGC_1',  # 2.2 line
        'IM1K3_PPM', 'BT1K3_AIR',  # TXI line
        'BT2K0_PLEG_TMO', 'LUSI'  # TMO line
    ]
    coor_id = [
        1, 1, 1,
        0, 0,
        -1, 0
    ]
    coor_val = [
        -0.8826040, -2.0921000, -0.7249275,
        1.0694435, 1.0480923,
        1.2500000, -1.2194000
    ]

    for n in range(len(name)):
        if coor_id[n] == -1:
            continue
        id_matches = N.index(name[n])
        coor[id_matches][coor_id[n]] = coor_val[n]
    return coor

coor = fix_dump_coords(N, P, coor)
coor = fix_aline_coords(N, P, coor)

# kicker/septum groups

KSname = [
    'BKRDG0', 'BLRDG0',
    'BKYSP0H', 'BKYSP1H', 'BKYSP2H', 'BKYSP3H', 'BKYSP4H', 'BKYSP5H', 'BLXSPH',
    'BKYSP0S', 'BKYSP1S', 'BKYSP2S', 'BKYSP3S', 'BKYSP4S', 'BKYSP5S', 'BLXSPS',
    'BKRDAS1', 'BKRDAS2', 'BKRDAS3', 'BKRDAS4', 'BKRDAS5', 'BKRDAS6', 'BLRDAS',
    'BKRCUS', 'BLRCUS'
]

# read FINT values for SBENs and undulator parameters from a special echo-file
# generated via MAD VALUE commands

C = []
for n in range(len(vfile)):
    fname = vfile[n]
    with open(fname, 'r') as f:
        C.extend(f.read().split())

P2 = np.zeros((Nelem, 2))

idb = [i for i,x in enumerate(K) if x == 'SBEN']
for m in range(0, len(idb), 2):
    na = idb[m]
    nb = idb[m+1]
    name = N[na].strip()
    name = name.split('.')[0]  # remove decoration, if any
    id_ = strmatch(name,C)[0]
    #id_ = [i for i, x in enumerate(C) if name in x][0]
    #if not id_:
    #    raise ValueError(f'No FINT for {name}')
    #elif len(id_) > 1:
    #    print('oops')
    fint = float(C[id_+6])
    P2[na][0] = fint
    P2[nb][0] = fint

idm = [i for i,x in enumerate(K) if x == 'MATR']
for m in range(0, len(idm), 2):
    n1 = idm[m]
    n2 = idm[m+1]
    name = N[n1].strip()
    Ktxt = f'"{name}_K"'
    Ltxt = f'"{name}_L"'
    idK = strmatch(Ktxt,C)[0]
    idL = strmatch(Ltxt,C)[0]
    undk = float(C[idK+2])
    undl = float(C[idL+2])
    P2[n1, :] = [undl, undk]
    P2[n2, :] = [undl, undk]

# make unique names

nfix = [['', 'MUQS', 'MPHS'],
        ['', 'MUQH', 'MPHH'],
        ['HOMCM', '', '']]
ioff = [[-3, -1, -1], [-2, -1, -1], [-3, -1, -1]]

for nr in range(len(nfix)):
    for nc in range(len(nfix[0])):
        if not nfix[nr][nc]:
            continue
        id = strmatch(nfix[nr][nc], N)
        if len(id) == 0:
            continue
        if nr < 2 and nc == 0:
            id = id[::2]
        lc = len(nfix[nr][nc])
        for m in range(len(id)):
            if nr < 2:
                cname = N[id[m] + ioff[nr][nc]].strip()
                N[id[m]] = N[id[m]][:lc] + '.' + cname[-2:]
            else:
                N[id[m]] = N[id[m]][:lc] + '.' + N[id[m]-2][2:4]

# Find indices of 'WOODDOOR' in N
jd = [i for i,x in enumerate(N) if x == 'WOODDOOR']
for j in jd:
    cname = 'WOODDOOR.{}'.format(area[ida[j]]['name'])
    N[j] = cname

# Shared devices (devices which see both kicked and unkicked beams)
aname_all = ['DIAG0', 'SPH', 'SPS', 'DASEL', 'CLTS']
name_all = ['BPMDG000', 'BPMSPH', 'BPMSPS', 'BPMDAS', 'BPMCUS']

for name, aname in zip(name_all,aname_all):
    jd = strmatch(name,N,True)
    for m in range(len(jd)):
        if aname == area[ida[jd[m]]]['name']:
            N[jd[m]] = name + '?'

# Change keyword for XTES mirrors from MULT to INST
name = ['MR1K1_BEND', 'SP1K1_MONO', 'MR3K1_GRATING',
        'MR1K3_TXI', 'MR2K3_TXI', 'MR1K4_SOMS',
        'MR1L0_HOMS_XTES', 'MR2L0_HOMS', 'MR1L0_HOMS_TXI', 'MR1L1_TXI']

for n in name:
    id_ = strmatch(n,N,True)
    for i in id_:
        K[i] = 'INST'

# Change keyword for GUN/GUNB solenoid correction quads from MULT to QUAD;
# copy T1 into TILT slot
name = ['CQ01', 'SQ01', 'CQ01B', 'SQ01B', 'SQ02B']

for n in name:
    id_ = strmatch(n,N,True)
    for i in id_:
        K[i] = 'QUAD'
        P[i][3] = P[i][5]  # T1 -> TILT

def assign_ucell(N, FDN, coor, idf):
    # NOTE: coordinates are assumed to be in MAD (not SYMBOLS) order

    debug = False
    filename = f'{script_dir}/sectors.xlsx'

    wb= pyxl.load_workbook(filename,data_only=True)

    # Read SXR sheet
    sxr_sht = wb.worksheets[2]
    sxr_data = sxr_sht['A2':'E36']
    ncell = len(sxr_data)
    ucell = []
    for row in sxr_data:
        ucell.append({'name':row[0].value,'froot':row[1].value,'Zbeg':row[3].value,'Zend':row[4].value})

    # Read HXR sheet
    hxr_sht = wb.worksheets[3]
    hxr_data = hxr_sht['A2':'E39']
    ncell = len(hxr_data)
    for row in hxr_data:
        ucell.append({'name':row[0].value,'froot':row[1].value,'Zbeg':row[3].value,'Zend':row[4].value})
    
    wb.close()
    # Combine data from both sheets
    
    coor3 = [c[2] for c in coor]
    for cell in ucell:
        id1 = [i for i,x in enumerate(idf) if x==cell['froot']]

        id2 = [i for i,x in enumerate(coor3) if x > cell['Zbeg']]
        inter1 = intersection(id1,id2)[0]

        id2 = [i for i,x in enumerate(coor3) if x < cell['Zend']]
        inter2 = [v for v in id1 if v in id2][-1]
        for jd in range(inter1,inter2+1):
            if not FDN[jd]:
                FDN[jd] = cell['name']

    return FDN

# Assign sector names (use empty FDN and FDN1 arrays)

# MAD SURVEY coordinates  [x,y,z,theta,phi   ,psi]
# correspond to SolidEdge [z,x,y,roll ,-pitch,yaw]

ic = [2, 0, 1, 5, 4, 3]
coor = np.array(coor)
coor = coor[:, ic]
coor[:, 4] = -coor[:, 4]

# generate ordered list of MAD keywords

tkeyw = sorted(list(set(K)))

MADK = [
    'LCAV', 'SBEN', 'QUAD', 'SEXT', 'SOLE', 'MATR', 'RCOL', 'ECOL', 'SROT',
    'HKIC', 'VKIC', 'MONI', 'WIRE', 'PROF', 'IMON', 'BLMO', 'INST',
    'MARK'
]
OUTK = [
    'LCAV', 'BEND', 'QUAD', 'SEXT', 'SOLE', 'USEG', 'COLL', 'PC  ', 'MARK',
    'XCOR', 'YCOR', 'BPM ', 'WIRE', 'PROF', 'IMON', 'BLMO', 'INST',
    'MARK'
]
XALK = [
    'BNCH', 'BEND', 'QUAD', 'SEXT', 'SOLE', 'USEG', 'COLL', 'ECOL', 'SROT',
    'XCOR', 'YCOR', 'BPM ', 'WIRE', 'PROF', 'TORO', 'BLMO', 'INST',
    'MARK'
]
idmisc = list(range(9, 16))  # non-magnet elements that might have nonzero lengths
keyw = []
keyo = []
keyx = []
jdmisc = []
for n in range(len(MADK)):
    id_ = strmatch(MADK[n], tkeyw)
    if len(id_) > 0:
        keyw.append(MADK[n])
        keyo.append(OUTK[n])
        keyx.append(XALK[n])
        if n in idmisc:
            jdmisc.append(len(keyw) - 1)
idmisc=jdmisc

# hard-wired list of bends that have energy polynomials in the database
Ebend = []  # ['BRB', 'BXSP', 'BYSP', 'BRSP', 'BYSP', 'BX3', 'BY1', 'BY2', 'BYD']

# Keyword Worksheets
# (NOTE: suml, energy, and coord are quoted at element's center)
# ==============================================================================
# LCAV (unsegmented):
# - idf,id,sequence,area,xkey,prim,name,type,dist(m),energy(GeV)
# - length(m),freq(MHz),ampl(MeV),phase(deg),grad(MeV/m),power(1)
# - sdsp(m),suml(m),coord(m,rad)
# - suml1(m),coord1(m,rad),suml2(m),coord2(m,rad)
# SBEN (unsplit):
# - idf,id,sequence,area,xkey,prim,name,type,dist(m),energy(GeV)
# - zleng(m),leng(m),gap(m),fint(1),tilt(deg),ang(deg),e1(deg),e2(deg),
#   BL(kG-m),B(T),k1(1/m^2),GL(kG),G(T/m),scale(name,value),polarity
# - sdsp(m),suml(m),coord(m,rad),mcoord(m,rad)
# - suml1(m),coord1(m,rad),mcoord1(m,rad),suml2(m),coord2(m,rad),mcoord2(m,rad)
# QUAD (unsplit):
# - idf,id,sequence,area,xkey,prim,name,type,dist(m),energy(GeV)
# - leng(m),bore(m),tilt(deg),k1(1/m^2),GL(kG),G(T/m),scale(name,value),polarity
# - sdsp(m),suml(m),coord(m,rad),mcoord(m,rad)
# - suml1(m),coord1(m,rad),suml2(m),coord2(m,rad)
# SEXT (unsplit):
# - idf,id,sequence,area,xkey,prim,name,type,dist(m),energy(GeV)
# - leng(m),bore(m),tilt(deg),k2(1/m^3),G'L(kG/m),G'(T/m^2),scale(name,value),
#   polarity
# - sdsp(m),suml(m),coord(m,rad)
# - suml1(m),coord1(m,rad),suml2(m),coord2(m,rad)
# SOLE (unsplit):
# - idf,id,sequence,area,xkey,prim,name,type,dist(m),energy(GeV)
# - leng(m),bore(m),ks(1/m),BL(kG-m),B(T),scale(name,value),polarity
# - sdsp(m),suml(m),coord(m,rad)
# - suml1(m),coord1(m,rad),suml2(m),coord2(m,rad)
# MATR (unsplit):
# - idf,id,sequence,area,xkey,prim,name,type,dist(m),energy(GeV)
# - leng(m),lambda(m),k(1)
# - sdsp(m),suml(m),coord(m,rad)
# - suml1(m),coord1(m,rad),suml2(m),coord2(m,rad)
# RCOL,ECOL:
# - idf,id,sequence,area,xkey,prim,name,type,dist(m),energy(GeV)
# - leng(m),xsize(m),ysize(m)
# - sdsp(m),suml(m),coord(m,rad)
# - suml1(m),coord1(m,rad),suml2(m),coord2(m,rad)
# SROT:
# - idf,id,sequence,area,xkey,prim,name,type,dist(m),energy(GeV)
# - leng(m),ang(deg)
# - sdsp(m),suml(m),coord(m,rad)
# - suml1(m),coord1(m,rad),suml2(m),coord2(m,rad)
# HKIC,VKIC,MONI,WIRE,PROF,IMON,BLMO,INST:
# - idf,id,sequence,area,xkey,prim,name,type,dist(m),energy(GeV)
# - leng(m)
# - sdsp(m),suml(m),coord(m,rad)
# - suml1(m),coord1(m,rad),suml2(m),coord2(m,rad)
# MARK:
# - idf,id,sequence,area,xkey,prim,name,type,dist(m),energy(GeV)
# - sdsp(m),suml(m),coord(m,rad)
# - suml1(m),coord1(m,rad),suml2(m),coord2(m,rad)

# Assumptions
# ==============================================================================
# - LCAVs may be segmented
# - the first 6 characters of an LCAV's name associate it with it's parent
# - all SBENs are split into two pieces
# - the last character of an SBEN's name differentiates it's pieces
# - all QUADs, SEXTs, and MATRs are split in half
# - all other keyword types are not (necessarily) split
# - SBENs with abs(ang)<amin will have ang set to zero
# - SBENs or QUADs with abs(k1)<kmin will have k1 set to zero
# - SOLEs with abs(ks)<kmin will have ks set to zero
# - SROTs with abs(ang)<amin will have ang set to zero

amin = 1e-9
kmin = 1e-6

# process by keyword

nHKIC = 0
nVKIC = 0
nMONI = 0
nWIRE = 0
nPROF = 0
nIMON = 0
nBLMO = 0

def slicer(N,ix):
    return [N[i] for i in ix]

#HEAD
Sd = np.array(Sd)
S = np.array(S)
L = np.array(L)
P = np.array(P)
P2 = np.array(P2)
E = np.array(E)
seq = np.array(seq)
A = np.array(A)
for kwn,kxn,kon in zip(keyw,keyx,keyo):
    id = strmatch(kwn,K)
    if kwn == 'LCAV':
        # create list of unique names that will allow unsplitting
        name = slicer(N,id)
        for i in range(len(name)):
            if name[i][0:4] in ['CAVL', 'CAVC']:  # unique in 7 characters
                name[i] = name[i][0:7]
            else:  # unique in 6 characters
                name[i] = name[i][0:6]
        name = list(dict.fromkeys(name))
        LCAV = []
        nLCAV = 0

        for mname in name:
            if mname.startswith('TCX'):
                id = strmatch(mname,N,True)
            else:
                id = strmatch(mname,N)
            id1 = id[0]  # first segment
            ide = [id1-1, id[-1]]  # [entrance, exit]
            sdsp = np.mean(Sd[ide])  # m (beam center)
            suml = np.mean(S[ide])  # m (beam center)
            dist = suml - seq[ids[id1]]['suml']  # m (sequence start to beam center)
            energy = np.mean(E[ide])  # GeV (beam center)
            leng = np.sum(L[id])  # m
            freq = P[id1, 4]  # MHz
            ampl = np.sum(P[id, 5])  # MeV
            phase = P[id1, 6]  # rad/2pi
            grad = ampl / leng  # MeV/m
            if re.match(r'K\d\d_\d[ABCD]', mname[0:6]):  # i.e. K27_3D
                id = strmatch(mname[0:5],N)
                grad0 = np.min(P[id, 5] / L[id])
                if grad0 == 0:
                    power = float("NaN")
                else:
                    power = 0.25 * round((grad / grad0) ** 2)  # KLYS power fraction (1)
            else:
                power = 1
            coorc = np.mean(coor[ide, :], axis=0)  # m,rad (beam center)
            nLCAV += 1
            LCAV.append({
                'idf': idf[id1],
                'id': idd[id1],
                'seq': seq[ids[id1]]['name'],
                'area': area[ida[id1]]['name'],
                'parent': area[ida[id1]]['parent'],
                'sector': FDN[id1].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon,
                'name': mname,
                'type': T[id1].strip(),
                'dist': dist,
                'energy': energy,
                'leng': leng,
                'freq': freq,
                'ampl': ampl,
                'phase': 360 * phase,  # deg
                'grad': grad,
                'power': power,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)}
            })

            # BSY coordinates

            LCAV[-1]['suml1'] = []
            for k in range(6):
                LCAV[-1][f'c1{k+1}'] = []
            # UND coordinates

            LCAV[-1]['suml2'] = []
            for k in range(6):
                LCAV[-1][f'c2{k+1}'] = []
    elif kwn == 'SBEN':
        name = [N[i] for i in id]
        SBEN = []
        Nelm = len(id) // 2
        for m in range(Nelm):
            mname1 = name[2 * m].strip()
            mname2 = name[2 * m + 1].strip()
            mname = mname1[:-1]  # remove last character from name
            id = [strmatch(mname1,N,True),strmatch(mname2,N,True)]
            id1 = id[0][0]  # first piece (beam center)
            idi = id1 - 1  # beam in
            ido = id[1][-1]  # beam out
            sdsp = Sd[id1]  # m
            suml = S[id1]  # m
            dist = suml - seq[ids[id1]]['suml']  # m
            energy = E[id1]  # GeV
            leng = np.sum(L[id])  # m
            gap = 2 * A[id1]  # m
            fint = P2[id1, 0]  # m
            tilt = P[id1, 3]  # rad
            ang = np.sum(P[id, 0])  # rad
            if abs(ang) < amin:
                ang = 0
                e1 = 0
                e2 = 0
            else:
                e1 = P[id1, 4]  # rad
                e2 = P[ido, 5]  # rad
            EeV = 1e9 * energy  # eV
            brho = np.sqrt(EeV ** 2 - Er ** 2) / clight  # T-m
            BL = brho * ang  # T-m
            B = BL / leng  # T
            k1 = P[id1, 1]  # 1/m^2
            if abs(k1) < kmin:
                k1 = 0
            G = brho * k1  # T/m
            GL = G * leng  # T
            if mname[:3] in Ebend:
                sname = 'GeV2T'
                sval = brho * abs(ang) / (leng * energy)
            else:
                sname = 'kG2T_Bdl2B'
                sval = 1 / (leng * T2kG)
            polarity = -np.sign(ang + np.finfo(float).eps)  # add eps so that sign=1 when ang=0
            coori = np.copy(coor[idi, :])  # m,rad
            coorc = np.copy(coor[id1, :])  # m,rad
            cooro = np.copy(coor[ido, :])  # m,rad
            coorm = np.zeros(coorc.shape)  # m,rad (magnet steel center)
            if mname in KSname:
                jd = strmatch(f'D{mname}',N)
                if len(jd) != 2:
                    raise ValueError(f'{mname} not split?')
                zleng = np.sum(L[jd])  # m
                coorm = np.copy(coor[jd[0], :])  
                coorm[3] = 0
            else:
                chicane1 = (e1 == 0) & (e2 != 0)
                chicane2 = (e1 != 0) & (e2 == 0)
                if chicane1 | chicane2:
                    zleng = leng * np.sinc(ang/np.pi)  # m
                    coorm[:3] = np.mean([coori[:3], cooro[:3]], axis=0)
                    if chicane1:
                        coorm[3:6] = np.copy(coori[3:6])
                    else:
                        coorm[3:6] = np.copy(cooro[3:6])
                else:
                    zleng = leng * np.sinc(ang / 2 / np.pi)  # m
                    coorm[:3] = (coori[:3] + cooro[:3] + 2 * coorc[:3]) / 4
                    coorm[3:6] = np.copy(coorc[3:6])
            pname = area[ida[id1]]['parent']
            if pname in ['DMPS', 'DMPH']:
                coorm[3] = coorc[3]  # dump line magnet coords set in FixDumpCoords
            elif pname == 'BSYA_2':
                coorm[3] = coorc[3]  # dump line magnet coords set in FixAlineCoords
            else:
                coorm[3] = tilt  # remove "creeping" rolls from non-rolled SBENs
            SBEN.append({
                'idf': idf[id1],
                'id': idd[id1],
                'seq': seq[ids[id1]]['name'],
                'area': area[ida[id1]]['name'],
                'parent': pname,
                'sector': FDN[id1].strip(),
                'ucell': [],
                'xkey': f'X{kxn}' if tilt == 0 else f'Y{kxn}' if abs(np.cos(tilt)) < 1e-9 else f'R{kxn}',
                'prim': 'KICK' if mname.startswith('BK') or 'KIK' in mname else 'SEPT' if mname.startswith('BL') else kon,
                'name': mname,
                'type': T[id1].strip(),
                'dist': dist,
                'energy': energy,
                'zleng': zleng,
                'leng': leng,
                'gap': 2 * A[id1],  # m
                'fint': fint,
                'tilt': np.rad2deg(tilt),  # deg
                'ang': np.rad2deg(ang),  # deg
                'e1': np.rad2deg(e1),  # deg
                'e2': np.rad2deg(e2),  # deg
                'BL': T2kG * BL,  # kG-m
                'B': charge * B,
                'k1': k1,
                'GL': T2kG * GL,  # kG
                'G': charge * G,
                'sname': sname,
                'sval': sval,
                'polarity': polarity,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)},
                **{f'm{k+1}': coorm[k] for k in range(6)}
            })
            # BSY coordinates

            SBEN[-1]['suml1'] = []
            for k in range(6):
                SBEN[-1][f'c1{k+1}'] = []
                SBEN[-1][f'm1{k+1}'] = []

            # UND coordinates

            SBEN[-1]['suml2'] = []
            for k in range(6):
                SBEN[-1][f'c2{k+1}'] = []
                SBEN[-1][f'm2{k+1}'] = []

    elif kwn == 'QUAD':
        name = list(dict.fromkeys([N[i] for i in id]))
        QUAD = []
        for mname in name:
            id = strmatch(mname,N,True)
            id1 = id[0]  # first segment (beam center)
            sdsp = Sd[id1]  # m
            suml = S[id1]  # m
            dist = suml - seq[ids[id1]]['suml']  # m
            energy = E[id1]  # GeV
            leng = np.sum(L[id])  # m
            bore = 2 * A[id1]  # m
            tilt = P[id1, 3]  # rad
            k1 = P[id1, 1]  # 1/m^2
            if abs(k1) < kmin:
                k1 = 0
            EeV = 1e9 * energy  # eV
            brho = np.sqrt(EeV ** 2 - Er ** 2) / clight  # T-m
            if leng == 0:
                G = 0  # T/m
                GL = brho * k1  # T
                sname = 'kG2T'
                sval = 1 / T2kG
            else:
                G = brho * k1  # T/m
                GL = G * leng  # T
                sname = 'kG2T_Gdl2G'
                sval = 1 / (leng * T2kG)
            polarity = -np.sign(k1 + np.finfo(float).eps)  # add eps so that sign=1 when k1=0
            coorc = np.copy(coor[id1, :])  # m,rad
            QUAD.append({
                'idf': idf[id1],
                'id': idd[id1],
                'seq': seq[ids[id1]]['name'],
                'area': area[ida[id1]]['name'],
                'parent': area[ida[id1]]['parent'],
                'sector': FDN[id1].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon,
                'name': mname,
                'type': T[id1].strip(),
                'dist': dist,
                'energy': energy,
                'leng': leng,
                'bore': bore,
                'tilt': np.rad2deg(tilt),  # deg
                'k1': k1,
                'GL': T2kG * GL,  # kG
                'G': charge * G,
                'sname': sname,
                'sval': sval,
                'polarity': polarity,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)},
                **{f'm{k+1}': [] for k in range(6)}
            })

            # BSY coordinates

            QUAD[-1]['suml1'] = []
            for k in range(6):
                QUAD[-1][f'c1{k+1}'] = []

            # UND coordinates

            QUAD[-1]['suml2'] = []
            for k in range(6):
                QUAD[-1][f'c2{k+1}'] = []

    elif kwn == 'SEXT':
        name = list(dict.fromkeys([N[i] for i in id]))
        SEXT = []
        for mname in name:
            id = strmatch(mname,N,True)
            id1 = id[0]  # first half (beam center)
            sdsp = Sd[id1]  # m
            suml = S[id1]  # m
            dist = suml - seq[ids[id1]]['suml']  # m
            energy = E[id1]  # GeV
            leng = np.sum(L[id])  # m
            bore = 2 * A[id1]  # m
            tilt = P[id1, 3]  # rad
            k2 = P[id1, 2]  # 1/m^3
            if abs(k2) < kmin:
                k2 = 0
            EeV = 1e9 * energy  # eV
            brho = np.sqrt(EeV ** 2 - Er ** 2) / clight  # T-m
            Gp = brho * k2  # T/m
            GpL = Gp * leng  # T
            if leng == 0:
                sname = 'kG2T'
                sval = 1 / T2kG
            else:
                sname = 'kG2T_Gpdl2Gp'
                sval = 1 / (leng * T2kG)
            polarity = -np.sign(k2 + np.finfo(float).eps)  # add eps so that sign=1 when k2=0
            coorc = np.copy(coor[id1, :])  # m,rad
            SEXT.append({
                'idf': idf[id1],
                'id': idd[id1],
                'seq': seq[ids[id1]]['name'],
                'area': area[ida[id1]]['name'],
                'parent': area[ida[id1]]['parent'],
                'sector': FDN[id1].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon,
                'name': mname,
                'type': T[id1].strip(),
                'dist': dist,
                'energy': energy,
                'leng': leng,
                'bore': bore,
                'tilt': np.rad2deg(tilt),  # deg
                'k2': k2,
                'GpL': T2kG * GpL,  # kG/m
                'Gp': charge * Gp,
                'sname': sname,
                'sval': sval,
                'polarity': polarity,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)}
            })

            # BSY coordinates

            SEXT[-1]['suml1'] = []
            for k in range(6):
                SEXT[-1][f'c1{k+1}'] = []

            # UND coordinates

            SEXT[-1]['suml2'] = []
            for k in range(6):
                SEXT[-1][f'c2{k+1}'] = []
    elif kwn == 'SOLE':
        name = list(dict.fromkeys([N[i] for i in id]))
        SOLE = []
        for mname in name:
            id = strmatch(mname,N,True)
            id1 = id[0]
            ide = [id1 - 1, id[-1]]  # [entrance, exit]
            sdsp = np.mean(Sd[ide])  # m
            suml = np.mean(S[ide])  # m
            dist = suml - seq[ids[id1]]['suml']  # m
            energy = np.mean(E[ide])  # GeV
            leng = np.sum(L[id])  # m
            bore = 2 * A[id1]  # m
            ks = P[id1, 4]  # 1/m
            if abs(ks) < kmin:
                ks = 0
            EeV = 1e9 * energy  # eV
            brho = np.sqrt(EeV**2 - Er**2) / clight  # T-m
            B = brho * ks  # T
            BL = B * leng  # T-m
            if leng == 0:
                sname = 'kG2T'
                sval = 1 / T2kG
            else:
                sname = 'kG2T_Bdl2B'
                sval = 1 / (leng * T2kG)
            polarity = -np.sign(ks + np.finfo(float).eps)  # add eps so that sign=1 when ks=0
            coorc = np.mean(coor[ide], axis=0)  # m, rad
            SOLE.append({
                'idf': idf[id1],
                'id': idd[id1],
                'seq': seq[ids[id1]]['name'],
                'area': area[ida[id1]]['name'],
                'parent': area[ida[id1]]['parent'],
                'sector': FDN[id1].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon,
                'name': mname,
                'type': T[id1].strip(),
                'dist': dist,
                'energy': energy,
                'leng': leng,
                'bore': bore,
                'ks': ks,
                'BL': T2kG * BL,  # kG-m
                'B': charge * B,
                'sname': sname,
                'sval': sval,
                'polarity': polarity,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)}
            })

            # BSY coordinates

            SOLE[-1]['suml1'] = []
            for k in range(6):
                SOLE[-1][f'c1{k+1}'] = []

            # UND coordinates

            SOLE[-1]['suml2'] = []
            for k in range(6):
                SOLE[-1][f'c2{k+1}'] = []

    elif kwn == 'MATR':
        name = list(dict.fromkeys([N[i] for i in id]))
        MATR = []
        for mname in name:
            id = strmatch(mname,N,True)
            id1 = id[0]  # first half (beam center)
            sdsp = Sd[id1]  # m
            suml = S[id1]  # m
            dist = suml - seq[ids[id1]]['suml']  # m
            energy = E[id1]  # GeV
            leng = np.sum(L[id])  # m
            undl = P2[id1, 0]  # m
            undk = P2[id1, 1]  # 1
            coorc = np.copy(coor[id1])  # m, rad
            MATR.append({
                'idf': idf[id1],
                'id': idd[id1],
                'seq': seq[ids[id1]]['name'],
                'area': area[ida[id1]]['name'],
                'parent': area[ida[id1]]['parent'],
                'sector': FDN[id1].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon,
                'name': mname,
                'type': T[id1].strip(),
                'dist': dist,
                'energy': energy,
                'leng': leng,
                'lambda': undl,
                'k': undk,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)}
            })

            # BSY coordinates

            MATR[-1]['suml1'] = []
            for k in range(6):
                MATR[-1][f'c1{k+1}'] = []
            # UND coordinates

            MATR[-1]['suml2'] = []
            for k in range(6):
                MATR[-1][f'c2{k+1}'] = []
    elif kwn == 'RCOL':
        name = [N[i] for i in id] # RCOLs are not split
        RCOL = []
        for mname in name:
            id = strmatch(mname,N,True)[0]
            ide = [id - 1, id]  # [entrance, exit]
            sdsp = np.mean(Sd[ide])  # m (beam center)
            suml = np.mean(S[ide])  # m (beam center)
            dist = suml - seq[ids[id]]['suml']  # m (sequence start to beam center)
            energy = E[id]  # GeV
            leng = L[id]  # m
            xgap = 2 * P[id, 3]  # m
            ygap = 2 * P[id, 4]  # m
            coorc = np.mean(coor[ide], axis=0)  # m, rad (beam center)
            RCOL.append({
                'idf': idf[id],
                'id': idd[id],
                'seq': seq[ids[id]]['name'],
                'area': area[ida[id]]['name'],
                'parent': area[ida[id]]['parent'],
                'sector': FDN[id].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon,
                'name': mname,
                'type': T[id].strip(),
                'dist': dist,
                'energy': energy,
                'leng': leng,
                'xgap': xgap,
                'ygap': ygap,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)}
            })

            # BSY coordinates

            RCOL[-1]['suml1'] = []
            for k in range(6):
                RCOL[-1][f'c1{k+1}'] = []
            # UND coordinates

            RCOL[-1]['suml2'] = []
            for k in range(6):
                RCOL[-1][f'c2{k+1}'] = []
    elif kwn == 'ECOL':
        name = [N[i] for i in id]  # ECOLs are not split
        ECOL = []
        for mname in name:
            id = strmatch(mname,N,True)[0]
            ide = [id - 1, id]  # [entrance, exit]
            sdsp = np.mean(Sd[ide])  # m (beam center)
            suml = np.mean(S[ide])  # m (beam center)
            dist = suml - seq[ids[id]]['suml']  # m (sequence start to beam center)
            energy = E[id]  # GeV
            leng = L[id]  # m
            xbore = 2 * P[id, 3]  # m
            ybore = 2 * P[id, 4]  # m
            coorc = np.mean(coor[ide], axis=0)  # m, rad (beam center)
            ECOL.append({
                'idf': idf[id],
                'id': idd[id],
                'seq': seq[ids[id]]['name'],
                'area': area[ida[id]]['name'],
                'parent': area[ida[id]]['parent'],
                'sector': FDN[id].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon,
                'name': mname,
                'type': T[id].strip(),
                'dist': dist,
                'energy': energy,
                'leng': leng,
                'xbore': xbore,
                'ybore': ybore,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)}
            })

            # BSY coordinates

            ECOL[-1]['suml1'] = []
            for k in range(6):
                ECOL[-1][f'c1{k+1}'] = []
            # UND coordinates

            ECOL[-1]['suml2'] = []
            for k in range(6):
                ECOL[-1][f'c2{k+1}'] = []
    elif kwn == 'SROT':
        name = [N[i] for i in id]  # SROTs are not split
        SROT = []
        for mname in name:
            id = strmatch(mname,N,True)[0]
            ide = [id - 1, id]  # [entrance, exit]
            sdsp = np.mean(Sd[ide])  # m (beam center)
            suml = np.mean(S[ide])  # m (beam center)
            dist = suml - seq[ids[id]]['suml']  # m (sequence start to beam center)
            energy = E[id]  # GeV
            leng = L[id]  # m
            ang = np.rad2deg(P[id, 4])  # deg
            if abs(ang) < amin:
                ang = 0
            coorc = np.mean(coor[ide], axis=0)  # m, rad (beam center)
            SROT.append({
                'idf': idf[id],
                'id': idd[id],
                'seq': seq[ids[id]]['name'],
                'area': area[ida[id]]['name'],
                'parent': area[ida[id]]['parent'],
                'sector': FDN[id].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon,
                'name': mname,
                'type': T[id].strip(),
                'dist': dist,
                'energy': energy,
                'leng': leng,
                'ang': ang,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)}
            })

            # BSY coordinates

            SROT[-1]['suml1'] = []
            for k in range(6):
                SROT[-1][f'c1{k+1}'] = []
            # UND coordinates

            SROT[-1]['suml2'] = []
            for k in range(6):
                SROT[-1][f'c2{k+1}'] = []
    elif kwn == 'INST':
        name = list(dict.fromkeys([N[i] for i in id]))
        INST = []
        for mname in name:
            id = strmatch(mname,N,True)
            if id[0] == 1:
                idi = 1
            else:
                idi = id[0] - 1  # beam in
            ide = [idi, id[-1]]  # [entrance, exit]
            sdsp = np.mean(Sd[ide])  # m (beam center)
            suml = np.mean(S[ide])  # m (beam center)
            dist = suml - seq[ids[id[0]]]['suml']  # m (sequence start to beam center)
            energy = E[id[0]]  # GeV
            leng = np.sum(L[id])  # m
            coorc = np.mean(coor[ide], axis=0)  # m, rad (beam center)
            t = T[id[0]].strip()
            if t == 'type-V':
                kon_use = 'EFC'
            elif t == '0.79K11.8':
                kon_use = 'BEND'
            else:
                kon_use = kon
            INST.append({
                'idf': idf[id[0]],
                'id': idd[id[0]],
                'seq': seq[ids[id[0]]]['name'],
                'area': area[ida[id[0]]]['name'],
                'parent': area[ida[id[0]]]['parent'],
                'sector': FDN[id[0]].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon_use,
                'name': mname,
                'type': T[id[0]].strip(),
                'dist': dist,
                'energy': energy,
                'leng': leng,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)},
                **{f'm{k+1}': [] for k in range(6)}
            })

            # BSY coordinates

            INST[-1]['suml1'] = []
            for k in range(6):
                INST[-1][f'c1{k+1}'] = []
                INST[-1][f'm1{k+1}'] = []
            # UND coordinates

            INST[-1]['suml2'] = []
            for k in range(6):
                INST[-1][f'c2{k+1}'] = []
                INST[-1][f'm2{k+1}'] = []  # for INST
    elif kwn in [keyw[i] for i in idmisc]: #MISC
        name = list(dict.fromkeys([N[i] for i in id]))
        for mname in name:
            id = strmatch(mname,N,True)
            if id[0] == 1:
                idi = 1
            else:
                idi = id[0] - 1  # beam in
            ide = [idi, id[-1]]  # [entrance, exit]
            sdsp = np.mean(Sd[ide])  # m (beam center)
            suml = np.mean(S[ide])  # m (beam center)
            dist = suml - seq[ids[id[0]]]['suml']  # m (sequence start to beam center)
            energy = E[id[0]]  # GeV
            leng = np.sum(L[id])  # m
            coorc = np.mean(coor[ide], axis=0)  # m, rad (beam center)
            t = T[id[0]].strip()
            if t == 'type-V':
                kon_use = 'EFC'
            else:
                kon_use = kon
            MISC = {
                'idf': idf[id[0]],
                'id': idd[id[0]],
                'seq': seq[ids[id[0]]]['name'],
                'area': area[ida[id[0]]]['name'],
                'parent': area[ida[id[0]]]['parent'],
                'sector': FDN[id[0]].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon_use,
                'name': mname,
                'type': T[id[0]].strip(),
                'dist': dist,
                'energy': energy,
                'leng': leng,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)}
            }

            # BSY coordinates

            MISC['suml1'] = []
            for k in range(6):
                MISC[f'c1{k+1}'] = []
            MISC['suml2'] = []
            for k in range(6):
                MISC[f'c2{k+1}'] = []
            globals()[f'n{kwn}'] += 1
            if f'{kwn}' in globals():
                globals()[f'{kwn}'].append(MISC)
            else:
                globals()[f'{kwn}'] = [MISC]

    elif kwn == 'MARK':
        name = [N[i] for i in id]
        MARK = []
        for mname in name:
            id = strmatch(mname,N,True)[0]
            sdsp = Sd[id]  # m
            suml = S[id]  # m
            dist = suml - seq[ids[id]]['suml']  # m (sequence start to beam center)
            energy = E[id]  # GeV
            coorc = np.copy(coor[id])
            MARK.append({
                'idf': idf[id],
                'id': idd[id],
                'seq': seq[ids[id]]['name'],
                'area': area[ida[id]]['name'],
                'parent': area[ida[id]]['parent'],
                'sector': FDN[id].strip(),
                'ucell': [],
                'xkey': kxn,
                'prim': kon,
                'name': mname,
                'type': T[id].strip(),
                'dist': dist,
                'energy': energy,
                'sdsp': sdsp,
                'suml': suml,
                **{f'c{k+1}': coorc[k] for k in range(6)}
            })

            # BSY coordinates

            MARK[-1]['suml1'] = []
            for k in range(6):
                MARK[-1][f'c1{k+1}'] = []
            # UND coordinates

            MARK[-1]['suml2'] = []
            for k in range(6):
                MARK[-1][f'c2{k+1}'] = []
import scipy.io as sio

def fix_power_fraction(lcav):
    #fname = r'V:\LCLS\Users\Woodley\AD_ACCEL\20190613_13JUN19\RDB\RDBdata'
    fname = f'{script_dir}/RDBdata.mat'
    old = sio.loadmat(fname)['LCAV']
    
    name = [x['name'] for x in old[0]]
    powr = [float(x['power'][0][0]) for x in old[0]]
    
    for cav in lcav:
        if np.isnan(cav['power']):
            id = [i for i, x in enumerate(name) if x == cav['name']]
            if len(id) != 1:
                raise ValueError(f"{cav['name']} not found!")
            cav['power'] = powr[id[0]]
    
    return lcav


# fix LCAV power fraction values (for deactivated klystrons)
LCAV = fix_power_fraction(LCAV)

def add_eic(inst):
    # Add EIC Faraday cup
    name = [x['name'] for x in inst]
    id = name.index('CATHODEB')
    temp = inst[id].copy()
    temp.id = 59
    temp.seq = 'CATHODE TO DIAG0'
    temp.area = 'EIC'
    temp.name = 'FC00EIC'
    temp.type = 'Faraday cup'
    temp.dist = 3.0
    temp.sdsp = -7.044667
    temp.suml = 3.0
    temp.c1 = -7.044667
    inst.append(temp)
    return inst

# deferred devices

DEPR = []
nDEPR = 0
pname = [x['parent'] for x in area]

for k_str in keyw:
    k = globals()[k_str]
    for m in range(len(k)):
        t = k[m]['type']
        if t and t[0] == '@':
            deplev = int(t[1])
            if len(t) > 2:
                t = t[3:]  # skip ","
            else:
                t = ''
            nDEPR += 1
            id = strmatch(k[m]['parent'],pname)
            if k == 'SBEN':
                z_use = k[m]['m1']
            else:
                z_use = k[m]['c1']
            DEPR.append({
                'id':k[m]['id'],
                'parent':k[m]['parent'],
                'ida': id[0]+1,
                'prim': k[m]['prim'],
                'name': k[m]['name'],
                'z': z_use,
                'type': t,
                'level': deplev,
                })

            k[m]['parent'] = '*' + k[m]['parent']
            k[m]['area'] = '*' + k[m]['area']
            k[m]['type'] = t

# ------------------------------------------------------------------------------
# Fix magnet coordinates ...
# ------------------------------------------------------------------------------

def FixMagnetCoords(SBEN, QUAD, INST, K, N, L, P, coor, cflag):
    # Set special magnet coordinates for:
    # - R56 compensation chicanes
    # - self-seeding chicanes and Cavity-Based-XFEL (CBXFEL) chicanes
    # - safety dump bends
    # - Lambertson septa
    # - rolled vertical bends, kickers, and septa
    # - spreader kickers
    # - CUSXR extraction magnets
    # - QDG001 and QDG003
    # - SXRSS optical components
    #
    #   cflag : []=linac, 1=BSY, 2=UND
    #
    # NOTE: see the MAD User's Reference Manual (v8.19), Section 1.3
    if cflag is None:
        n1 = 0  # all chicanes have linac coordinates
    else:
        n1 = 2  # CCDLU and CCDLD chicanes do not have BSY or UND coordinates
    Xgun = 0.28
    Ygun = -0.99
    Bname = [i['name'] for i in SBEN]
    Qname = [i['name'] for i in QUAD]
    Iname = [i['name'] for i in INST]

    # R56 compensation chicanes
    # coor=[z,x,y,roll,-pitch,yaw] (SYMBOLS coordinates)

    R56name = ['CCDLU', 'CCDLD', 'CC31B', 'CC32B', 'CC31', 'CC32', 'CC35', 'CC36']
    off = 0.005  # 5 mm offset
    for name in R56name[n1:]:
        id1 = strmatch(f"{name}BEG",N)[0]
        id2 = strmatch(f"{name}END",N)[0]
        id = range(id1, id2 + 1)
        jd1 = [i for i,x in enumerate(K) if x == 'SBEN']
        idb = intersection(id,jd1)[::2]
        ang = P[idb[0], 0]
        X, Y, Z, yaw, pitch, roll = coor[id, 1], coor[id, 2], coor[id, 0], coor[id1, 5], -coor[id1, 4], coor[id1, 3]
        O1 = np.array([[np.cos(yaw), 0, np.sin(yaw)], [0, 1, 0], [-np.sin(yaw), 0, np.cos(yaw)]])
        O2 = np.array([[1, 0, 0], [0, np.cos(pitch), np.sin(pitch)], [0, -np.sin(pitch), np.cos(pitch)]])
        O3 = np.array([[np.cos(roll), -np.sin(roll), 0], [np.sin(roll), np.cos(roll), 0], [0, 0, 1]])
        O = O1 @ O2 @ O3
        t = np.linalg.solve(O, np.array([X, Y, Z]))  # remove yaw, pitch, and roll
        Xr, Yr, Zr = t[0], t[1], t[2]
        dX = -off * np.sign(ang + np.finfo(float).eps)  # offset in chicane direction
        Xr = Xr[0] + dX * np.ones_like(Xr)
        t = O @ np.array([Xr, Yr, Zr])  # restore roll, pitch, and yaw
        Xm, Ym, Zm = t[0], t[1], t[2]
        for m in range(len(idb)):
            name = N[idb[m]].strip()[:-1] #remove last character
            jdb = strmatch(name,[N[i] for i in id])[0]
            jd = strmatch(name,Bname)[0]
            SBEN[jd][f'm{cflag or ""}2'] = Xm[jdb]
            SBEN[jd][f'm{cflag or ""}3'] = Ym[jdb]
            SBEN[jd][f'm{cflag or ""}1'] = Zm[jdb]

    # self-seeding chicane bends and Cavity-Based-XFEL bends

    name = ['BCXHS1', 'BCXHS2', 'BCXHS3', 'BCXHS4',  # HXRSS self-seeding chicane
            'BCXSS1', 'BCXSS2', 'BCXSS3', 'BCXSS4',  # SXRSS self-seeding chicane
            'BCXXL1', 'BCXXL2', 'BCXXL3', 'BCXXL4',  # XLEAP-II self-seeding chicane
            'BCXCBX11', 'BCXCBX12', 'BCXCBX13', 'BCXCBX14',  # CBXFEL chicane #1
            'BCXCBX21', 'BCXCBX22', 'BCXCBX23', 'BCXCBX24']  # CBXFEL chicane #2
    dX = 1e-3 * np.array([0, -2.39, -2.39, 0,  # HXRSS
                          +1, +9.7, +9.7, +1,  # SXRSS
                          -5, -12, -12, -5,  # XLEAP-II
                          +1, +9.7, +9.7, +1,  # CBXFEL #1
                          +1, +9.7, +9.7, +1])  # CBXFEL #2
    for n in range(len(name)):
        id = strmatch(name[n],Bname,True)[0]
        X0 = SBEN[id][f'm{cflag or ""}2']
        X = X0 + dX[n]
        SBEN[id][f'm{cflag or ""}2'] = X

    # safety dump bends (permanent magnet dipoles)

    name = ['BXPM1B', 'BXPM1', 'BXPM2']
    for n in range(len(name)):
        if n == 0:  # SXR
            Xm = 1.25
        else:  # HXR
            Xm = -1.215
        id1 = strmatch(f"{name[n]}1",N)[0] #center
        id0 = id1 - 1  # entrance
        if n != 2:
            pitch = -coor[id0, 4]
            z0 = coor[id0, 0]
            y0 = coor[id0, 2]
        z1 = coor[id1, 0]
        Ym = y0 + np.tan(pitch) * (z1 - z0)
        yaw = 0
        id = strmatch(name[n],Bname,True)[0]
        SBEN[id][f'm{cflag or ""}2'] = Xm
        SBEN[id][f'm{cflag or ""}3'] = Ym
        SBEN[id][f'm{cflag or ""}6'] = yaw
        SBEN[id][f'm{cflag or ""}5'] = -pitch

    # Lambertson septa
    # coor=[z,x,y,roll,-pitch,yaw] (SYMBOLS coordinates)

    name = ['BLRDG0', 'BLXSPS', 'BLXSPH', 'BLRDAS', 'BLRCUS']
    r = 0.010  # radius of field-free channel
    off = -0.004  # beam is 6 mm from top of field-free channel
    for n in range(len(name)):
        if n <= 1 and cflag is not None:
            continue  # no BSY or UND coords for BLRDG0 or BLRL3X
        id = strmatch(name[n],Bname,True)[0]
        Xm0 = SBEN[id][f'm{cflag or ""}2']
        Ym0 = SBEN[id][f'm{cflag or ""}3']
        Zm0 = SBEN[id][f'm{cflag or ""}1']
        yaw = SBEN[id][f'm{cflag or ""}6']
        pitch = -SBEN[id][f'm{cflag or ""}5']
        roll = (np.pi / 180) * SBEN[id]['tilt']
        O1 = np.array([[np.cos(yaw), 0, np.sin(yaw)], [0, 1, 0], [-np.sin(yaw), 0, np.cos(yaw)]])
        O2 = np.array([[1, 0, 0], [0, np.cos(pitch), np.sin(pitch)], [0, -np.sin(pitch), np.cos(pitch)]])
        O3 = np.array([[np.cos(roll), -np.sin(roll), 0], [np.sin(roll), np.cos(roll), 0], [0, 0, 1]])
        O = O1 @ O2 @ O3
        t = np.linalg.solve(O, np.array([Xm0, Ym0, Zm0]))  # remove roll, pitch, and yaw
        Xr, Yr, Zr = t[0], t[1], t[2]
        Yr = Yr + off  # apply vertical offset
        t = O @ np.array([Xr, Yr, Zr])  # restore roll, pitch, and yaw
        Xm, Ym, Zm = t[0], t[1], t[2]
        # CheckMagnetCoords
        SBEN[id][f'm{cflag or ""}2'] = Xm
        SBEN[id][f'm{cflag or ""}3'] = Ym
        SBEN[id][f'm{cflag or ""}1'] = Zm

    # set magnet roll for bends, kickers, and septa
    # coor=[z,x,y,roll,-pitch,yaw] (SYMBOLS coordinates)

    for n in range(len(Bname)):
        mname = Bname[n].strip()
        if mname in ['BKRDG0',  # DIAG0 kicker
                     'BKRCUS',  # cuS kicker
                     'BKRDAS1', 'BKRDAS2', 'BKRDAS3', 'BKRDAS4', 'BKRDAS5', 'BKRDAS6',  # DASEL kickers
                     'WIG1S', 'WIG2S', 'WIG3S', 'WIG1H', 'WIG2H', 'WIG3H']:  # SLC-style wigglers
            roll = (np.pi / 180) * (SBEN[n]['tilt'] - 90)  # rolled vertical bends and kickers
        elif mname in ['BYDSS', 'BYD1B', 'BYD2B', 'BYD3B',  # SXR dump line
                       'BYDSH', 'BYD1', 'BYD2', 'BYD3']:  # HXR dump line
            continue  # see FixDumpCoords
        elif mname in ['B11', 'B12', 'B13', 'B14', 'B15', 'B16',
                       'B21', 'B22', 'B23', 'B24', 'B25', 'B26']:  # A-line
            continue  # see FixAlineCoords
        else:
            if 'BKY' in mname:
                roll = (np.pi / 180) * (SBEN[n]['tilt'] - 90)  # non-rolled vertical kickers
            else:
                roll = (np.pi / 180) * SBEN[n]['tilt']  # other bends, kickers, or septa
        SBEN[n][f'm{cflag or ""}4'] = roll

    # spreader kickers
    # coor=[z,x,y,roll,-pitch,yaw] (SYMBOLS coordinates)

    name = ['BKYSP0H', 'BKYSP1H', 'BKYSP2H', 'BKYSP3H', 'BKYSP4H', 'BKYSP5H',
            'BKYSP0S', 'BKYSP1S', 'BKYSP2S', 'BKYSP3S', 'BKYSP4S', 'BKYSP5S']
    off = [0.0002, 0.0004]  # to mitigate resistive wall wakefield effects
    for n in range(len(name)):
        if n <= 5:
            yoff = off[0]
        else:
            yoff = off[1]
        id = strmatch(name[n],Bname,True)[0]
        Ym0 = SBEN[id][f'm{cflag or ""}3']
        SBEN[id][f'm{cflag or ""}3'] = Ym0 + yoff

    # CUSXR extraction magnets
    # coor=[z,x,y,roll,-pitch,yaw] (SYMBOLS coordinates)

    name = ['BRCUSDC1', 'BKRCUS', 'BRCUSDC2']
    for namen in name:
        dname = f'D{namen}A'
        idd = strmatch(dname,N,True)[0]
        id = strmatch(namen,Bname,True)[0]
        for m in [1,2,3,5,6]:
            SBEN[id][f'm{cflag or ""}{m}'] = coor[idd, m-1]

    # QDG001 and QDG003

    if cflag is None:  # only linac coordinates
        name = ['QDG001', 'QDG003']
        for n in range(len(name)):
            id = strmatch(name[n], N)
            KL = np.sum(P[id, 1] * L[id])
            id = strmatch(f'DY{name[n]}', N)[0]
            kick = P[id, 0]
            off = kick / KL
            id = strmatch(name[n],Qname,True)[0]
            Xm0 = QUAD[id]['c2']
            Ym0 = QUAD[id]['c3']
            Zm0 = QUAD[id]['c1']
            yaw = QUAD[id]['c6']
            pitch = -QUAD[id]['c5']
            roll = QUAD[id]['c4']
            O1 = np.array([[np.cos(yaw), 0, np.sin(yaw)], [0, 1, 0], [-np.sin(yaw), 0, np.cos(yaw)]])
            O2 = np.array([[1, 0, 0], [0, np.cos(pitch), np.sin(pitch)], [0, -np.sin(pitch), np.cos(pitch)]])
            O3 = np.array([[np.cos(roll), -np.sin(roll), 0], [np.sin(roll), np.cos(roll), 0], [0, 0, 1]])
            O = O1 @ O2 @ O3
            t = np.linalg.solve(O, np.array([Xm0, Ym0, Zm0]))  # remove roll, pitch, and yaw
            Xr, Yr, Zr = t[0], t[1], t[2]
            Yr = Yr + off  # apply vertical offset
            t = O @ np.array([Xr, Yr, Zr])  # restore roll, pitch, and yaw
            Xm, Ym, Zm = t[0], t[1], t[2]
            # CheckMagnetCoords
            QUAD[id]['m2'] = Xm
            QUAD[id]['m3'] = Ym
            QUAD[id]['m1'] = Zm
            QUAD[id]['m6'] = yaw
            QUAD[id]['m5'] = -pitch
            QUAD[id]['m4'] = roll

    # SXRSS optical components

    name = ['GSXS1', 'MSXS1', 'SLSXS1', 'MSXS2', 'MSXS3']
    dX = 1e-3 * np.array([0, -1.93, -3.85, -3.85, 0])
    for n in range(len(name)):
        id = strmatch(name[n],Iname,True)[0]
        X0 = INST[id][f'c{cflag or ""}2']
        Y0 = INST[id][f'c{cflag or ""}3']
        Z0 = INST[id][f'c{cflag or ""}1']
        yaw = INST[id][f'c{cflag or ""}6']
        pitch = -INST[id][f'c{cflag or ""}5']
        roll = INST[id][f'c{cflag or ""}4']
        O1 = np.array([[np.cos(yaw), 0, np.sin(yaw)], [0, 1, 0], [-np.sin(yaw), 0, np.cos(yaw)]])
        O2 = np.array([[1, 0, 0], [0, np.cos(pitch), np.sin(pitch)], [0, -np.sin(pitch), np.cos(pitch)]])
        O3 = np.array([[np.cos(roll), -np.sin(roll), 0], [np.sin(roll), np.cos(roll), 0], [0, 0, 1]])
        O = O1 @ O2 @ O3
        t = np.linalg.solve(O, np.array([X0, Y0, Z0]))  # remove roll, pitch, and yaw
        Xr, Yr, Zr = t[0], t[1], t[2]
        Xr = Xr + dX[n]  # apply horizontal offset
        t = O @ np.array([Xr, Yr, Zr])  # restore roll, pitch, and yaw
        Xm, Ym, Zm = t[0], t[1], t[2]
        INST[id][f'm{cflag or ""}2'] = Xm
        INST[id][f'm{cflag or ""}3'] = Ym
        INST[id][f'm{cflag or ""}1'] = Zm
        INST[id][f'm{cflag or ""}6'] = yaw
        INST[id][f'm{cflag or ""}5'] = -pitch
        INST[id][f'm{cflag or ""}4'] = roll

    return SBEN, QUAD,INST

SBEN, QUAD, INST = FixMagnetCoords(SBEN, QUAD, INST, K, N, L, P, coor, None)

# ------------------------------------------------------------------------------

def fix_mark_prim(mark):
    # change MARK().prim from MARK to INST for MARKERs that define the beginning
    # and end of Areas, Sectors, and Bypass Line Sectors ... per K. Luchini
    mname = [m['name'] for m in mark]
    idabeg = [i for i, name in enumerate(mname) if name.startswith('BEG')]
    idaend = [i for i, name in enumerate(mname) if name.startswith('END')]
    idsbeg = [i for i, name in enumerate(mname) if name.startswith('LI') and name.endswith('BEG')]
    idsend = [i for i, name in enumerate(mname) if name.startswith('LI') and name.endswith('END')]
    idbbeg = [i for i, name in enumerate(mname) if name.startswith('BPN') and name.endswith('BEG')]
    idbend = [i for i, name in enumerate(mname) if name.startswith('BPN') and name.endswith('END')]
    id = idabeg + idaend + idsbeg + idsend + idbbeg + idbend

    for n in id:
        prim = mark[n]['prim']
        mark[n]['prim'] = 'INST'

    return mark


# change MARK().prim='MARK' to MARK().prim='INST' for selected MARKer elements
MARK = fix_mark_prim(MARK)

# ------------------------------------------------------------------------------

# Worksheets
# ==========
# one sheet per keyword
# one sheet for segments
# one sheet for sequences
# one sheet for deferred devices

# common worksheet header

head1 = ['MAD #', 'Sequence', 'Area', 'Sector', 'Undulator Cell', 'XAL Keyword',
         'DB Keyword', 'MAD Name', 'Engineering Type', 'SeqDist', 'Energy']
head2 = ['Display S', 'SumL (linac)', 'X Coor (linac)', 'Y Coor (linac)',
         'Z Coor (linac)', 'X Angle (linac)', 'Y Angle (linac)', 'Z Angle (linac)']
head3 = ['SumL (BSY)', 'X Coor (BSY)', 'Y Coor (BSY)', 'Z Coor (BSY)',
         'X Angle (BSY)', 'Y Angle (BSY)', 'Z Angle (BSY)']
head4 = ['SumL (UND)', 'X Coor (UND)', 'Y Coor (UND)', 'Z Coor (UND)',
         'X Angle (UND)', 'Y Angle (UND)', 'Z Angle (UND)']
foot1 = head1
foot2 = ['Display S', 'SumL (linac)', 'MAD Z (linac)', 'MAD X (linac)',
         'MAD Y (linac)', 'MAD Psi (linac)', 'MAD Phi (linac)', 'MAD Theta (linac)']
foot3 = ['SumL (BSY)', 'MAD Z (BSY)', 'MAD X (BSY)', 'MAD Y (BSY)',
         'MAD Psi (BSY)', 'MAD Phi (BSY)', 'MAD Theta (BSY)']
foot4 = ['SumL (UND)', 'MAD Z (UND)', 'MAD X (UND)', 'MAD Y (UND)',
         'MAD Psi (UND)', 'MAD Phi (UND)', 'MAD Theta (UND)']
unit1 = ['', '', '', '', '', '',
         '', '', '', 'm', 'GeV']
unit2 = ['m', 'm', 'm', 'm', 'm', 'rad', 'rad', 'rad']
unit3 = ['m', 'm', 'm', 'm', 'rad', 'rad', 'rad']
unit4 = ['m', 'm', 'm', 'm', 'rad', 'rad', 'rad']


# Precision for coordinate output
prec = 1e-6

wb = pyxl.Workbook()
for keywn in keyw:
    thead2, tfoot2, tunit2 = head2.copy(), foot2.copy(), unit2.copy()
    thead3, tfoot3, tunit3 = head3.copy(), foot3.copy(), unit3.copy()
    thead4, tfoot4, tunit4 = head4.copy(), foot4.copy(), unit4.copy()
    
    if keywn == 'LCAV':
        khead = ['Length', 'Frequency', 'Amplitude', 'Phase', 'Gradient', 'Power']
        kunit = ['m', 'MHz', 'MeV', 'degree', 'MeV/m', 'fraction']
    elif keywn == 'SBEN':
        khead = ['Z Length', 'Effective Length', 'Gap', 'Field Integral', 'Tilt',
                 'Angle', 'E1', 'E2', 'BL', 'B', 'K1', 'GL', 'G',
                 'XAL Scale Name', 'XAL Scale Value', 'XAL Polarity']
        kunit = ['m', 'm', 'm', '', 'degree',
                 'degree', 'degree', 'degree', 'kG-m', 'T', '1/m^2', 'kG', 'T/m',
                 '', '', '']
        thead2.extend(['Magnet X Coor (linac)', 'Magnet Y Coor (linac)',
                       'Magnet Z Coor (linac)', 'Magnet X Angle (linac)',
                       'Magnet Y Angle (linac)', 'Magnet Z Angle (linac)'])
        tfoot2.extend(['Magnet MAD Z (linac)', 'Magnet MAD X (linac)',
                       'Magnet MAD Y (linac)', 'Magnet MAD Psi (linac)',
                       'Magnet MAD Phi (linac)', 'Magnet MAD Theta (linac)'])
        tunit2.extend(['m', 'm', 'm', 'rad', 'rad', 'rad'])
        # Similar extensions for thead3, tfoot3, tunit3, thead4, tfoot4, tunit4
    # ... (similar cases for other keywords)
    
    kfoot = khead.copy()
    
    # Write to Excel file
    ws = wb.create_sheet(keywn)
    data = head1+khead+thead2+thead3+thead4
    for col, value in enumerate(data,start=1):
        ws.cell(row=1,column=col,value=value)
    data = unit1+kunit+tunit2+tunit3+tunit4
    for col, value in enumerate(data,start=1):
        ws.cell(row=2,column=col,value=value)
    
    # Process TEMP data
    M = []
    TEMP=globals()[keywn]
    for TEMPm in TEMP:
        for j in range(1, 7):
            TEMPm[f'c{j}'] = np.round(TEMPm[f'c{j}'], decimals=int(-np.log10(prec)))
            TEMPm[f'c1{j}'] = np.round(TEMPm[f'c1{j}'], decimals=int(-np.log10(prec)))
            TEMPm[f'c2{j}'] = np.round(TEMPm[f'c2{j}'], decimals=int(-np.log10(prec)))
            if TEMPm['prim'] == 'BEND' and TEMPm['type'] != '0.79K11.8':
                TEMPm[f'm{j}'] = np.round(TEMPm[f'm{j}'], decimals=int(-np.log10(prec)))
                TEMPm[f'm1{j}'] = np.round(TEMPm[f'm1{j}'], decimals=int(-np.log10(prec)))
                TEMPm[f'm2{j}'] = np.round(TEMPm[f'm2{j}'], decimals=int(-np.log10(prec)))
            elif TEMPm['prim'] == 'QUAD':
                TEMPm[f'm{j}'] = np.round(TEMPm[f'm{j}'], decimals=int(-np.log10(prec)))
            elif TEMPm['prim'] == 'INST':
                if TEMPm['xkey'] != 'MARK':
                    TEMPm[f'm{j}'] = np.round(TEMPm[f'm{j}'], decimals=int(-np.log10(prec)))
                    TEMPm[f'm1{j}'] = np.round(TEMPm[f'm1{j}'], decimals=int(-np.log10(prec)))
                    TEMPm[f'm2{j}'] = np.round(TEMPm[f'm2{j}'], decimals=int(-np.log10(prec)))
        
        if not (TEMPm['idf'] in [3, 4, 5, 12, 13] and TEMPm['name'].startswith('TEMP')):
            TEMPwr = {k: v for k, v in TEMPm.items() if k not in ['idf', 'area']}
            M.append(list(TEMPwr.values()))
    
    M.append(foot1 + kfoot + tfoot2 + tfoot3 + tfoot4)
    M.append(unit1 + kunit + tunit2 + tunit3 + tunit4)

    for i,datarow in enumerate(M,3):
        for j,data in enumerate(datarow,1):
            if j == 1 and isinstance(data,int):
                data = data + 1
            if isinstance(data,list):
                if data != []:
                    ws.cell(row=i,column=j,value=data)
            elif isinstance(data,np.ndarray):
                if data.size > 0:
                    ws.cell(row=i,column=j,value=data[0])
            else:
                ws.cell(row=i,column=j,value=data)



# ------------------------------------------------------------------------------
# Write Sequences Worksheet ...

# Sequences header

head = ['Sequence', 'Previous Sequence', 'Begin Element', 'End Element', 'Length']
unit = ['', '', '', '', 'm']

ws_seq = wb.create_sheet('Sequences')

for i,data in enumerate(head,1):
    ws_seq.cell(row=1,column=i,value=data)
for i,data in enumerate(unit,1):
    ws_seq.cell(row=2,column=i,value=data)

# Sequences Worksheet

M = []
for s in seq:
    temp = {
        'name': s['name'],
        'prev': '',
        'beg': s['beg'],
        'end': s['end'],
        'leng': s['length']
    }
    m = s['prev']
    if m > 0:
        temp['prev'] = seq[m-1]['name']
    M.append(list(temp.values()))

for i,row in enumerate(M,3):
    for j,data in enumerate(row,1):
        ws_seq.cell(row=i,column=j,value=data)

# ------------------------------------------------------------------------------
# Deferred devices header

head = ['MAD #', 'Area Name', 'Area Id', 'DB Keyword', 'MAD Name', 'Z',
        'Engineering Type', 'Level']

ws_def = wb.create_sheet('Deferred')
for i,data in enumerate(head,1):
    ws_def.cell(row=1,column=i,value=data)

# Deferred devices Worksheet

M = [list(DEPR[n].values()) for n in range(nDEPR)]
for i,row in enumerate(M,2):
    for j,data in enumerate(row,1):
        if j == 1 and isinstance(data,int):
            data = data + 1
        ws_def.cell(row=i,column=j,value=data)

Path(outdir).mkdir(parents=True,exist_ok=True)
wb.save(outdir+'/'+xfile)
# ------------------------------------------------------------------------------
# Write SYMBOLS txt-files ...

# set up pointers

seqname = [x['name'] for x in seq]
ip = []
for n in range(len(keyw)):
    TEMP = globals()[keyw[n]]
    for m in range(len(TEMP)):
        id = strmatch(TEMP[m]['seq'],seqname,True)[0]
        ip.append([seq[id]['froot'], n, m, TEMP[m]['id']])
ip = sorted(ip, key=lambda x: (x[0], x[3]))


# SYMBOLS text-file headers and footers

head = ('Solid Edge,AREA,KeyW,ELEMENT,Eng_Name,L_EFF,'
        'APER,ANGLE,K1,K2,TILT,E1,E2,H1,H2,ENERGY,'
        'SUML,X Coor,Y Coor,Z Coor,X Angle,Y Angle,Z Angle,'
        'RF_Frequency,RF_Amplitude,RF_Phase,RF_Gradient,RF_Power_Fraction,'
        'Z_Length,Fringe_Field_Integral,Integrated_Field_BL,Field_B,'
        'Integrated_Field_Gradient_GL,Field_Gradient_G,'
        'XAL_Scale_Name,XAL_Scale_Value,XAL_Polarity,'
        'Magnet_X_Coor,Magnet_Y_Coor,Magnet_Z_Coor,'
        'Magnet_X_Angle,Magnet_Y_Angle,Magnet_Z_Angle,'
        'Solenoid_Strength_KS,Undulator_Period_Length,Undulator_Strength_K,'
        'X_Size,Y_Size,'
        'Section,Distance_From_Section_Start,XAL_Keyword,S_Display')

foot = ('MAD #,AREA,KeyW,ELEMENT,Eng_Name,L_EFF,'
        'APER,ANGLE,K1,K2,TILT,E1,E2,H1,H2,ENERGY,'
        'SUML,MAD Z,MAD X,MAD Y,MAD Psi,MAD Phi,MAD Theta,'
        'RF_Frequency,RF_Amplitude,RF_Phase,RF_Gradient,RF_Power_Fraction,'
        'Z_Length,Fringe_Field_Integral,Integrated_Field_BL,Field_B,'
        'Integrated_Field_Gradient_GL,Field_Gradient_G,'
        'XAL_Scale_Name,XAL_Scale_Value,XAL_Polarity,'
        'Magnet_MAD_Z,Magnet_MAD_X,Magnet_MAD_Y,'
        'Magnet_MAD_Psi,Magnet_MAD_Phi,Magnet_MAD_Theta,'
        'Solenoid_Strength_KS,Undulator_Period_Length,Undulator_Strength_K,'
        'X_Size,Y_Size,'
        'Section,Distance_From_Section_Start,XAL_Keyword,S_Display')

unit = (',,,,,m,'
        'm,deg,1/m^2,1/m^3,deg,deg,deg,1/m,1/m,GeV,'
        'm,m,m,m,rad,rad,rad,'
        'MHz,MeV,deg,MeV/m,1,'
        'm,1,kG-m,T,'
        'kG,T/m,'
        ',,,'
        'm,m,m,'
        'rad,rad,rad,'
        '1/m,m,1,'
        'm,m,'
        ',m,,m')

Ncol = head.count(',') + 1

# SYMBOLS text-file (linac coordinates)

def madval(rval):
    if rval == 0:
        return '0.0'
    else:
        aval = abs(rval)
        if 0.01 < aval < 100:
            iexp = 0.0
        else:
            iexp = int(math.log10(aval))
            rval = rval * 10**(-iexp)
        
        s = f'{rval:.12f}'
        s = s.rstrip('0')
        
        if s.endswith('.'):
            s = s + '0'
        
        if iexp != 0:
            s += f'E{iexp}'
        
        return s.strip()

def roundoff(val, prec=None):
    if isinstance(val,list):
        return None
    if prec is None:
        return val
    else:
        return prec * np.round(val / prec)

fname = f'AD_ACCEL-{optics}.txt'
with open(outdir+'/'+fname, 'wt') as fid:
    fid.write(f'{head}\n')
    fid.write(f'{unit}\n')

    for nf in range(1, len(froot) + 1):
        id = [i for i,x in enumerate(ip) if x[0]==nf]
        for n in id:
            idk = ip[n][1]
            idn = ip[n][2]
            TEMP = globals()[keyw[idk]][idn]
            s = [None] * Ncol

            # common data
            s[0] = TEMP['id']
            s[1] = TEMP['parent']
            if TEMP['prim'] == 'EFC':
                s[2] = 'EFC'
            elif TEMP['type'] == '0.79K11.8':
                s[2] = 'BEND'
            elif TEMP['xkey'] == 'MARK' and TEMP['prim'] == 'INST':
                s[2] = 'INST'
            else:
                s[2] = keyo[idk]
            s[3] = TEMP['name']
            s[4] = TEMP['type']
            s[15] = TEMP['energy']
            s[16] = TEMP['suml']
            s[17] = roundoff(TEMP['c1'], prec)
            s[18] = roundoff(TEMP['c2'], prec)
            s[19] = roundoff(TEMP['c3'], prec)
            s[20] = roundoff(TEMP['c4'], prec)
            s[21] = roundoff(TEMP['c5'], prec)
            s[22] = roundoff(TEMP['c6'], prec)
            s[48] = TEMP['seq']
            s[49] = TEMP['dist']
            s[50] = TEMP['xkey']
            s[51] = TEMP['sdsp']

            # keyword data
            if keyw[idk] == 'LCAV':
                s[5] = TEMP['leng']
                s[23] = TEMP['freq']
                s[24] = TEMP['ampl']
                s[25] = TEMP['phase']
                s[26] = TEMP['grad']
                s[27] = TEMP['power']
            elif keyw[idk] == 'SBEN':
                s[5] = TEMP['leng']
                s[6] = TEMP['gap']
                s[7] = TEMP['ang']
                s[8] = TEMP['k1']
                s[10] = TEMP['tilt']
                s[11] = TEMP['e1']
                s[12] = TEMP['e2']
                s[28] = TEMP['zleng']
                s[29] = TEMP['fint']
                s[30] = TEMP['BL']
                s[31] = TEMP['B']
                s[32] = TEMP['GL']
                s[33] = TEMP['G']
                s[34] = TEMP['sname']
                s[35] = TEMP['sval']
                s[36] = TEMP['polarity']
                s[37] = roundoff(TEMP['m1'], prec)
                s[38] = roundoff(TEMP['m2'], prec)
                s[39] = roundoff(TEMP['m3'], prec)
                s[40] = roundoff(TEMP['m4'], prec)
                s[41] = roundoff(TEMP['m5'], prec)
                s[42] = roundoff(TEMP['m6'], prec)
            elif keyw[idk] == 'QUAD':
                s[5] = TEMP['leng']
                s[6] = TEMP['bore']
                s[8] = TEMP['k1']
                s[10] = TEMP['tilt']
                s[32] = TEMP['GL']
                s[33] = TEMP['G']
                s[34] = TEMP['sname']
                s[35] = TEMP['sval']
                s[36] = TEMP['polarity']
                s[37] = roundoff(TEMP['m1'], prec)
                s[38] = roundoff(TEMP['m2'], prec)
                s[39] = roundoff(TEMP['m3'], prec)
                s[40] = roundoff(TEMP['m4'], prec)
                s[41] = roundoff(TEMP['m5'], prec)
                s[42] = roundoff(TEMP['m6'], prec)
            elif keyw[idk] == 'SEXT':
                s[5] = TEMP['leng']
                s[6] = TEMP['bore']
                s[9] = TEMP['k2']
                s[10] = TEMP['tilt']
                s[32] = TEMP['GpL']
                s[33] = TEMP['Gp']
                s[34] = TEMP['sname']
                s[35] = TEMP['sval']
                s[36] = TEMP['polarity']
            elif keyw[idk] == 'SOLE':
                s[5] = TEMP['leng']
                s[6] = TEMP['bore']
                s[30] = TEMP['BL']
                s[31] = TEMP['B']
                s[34] = TEMP['sname']
                s[35] = TEMP['sval']
                s[36] = TEMP['polarity']
                s[43] = TEMP['ks']
            elif keyw[idk] == 'MATR':
                s[5] = TEMP['leng']
                s[44] = TEMP['lambda']
                s[45] = TEMP['k']
            elif keyw[idk] == 'RCOL':
                s[5] = TEMP['leng']
                s[46] = TEMP['xgap']
                s[47] = TEMP['ygap']
            elif keyw[idk] == 'ECOL':
                s[5] = TEMP['leng']
                s[46] = TEMP['xbore']
                s[47] = TEMP['ybore']
            elif keyw[idk] == 'SROT':
                s[5] = TEMP['leng']
                s[7] = TEMP['ang']
            elif keyw[idk] == 'INST':
                s[5] = TEMP['leng']
                if TEMP['xkey'] != 'MARK':
                    s[37] = roundoff(TEMP['m1'], prec)
                    s[38] = roundoff(TEMP['m2'], prec)
                    s[39] = roundoff(TEMP['m3'], prec)
                    s[40] = roundoff(TEMP['m4'], prec)
                    s[41] = roundoff(TEMP['m5'], prec)
                    s[42] = roundoff(TEMP['m6'], prec)
            elif keyw[idk] in [keyw[i] for i in idmisc]:
                s[5] = TEMP['leng']

            fid.write(f"{s[0]+1},")
            for k in range(1, Ncol - 1):
                if s[k] is None:
                    fid.write(",")
                elif isinstance(s[k], np.ndarray):
                    fid.write(",")
                elif isinstance(s[k], str):
                    fid.write(f"{s[k]},")
                else:
                    fid.write(f"{madval(s[k])},")
            if isinstance(s[Ncol - 1], str):
                fid.write(f"{s[Ncol - 1]}\n")
            else:
                fid.write(f"{madval(s[Ncol - 1])}\n")
    
    fid.write(f'{foot}\n')
    fid.write(f'{unit}\n')

# ------------------------------------------------------------------------------
# Write extra SYMBOLS txt-file ...

# Element name, area name, undulator cell, sector

fname = f'AD_ACCEL-extra-{optics}.txt'
with open(outdir+'/'+fname, 'wt') as fid:
    fid.write('ELEMENT,Area2,Undulator Cell,Sector\n')
    for nf in range(1,len(froot)+1):
        id = [i for i,x in enumerate(ip) if x[0]==nf]
        for n in id:
            idk = ip[n][1]
            if keyw[idk] == 'MARK' or keyw[idk] == 'SROT':
                continue
            idn = ip[n][2]
            TEMP = globals()[keyw[idk]][idn]
            TEMPucell = TEMP['ucell']
            TEMPucell = '' if isinstance(TEMPucell,list) else TEMPucell
            fid.write(f"{TEMP['name']},{TEMP['area']},{TEMPucell},{TEMP['sector']}\n")
    fid.write('ELEMENT,Area2,Undulator Cell,Sector\n')

# ------------------------------------------------------------------------------

# save RDBdata

print(f'Be sure to add FACET2 elements to {fname}!\n')


