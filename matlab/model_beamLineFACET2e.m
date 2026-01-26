function beamLine=model_beamLineFACET2e()
%
% -----------------------------------------------------------------------------
% *** OPTICS=28JAN2026 ***
% -----------------------------------------------------------------------------
%
% Returns Matlab model beam lines that correspond to defined FACET-II
% electron beampaths:
%
%  beamLine.F2_ELEC     = gun to main electron beam dump (PAX chicane bypassed)
%  beamLine.F2_ELEC_PAX = gun to main electron beam dump (PAX chicane included)
%  beamLine.F2_SCAV     = gun to positron production target
%
% Additional beam lines used for comparison with MAD (the start point is at
% element BEGDL10, at 135 MeV):
%
%  beamLine.F2_ELECI     = BEGDL10 to main electron beam dump (PAX chicane bypassed)
%  beamLine.F2_ELECI_PAX = BEGDL10 to main electron beam dump (PAX chicane included)
%  beamLine.F2_SCAVI     = BEGDL10 to positron production target
%
% -----------------------------------------------------------------------------
 
% check for mat-file version ... load and return beamLine if found
if (exist('model_beamLineFACET2e.mat')==2)
  load model_beamLineFACET2e.mat
  return
end

[FACET2E,FACET2E_PAX,FACET2S,INJ]=bl();
beamLine.F2_ELEC=[INJ,FACET2E]';


[FACET2E,FACET2E_PAX,FACET2S,INJ]=bl();
beamLine.F2_ELEC_PAX=[INJ,FACET2E_PAX]';


[FACET2E,FACET2E_PAX,FACET2S,INJ]=bl();
beamLine.F2_SCAV=[INJ,FACET2S]';


[FACET2E,FACET2E_PAX,FACET2S,INJ]=bl();
beamLine.F2_ELECI=[FACET2E]';


[FACET2E,FACET2E_PAX,FACET2S,INJ]=bl();
beamLine.F2_ELECI_PAX=[FACET2E_PAX]';


[FACET2E,FACET2E_PAX,FACET2S,INJ]=bl();
beamLine.F2_SCAVI=[FACET2S]';



function [FACET2E,FACET2E_PAX,FACET2S,INJ]=bl()


PI     = pi;
TWOPI  = 2*pi;
DEGRAD = 180/pi;
RADDEG = pi/180;
E      = exp(1);
EMASS  = 0.510998902e-3; % electron rest mass [GeV]
PMASS  = 0.938271998;    % proton rest mass [GeV]
CLIGHT = 2.99792458e8;   % speed of light [m/s]

% *** OPTICS=FACET2-28JAN2026 ***



% ==============================================================================
% Modification History
% ------------------------------------------------------------------------------
% 06-SEP-2023, M. Woodley
%  * fix BC14E/P drifts
%  * set all QUAD's and SEXT's to Lucretia values
% ------------------------------------------------------------------------------
% 18-MAY-2021, G. White
%  * fixed LI19 corrector locations, removed 19-8a and updated S20 experiment
%    object locations @ Q0-2D quad locations
%  * added Q0FF, Q1FF, Q2FF, removed QFF4, locations according to metrology
%    measurements
% ------------------------------------------------------------------------------
% 13-NOV-2017, M. Woodley
%  * move CATHODE closer to L0; remove BPM3F
% ------------------------------------------------------------------------------
% 05-SEP-2017, M. Woodley
%  * tweak XLL, ZLL, and LINJ to get correct SURVEY coordinates
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% element and line definitions
% ------------------------------------------------------------------------------
% *** OPTICS=FACET2-28JAN2026 ***
% ------------------------------------------------------------------------------
% constants and global parameters (from FACET2e_baseline.mat)
% ------------------------------------------------------------------------------
% constants
CB=1.0E10/CLIGHT;%energy (GeV) to magnetic rigidity (kG-m)
IN2M=0.0254;%inches to meters
QELEC=1.602176462E-19;%electron charge (C)
SBANDF=2856;%S-band rf frequency (MHz)
XBANDF=11424;%X-band rf frequency (MHz)
DLWL10=3.0441;%"10 ft" structure length (m)
DLWL9=2.8692;%"9.4 ft" structure length (m)
DLWL7=2.1694;%"7 ft" structure length (m)
DLWLX=0.5948;%Xband structure length (m)
P25=1;%25% power factor
P50=sqrt(2);%50% power factor
% global parameters
Q0 =  2.0E-9 ;%C
R56_HTR =   0.007892 ;%m
R56_DL10 =  -0.006286 ;%m
R56_BC11 =   0.045898 ;%m
R56_BC14 =   0.036021 ;%m
R56_BC20 =   0;
QSIGN =  +1 ;%electron=+1; positron=-1
% energy profile (treaty values except for E19)
E0 =   0.006;
E0A =   0.064;
EI =   0.125;
E11 =   0.335;
E14 =   4.5;
E19 =   9.778986367937 ;%at MSCAVEXT
E20 =  10.0;
% BC14 parameters
R11E =   0.932028241295;
R12E =  14.0;
R21E =  -0.01;
R33E =   0.932028241295;
R34E =  14.0;
R43E =  -0.01;
R11P =  R33E;
R12P =  R34E;
R21P =  R43E;
R33P =  R11E;
R34P =  R12E;
R43P =  R21E;
% ------------------------------------------------------------------------------
% Twiss (from Lucretia/FACET2e.mat unless otherwise noted)
% ------------------------------------------------------------------------------
% at BEGDL10
BXI =   0.137761791898;
AXI =   0.620280308601;
BYI =   7.063979455311;
AYI =  -5.750562653636;
% at CATHODEF (matched to BXi/AXi/BYi/AYi)
BX0 =  0.132923615689 ;%  0.132923615646
AX0 =  0.896770130014 ;%  0.896770129714
BY0 =  0.365931285843 ;%  0.36593128579
AY0 =  1.863967282833 ;%  1.863967282674
% at MRK0F
BX10 =   5.285886040780;
AX10 =  -2.039021010174;
BY10 =   2.581889827366;
AY10 =   0.052047744707;
% at BC11CEND (treaty values)
BX11 =  3.0;
AX11 =  0.0;
BY11 =  3.0;
AY11 =  0.0;
% at BEGBC14E
BX14I =  70.22929918739;
AX14I =   2.506815398918;
BY14I =  65.681299785911;
AY14I =   2.363169950675;
% at ENDBC14E
BX14 =   8.400776344096;
AX14 =  -0.004252878348;
BY14 =   8.671895574182;
AY14 =   0.027909144552;
% at MSCAVEXT
BX19 =  13.114920013535;
AX19 =   0.678453219664;
BY19 =  41.689116108638;
AY19 =  -1.989996071505;
% at BEGBC20
BX20 =  11.502858236938;
AX20 =   0.704134099969;
BY20 =  27.275402080101;
AY20 =   1.224927250207;
% at MIP (treaty values)
BXIP =  0.5;
AXIP =  0.0;
BYIP =  0.5;
AYIP =  0.0;
% at MAINDUMP
BXD =  40.758986780294;
AXD =  -3.969312435200;
BYD =   3.127601194552;
AYD =   0.128435821760;
% ------------------------------------------------------------------------------
% misc
% ------------------------------------------------------------------------------
BMAXL2 =  48.250          ;%45 degree cells
BMAXL3 =  41.845226568382 ;%65 degree cells
BMAX19 =  70;
% ------------------------------------------------------------------------------
% load lattice definitions
% ------------------------------------------------------------------------------
% *** OPTICS=FACET2-28JAN2026 ***
% FACET2 common parameters
% ==============================================================================
% Modification History
% ------------------------------------------------------------------------------
% 31-JAN-2017, M. Woodley
%  * from FACET2e.xsif
% ------------------------------------------------------------------------------
% parameters
ADL1 =  -35.0*RADDEG ;%injection line angle [rad]
% FACET2 quadrupoles
% Qx  : Everson-Tesla quadrupole (1.26Q3.5)
% QE  : linac QE4 (1.085Q4.31)
% Qc  : SigmaPhi "tweaker" quadrupole (1.69Q3.4)
% Qc2 : BC14 "tweaker" quadrupole (2.1Q5.87)
LQX =  0.108  ;
RQX =  0.016        ;%GLmax=  20 kG @  12 A
LQE =  0.1068 ;
RQE =  1.085*IN2M/2 ;%GLmax= 106 kG @ 220 A
LQC =  0.108  ;
RQC =  0.043/2      ;%GLmax= 2.1 kG @ 12 A
LQC2 =  0.1759 ;
RQC2 =  2.1*IN2M/2   ;%GLmax=   4 kG @ 8.67 A
LSQ =  0.3813  ;
ASQ =  2.026*IN2M/2;
LC260 =  0.260 ;%3D8.8 corrector effective length
% LI19 scavenger kicker
ZKS =  0.5;
% common drift lengths
LDAQ1 =  0.0342;
LDAQ2 =  0.027;
LDAQ3 =  0.3533;
LDAQ4 =  2.5527;
DAQ1={'dr' '' LDAQ1 []}';
DAQ2={'dr' '' LDAQ2 []}';
DAQ3={'dr' '' LDAQ3 []}';
DAQ4={'dr' '' LDAQ4 []}';
% ------------------------------------------------------------------------------

% *** OPTICS=FACET2-28JAN2026 ***
% ==============================================================================
% Modification History
% ------------------------------------------------------------------------------
% 17-DEC-2019, M. Woodley
%  * define unsplit L0AF and L0BF for Bmad
% ------------------------------------------------------------------------------
% 23-AUG-2018, M. Woodley
%  * gradL0A, gradL0B, and L0phase from FACET2e_baseline.mat
%  * remove old LCLS-derived common names
% ------------------------------------------------------------------------------
% 08-NOV-2017, M. Woodley
%  * move CATHODE 442.5443 mm closer to L0 than LCLS-II Phase I per G. Bouchard
%  * remove BPM3F (BPM10235)
% ------------------------------------------------------------------------------
% 07-APR-2017, M. Woodley
%  * from LCLS2 Phase 1 (03MAR2013)
%  * remove nonexistent components
%  * rename components (use LISTs until we're sure ... )
%  * switch from appending "B" to LCLS element names to appending "F"
% ------------------------------------------------------------------------------
% ==============================================================================
% L0AF/L0BF (analytic gradients)
% ------------------------------------------------------------------------------
DEL0A =  E0A-E0;
DEL0B =  EI-E0A;
LL0ACC =  3.095244 ;%length of L0A and L0B accelerating structures (m)
% L0phase := L0A/L0B S-band rf phase (deg)
% PhiL0   := L0A/L0B S-band rf phase (radians/2pi)
% gradL0A := L0A accelerating gradient (MeV/m)
% gradL0B := L0B accelerating gradient (MeV/m)
L0PHASE =  -2.5;
PHIL0 =  L0PHASE/360;
GRADL0A =  1.0E3*DEL0A/(LL0ACC*cos(PHIL0*TWOPI)) ;%MeV/m
GRADL0B =  1.0E3*DEL0B/(LL0ACC*cos(PHIL0*TWOPI)) ;%MeV/m
L0AF__1={'lc' 'L0AF' 0.0586460 [SBANDF GRADL0A*0.0586460 PHIL0*TWOPI]}';
L0AF__2={'lc' 'L0AF' 0.1993540 [SBANDF GRADL0A*0.1993540 PHIL0*TWOPI]}';
L0AF__3={'lc' 'L0AF' 0.6493198 [SBANDF GRADL0A*0.6493198 PHIL0*TWOPI]}';
L0AF__4={'lc' 'L0AF' 0.6403022 [SBANDF GRADL0A*0.6403022 PHIL0*TWOPI]}';
L0AF__5={'lc' 'L0AF' 1.1518464 [SBANDF GRADL0A*1.1518464 PHIL0*TWOPI]}';
L0AF__6={'lc' 'L0AF' 0.3348566 [SBANDF GRADL0A*0.3348566 PHIL0*TWOPI]}';
L0AF__7={'lc' 'L0AF' 0.0609190 [SBANDF GRADL0A*0.0609190 PHIL0*TWOPI]}';
L0BF__1={'lc' 'L0BF' 0.0586460 [SBANDF GRADL0B*0.0586460 PHIL0*TWOPI]}';
L0BF__2={'lc' 'L0BF' 0.3371281 [SBANDF GRADL0B*0.3371281 PHIL0*TWOPI]}';
L0BF__3={'lc' 'L0BF' 1.1518479 [SBANDF GRADL0B*1.1518479 PHIL0*TWOPI]}';
L0BF__4={'lc' 'L0BF' 1.1515630 [SBANDF GRADL0B*1.1515630 PHIL0*TWOPI]}';
L0BF__5={'lc' 'L0BF' 0.3351400 [SBANDF GRADL0B*0.3351400 PHIL0*TWOPI]}';
L0BF__6={'lc' 'L0BF' 0.0609190 [SBANDF GRADL0B*0.0609190 PHIL0*TWOPI]}';
% define unsplit LCAVs for BMAD ... not used by MAD
L0AF={'lc' 'L0AF' LL0ACC [SBANDF GRADL0A*LL0ACC PHIL0*TWOPI]}';
L0BF={'lc' 'L0BF' LL0ACC [SBANDF GRADL0B*LL0ACC PHIL0*TWOPI]}';
% ==============================================================================
% QUAD
% ------------------------------------------------------------------------------
% CQ10121 = correction quad in 1st solenoid at gun (nominally set to 0)
% SQ10122 = correction skew-quad in 1st solenoid at gun (nominally set to 0)
% ------------------------------------------------------------------------------
CQ10121={'mu' 'CQ10121' 0 [0 0 0 0]}';%solenoid trim
SQ10122={'mu' 'SQ10122' 0 [0 0 0 pi/4]}';%solenoid trim
KQA10361 =  -11.288886255557;
KQA10371 =   11.405472211446;
QA10361={'qu' 'QA10361' LQX/2 [KQA10361 0]}';
QA10371={'qu' 'QA10371' LQX/2 [KQA10371 0]}';
% ==============================================================================
% SOLE
% ------------------------------------------------------------------------------
% - SOL10111 = gun bucking-solenoid (set to ~zero length and strength, with
%              longitudinal unknown for now)
% - SOL10121 = gun solenoid
% ------------------------------------------------------------------------------
LSOL1 =  0.2;
SOL10111={'so' 'SOL10111' 0 [0]}';
SOL10121={'so' 'SOL10121' LSOL1/2 [0]}';%design= 0.38 kG-m
% ==============================================================================
% DRIF
% ------------------------------------------------------------------------------
LGGUN =  7.51*0.3048;
LOADLOCKF={'dr' '' LGGUN-1.42+1.E-9 []}';
DL00={'dr' '' -(LOADLOCKF{3}+SOL10111{3}) []}';%from cathode back to u/s end of loadlock
DL01={'dr' '' 1.0322037 []}';%1.474748003 (dS= -442.544303 mm)
DL02={'dr' '' 0.2309416 []}';
DL03={'dr' '' 0.220376 []}';
DL04={'dr' '' 0.065888 []}';
DL01A={'dr' '' 0.097294699+1.E-9 []}';
DL01B={'dr' '' 0.078510099+1.E-9 []}';
DL01C={'dr' '' 0.1160862 []}';
DL01D={'dr' '' 0.0811358 []}';
DL01E={'dr' '' 0.1263224 []}';
DL01F={'dr' '' 0.8509E-2 []}';
DL01G={'dr' '' 0.2342261 []}';
DL01H={'dr' '' 0.0901194 []}';
DL03A={'dr' '' 0.094484 []}';
DL03B={'dr' '' 0.125892 []}';
%VALUE, DL01a[L]+LSOL1+2*CQ10121[L]+2*SQ10122[L]+DL01b[L]+DL01c[L]+DL01d[L]+ &
%       DL01e[L]+DL01f[L]+DL01g[L]+DL01h[L]
% Injector AIP
DL01D2={'dr' '' 0.219 []}';
DL01E2={'dr' '' 0.15 []}';
DL01F2={'dr' '' 0.11 []}';
% ==============================================================================
% XCORs and YCORs
% ------------------------------------------------------------------------------
XC10121={'mo' 'XC10121' 0 []}';%solenoid dipole trim
XC10221={'mo' 'XC10221' 0 []}';
XC10311={'mo' 'XC10311' 0 []}';
XC10381={'mo' 'XC10381' 0 []}';%fast-feedback (loop-1)
XC10411={'mo' 'XC10411' 0 []}';%calibrated to <1%
YC10122={'mo' 'YC10122' 0 []}';%solenoid dipole trim
YC10222={'mo' 'YC10222' 0 []}';
YC10312={'mo' 'YC10312' 0 []}';
YC10382={'mo' 'YC10382' 0 []}';%fast-feedback (loop-1)
YC10412={'mo' 'YC10412' 0 []}';%calibrated to <1%
% ==============================================================================
% diagnostics
% ------------------------------------------------------------------------------
% BPMs
BPM10221={'mo' 'BPM10221' 0 []}';
BPM10371={'mo' 'BPM10371' 0 []}';
% misc
CATHODEF={'mo' 'CATHODEF' 0 []}';
VV10155={'mo' 'VV10155' 0 []}';%vacuum valve near gun
MIR10181={'mo' 'MIR10181' 0 []}';%gun laser normal incidence mirror
VV10215={'mo' 'VV10215' 0 []}';%vacuum valve near gun
PR10241={'mo' 'PR10241' 0 []}';%gun
PH10365={'mo' 'PH10365' 0 []}';%phase measurement cavity between L0a and L0b
% Injector AIP
FCUP={'mo' 'FCUP' 0 []}';%Faraday cup used as injector AIP test beam dump
% ==============================================================================
% MARK
% ------------------------------------------------------------------------------
BEGINJ={'mo' 'BEGINJ' 0 []}';
L0AFBEG={'mo' 'L0AFBEG' 0 []}';
FLNGAF1={'mo' 'FLNGAF1' 0 []}';%upstream face of L0a entrance flange
DLFDAF={'mo' 'DLFDAF' 0 []}';%dual-feed input coupler location for L0a structure
L0AFMID={'mo' 'L0AFMID' 0 []}';
OUTCPAF={'mo' 'OUTCPAF' 0 []}';%output coupler location for L0a structure
FLNGAF2={'mo' 'FLNGAF2' 0 []}';%downstream face of L0a exit flange
L0AFWAKE={'mo' 'L0AFWAKE' 0 []}';
L0AFEND={'mo' 'L0AFEND' 0 []}';
L0BFBEG={'mo' 'L0BFBEG' 0 []}';
FLNGBF1={'mo' 'FLNGBF1' 0 []}';%upstream face of L0b entrance flange
DLFDBF={'mo' 'DLFDBF' 0 []}';%dual-feed input coupler location for L0b structure
L0BFMID={'mo' 'L0BFMID' 0 []}';
OUTCPBF={'mo' 'OUTCPBF' 0 []}';%output coupler location for L0b structure
FLNGBF2={'mo' 'FLNGBF2' 0 []}';%downstream face of L0b exit flange
L0BFWAKE={'mo' 'L0BFWAKE' 0 []}';
L0BFEND={'mo' 'L0BFEND' 0 []}';
ENDINJ={'mo' 'ENDINJ' 0 []}';
% ==============================================================================
% BEAMLINEs
% ------------------------------------------------------------------------------
SC1F=[XC10221,YC10222];
SC2F=[XC10311,YC10312];
SC4F=[XC10381,YC10382];
SC5F=[XC10411,YC10412];
L0AF_FULL=[L0AFBEG,FLNGAF1,L0AF__1,DLFDAF,L0AF__2,L0AF__3,SC2F,L0AF__4,L0AFMID,L0AF__5,L0AF__6,OUTCPAF,L0AF__7,FLNGAF2,L0AFWAKE,L0AFEND];
L0BF_FULL=[L0BFBEG,FLNGBF1,L0BF__1,DLFDBF,L0BF__2,SC4F,L0BF__3,L0BFMID,L0BF__4,SC5F,L0BF__5,OUTCPBF,L0BF__6,FLNGBF2,L0BFWAKE,L0BFEND];
% FACET-II configuration with L0A 0.97m from the cathode
QA10361_FULL=[QA10361,QA10361];
QA10371_FULL=[QA10371,BPM10371,QA10371];
SOL10121_FULL=[SOL10121,CQ10121,XC10121,YC10122,SQ10122,SOL10121];
INJ=[DL00,LOADLOCKF,BEGINJ,SOL10111,CATHODEF,DL01A,SOL10121_FULL,DL01B,VV10155,DL01C,MIR10181,DL01D,VV10215,DL01E,BPM10221,DL01F,SC1F,DL01G,PR10241,DL01H,L0AF_FULL,DL02,QA10361_FULL,DL03A,PH10365,DL03B,QA10371_FULL,DL04,L0BF_FULL,ENDINJ];
INJS10AIP=[DL00,LOADLOCKF,BEGINJ,SOL10111,CATHODEF,DL01A,SOL10121_FULL,DL01B,VV10155,DL01C,MIR10181,DL01D2,SC1F,DL01E2,PR10241,DL01F2,FCUP,ENDINJ];

% ==============================================================================
LL0A =  L0AF__1{3}+L0AF__2{3}+L0AF__3{3}+XC10311{3}+YC10312{3}+L0AF__4{3}+ L0AF__5{3}+L0AF__6{3}+L0AF__7{3};
LL0B =  L0BF__1{3}+L0BF__2{3}+XC10381{3}+YC10382{3}+L0BF__3{3}+L0BF__4{3}+ XC10411{3}+YC10412{3}+L0BF__5{3}+L0BF__6{3};

% ------------------------------------------------------------------------------

% *** OPTICS=FACET2-28JAN2026 ***
% ==============================================================================
% Modification History
% ------------------------------------------------------------------------------
% 05-SEP-2023, M. Woodley
%  * undefer laser heater undulator
% ------------------------------------------------------------------------------
% 20-OCT-2022, M. Woodley
%  * undefer laser heater chicane dipoles ... rematch (undulator still deferred)
% ------------------------------------------------------------------------------
% 19-FEB-2021, M. Woodley
%  * remove IN10 ACMs (IM10607 and IM10618)
% ------------------------------------------------------------------------------
% 25-JUL-2019, M. Woodley
%  * laser heater chicane, undulator, and OTR's won't be installed for startup
%    per N. Lipkowitz ... defer (level 0)
% ------------------------------------------------------------------------------
% 05-MAR-2019, M. Woodley
%  * add WALLBEG, WALLEND (INSTs), IM10607, IM10618 (IMONs) per N. Lipkowitz
% ------------------------------------------------------------------------------
% 23-AUG-2018, M. Woodley
%  * correct TYPE designation of dogleg bends
%  * quadrupole K1 values from FACET2e_baseline.mat
% ------------------------------------------------------------------------------
% 16-MAY-2017, M. Woodley
%  * add bunch length monitor (ceramic gap) BZ10596
% 06-APR-2017, M. Woodley
%  * use LCLS unit numbers in element names
%  * undefer TCY10490
%  * remove OTR1B
%  * change keyword for vacuum valves (VV*) to INST
% 31-JAN-2017, M. Woodley
%  * from FACET2e.xsif
% ------------------------------------------------------------------------------
% ==============================================================================
% LCAV
% ------------------------------------------------------------------------------
% vertically deflecting transverse cavity
TCY10490={'tc' 'TCY10490' 0.6680236/2 [SBANDF 0 0*TWOPI]}';
% ==============================================================================
% SBEN
% ------------------------------------------------------------------------------
% laser heater chicane (1.18D3.17 dipoles)
% - use series approximation for sinc(x)=sin(x)/x to allow BBh=0
% GBH   : gap height (m)
% ZBH   : "Z" length (m)
% FBH   : measured fringe field integral (1)
% ABH   : chicane bend angle (rad)
% LBH   : chicane bend path length (m)
% ABHs  : "short" half chicane bend angle (rad)
% LBHs  : "short" half chicane bend path length (m)
% ABHl  : "long" half chicane bend angle (rad)
% LBHl  : "long" half chicane bend path length (m)
GBH =  0.03;
ZBH =  0.1244;
FBH =  0.3997;
ABH =  0.1316410831;
ABH_2 =  ABH*ABH;
ABH_4 =  ABH_2*ABH_2;
ABH_6 =  ABH_4*ABH_2;
SINCABH =  1-ABH_2/6+ABH_4/120-ABH_6/5040 ;%~sinc(ABH)=sin(ABH)/ABH
LBH =  ZBH/SINCABH;
ABHS =  asin(sin(ABH)/2);
ABHS_2 =  ABHS*ABHS;
ABHS_4 =  ABHS_2*ABHS_2;
ABHS_6 =  ABHS_4*ABHS_2;
SINCABHS =  1-ABHS_2/6+ABHS_4/120-ABHS_6/5040 ;%~sinc(ABHs)=sin(ABHs)/ABHs
LBHS =  (ZBH/2)/SINCABHS;
ABHL =  ABH-ABHS;
LBHL =  LBH-LBHS;
BCX10451A={'be' 'BCX10451' LBHS [-ABHS GBH/2 0 0 FBH 0 0]}';
BCX10451B={'be' 'BCX10451' LBHL [-ABHL GBH/2 0 -ABH 0 FBH 0]}';
BCX10461A={'be' 'BCX10461' LBHL [+ABHL GBH/2 +ABH 0 FBH 0 0]}';
BCX10461B={'be' 'BCX10461' LBHS [+ABHS GBH/2 0 0 0 FBH 0]}';
BCX10475A={'be' 'BCX10475' LBHS [+ABHS GBH/2 0 0 FBH 0 0]}';
BCX10475B={'be' 'BCX10475' LBHL [+ABHL GBH/2 0 +ABH 0 FBH 0]}';
BCX10481A={'be' 'BCX10481' LBHL [-ABHL GBH/2 -ABH 0 FBH 0 0]}';
BCX10481B={'be' 'BCX10481' LBHS [-ABHS GBH/2 0 0 0 FBH 0]}';
% define unsplit SBENs for BMAD ... not used by MAD
BCX10451={'be' 'BCX1045' LBH [-ABH GBH/2 0 -ABH FBH FBH 0]}';
BCX10461={'be' 'BCX1046' LBH [+ABH GBH/2 +ABH 0 FBH FBH 0]}';
BCX10475={'be' 'BCX10475' LBH [+ABH GBH/2 0 +ABH FBH FBH 0]}';
BCX10481={'be' 'BCX1048' LBH [-ABH GBH/2 -ABH 0 FBH FBH 0]}';
% dogleg (1.182D6.82T dipoles)
GB0 =  0.03                   ;%gap height (m)
ZB0 =  0.2032                 ;%full "Z" length (m)
FB0 =  0.45                   ;%measured fringe field integral (1)
AB0 =  ADL1/2                 ;%full bend angle (rad)
LB0 =  ZB0*AB0/(2*sin(AB0/2)) ;%full bend path length (m)
BX10661A={'be' 'BX10661' LB0/2 [AB0/2 GB0/2 AB0/2 0 FB0 0 0]}';
BX10661B={'be' 'BX10661' LB0/2 [AB0/2 GB0/2 0 AB0/2 0 FB0 0]}';
BX10751A={'be' 'BX10751' LB0/2 [AB0/2 GB0/2 AB0/2 0 FB0 0 0]}';
BX10751B={'be' 'BX10751' LB0/2 [AB0/2 GB0/2 0 AB0/2 0 FB0 0]}';
% define unsplit SBENs for BMAD ... not used by MAD
BX10661={'be' 'BX1066' LB0 [AB0 GB0/2 AB0/2 AB0/2 FB0 FB0 0]}';
BX10751={'be' 'BX1075' LB0 [AB0 GB0/2 AB0/2 AB0/2 FB0 FB0 0]}';
% ==============================================================================
% MATR
% ------------------------------------------------------------------------------
% laser heater undulator
% - half-undulator modeled as R-matrix to include vertical natural focusing
% lam   : laser-heater undulator period [m]
% lamr  : heater laser wavelength [m]
% gami  : Lorentz energy factor in laser-heater undulator [1]
% K_und : undulator K for laser heater undulator
% Lhun  : half-length of laser-heater undulator (5 periods) [m]
% kqlh  : natural undulator focusing "k" in y-plane [1/m2]
LAM =  0.053855 ;% 0.054 -> changed to reflect K=1.17 measurememnt data
LAMR =  758E-9;
GAMI =  EI/EMASS;
K_UND =  sqrt(2*(LAMR*2*GAMI^2/LAM-1));
LHUN =  0.506263/2;
KQLH =  (K_UND*2*PI/LAM/sqrt(2)/GAMI)^2;
% handle K_und->0 by expressing R34 as an approximate SINC function
ARGU =  LHUN*sqrt(KQLH);
ARGU2 =  ARGU*ARGU;
ARGU4 =  ARGU2*ARGU2;
ARGU6 =  ARGU4*ARGU2;
SINCARGU =  1-ARGU2/6+ARGU4/120-ARGU6/5040 ;%~sinc(ARGu)=sin(ARGu)/ARGu
R34U =  LHUN*SINCARGU;
%comment
UM10466={'un' 'UM10466' LHUN [KQLH LAM 1]}';%,                            &
%RM(5,6)=Lhun/(gami^2)*(1+(K_und^2)/2)
%endcomment
%UM10466 : DRIF, TYPE="LHund", L=Lhun
% ==============================================================================
% QUAD
% ------------------------------------------------------------------------------
% match around laser heater (installed)
%                heater ON      
%             ----------------  
KQE10425 =  -19.483733913907;
KQE10441 =   22.738480699594;
KQE10511 =   12.1393;
KQE10525 =  -10.9134;
% match to BC11CEND (X-band not installed)
KQM10631 =   11.3442;
KQM10651 =  -11.5137;
KQB10731 =   22.169701529671;
KQM10771 =  -11.071;
KQM10781 =   12.275640217878;
QE10425={'qu' 'QE10425' LQX/2 [KQE10425 0]}';
QE10441={'qu' 'QE10441' LQX/2 [KQE10441 0]}';
QE10511={'qu' 'QE10511' LQX/2 [KQE10511 0]}';
QE10525={'qu' 'QE10525' LQX/2 [KQE10525 0]}';
QM10631={'qu' 'QM10631' LQX/2 [KQM10631 0]}';
QM10651={'qu' 'QM10651' LQX/2 [KQM10651 0]}';
QB10731={'qu' 'QB10731' LQE/2 [KQB10731 0]}';
QM10771={'qu' 'QM10771' LQX/2 [KQM10771 0]}';
QM10781={'qu' 'QM10781' LQX/2 [KQM10781 0]}';
% ==============================================================================
% DRIF
% ------------------------------------------------------------------------------
DE00={'dr' '' 0.041607 []}';
DE01={'dr' '' 0.293513 []}';
DH00={'dr' '' 0.12733045 []}';
DH01={'dr' '' 0.1406/cos(ABH) []}';
DH02={'dr' '' 0.165597479262 []}';
DH03={'dr' '' 0.159129470738 []}';
DH04={'dr' '' 0.1406/cos(ABH) []}';
DH05={'dr' '' 0.8430976 []}';
DE03={'dr' '' 0.2819662 []}';
DE04={'dr' '' 0.411914814698 []}';
IN10WALL={'dr' '' 1.4224 []}';
DE05={'dr' '' 2.439604885302 []}';
DE06={'dr' '' 0.4478221 []}';
DE07={'dr' '' 0.22345 []}';
DB00={'dr' '' 0.7264 []}';
DB01={'dr' '' 0.7264 []}';
DM00={'dr' '' 0.361983 []}';
DM01={'dr' '' 0.297167 []}';
DM02={'dr' '' 0.297648 []}';
DE01A={'dr' '' 0.076373 []}';
DE01B={'dr' '' 0.122359 []}';
DE01C={'dr' '' 0.094781 []}';
DH02A={'dr' '' 0.0809678529 []}';
DH02B={'dr' '' 0.0846296264 []}';
DH03A={'dr' '' 0.0845477112 []}';
DH03B={'dr' '' 0.0745817596 []}';
DH05A={'dr' '' 0.1229007 []}';
DH05B={'dr' '' 0.0521733 []}';
DE03A={'dr' '' 0.1454672 []}';
DE03B={'dr' '' 0.136499 []}';
DE04A={'dr' '' 0.290414154052 []}';
DE04B={'dr' '' 0.121500660646 []}';
DE05A={'dr' '' 0.280100320492 []}';
DE05B={'dr' '' 0.317702877106 []}';
DE05C={'dr' '' 0.215468253984 []}';
DE05D={'dr' '' 0.16498649439 []}';
DE05E={'dr' '' 0.391744624197 []}';
DE05F={'dr' '' 0.15 []}';
DE05G={'dr' '' 1.069602315133-DE05F{3} []}';
DE06A={'dr' '' 0.1478721 []}';
DE06B={'dr' '' 0.29995 []}';
DB00A={'dr' '' 0.3997 []}';
DB00B={'dr' '' 0.161 []}';
DB00C={'dr' '' 0.1657 []}';
DM00A={'dr' '' 0.2213 []}';
DM00B={'dr' '' 0.140683 []}';
DM02A={'dr' '' 0.1402 []}';
DM02B={'dr' '' 0.057448 []}';
DM02C={'dr' '' 0.1 []}';
DE05C1={'dr' '' DE05C{3}/2 []}';
DE05C2={'dr' '' DE05C{3}/2 []}';
DE05G1={'dr' '' 0.3 []}';
DE05G2={'dr' '' 0.3 []}';
DE05G3={'dr' '' DE05G{3}-DE05G1{3}-DE05G2{3} []}';
% ==============================================================================
% XCORs and YCORs
% ------------------------------------------------------------------------------
XC10491={'mo' 'XC10491' 0 []}';
XC10521={'mo' 'XC10521' 0 []}';
XC10641={'mo' 'XC10641' 0 []}';
XC10721={'mo' 'XC10721' 0 []}';
XC10761={'mo' 'XC10761' 0 []}';
YC10492={'mo' 'YC10492' 0 []}';
YC10522={'mo' 'YC10522' 0 []}';
YC10642={'mo' 'YC10642' 0 []}';
YC10722={'mo' 'YC10722' 0 []}';
YC10762={'mo' 'YC10762' 0 []}';
% ==============================================================================
% diagnostics
% ------------------------------------------------------------------------------
% BPMs (LCLS-II type designations)
BPM10425={'mo' 'BPM10425' 0 []}';
BPM10511={'mo' 'BPM10511' 0 []}';
BPM10525={'mo' 'BPM10525' 0 []}';
BPM10581={'mo' 'BPM10581' 0 []}';
BPM10631={'mo' 'BPM10631' 0 []}';
BPM10651={'mo' 'BPM10651' 0 []}';
BPM10731={'mo' 'BPM10731' 0 []}';
BPM10771={'mo' 'BPM10771' 0 []}';
BPM10781={'mo' 'BPM10781' 0 []}';
% misc
IM10431={'mo' 'IM10431' 0 []}';
VV10435={'mo' 'VV10435' 0 []}';
PR10465={'mo' 'PR10465' 0 []}';
PR10471={'mo' 'PR10471' 0 []}';
VV10545={'mo' 'VV10545' 0 []}';
WALLBEG={'mo' 'WALLBEG' 0 []}';
WALLEND={'mo' 'WALLEND' 0 []}';
RST10551={'mo' 'RST10551' 0 []}';
WS10561={'mo' 'WS10561' 0 []}';
PR10571={'mo' 'PR10571' 0 []}';
IM10591={'mo' 'IM10591' 0 []}';
BZ10596={'mo' 'BZ10596' 0 []}';
PR10711={'mo' 'PR10711' 0 []}';
IM10791={'mo' 'IM10791' 0 []}';
VV10795={'mo' 'VV10795' 0 []}';
% ==============================================================================
% MARK
% ------------------------------------------------------------------------------
BEGDL10={'mo' 'BEGDL10' 0 []}';
LH10BEG={'mo' 'LH10BEG' 0 []}';
HTRUNDF={'mo' 'HTRUNDF' 0 []}';
LH10END={'mo' 'LH10END' 0 []}';
MRK0F={'mo' 'MRK0F' 0 []}';%beam waist location
BX0FBEG={'mo' 'BX0FBEG' 0 []}';
BX0FEND={'mo' 'BX0FEND' 0 []}';
CNT0F={'mo' 'CNT0F' 0 []}';
ENDDL10={'mo' 'ENDDL10' 0 []}';
% ==============================================================================
% BEAMLINEs
% ------------------------------------------------------------------------------
BCX10451_FULL=[BCX10451A,BCX10451B];
BCX10461_FULL=[BCX10461A,BCX10461B];
UM10466_FULL=[UM10466,HTRUNDF,UM10466];
BCX10475_FULL=[BCX10475A,BCX10475B];
BCX10481_FULL=[BCX10481A,BCX10481B];
LH10=[LH10BEG,BCX10451_FULL,DH01,BCX10461_FULL,DH02A,PR10465,DH02B,UM10466_FULL,DH03A,PR10471,DH03B,BCX10475_FULL,DH04,BCX10481_FULL,LH10END];
BX10661_FULL=[BX10661A,BX10661B];
QB10731_FULL=[QB10731,BPM10731,QB10731];
BX10751_FULL=[BX10751A,BX10751B];
BX0F=[BX0FBEG,BX10661_FULL,DB00A,PR10711,DB00B,XC10721,YC10722,DB00C,QB10731_FULL,DB01,BX10751_FULL,CNT0F,BX0FEND];
TCY10490_FULL=[TCY10490,XC10491,YC10492,TCY10490];
QE10425_FULL=[QE10425,BPM10425,QE10425];
QE10441_FULL=[QE10441,QE10441];
QE10511_FULL=[QE10511,BPM10511,QE10511];
QE10525_FULL=[QE10525,BPM10525,QE10525];
QM10631_FULL=[QM10631,BPM10631,QM10631];
QM10651_FULL=[QM10651,BPM10651,QM10651];
QM10771_FULL=[QM10771,BPM10771,QM10771];
QM10781_FULL=[QM10781,BPM10781,QM10781];
DL10=[BEGDL10,DE00,QE10425_FULL,DE01A,IM10431,DE01B,VV10435,DE01C,QE10441_FULL,DH00,LH10,DH05A,TCY10490_FULL,DH05B,QE10511_FULL,DE03A,XC10521,YC10522,DE03B,QE10525_FULL,DE04A,VV10545,DE04B,WALLBEG,IN10WALL,WALLEND,DE05A,RST10551,DE05B,WS10561,DE05C1,MRK0F,DE05C2,PR10571,DE05D,BPM10581,DE05E,IM10591,DE05F,BZ10596,DE05G,QM10631_FULL,DE06A,XC10641,YC10642,DE06B,QM10651_FULL,DE07,BX0F,DM00A,XC10761,YC10762,DM00B,QM10771_FULL,DM01,QM10781_FULL,DM02A,IM10791,DM02B,VV10795,DM02C,ENDDL10];
% ==============================================================================

% *** OPTICS=FACET2-28JAN2026 ***
% ==============================================================================
% Modification History
% ------------------------------------------------------------------------------
% 20-NOV-2022, M. Woodley
%  * increase fL1 (gradL1) value for 125 MeV out of L0 per G. White
% ------------------------------------------------------------------------------
% 23-AUG-2018, M. Woodley
%  * gradL1 and phiL1 values from FACET2e_baseline.mat
%  * quadrupole K1 values from FACET2e_baseline.mat
% ------------------------------------------------------------------------------
% 31-JAN-2017, M. Woodley
%  * from FACET2e.xsif
% ------------------------------------------------------------------------------
% ==============================================================================
% accelerating structures
% ------------------------------------------------------------------------------
% the L1 S-band linac consists of: 3 x  10' structure @ 25% power
%                                  1 x  10' structure @ 50% power
%                                  1 x 9.4' structure @ 25% power
%                                  1 x 9.4' structure @ 50% power
% ------------------------------------------------------------------------------
FL1 =    1.035858624668     ;%0.986532023493
GRADL1 =   10.628316745522*FL1 ;%MeV/m
PHIL1 =  -20.5/360            ;%rad/2pi
KLOSSL1 =    0                  ;%V/C/m
K11_1B1={'lc' 'K11_1B' 0.453432 [SBANDF P50*GRADL1*0.453432 PHIL1*TWOPI]}';
K11_1B2={'lc' 'K11_1B' 2.415768 [SBANDF P50*GRADL1*2.415768 PHIL1*TWOPI]}';
K11_1C1={'lc' 'K11_1C' 0.6689 [SBANDF P25*GRADL1*0.6689 PHIL1*TWOPI]}';
K11_1C2={'lc' 'K11_1C' 2.3752 [SBANDF P25*GRADL1*2.3752 PHIL1*TWOPI]}';
K11_1D={'lc' 'K11_1D' DLWL10 [SBANDF P25*GRADL1*DLWL10 PHIL1*TWOPI]}';
K11_2A1={'lc' 'K11_2A' 0.3256 [SBANDF P25*GRADL1*0.3256 PHIL1*TWOPI]}';
K11_2A2={'lc' 'K11_2A' 0.2870 [SBANDF P25*GRADL1*0.2870 PHIL1*TWOPI]}';
K11_2A3={'lc' 'K11_2A' 2.4315 [SBANDF P25*GRADL1*2.4315 PHIL1*TWOPI]}';
K11_2B={'lc' 'K11_2B' DLWL9 [SBANDF P25*GRADL1*DLWL9 PHIL1*TWOPI]}';
K11_2C1={'lc' 'K11_2C' 0.5000 [SBANDF P50*GRADL1*0.5000 PHIL1*TWOPI]}';
K11_2C2={'lc' 'K11_2C' 2.5441 [SBANDF P50*GRADL1*2.5441 PHIL1*TWOPI]}';
% define unsplit LCAVs for BMAD ... not used by MAD
K11_1B={'lc' 'K11_1B' DLWL9 [SBANDF P50*GRADL1*DLWL9 PHIL1*TWOPI]}';
K11_1C={'lc' 'K11_1C' DLWL10 [SBANDF P25*GRADL1*DLWL10 PHIL1*TWOPI]}';
K11_2A={'lc' 'K11_2A' DLWL10 [SBANDF P25*GRADL1*DLWL10 PHIL1*TWOPI]}';
K11_2C={'lc' 'K11_2C' DLWL10 [SBANDF P50*GRADL1*DLWL10 PHIL1*TWOPI]}';
% ------------------------------------------------------------------------------
% L1 X-band
% ------------------------------------------------------------------------------
L1XF__1={'lc' 'L1XF' DLWLX/2 [XBANDF 0 0*TWOPI]}';
L1XF__2={'lc' 'L1XF' DLWLX/2 [XBANDF 0 0*TWOPI]}';
% define unsplit LCAVs for BMAD ... not used by MAD
L1XF={'lc' 'L1XF' DLWLX [XBANDF 0 0*TWOPI]}';
% ==============================================================================
% QUADs
% ------------------------------------------------------------------------------
KQL1 =  3.789198342593;
QDL1={'qu' 'QDL1' LQE/2 [-KQL1 0]}';
QFL1={'qu' 'QFL1' LQE/2 [+KQL1 0]}';
KQA11132 =  -3.031820754245;
KQ11201 =   1.772655565069;
KQA11265 =   0.0;
KQ11301 =  -3.164932554748;
QA11132={'qu' 'QA11132' LQE/2 [KQA11132 0]}';
Q11201={'qu' 'Q11201' LQE/2 [KQ11201 0]}';
QA11265={'qu' 'QA11265' LQE/2 [KQA11265 0]}';
Q11301={'qu' 'Q11301' LQE/2 [KQ11301 0]}';
% ==============================================================================
% drifts
% ------------------------------------------------------------------------------
D9={'dr' '' DLWL9 []}';
DAQA={'dr' '' 0.03345 []}';
DAQA1={'dr' '' 0.03747 []}';
DAQA2={'dr' '' 0.03063 []}';
DAQA3={'dr' '' 0.03405 []}';
DAQA4={'dr' '' 0.03405 []}';
DAQA5={'dr' '' 0.055970350578 []}';%0.0559696203
DAQA6={'dr' '' 0.15733138 []}';
DL1X={'dr' '' 0.14146862 []}';
DM10A={'dr' '' 0.26360938 []}';
% ==============================================================================
% XCORs and YCORs
% ------------------------------------------------------------------------------
XC11104={'mo' 'XC11104' 0 []}';
XC11140={'mo' 'XC11140' 0 []}';
XC11202={'mo' 'XC11202' 0 []}';
XC11272={'mo' 'XC11272' 0 []}';
XC11304={'mo' 'XC11304' 0 []}';
YC11105={'mo' 'YC11105' 0 []}';
YC11141={'mo' 'YC11141' 0 []}';
YC11203={'mo' 'YC11203' 0 []}';
YC11273={'mo' 'YC11273' 0 []}';
YC11305={'mo' 'YC11305' 0 []}';
% ==============================================================================
% diagnostics
% ------------------------------------------------------------------------------
% BPMs (LCLS-II type designations)
BPM11132={'mo' 'BPM11132' 0 []}';
BPM11201={'mo' 'BPM11201' 0 []}';
BPM11265={'mo' 'BPM11265' 0 []}';
BPM11301={'mo' 'BPM11301' 0 []}';
% ==============================================================================
% MARK
% ------------------------------------------------------------------------------
BEGL1F={'mo' 'BEGL1F' 0 []}';
VV11302={'mo' 'VV11302' 0 []}';%vacuum valve                          
L1XFBEG={'mo' 'L1XFBEG' 0 []}';
L1XFEND={'mo' 'L1XFEND' 0 []}';
VV11308={'mo' 'VV11308' 0 []}';%vacuum valve                           
ENDL1F={'mo' 'ENDL1F' 0 []}';
% ==============================================================================
% BEAMLINEs
% ------------------------------------------------------------------------------
L1C=[D9,DAQA,QFL1,QFL1,DAQA,D9,DAQA,QDL1,QDL1,DAQA];
K11_1B_FULL=[K11_1B1,XC11104,YC11105,K11_1B2];
K11_1C_FULL=[K11_1C1,XC11140,YC11141,K11_1C2];
K11_1D_FULL=[K11_1D];
K11_2A_FULL=[K11_2A1,XC11202,K11_2A2,YC11203,K11_2A3];
K11_2B_FULL=[K11_2B];
K11_2C_FULL=[K11_2C1,XC11272,YC11273,K11_2C2];
L1XF_FULL=[L1XFBEG,L1XF__1,XC11304,YC11305,L1XF__2,L1XFEND];
QA11132_FULL=[QA11132,BPM11132,QA11132];
Q11201_FULL=[Q11201,BPM11201,Q11201];
QA11265_FULL=[QA11265,BPM11265,QA11265];
Q11301_FULL=[Q11301,BPM11301,Q11301];
L1F=[BEGL1F,K11_1B_FULL,DAQA1,QA11132_FULL,DAQA2,K11_1C_FULL,K11_1D_FULL,DAQ1,Q11201_FULL,DAQ2,K11_2A_FULL,K11_2B_FULL,DAQA3,QA11265_FULL,DAQA4,K11_2C_FULL,DAQA5,Q11301_FULL,DAQA6,VV11302,DL1X,L1XF_FULL,DM10A,VV11308,ENDL1F];
% ==============================================================================

% *** OPTICS=FACET2-28JAN2026 ***
% ==============================================================================
% Modification History
% ------------------------------------------------------------------------------
% 26-JUL-2023, M. Woodley
%  * clean up drift length definitions to avoid negative drift lengths
% 17-OCT-2022, G. White - Survey co-ordinates for BC11 from Georg
%  * All previous changes verified, other than the following:
%    BPM11333 z=1042.6467 -> 1042.549
%    SQ11340 z=1043.7724 -> 1044.804
%  * XLS sheet with data stored in doc directory
% 08-JUL-2022, G. White - Position and name changes from B. O'Shea
%  * QM11312 z = 1039.1048 -> 1039.1590 (dz = +54.2 mm)
%  * CQ11317 z = 1040.5585 -> 1040.6694 (dz = +110.9 mm)
%  * YC11321 z = 1040.6694 -> 1041.0487 (dz = +373.9 mm)
%  * CE11345 -> CE11334 z = 1043.6915 -> 1044.2771 (dz = +585.6 mm)
%  * CQ11352 z = 1045.0678 -> 1045.0798 (dz = +12.0 mm)
%  * QM11358 z = 1046.547 -> 1046.4623 (dz = -84.7 mm)
%  * BL11359 -> BL11356 z = 1046.8851 -> 1046.7396 (dz = -145.5 mm)
% 08-MAR-2022, G. White
%  * put XC11398, YC11399, and SQ11340 back in
% ------------------------------------------------------------------------------
% 22-JAN-2021, G. White - Changes per B. O'Shea for edge radiation equipment
%  * Move QM11312 27cm d/s
%  * Move CQ11317 50cm d/s
%  * Move CQ11352 50cm u/s
%  * Move QM11358 10cm d/s
%  * Move CE11345 13cm d/s
%  * Rematched into L2 with Q358/Q362/Q393/Q401
% ------------------------------------------------------------------------------
% 28-FEB-2020, G. White - changes after visual inspection of beamline
%  * XC11398, YC11399 & SQ11340 devices removed
%  * BL11357 moved downstream of QM11358 and changed unit # to 359
%  * Changed CE11334 unit # to 345
% 23-AUG-2018, M. Woodley
%  * quadrupole K1 values from FACET2e_baseline.mat
% ------------------------------------------------------------------------------
% 11-APR-2017, M. Woodley
%  * use LCLS-II Phase 1 BC1B definitions for BC11 chicane
%  * remove BPM11374
%  * add PR11375 and TD11390
% 31-JAN-2017, M. Woodley
%  * from FACET2e.xsif
% ------------------------------------------------------------------------------
% ==============================================================================
% SBEN
% ------------------------------------------------------------------------------
% BC11
% - use series approximation for sinc(x)=sin(x)/x to allow AB11=0
% GB11  : chicane bend gap height (m)
% ZB11  : chicane bend "Z" length (m)
% AB11  : chicane bend angle (rad)
% LB11  : chicane bend path length (m)
% AB11s : "short" half chicane bend angle (rad)
% LB11s : "short" half chicane bend path length (m)
% AB11l : "long" half chicane bend angle (rad)
% LB11l : "long" half chicane bend path length (m)
GB11 =  0.043;
ZB11 =  0.2032;
AB11 =  0.09410384256 ;%0.094
AB11_2 =  AB11*AB11;
AB11_4 =  AB11_2*AB11_2;
AB11_6 =  AB11_4*AB11_2;
SINCAB11 =  1-AB11_2/6+AB11_4/120-AB11_6/5040 ;%~sinc(AB11)=sin(AB11)/AB11
LB11 =  ZB11/SINCAB11;
AB11S =  asin(sin(AB11)/2);
AB11S_2 =  AB11S*AB11S;
AB11S_4 =  AB11S_2*AB11S_2;
AB11S_6 =  AB11S_4*AB11S_2;
SINCAB11S =  1-AB11S_2/6+AB11S_4/120-AB11S_6/5040 ;%~sinc(AB11s)=sin(AB11s)/AB11s
LB11S =  ZB11/2/SINCAB11S;
AB11L =  AB11-AB11S;
LB11L =  LB11-LB11S;
BCX11314A={'be' 'BCX11314' LB11S [-AB11S GB11/2 0 0 0.387 0 0]}';
BCX11314B={'be' 'BCX11314' LB11L [-AB11L GB11/2 0 -AB11 0 0.387 0]}';
BCX11331A={'be' 'BCX11331' LB11L [+AB11L GB11/2 +AB11 0 0.387 0 0]}';
BCX11331B={'be' 'BCX11331' LB11S [+AB11S GB11/2 0 0 0 0.387 0]}';
BCX11338A={'be' 'BCX11338' LB11S [+AB11S GB11/2 0 0 0.387 0 0]}';
BCX11338B={'be' 'BCX11338' LB11L [+AB11L GB11/2 0 +AB11 0 0.387 0]}';
BCX11355A={'be' 'BCX11355' LB11L [-AB11L GB11/2 -AB11 0 0.387 0 0]}';
BCX11355B={'be' 'BCX11355' LB11S [-AB11S GB11/2 0 0 0 0.387 0]}';
% define unsplit SBENs for BMAD ... not used by MAD
BCX11314={'be' 'BCX11314' LB11 [-AB11 GB11/2 0 -AB11 0.387 0.387 0]}';
BCX11331={'be' 'BCX1133' LB11 [+AB11 GB11/2 +AB11 0 0.387 0.387 0]}';
BCX11338={'be' 'BCX11338' LB11 [+AB11 GB11/2 0 +AB11 0.387 0.387 0]}';
BCX11355={'be' 'BCX11355' LB11 [-AB11 GB11/2 -AB11 0 0.387 0.387 0]}';
% ==============================================================================
% QUAD
% ------------------------------------------------------------------------------
% electron
KQM11312 =  3.409492877395;
KCQ11317 =  0.0;
KSQ11340 =  0.0;
KCQ11352 =  0.0;
QM11312={'qu' 'QM11312' LQE/2 [KQM11312 0]}';
CQ11317={'qu' 'CQ11317' LQC/2 [KCQ11317 0]}';
SQ11340={'qu' 'SQ11340' LQC/2 [KSQ11340 pi/4]}';
CQ11352={'qu' 'CQ11352' LQC/2 [KCQ11352 0]}';
% common
KQM11358 =  -6.1931968;
KQM11362 =   9.35248294;
KQM11393 =  -6.04847226;
KQ11401 =   5.40858133;
QM11358={'qu' 'QM11358' LQE/2 [QSIGN*(KQM11358) 0]}';
QM11362={'qu' 'QM11362' LQE/2 [QSIGN*(KQM11362) 0]}';
QM11393={'qu' 'QM11393' LQE/2 [QSIGN*(KQM11393) 0]}';
Q11401={'qu' 'Q11401' LQE/2 [QSIGN*(KQ11401) 0]}';
% ==============================================================================
% DRIF
% ------------------------------------------------------------------------------
% clean up drift length definitions to avoid negative drift lengths
DM11={'dr' '' 0.317889 []}';
DM12={'dr' '' 0.344401 []}';
D11OA1={'dr' '' 0.589134/cos(AB11) []}';
D11OA2={'dr' '' 0.266547/cos(AB11) []}';
D11OB={'dr' '' 0.325510/cos(AB11) []}';
D11OC={'dr' '' 1.146186/cos(AB11) []}';
DDG11={'dr' '' 0.150940 []}';
DDG12={'dr' '' 0.147487 []}';
DDG13={'dr' '' 0.297279 []}';
DDG14={'dr' '' 0.234494 []}';
D11OD1A={'dr' '' 0.597456/cos(AB11) []}';
D11OD1B={'dr' '' 0.248217/cos(AB11) []}';
D11OD2={'dr' '' 0.473113/cos(AB11) []}';
D11OE={'dr' '' 0.168243/cos(AB11) []}';
D11OF={'dr' '' 0.732827/cos(AB11) []}';
DM13A={'dr' '' 0.156951 []}';
DM13B={'dr' '' 0.182399 []}';
DM14A={'dr' '' 0.223924 []}';
DM14B={'dr' '' 0.17946475+145.5E-3 []}';%0.17886475
DM15A={'dr' '' 0.352146 []}';%0.351546
DM15B={'dr' '' 2.155039035421 []}';
DM15C={'dr' '' 1.81845801814 []}';
DM15D={'dr' '' 0.488187 []}';%0.487587
DM16A={'dr' '' 0.5760517 []}';%0.5754517
DM16B={'dr' '' 0.23128+2.0E-07 []}';%tweak to make sector boundaries line up
DM17={'dr' '' 0.027 []}';
% ==============================================================================
% XCORs and YCORs
% ------------------------------------------------------------------------------
YC11321={'mo' 'YC11321' 0 []}';
YC11365={'mo' 'YC11365' 0 []}';
XC11398={'mo' 'XC11398' 0 []}';
YC11399={'mo' 'YC11399' 0 []}';
% ==============================================================================
% diagnostics
% ------------------------------------------------------------------------------
% BPMs (LCLS-II type designations)
BPM11312={'mo' 'BPM11312' 0 []}';
BPM11333={'mo' 'BPM11333' 0 []}';
BPM11358={'mo' 'BPM11358' 0 []}';%was BPM11357
BPM11362={'mo' 'BPM11362' 0 []}';%was BPM11363
BPM11393={'mo' 'BPM11393' 0 []}';
BPM11401={'mo' 'BPM11401' 0 []}';
% misc
PR11316={'mo' 'PR11316' 0 []}';
PR11334={'mo' 'PR11334' 0 []}';
PR11335={'mo' 'PR11335' 0 []}';
PR11342={'mo' 'PR11342' 0 []}';
PR11357={'mo' 'PR11357' 0 []}';
BL11356={'mo' 'BL11356' 0 []}';
IM11360={'mo' 'IM11360' 0 []}';
PR11375={'mo' 'PR11375' 0 []}';
TD11390={'mo' 'TD11390' 0 []}';
% ==============================================================================
% collimators
% ------------------------------------------------------------------------------
CE11334={'dr' 'CE11334' 0 []}';%energy collimator
% ==============================================================================
% MARK
% ------------------------------------------------------------------------------
BEGBC11_1={'mo' 'BEGBC11_1' 0 []}';
BC11CBEG={'mo' 'BC11CBEG' 0 []}';
CNT1B={'mo' 'CNT1B' 0 []}';
BC11CEND={'mo' 'BC11CEND' 0 []}';
ENDBC11_1={'mo' 'ENDBC11_1' 0 []}';
BEGBC11_2={'mo' 'BEGBC11_2' 0 []}';
ENDBC11_2={'mo' 'ENDBC11_2' 0 []}';
% ==============================================================================
% BEAMLINEs
% ------------------------------------------------------------------------------
BCX11314_FULL=[BCX11314A,BCX11314B];
CQ11317_FULL=[CQ11317,CQ11317];
BCX11331_FULL=[BCX11331A,BCX11331B];
BCX11338_FULL=[BCX11338A,BCX11338B];
SQ11340_FULL=[SQ11340,SQ11340];
CQ11352_FULL=[CQ11352,CQ11352];
BCX11355_FULL=[BCX11355A,BCX11355B];
BC11C=[BC11CBEG,BCX11314_FULL,D11OA1,PR11316,D11OA2,CQ11317_FULL,D11OB,YC11321,D11OC,BCX11331_FULL,DDG11,BPM11333,DDG12,PR11334,DDG13,PR11335,DDG14,BCX11338_FULL,D11OD1A,PR11342,D11OD1B,CE11334,D11OD2,SQ11340_FULL,D11OE,CQ11352_FULL,D11OF,BCX11355_FULL,CNT1B,BC11CEND];
QM11312_FULL=[QM11312,BPM11312,QM11312];
BC11_1=[BEGBC11_1,DM11,QM11312_FULL,DM12,BC11C,ENDBC11_1];
QM11358_FULL=[QM11358,BPM11358,QM11358];
QM11362_FULL=[QM11362,BPM11362,QM11362];
QM11393_FULL=[QM11393,BPM11393,QM11393];
Q11401_FULL=[Q11401,BPM11401,Q11401];
BC11_2=[BEGBC11_2,DM13A,PR11357,DM13B,QM11358_FULL,DM14A,BL11356,IM11360,DM14B,QM11362_FULL,DM15A,YC11365,DM15B,PR11375,DM15C,TD11390,DM15D,QM11393_FULL,DM16A,XC11398,YC11399,DM16B,Q11401_FULL,DM17,ENDBC11_2];
BC11=[BC11_1,BC11_2];
% ==============================================================================

% *** OPTICS=FACET2-28JAN2026 ***
% ==============================================================================
% Modification History
% ------------------------------------------------------------------------------
% 26-APR-2023, G. White
%  * XC14702 dz -0.1352m as per Georg measurements, also:
%  * YC14703 dz -0.0158m 
% 23-AUG-2018, M. Woodley
%  * gradL2, phiL1, gradFB2, and phiFB2 values from FACET2e_baseline.mat;
%    set KlossL2=0
%  * quadrupole K1 values from FACET2e_baseline.mat
% ------------------------------------------------------------------------------
% 31-JAN-2017, M. Woodley
%  * from FACET2e.xsif
% ------------------------------------------------------------------------------
% ==============================================================================
% accelerating structures
% ------------------------------------------------------------------------------
% the L2 S-band linac consists of: 100 x 10' structure @ 25% power
%                                    8 x  7' structure @ 25% power
% ------------------------------------------------------------------------------
FL2 =    0.987979678125;
GRADL2 =   17.995282755779*FL2 ;%MeV/m
PHIL2 =  -38.25/360           ;%rad/2pi
GRADFB2 =   15.170209548027*FL2 ;%MeV/m
PHIFB2 =   60.0/360            ;%rad/2pi
KLOSSL2 =    0                  ;%V/C/m
G11_4 =  GRADL2;
G11_5 =  GRADL2  ;
G11_6 =  GRADL2 ;
G11_7 =  GRADL2 ;
G11_8 =  GRADL2;
G12_1 =  GRADL2  ;
G12_2 =  GRADL2 ;
G12_3 =  GRADL2 ;
G12_4 =  GRADL2;
G12_5 =  GRADL2  ;
G12_6 =  GRADL2 ;
G12_7 =  GRADL2 ;
G12_8 =  GRADL2;
G13_1 =  GRADL2  ;
G13_2 =  GRADL2 ;
G13_3 =  GRADL2 ;
G13_4 =  GRADL2;
G13_5 =  GRADL2  ;
G13_6 =  GRADL2 ;
G13_7 =  GRADL2 ;
G13_8 =  GRADL2;
G14_1 =  GRADL2  ;
G14_2 =  GRADL2 ;
G14_3 =  GRADL2 ;
G14_4 =  GRADFB2;
G14_5 =  GRADFB2 ;
G14_6 =  0;
K11_4A1={'lc' 'K11_4A' 0.3268 [SBANDF P25*G11_4*0.3268 PHIL2*TWOPI]}';
K11_4A2={'lc' 'K11_4A' 0.2500 [SBANDF P25*G11_4*0.2500 PHIL2*TWOPI]}';
K11_4A3={'lc' 'K11_4A' 1.5926 [SBANDF P25*G11_4*1.5926 PHIL2*TWOPI]}';
K11_4B={'lc' 'K11_4B' DLWL7 [SBANDF P25*G11_4*DLWL7 PHIL2*TWOPI]}';
K11_4C={'lc' 'K11_4C' DLWL10 [SBANDF P25*G11_4*DLWL10 PHIL2*TWOPI]}';
K11_4D={'lc' 'K11_4D' DLWL7 [SBANDF P25*G11_4*DLWL7 PHIL2*TWOPI]}';
K11_5A1={'lc' 'K11_5A' 0.6689 [SBANDF P25*G11_5*0.6689 PHIL2*TWOPI]}';
K11_5A2={'lc' 'K11_5A' 2.3752 [SBANDF P25*G11_5*2.3752 PHIL2*TWOPI]}';
K11_5B={'lc' 'K11_5B' DLWL7 [SBANDF P25*G11_5*DLWL7 PHIL2*TWOPI]}';
K11_5C={'lc' 'K11_5C' DLWL10 [SBANDF P25*G11_5*DLWL10 PHIL2*TWOPI]}';
K11_5D={'lc' 'K11_5D' DLWL10 [SBANDF P25*G11_5*DLWL10 PHIL2*TWOPI]}';
K11_6A1={'lc' 'K11_6A' 0.3280 [SBANDF P25*G11_6*0.3280 PHIL2*TWOPI]}';
K11_6A2={'lc' 'K11_6A' 0.2500 [SBANDF P25*G11_6*0.2500 PHIL2*TWOPI]}';
K11_6A3={'lc' 'K11_6A' 1.5914 [SBANDF P25*G11_6*1.5914 PHIL2*TWOPI]}';
K11_6B={'lc' 'K11_6B' DLWL10 [SBANDF P25*G11_6*DLWL10 PHIL2*TWOPI]}';
K11_6C={'lc' 'K11_6C' DLWL10 [SBANDF P25*G11_6*DLWL10 PHIL2*TWOPI]}';
K11_6D={'lc' 'K11_6D' DLWL7 [SBANDF P25*G11_6*DLWL7 PHIL2*TWOPI]}';
K11_7A={'lc' 'K11_7A' DLWL10 [SBANDF P25*G11_7*DLWL10 PHIL2*TWOPI]}';
K11_7B={'lc' 'K11_7B' DLWL10 [SBANDF P25*G11_7*DLWL10 PHIL2*TWOPI]}';
K11_7C={'lc' 'K11_7C' DLWL10 [SBANDF P25*G11_7*DLWL10 PHIL2*TWOPI]}';
K11_7D={'lc' 'K11_7D' DLWL7 [SBANDF P25*G11_7*DLWL7 PHIL2*TWOPI]}';
K11_8A1={'lc' 'K11_8A' 0.38633 [SBANDF P25*G11_8*0.38633 PHIL2*TWOPI]}';
K11_8A2={'lc' 'K11_8A' 2.65777 [SBANDF P25*G11_8*2.65777 PHIL2*TWOPI]}';
K11_8B={'lc' 'K11_8B' DLWL10 [SBANDF P25*G11_8*DLWL10 PHIL2*TWOPI]}';
K11_8C={'lc' 'K11_8C' DLWL10 [SBANDF P25*G11_8*DLWL10 PHIL2*TWOPI]}';
K11_8D={'lc' 'K11_8D' DLWL10 [SBANDF P25*G11_8*DLWL10 PHIL2*TWOPI]}';
K12_1A={'lc' 'K12_1A' DLWL10 [SBANDF P25*G12_1*DLWL10 PHIL2*TWOPI]}';
K12_1B={'lc' 'K12_1B' DLWL10 [SBANDF P25*G12_1*DLWL10 PHIL2*TWOPI]}';
K12_1C={'lc' 'K12_1C' DLWL10 [SBANDF P25*G12_1*DLWL10 PHIL2*TWOPI]}';
K12_1D={'lc' 'K12_1D' DLWL10 [SBANDF P25*G12_1*DLWL10 PHIL2*TWOPI]}';
K12_2A1={'lc' 'K12_2A' 0.3280 [SBANDF P25*G12_2*0.3280 PHIL2*TWOPI]}';
K12_2A2={'lc' 'K12_2A' 0.2500 [SBANDF P25*G12_2*0.2500 PHIL2*TWOPI]}';
K12_2A3={'lc' 'K12_2A' 1.5914 [SBANDF P25*G12_2*1.5914 PHIL2*TWOPI]}';
K12_2B={'lc' 'K12_2B' DLWL10 [SBANDF P25*G12_2*DLWL10 PHIL2*TWOPI]}';
K12_2C={'lc' 'K12_2C' DLWL10 [SBANDF P25*G12_2*DLWL10 PHIL2*TWOPI]}';
K12_2D={'lc' 'K12_2D' DLWL10 [SBANDF P25*G12_2*DLWL10 PHIL2*TWOPI]}';
K12_3A1={'lc' 'K12_3A' 0.3312 [SBANDF P25*G12_3*0.3312 PHIL2*TWOPI]}';
K12_3A2={'lc' 'K12_3A' 0.2500 [SBANDF P25*G12_3*0.2500 PHIL2*TWOPI]}';
K12_3A3={'lc' 'K12_3A' 2.4629 [SBANDF P25*G12_3*2.4629 PHIL2*TWOPI]}';
K12_3B={'lc' 'K12_3B' DLWL10 [SBANDF P25*G12_3*DLWL10 PHIL2*TWOPI]}';
K12_3C={'lc' 'K12_3C' DLWL10 [SBANDF P25*G12_3*DLWL10 PHIL2*TWOPI]}';
K12_3D={'lc' 'K12_3D' DLWL10 [SBANDF P25*G12_3*DLWL10 PHIL2*TWOPI]}';
K12_4A1={'lc' 'K12_4A' 0.3268 [SBANDF P25*G12_4*0.3268 PHIL2*TWOPI]}';
K12_4A2={'lc' 'K12_4A' 0.2500 [SBANDF P25*G12_4*0.2500 PHIL2*TWOPI]}';
K12_4A3={'lc' 'K12_4A' 2.4673 [SBANDF P25*G12_4*2.4673 PHIL2*TWOPI]}';
K12_4B={'lc' 'K12_4B' DLWL10 [SBANDF P25*G12_4*DLWL10 PHIL2*TWOPI]}';
K12_4C={'lc' 'K12_4C' DLWL10 [SBANDF P25*G12_4*DLWL10 PHIL2*TWOPI]}';
K12_4D={'lc' 'K12_4D' DLWL10 [SBANDF P25*G12_4*DLWL10 PHIL2*TWOPI]}';
K12_5A1={'lc' 'K12_5A' 0.3324 [SBANDF P25*G12_5*0.3324 PHIL2*TWOPI]}';
K12_5A2={'lc' 'K12_5A' 0.2500 [SBANDF P25*G12_5*0.2500 PHIL2*TWOPI]}';
K12_5A3={'lc' 'K12_5A' 2.4617 [SBANDF P25*G12_5*2.4617 PHIL2*TWOPI]}';
K12_5B={'lc' 'K12_5B' DLWL10 [SBANDF P25*G12_5*DLWL10 PHIL2*TWOPI]}';
K12_5C={'lc' 'K12_5C' DLWL10 [SBANDF P25*G12_5*DLWL10 PHIL2*TWOPI]}';
K12_5D={'lc' 'K12_5D' DLWL10 [SBANDF P25*G12_5*DLWL10 PHIL2*TWOPI]}';
K12_6A1={'lc' 'K12_6A' 0.3280 [SBANDF P25*G12_6*0.3280 PHIL2*TWOPI]}';
K12_6A2={'lc' 'K12_6A' 0.4330 [SBANDF P25*G12_6*0.4330 PHIL2*TWOPI]}';
K12_6A3={'lc' 'K12_6A' 2.2831 [SBANDF P25*G12_6*2.2831 PHIL2*TWOPI]}';
K12_6B={'lc' 'K12_6B' DLWL10 [SBANDF P25*G12_6*DLWL10 PHIL2*TWOPI]}';
K12_6C={'lc' 'K12_6C' DLWL10 [SBANDF P25*G12_6*DLWL10 PHIL2*TWOPI]}';
K12_6D={'lc' 'K12_6D' DLWL10 [SBANDF P25*G12_6*DLWL10 PHIL2*TWOPI]}';
K12_7A1={'lc' 'K12_7A' 0.3336 [SBANDF P25*G12_7*0.3336 PHIL2*TWOPI]}';
K12_7A2={'lc' 'K12_7A' 0.3575 [SBANDF P25*G12_7*0.3575 PHIL2*TWOPI]}';
K12_7A3={'lc' 'K12_7A' 2.3530 [SBANDF P25*G12_7*2.3530 PHIL2*TWOPI]}';
K12_7B={'lc' 'K12_7B' DLWL10 [SBANDF P25*G12_7*DLWL10 PHIL2*TWOPI]}';
K12_7C={'lc' 'K12_7C' DLWL10 [SBANDF P25*G12_7*DLWL10 PHIL2*TWOPI]}';
K12_7D={'lc' 'K12_7D' DLWL10 [SBANDF P25*G12_7*DLWL10 PHIL2*TWOPI]}';
K12_8A1={'lc' 'K12_8A' 0.3292 [SBANDF P25*G12_8*0.3292 PHIL2*TWOPI]}';
K12_8A2={'lc' 'K12_8A' 0.4032 [SBANDF P25*G12_8*0.4032 PHIL2*TWOPI]}';
K12_8A3={'lc' 'K12_8A' 2.3117 [SBANDF P25*G12_8*2.3117 PHIL2*TWOPI]}';
K12_8B={'lc' 'K12_8B' DLWL10 [SBANDF P25*G12_8*DLWL10 PHIL2*TWOPI]}';
K12_8C={'lc' 'K12_8C' DLWL10 [SBANDF P25*G12_8*DLWL10 PHIL2*TWOPI]}';
K12_8D1={'lc' 'K12_8D' 2.3869 [SBANDF P25*G12_8*2.3869 PHIL2*TWOPI]}';
K12_8D2={'lc' 'K12_8D' 0.2500 [SBANDF P25*G12_8*0.2500 PHIL2*TWOPI]}';
K12_8D3={'lc' 'K12_8D' 0.4072 [SBANDF P25*G12_8*0.4072 PHIL2*TWOPI]}';
K13_1A={'lc' 'K13_1A' DLWL10 [SBANDF P25*G13_1*DLWL10 PHIL2*TWOPI]}';
K13_1B={'lc' 'K13_1B' DLWL10 [SBANDF P25*G13_1*DLWL10 PHIL2*TWOPI]}';
K13_1C={'lc' 'K13_1C' DLWL10 [SBANDF P25*G13_1*DLWL10 PHIL2*TWOPI]}';
K13_1D={'lc' 'K13_1D' DLWL10 [SBANDF P25*G13_1*DLWL10 PHIL2*TWOPI]}';
K13_2A1={'lc' 'K13_2A' 0.3256 [SBANDF P25*G13_2*0.3256 PHIL2*TWOPI]}';
K13_2A2={'lc' 'K13_2A' 0.3846 [SBANDF P25*G13_2*0.3846 PHIL2*TWOPI]}';
K13_2A3={'lc' 'K13_2A' 2.3339 [SBANDF P25*G13_2*2.3339 PHIL2*TWOPI]}';
K13_2B={'lc' 'K13_2B' DLWL10 [SBANDF P25*G13_2*DLWL10 PHIL2*TWOPI]}';
K13_2C={'lc' 'K13_2C' DLWL10 [SBANDF P25*G13_2*DLWL10 PHIL2*TWOPI]}';
K13_2D={'lc' 'K13_2D' DLWL10 [SBANDF P25*G13_2*DLWL10 PHIL2*TWOPI]}';
K13_3A1={'lc' 'K13_3A' 0.3312 [SBANDF P25*G13_3*0.3312 PHIL2*TWOPI]}';
K13_3A2={'lc' 'K13_3A' 0.3822 [SBANDF P25*G13_3*0.3822 PHIL2*TWOPI]}';
K13_3A3={'lc' 'K13_3A' 2.3307 [SBANDF P25*G13_3*2.3307 PHIL2*TWOPI]}';
K13_3B={'lc' 'K13_3B' DLWL10 [SBANDF P25*G13_3*DLWL10 PHIL2*TWOPI]}';
K13_3C={'lc' 'K13_3C' DLWL10 [SBANDF P25*G13_3*DLWL10 PHIL2*TWOPI]}';
K13_3D={'lc' 'K13_3D' DLWL10 [SBANDF P25*G13_3*DLWL10 PHIL2*TWOPI]}';
K13_4A1={'lc' 'K13_4A' 0.3268 [SBANDF P25*G13_4*0.3268 PHIL2*TWOPI]}';
K13_4A2={'lc' 'K13_4A' 0.4501 [SBANDF P25*G13_4*0.4501 PHIL2*TWOPI]}';
K13_4A3={'lc' 'K13_4A' 2.2672 [SBANDF P25*G13_4*2.2672 PHIL2*TWOPI]}';
K13_4B={'lc' 'K13_4B' DLWL10 [SBANDF P25*G13_4*DLWL10 PHIL2*TWOPI]}';
K13_4C={'lc' 'K13_4C' DLWL10 [SBANDF P25*G13_4*DLWL10 PHIL2*TWOPI]}';
K13_4D={'lc' 'K13_4D' DLWL10 [SBANDF P25*G13_4*DLWL10 PHIL2*TWOPI]}';
K13_5A1={'lc' 'K13_5A' 0.3324 [SBANDF P25*G13_5*0.3324 PHIL2*TWOPI]}';
K13_5A2={'lc' 'K13_5A' 0.2500 [SBANDF P25*G13_5*0.2500 PHIL2*TWOPI]}';
K13_5A3={'lc' 'K13_5A' 2.4617 [SBANDF P25*G13_5*2.4617 PHIL2*TWOPI]}';
K13_5B={'lc' 'K13_5B' DLWL10 [SBANDF P25*G13_5*DLWL10 PHIL2*TWOPI]}';
K13_5C={'lc' 'K13_5C' DLWL10 [SBANDF P25*G13_5*DLWL10 PHIL2*TWOPI]}';
K13_5D={'lc' 'K13_5D' DLWL10 [SBANDF P25*G13_5*DLWL10 PHIL2*TWOPI]}';
K13_6A1={'lc' 'K13_6A' 0.3280 [SBANDF P25*G13_6*0.3280 PHIL2*TWOPI]}';
K13_6A2={'lc' 'K13_6A' 0.2500 [SBANDF P25*G13_6*0.2500 PHIL2*TWOPI]}';
K13_6A3={'lc' 'K13_6A' 2.4661 [SBANDF P25*G13_6*2.4661 PHIL2*TWOPI]}';
K13_6B={'lc' 'K13_6B' DLWL10 [SBANDF P25*G13_6*DLWL10 PHIL2*TWOPI]}';
K13_6C={'lc' 'K13_6C' DLWL10 [SBANDF P25*G13_6*DLWL10 PHIL2*TWOPI]}';
K13_6D={'lc' 'K13_6D' DLWL10 [SBANDF P25*G13_6*DLWL10 PHIL2*TWOPI]}';
K13_7A1={'lc' 'K13_7A' 0.3336 [SBANDF P25*G13_7*0.3336 PHIL2*TWOPI]}';
K13_7A2={'lc' 'K13_7A' 0.2500 [SBANDF P25*G13_7*0.2500 PHIL2*TWOPI]}';
K13_7A3={'lc' 'K13_7A' 2.4605 [SBANDF P25*G13_7*2.4605 PHIL2*TWOPI]}';
K13_7B={'lc' 'K13_7B' DLWL10 [SBANDF P25*G13_7*DLWL10 PHIL2*TWOPI]}';
K13_7C={'lc' 'K13_7C' DLWL10 [SBANDF P25*G13_7*DLWL10 PHIL2*TWOPI]}';
K13_7D={'lc' 'K13_7D' DLWL10 [SBANDF P25*G13_7*DLWL10 PHIL2*TWOPI]}';
K13_8A1={'lc' 'K13_8A' 0.3292 [SBANDF P25*G13_8*0.3292 PHIL2*TWOPI]}';
K13_8A2={'lc' 'K13_8A' 0.4064 [SBANDF P25*G13_8*0.4064 PHIL2*TWOPI]}';
K13_8A3={'lc' 'K13_8A' 2.3085 [SBANDF P25*G13_8*2.3085 PHIL2*TWOPI]}';
K13_8B={'lc' 'K13_8B' DLWL10 [SBANDF P25*G13_8*DLWL10 PHIL2*TWOPI]}';
K13_8C={'lc' 'K13_8C' DLWL10 [SBANDF P25*G13_8*DLWL10 PHIL2*TWOPI]}';
K13_8D1={'lc' 'K13_8D' 2.3869 [SBANDF P25*G13_8*2.3869 PHIL2*TWOPI]}';
K13_8D2={'lc' 'K13_8D' 0.2500 [SBANDF P25*G13_8*0.2500 PHIL2*TWOPI]}';
K13_8D3={'lc' 'K13_8D' 0.4072 [SBANDF P25*G13_8*0.4072 PHIL2*TWOPI]}';
K14_1A={'lc' 'K14_1A' DLWL10 [SBANDF P25*G14_1*DLWL10 PHIL2*TWOPI]}';
K14_1B={'lc' 'K14_1B' DLWL10 [SBANDF P25*G14_1*DLWL10 PHIL2*TWOPI]}';
K14_1C={'lc' 'K14_1C' DLWL10 [SBANDF P25*G14_1*DLWL10 PHIL2*TWOPI]}';
K14_1D={'lc' 'K14_1D' DLWL10 [SBANDF P25*G14_1*DLWL10 PHIL2*TWOPI]}';
K14_2A1={'lc' 'K14_2A' 0.3256 [SBANDF P25*G14_2*0.3256 PHIL2*TWOPI]}';
K14_2A2={'lc' 'K14_2A' 0.2500 [SBANDF P25*G14_2*0.2500 PHIL2*TWOPI]}';
K14_2A3={'lc' 'K14_2A' 2.4685 [SBANDF P25*G14_2*2.4685 PHIL2*TWOPI]}';
K14_2B={'lc' 'K14_2B' DLWL10 [SBANDF P25*G14_2*DLWL10 PHIL2*TWOPI]}';
K14_2C={'lc' 'K14_2C' DLWL10 [SBANDF P25*G14_2*DLWL10 PHIL2*TWOPI]}';
K14_2D={'lc' 'K14_2D' DLWL10 [SBANDF P25*G14_2*DLWL10 PHIL2*TWOPI]}';
K14_3A1={'lc' 'K14_3A' 0.3312 [SBANDF P25*G14_3*0.3312 PHIL2*TWOPI]}';
K14_3A2={'lc' 'K14_3A' 0.4012 [SBANDF P25*G14_3*0.4012 PHIL2*TWOPI]}';
K14_3A3={'lc' 'K14_3A' 2.3117 [SBANDF P25*G14_3*2.3117 PHIL2*TWOPI]}';
K14_3B={'lc' 'K14_3B' DLWL10 [SBANDF P25*G14_3*DLWL10 PHIL2*TWOPI]}';
K14_3C={'lc' 'K14_3C' DLWL10 [SBANDF P25*G14_3*DLWL10 PHIL2*TWOPI]}';
K14_3D={'lc' 'K14_3D' DLWL10 [SBANDF P25*G14_3*DLWL10 PHIL2*TWOPI]}';
K14_4A1={'lc' 'K14_4A' 0.3268 [SBANDF P25*G14_4*0.3268 -PHIFB2*TWOPI]}';
K14_4A2={'lc' 'K14_4A' 0.2500 [SBANDF P25*G14_4*0.2500 -PHIFB2*TWOPI]}';
K14_4A3={'lc' 'K14_4A' 2.4673 [SBANDF P25*G14_4*2.4673 -PHIFB2*TWOPI]}';
K14_4B={'lc' 'K14_4B' DLWL10 [SBANDF P25*G14_4*DLWL10 -PHIFB2*TWOPI]}';
K14_4C={'lc' 'K14_4C' DLWL10 [SBANDF P25*G14_4*DLWL10 -PHIFB2*TWOPI]}';
K14_4D={'lc' 'K14_4D' DLWL10 [SBANDF P25*G14_4*DLWL10 -PHIFB2*TWOPI]}';
K14_5A1={'lc' 'K14_5A' 0.3324 [SBANDF P25*G14_5*0.3324 +PHIFB2*TWOPI]}';
K14_5A2={'lc' 'K14_5A' 0.2500 [SBANDF P25*G14_5*0.2500 +PHIFB2*TWOPI]}';
K14_5A3={'lc' 'K14_5A' 2.4617 [SBANDF P25*G14_5*2.4617 +PHIFB2*TWOPI]}';
K14_5B={'lc' 'K14_5B' DLWL10 [SBANDF P25*G14_5*DLWL10 +PHIFB2*TWOPI]}';
K14_5C={'lc' 'K14_5C' DLWL10 [SBANDF P25*G14_5*DLWL10 +PHIFB2*TWOPI]}';
K14_5D={'lc' 'K14_5D' DLWL10 [SBANDF P25*G14_5*DLWL10 +PHIFB2*TWOPI]}';
K14_6A1={'lc' 'K14_6A' 0.3280 [SBANDF P25*G14_6*0.3280 PHIL2*TWOPI]}';
K14_6A2={'lc' 'K14_6A' 0.4044 [SBANDF P25*G14_6*0.4044 PHIL2*TWOPI]}';
K14_6A3={'lc' 'K14_6A' 2.3117 [SBANDF P25*G14_6*2.3117 PHIL2*TWOPI]}';
K14_6B={'lc' 'K14_6B' DLWL10 [SBANDF P25*G14_6*DLWL10 PHIL2*TWOPI]}';
K14_6C={'lc' 'K14_6C' DLWL10 [SBANDF P25*G14_6*DLWL10 PHIL2*TWOPI]}';
K14_6D1={'lc' 'K14_6D' 2.1969 [SBANDF P25*G14_6*2.1969 PHIL2*TWOPI]}';
K14_6D2={'lc' 'K14_6D' 0.4242 [SBANDF P25*G14_6*0.4242 PHIL2*TWOPI]}';
K14_6D3={'lc' 'K14_6D' 0.423 [SBANDF P25*G14_6*0.423 PHIL2*TWOPI]}';
% define unsplit LCAVs for BMAD ... not used by MAD
K11_4A={'lc' 'K11_4A' DLWL7 [SBANDF P25*G11_4*DLWL7 PHIL2*TWOPI]}';
K11_5A={'lc' 'K11_5A' DLWL10 [SBANDF P25*G11_5*DLWL10 PHIL2*TWOPI]}';
K11_6A={'lc' 'K11_6A' DLWL7 [SBANDF P25*G11_6*DLWL7 PHIL2*TWOPI]}';
K11_8A={'lc' 'K11_8A' DLWL10 [SBANDF P25*G11_8*DLWL10 PHIL2*TWOPI]}';
K12_2A={'lc' 'K12_2A' DLWL7 [SBANDF P25*G12_2*DLWL7 PHIL2*TWOPI]}';
K12_3A={'lc' 'K12_3A' DLWL10 [SBANDF P25*G12_3*DLWL10 PHIL2*TWOPI]}';
K12_4A={'lc' 'K12_4A' DLWL10 [SBANDF P25*G12_4*DLWL10 PHIL2*TWOPI]}';
K12_5A={'lc' 'K12_5A' DLWL10 [SBANDF P25*G12_5*DLWL10 PHIL2*TWOPI]}';
K12_6A={'lc' 'K12_6A' DLWL10 [SBANDF P25*G12_6*DLWL10 PHIL2*TWOPI]}';
K12_7A={'lc' 'K12_7A' DLWL10 [SBANDF P25*G12_7*DLWL10 PHIL2*TWOPI]}';
K12_8A={'lc' 'K12_8A' DLWL10 [SBANDF P25*G12_8*DLWL10 PHIL2*TWOPI]}';
K12_8D={'lc' 'K12_8D' DLWL10 [SBANDF P25*G12_8*DLWL10 PHIL2*TWOPI]}';
K13_2A={'lc' 'K13_2A' DLWL10 [SBANDF P25*G13_2*DLWL10 PHIL2*TWOPI]}';
K13_3A={'lc' 'K13_3A' DLWL10 [SBANDF P25*G13_3*DLWL10 PHIL2*TWOPI]}';
K13_4A={'lc' 'K13_4A' DLWL10 [SBANDF P25*G13_4*DLWL10 PHIL2*TWOPI]}';
K13_5A={'lc' 'K13_5A' DLWL10 [SBANDF P25*G13_5*DLWL10 PHIL2*TWOPI]}';
K13_6A={'lc' 'K13_6A' DLWL10 [SBANDF P25*G13_6*DLWL10 PHIL2*TWOPI]}';
K13_7A={'lc' 'K13_7A' DLWL10 [SBANDF P25*G13_7*DLWL10 PHIL2*TWOPI]}';
K13_8A={'lc' 'K13_8A' DLWL10 [SBANDF P25*G13_8*DLWL10 PHIL2*TWOPI]}';
K13_8D={'lc' 'K13_8D' DLWL10 [SBANDF P25*G13_8*DLWL10 PHIL2*TWOPI]}';
K14_2A={'lc' 'K14_2A' DLWL10 [SBANDF P25*G14_2*DLWL10 PHIL2*TWOPI]}';
K14_3A={'lc' 'K14_3A' DLWL10 [SBANDF P25*G14_3*DLWL10 PHIL2*TWOPI]}';
K14_4A={'lc' 'K14_4A' DLWL10 [SBANDF P25*G14_4*DLWL10 -PHIFB2*TWOPI]}';
K14_5A={'lc' 'K14_5A' DLWL10 [SBANDF P25*G14_5*DLWL10 +PHIFB2*TWOPI]}';
K14_6A={'lc' 'K14_6A' DLWL10 [SBANDF P25*G14_6*DLWL10 PHIL2*TWOPI]}';
K14_6D={'lc' 'K14_6D' DLWL10 [SBANDF P25*G14_6*DLWL10 PHIL2*TWOPI]}';
% ==============================================================================
% QUADs
% ------------------------------------------------------------------------------
KQL2 =  0.582216636003;
QFL2={'qu' 'QFL2' LQE/2 [+KQL2 0]}';
QDL2={'qu' 'QDL2' LQE/2 [-KQL2 0]}';
KQ11501 =  -1.425177997836;
KQ11601 =   0.83103987042;
KQ11701 =  -0.626902111414;
KQ11801 =   0.547907109982;
KQ11901 =  -0.568683821067;
KQ12201 =   0.550253669725;
KQ12301 =  -0.502994501998;
KQ12401 =   0.477741317823;
KQ12501 =  -0.538885575405;
KQ12601 =   0.567408017941;
KQ12701 =  -0.588457931187;
KQ12801 =   0.611779052226;
KQ12901 =  -0.58238706032;
KQ13201 =   0.588339804379;
KQ13301 =  -0.617108363677;
KQ13401 =   0.593846036961;
KQ13501 =  -0.587883009947;
KQ13601 =   0.591692793351;
KQ13701 =  -0.609602245911;
KQ13801 =   0.650373670386;
KQ13901 =  -0.671099686999;
KQ14201 =   0.841694945304;
KQ14301 =  -0.838009749369;
KQ14401 =   0.916835066127;
KQ14501 =  -0.427478762325;
KQ14601 =   0.455009885362;
Q11501={'qu' 'Q11501' LQE/2 [QSIGN*(KQ11501) 0]}';
Q11601={'qu' 'Q11601' LQE/2 [QSIGN*(KQ11601) 0]}';
Q11701={'qu' 'Q11701' LQE/2 [QSIGN*(KQ11701) 0]}';
Q11801={'qu' 'Q11801' LQE/2 [QSIGN*(KQ11801) 0]}';
Q11901={'qu' 'Q11901' LQE/2 [QSIGN*(KQ11901) 0]}';
Q12201={'qu' 'Q12201' LQE/2 [QSIGN*(KQ12201) 0]}';
Q12301={'qu' 'Q12301' LQE/2 [QSIGN*(KQ12301) 0]}';
Q12401={'qu' 'Q12401' LQE/2 [QSIGN*(KQ12401) 0]}';
Q12501={'qu' 'Q12501' LQE/2 [QSIGN*(KQ12501) 0]}';
Q12601={'qu' 'Q12601' LQE/2 [QSIGN*(KQ12601) 0]}';
Q12701={'qu' 'Q12701' LQE/2 [QSIGN*(KQ12701) 0]}';
Q12801={'qu' 'Q12801' LQE/2 [QSIGN*(KQ12801) 0]}';
Q12901={'qu' 'Q12901' LQE/2 [QSIGN*(KQ12901) 0]}';
Q13201={'qu' 'Q13201' LQE/2 [QSIGN*(KQ13201) 0]}';
Q13301={'qu' 'Q13301' LQE/2 [QSIGN*(KQ13301) 0]}';
Q13401={'qu' 'Q13401' LQE/2 [QSIGN*(KQ13401) 0]}';
Q13501={'qu' 'Q13501' LQE/2 [QSIGN*(KQ13501) 0]}';
Q13601={'qu' 'Q13601' LQE/2 [QSIGN*(KQ13601) 0]}';
Q13701={'qu' 'Q13701' LQE/2 [QSIGN*(KQ13701) 0]}';
Q13801={'qu' 'Q13801' LQE/2 [QSIGN*(KQ13801) 0]}';
Q13901={'qu' 'Q13901' LQE/2 [QSIGN*(KQ13901) 0]}';
Q14201={'qu' 'Q14201' LQE/2 [QSIGN*(KQ14201) 0]}';
Q14301={'qu' 'Q14301' LQE/2 [QSIGN*(KQ14301) 0]}';
Q14401={'qu' 'Q14401' LQE/2 [QSIGN*(KQ14401) 0]}';
Q14501={'qu' 'Q14501' LQE/2 [QSIGN*(KQ14501) 0]}';
Q14601={'qu' 'Q14601' LQE/2 [QSIGN*(KQ14601) 0]}';
% ==============================================================================
% drifts
% ------------------------------------------------------------------------------
DL2FODO={'dr' '' 12.2376 []}';
LDAA7 =  DLWL10-DLWL7;
DAA7={'dr' '' LDAA7 []}';
DAA7A={'dr' '' 0.2511 []}';
DAA7B={'dr' '' 0.4238 []}';
DAA7C={'dr' '' LDAA7-(DAA7A{3}+DAA7B{3}) []}';
DAA7D={'dr' '' 0.36121 []}';
DAA7E={'dr' '' LDAA7-DAA7D{3} []}';
DAA7F={'dr' '' 0.34 []}';
DAA7G={'dr' '' 0.2921 []}';
DAA7H={'dr' '' LDAA7-(DAA7F{3}+DAA7G{3}) []}';
DAQ4A={'dr' '' 0.84513 []}';
DAQ4B={'dr' '' 0.27622 []}';
DAQ4C={'dr' '' LDAQ4-(DAQ4A{3}+DAQ4B{3}) []}';
DAA7I={'dr' '' 0.32195 []}';
DAA7J={'dr' '' 0.2794 []}';
DAA7K={'dr' '' LDAA7-(DAA7I{3}+DAA7J{3}) []}';
% ==============================================================================
% XCORs and YCORs
% ------------------------------------------------------------------------------
XC11402={'mo' 'XC11402' 0 []}';
XC11502={'mo' 'XC11502' 0 []}';
XC11602={'mo' 'XC11602' 0 []}';
XC11702={'mo' 'XC11702' 0 []}';
XC11802={'mo' 'XC11802' 0 []}';
XC11900={'mo' 'XC11900' 0 []}';
XC12202={'mo' 'XC12202' 0 []}';
XC12302={'mo' 'XC12302' 0 []}';
XC12402={'mo' 'XC12402' 0 []}';
XC12502={'mo' 'XC12502' 0 []}';
XC12602={'mo' 'XC12602' 0 []}';
XC12702={'mo' 'XC12702' 0 []}';
XC12802={'mo' 'XC12802' 0 []}';
XC12900={'mo' 'XC12900' 0 []}';
XC13202={'mo' 'XC13202' 0 []}';
XC13302={'mo' 'XC13302' 0 []}';
XC13402={'mo' 'XC13402' 0 []}';
XC13502={'mo' 'XC13502' 0 []}';
XC13602={'mo' 'XC13602' 0 []}';
XC13702={'mo' 'XC13702' 0 []}';
XC13802={'mo' 'XC13802' 0 []}';
XC13900={'mo' 'XC13900' 0 []}';
XC14202={'mo' 'XC14202' 0 []}';
XC14302={'mo' 'XC14302' 0 []}';
XC14402={'mo' 'XC14402' 0 []}';
XC14502={'mo' 'XC14502' 0 []}';
XC14602={'mo' 'XC14602' 0 []}';
XC14702={'mo' 'XC14702' 0 []}';
YC11403={'mo' 'YC11403' 0 []}';
YC11503={'mo' 'YC11503' 0 []}';
YC11603={'mo' 'YC11603' 0 []}';
YC11703={'mo' 'YC11703' 0 []}';
YC11803={'mo' 'YC11803' 0 []}';
YC11900={'mo' 'YC11900' 0 []}';
YC12203={'mo' 'YC12203' 0 []}';
YC12303={'mo' 'YC12303' 0 []}';
YC12403={'mo' 'YC12403' 0 []}';
YC12503={'mo' 'YC12503' 0 []}';
YC12603={'mo' 'YC12603' 0 []}';
YC12703={'mo' 'YC12703' 0 []}';
YC12803={'mo' 'YC12803' 0 []}';
YC12900={'mo' 'YC12900' 0 []}';
YC13203={'mo' 'YC13203' 0 []}';
YC13303={'mo' 'YC13303' 0 []}';
YC13403={'mo' 'YC13403' 0 []}';
YC13503={'mo' 'YC13503' 0 []}';
YC13603={'mo' 'YC13603' 0 []}';
YC13703={'mo' 'YC13703' 0 []}';
YC13803={'mo' 'YC13803' 0 []}';
YC13900={'mo' 'YC13900' 0 []}';
YC14203={'mo' 'YC14203' 0 []}';
YC14303={'mo' 'YC14303' 0 []}';
YC14403={'mo' 'YC14403' 0 []}';
YC14503={'mo' 'YC14503' 0 []}';
YC14603={'mo' 'YC14603' 0 []}';
YC14703={'mo' 'YC14703' 0 []}';
% ==============================================================================
% diagnostics
% ------------------------------------------------------------------------------
% BPMs (LCLS-II type designations)
BPM11501={'mo' 'BPM11501' 0 []}';
BPM11601={'mo' 'BPM11601' 0 []}';
BPM11701={'mo' 'BPM11701' 0 []}';
BPM11801={'mo' 'BPM11801' 0 []}';
BPM11901={'mo' 'BPM11901' 0 []}';
BPM12201={'mo' 'BPM12201' 0 []}';
BPM12301={'mo' 'BPM12301' 0 []}';
BPM12401={'mo' 'BPM12401' 0 []}';
BPM12501={'mo' 'BPM12501' 0 []}';
BPM12601={'mo' 'BPM12601' 0 []}';
BPM12701={'mo' 'BPM12701' 0 []}';
BPM12801={'mo' 'BPM12801' 0 []}';
BPM12901={'mo' 'BPM12901' 0 []}';
BPM13201={'mo' 'BPM13201' 0 []}';
BPM13301={'mo' 'BPM13301' 0 []}';
BPM13401={'mo' 'BPM13401' 0 []}';
BPM13501={'mo' 'BPM13501' 0 []}';
BPM13601={'mo' 'BPM13601' 0 []}';
BPM13701={'mo' 'BPM13701' 0 []}';
BPM13801={'mo' 'BPM13801' 0 []}';
BPM13901={'mo' 'BPM13901' 0 []}';
BPM14201={'mo' 'BPM14201' 0 []}';
BPM14301={'mo' 'BPM14301' 0 []}';
BPM14401={'mo' 'BPM14401' 0 []}';
BPM14501={'mo' 'BPM14501' 0 []}';
BPM14601={'mo' 'BPM14601' 0 []}';
% misc
WS11444={'mo' 'WS11444' 0 []}';
WS11614={'mo' 'WS11614' 0 []}';
WS11744={'mo' 'WS11744' 0 []}';
WS12214={'mo' 'WS12214' 0 []}';
% ==============================================================================
% miscellaneous diagnostics, collimators, MARKERs, etc.
% ------------------------------------------------------------------------------
BEGL2F={'mo' 'BEGL2F' 0 []}';
LI11STRT={'mo' 'LI11STRT' 0 []}';
LI11END={'mo' 'LI11END' 0 []}';
LI12BEG={'mo' 'LI12BEG' 0 []}';
LI12END={'mo' 'LI12END' 0 []}';
LI13BEG={'mo' 'LI13BEG' 0 []}';
LI13END={'mo' 'LI13END' 0 []}';
LI14BEG={'mo' 'LI14BEG' 0 []}';
LI14TERM={'mo' 'LI14TERM' 0 []}';
ENDL2F={'mo' 'ENDL2F' 0 []}';
% ==============================================================================
% BEAMLINEs
% ------------------------------------------------------------------------------
L2C=[QFL2,DL2FODO,QDL2,QDL2,DL2FODO,QFL2];
% ------------------------------------------------------------------------------
K11_4A_FULL=[K11_4A1,XC11402,K11_4A2,YC11403,K11_4A3];
K11_4B_FULL=[K11_4B];
K11_4C_FULL=[K11_4C];
K11_4D_FULL=[K11_4D];
K11_5A_FULL=[K11_5A1,YC11503,K11_5A2];
K11_5B_FULL=[K11_5B];
K11_5C_FULL=[K11_5C];
K11_5D_FULL=[K11_5D];
K11_6A_FULL=[K11_6A1,XC11602,K11_6A2,YC11603,K11_6A3];
K11_6B_FULL=[K11_6B];
K11_6C_FULL=[K11_6C];
K11_6D_FULL=[K11_6D];
K11_7A_FULL=[K11_7A];
K11_7B_FULL=[K11_7B];
K11_7C_FULL=[K11_7C];
K11_7D_FULL=[K11_7D];
K11_8A_FULL=[K11_8A1,YC11803,K11_8A2];
K11_8B_FULL=[K11_8B];
K11_8C_FULL=[K11_8C];
K11_8D_FULL=[K11_8D];
K11_4=[K11_4A_FULL,DAA7,K11_4B_FULL,DAA7,K11_4C_FULL,K11_4D_FULL];
K11_5=[K11_5A_FULL,K11_5B_FULL,DAA7,K11_5C_FULL,K11_5D_FULL];
K11_6=[K11_6A_FULL,DAA7D,WS11614,DAA7E,K11_6B_FULL,K11_6C_FULL,K11_6D_FULL];
K11_7=[K11_7A_FULL,K11_7B_FULL,K11_7C_FULL,K11_7D_FULL];
K11_8=[K11_8A_FULL,K11_8B_FULL,K11_8C_FULL,K11_8D_FULL];
Q11501_FULL=[Q11501,BPM11501,Q11501];
Q11601_FULL=[Q11601,BPM11601,Q11601];
Q11701_FULL=[Q11701,BPM11701,Q11701];
Q11801_FULL=[Q11801,BPM11801,Q11801];
Q11901_FULL=[Q11901,BPM11901,Q11901];
LI11=[LI11STRT,K11_4,DAA7A,XC11502,DAA7B,WS11444,DAA7C,DAQ1,Q11501_FULL,DAQ2,K11_5,DAQ1,Q11601_FULL,DAQ2,K11_6,DAA7F,XC11702,DAA7G,YC11703,DAA7H,DAQ1,Q11701_FULL,DAQ2,K11_7,DAA7A,XC11802,DAA7B,WS11744,DAA7C,DAQ1,Q11801_FULL,DAQ2,K11_8,DAQ3,Q11901_FULL,DAQ4A,XC11900,DAQ4B,YC11900,DAQ4C,LI11END];
% ------------------------------------------------------------------------------
K12_1A_FULL=[K12_1A];
K12_1B_FULL=[K12_1B];
K12_1C_FULL=[K12_1C];
K12_1D_FULL=[K12_1D];
K12_2A_FULL=[K12_2A1,XC12202,K12_2A2,YC12203,K12_2A3];
K12_2B_FULL=[K12_2B];
K12_2C_FULL=[K12_2C];
K12_2D_FULL=[K12_2D];
K12_3A_FULL=[K12_3A1,XC12302,K12_3A2,YC12303,K12_3A3];
K12_3B_FULL=[K12_3B];
K12_3C_FULL=[K12_3C];
K12_3D_FULL=[K12_3D];
K12_4A_FULL=[K12_4A1,XC12402,K12_4A2,YC12403,K12_4A3];
K12_4B_FULL=[K12_4B];
K12_4C_FULL=[K12_4C];
K12_4D_FULL=[K12_4D];
K12_5A_FULL=[K12_5A1,XC12502,K12_5A2,YC12503,K12_5A3];
K12_5B_FULL=[K12_5B];
K12_5C_FULL=[K12_5C];
K12_5D_FULL=[K12_5D];
K12_6A_FULL=[K12_6A1,XC12602,K12_6A2,YC12603,K12_6A3];
K12_6B_FULL=[K12_6B];
K12_6C_FULL=[K12_6C];
K12_6D_FULL=[K12_6D];
K12_7A_FULL=[K12_7A1,XC12702,K12_7A2,YC12703,K12_7A3];
K12_7B_FULL=[K12_7B];
K12_7C_FULL=[K12_7C];
K12_7D_FULL=[K12_7D];
K12_8A_FULL=[K12_8A1,XC12802,K12_8A2,YC12803,K12_8A3];
K12_8B_FULL=[K12_8B];
K12_8C_FULL=[K12_8C];
K12_8D_FULL=[K12_8D1,XC12900,K12_8D2,YC12900,K12_8D3];
K12_1=[K12_1A_FULL,K12_1B_FULL,K12_1C_FULL,K12_1D_FULL];
K12_2=[K12_2A_FULL,DAA7D,WS12214,DAA7E,K12_2B_FULL,K12_2C_FULL,K12_2D_FULL];
K12_3=[K12_3A_FULL,K12_3B_FULL,K12_3C_FULL,K12_3D_FULL];
K12_4=[K12_4A_FULL,K12_4B_FULL,K12_4C_FULL,K12_4D_FULL];
K12_5=[K12_5A_FULL,K12_5B_FULL,K12_5C_FULL,K12_5D_FULL];
K12_6=[K12_6A_FULL,K12_6B_FULL,K12_6C_FULL,K12_6D_FULL];
K12_7=[K12_7A_FULL,K12_7B_FULL,K12_7C_FULL,K12_7D_FULL];
K12_8=[K12_8A_FULL,K12_8B_FULL,K12_8C_FULL,K12_8D_FULL];
Q12201_FULL=[Q12201,BPM12201,Q12201];
Q12301_FULL=[Q12301,BPM12301,Q12301];
Q12401_FULL=[Q12401,BPM12401,Q12401];
Q12501_FULL=[Q12501,BPM12501,Q12501];
Q12601_FULL=[Q12601,BPM12601,Q12601];
Q12701_FULL=[Q12701,BPM12701,Q12701];
Q12801_FULL=[Q12801,BPM12801,Q12801];
Q12901_FULL=[Q12901,BPM12901,Q12901];
LI12=[LI12BEG,K12_1,DAQ1,Q12201_FULL,DAQ2,K12_2,DAQ1,Q12301_FULL,DAQ2,K12_3,DAQ1,Q12401_FULL,DAQ2,K12_4,DAQ1,Q12501_FULL,DAQ2,K12_5,DAQ1,Q12601_FULL,DAQ2,K12_6,DAQ1,Q12701_FULL,DAQ2,K12_7,DAQ1,Q12801_FULL,DAQ2,K12_8,DAQ3,Q12901_FULL,DAQ4,LI12END];
% ------------------------------------------------------------------------------
K13_1A_FULL=[K13_1A];
K13_1B_FULL=[K13_1B];
K13_1C_FULL=[K13_1C];
K13_1D_FULL=[K13_1D];
K13_2A_FULL=[K13_2A1,XC13202,K13_2A2,YC13203,K13_2A3];
K13_2B_FULL=[K13_2B];
K13_2C_FULL=[K13_2C];
K13_2D_FULL=[K13_2D];
K13_3A_FULL=[K13_3A1,XC13302,K13_3A2,YC13303,K13_3A3];
K13_3B_FULL=[K13_3B];
K13_3C_FULL=[K13_3C];
K13_3D_FULL=[K13_3D];
K13_4A_FULL=[K13_4A1,XC13402,K13_4A2,YC13403,K13_4A3];
K13_4B_FULL=[K13_4B];
K13_4C_FULL=[K13_4C];
K13_4D_FULL=[K13_4D];
K13_5A_FULL=[K13_5A1,XC13502,K13_5A2,YC13503,K13_5A3];
K13_5B_FULL=[K13_5B];
K13_5C_FULL=[K13_5C];
K13_5D_FULL=[K13_5D];
K13_6A_FULL=[K13_6A1,XC13602,K13_6A2,YC13603,K13_6A3];
K13_6B_FULL=[K13_6B];
K13_6C_FULL=[K13_6C];
K13_6D_FULL=[K13_6D];
K13_7A_FULL=[K13_7A1,XC13702,K13_7A2,YC13703,K13_7A3];
K13_7B_FULL=[K13_7B];
K13_7C_FULL=[K13_7C];
K13_7D_FULL=[K13_7D];
K13_8A_FULL=[K13_8A1,XC13802,K13_8A2,YC13803,K13_8A3];
K13_8B_FULL=[K13_8B];
K13_8C_FULL=[K13_8C];
K13_8D_FULL=[K13_8D1,XC13900,K13_8D2,YC13900,K13_8D3];
K13_1=[K13_1A_FULL,K13_1B_FULL,K13_1C_FULL,K13_1D_FULL];
K13_2=[K13_2A_FULL,K13_2B_FULL,K13_2C_FULL,K13_2D_FULL];
K13_3=[K13_3A_FULL,K13_3B_FULL,K13_3C_FULL,K13_3D_FULL];
K13_4=[K13_4A_FULL,K13_4B_FULL,K13_4C_FULL,K13_4D_FULL];
K13_5=[K13_5A_FULL,K13_5B_FULL,K13_5C_FULL,K13_5D_FULL];
K13_6=[K13_6A_FULL,K13_6B_FULL,K13_6C_FULL,K13_6D_FULL];
K13_7=[K13_7A_FULL,K13_7B_FULL,K13_7C_FULL,K13_7D_FULL];
K13_8=[K13_8A_FULL,K13_8B_FULL,K13_8C_FULL,K13_8D_FULL];
Q13201_FULL=[Q13201,BPM13201,Q13201];
Q13301_FULL=[Q13301,BPM13301,Q13301];
Q13401_FULL=[Q13401,BPM13401,Q13401];
Q13501_FULL=[Q13501,BPM13501,Q13501];
Q13601_FULL=[Q13601,BPM13601,Q13601];
Q13701_FULL=[Q13701,BPM13701,Q13701];
Q13801_FULL=[Q13801,BPM13801,Q13801];
Q13901_FULL=[Q13901,BPM13901,Q13901];
LI13=[LI13BEG,K13_1,DAQ1,Q13201_FULL,DAQ2,K13_2,DAQ1,Q13301_FULL,DAQ2,K13_3,DAQ1,Q13401_FULL,DAQ2,K13_4,DAQ1,Q13501_FULL,DAQ2,K13_5,DAQ1,Q13601_FULL,DAQ2,K13_6,DAQ1,Q13701_FULL,DAQ2,K13_7,DAQ1,Q13801_FULL,DAQ2,K13_8,DAQ3,Q13901_FULL,DAQ4,LI13END];
% ------------------------------------------------------------------------------
K14_1A_FULL=[K14_1A];
K14_1B_FULL=[K14_1B];
K14_1C_FULL=[K14_1C];
K14_1D_FULL=[K14_1D];
K14_2A_FULL=[K14_2A1,XC14202,K14_2A2,YC14203,K14_2A3];
K14_2B_FULL=[K14_2B];
K14_2C_FULL=[K14_2C];
K14_2D_FULL=[K14_2D];
K14_3A_FULL=[K14_3A1,XC14302,K14_3A2,YC14303,K14_3A3];
K14_3B_FULL=[K14_3B];
K14_3C_FULL=[K14_3C];
K14_3D_FULL=[K14_3D];
K14_4A_FULL=[K14_4A1,XC14402,K14_4A2,YC14403,K14_4A3];
K14_4B_FULL=[K14_4B];
K14_4C_FULL=[K14_4C];
K14_4D_FULL=[K14_4D];
K14_5A_FULL=[K14_5A1,XC14502,K14_5A2,YC14503,K14_5A3];
K14_5B_FULL=[K14_5B];
K14_5C_FULL=[K14_5C];
K14_5D_FULL=[K14_5D];
K14_6A_FULL=[K14_6A1,XC14602,K14_6A2,YC14603,K14_6A3];
K14_6B_FULL=[K14_6B];
K14_6C_FULL=[K14_6C];
K14_6D_FULL=[K14_6D1,XC14702,K14_6D2,YC14703,K14_6D3];
K14_1=[K14_1A_FULL,K14_1B_FULL,K14_1C_FULL,K14_1D_FULL];
K14_2=[K14_2A_FULL,K14_2B_FULL,K14_2C_FULL,K14_2D_FULL];
K14_3=[K14_3A_FULL,K14_3B_FULL,K14_3C_FULL,K14_3D_FULL];
K14_4=[K14_4A_FULL,K14_4B_FULL,K14_4C_FULL,K14_4D_FULL];
K14_5=[K14_5A_FULL,K14_5B_FULL,K14_5C_FULL,K14_5D_FULL];
K14_6=[K14_6A_FULL,K14_6B_FULL,K14_6C_FULL,K14_6D_FULL];
Q14201_FULL=[Q14201,BPM14201,Q14201];
Q14301_FULL=[Q14301,BPM14301,Q14301];
Q14401_FULL=[Q14401,BPM14401,Q14401];
Q14501_FULL=[Q14501,BPM14501,Q14501];
Q14601_FULL=[Q14601,BPM14601,Q14601];
LI14=[LI14BEG,K14_1,DAQ1,Q14201_FULL,DAQ2,K14_2,DAQ1,Q14301_FULL,DAQ2,K14_3,DAQ1,Q14401_FULL,DAQ2,K14_4,DAQ1,Q14501_FULL,DAQ2,K14_5,DAQ1,Q14601_FULL,DAQ2,K14_6,LI14TERM];
% ------------------------------------------------------------------------------
L2F=[BEGL2F,LI11,LI12,LI13,LI14,ENDL2F];
% ==============================================================================

% *** OPTICS=FACET2-28JAN2026 ***
% ==============================================================================
% Modification History
% ------------------------------------------------------------------------------
% 14-JAN-2024, M. Woodley
%  * restore changed element names per K. Luchini
%    > CE141802 -> CE141815
%    > YC141820 -> YC141780
%    > VV14885  -> VV14887 
%    > IM14895  -> IM14890 
% ------------------------------------------------------------------------------
% 19-OCT-2023, M. Woodley
%  * restore element names that were changed per K. Luchini
%   * CE141815 -> CE141802
%   * YC141780 -> YC141820
%   * VV14887  -> VV14885
%   * IM14890  -> IM14895
% 01-SEP-2023, M. Woodley
%  * rename/redefine BC14E drift lengths (see G. White comments below)
%  * make BC14P equivalent to BC14E (same complement/Z-locations of devices)
% ------------------------------------------------------------------------------
% 26-APR-2023, G. White
%  * From measurements by Georg:
%   * YC14780 dz +0.0841m
%   * BPM14801 dz -0.0188m
%   * PR14803 dz -0.2212m
%   * CE14815 dz +0.863m
%   * CQ14866 dz +0.0651m
% 13-APR-2023, G. White
%  * Changes according to tunnel inspection by L. Alsburg & G. White:
%   * VV14885 -> VV14887 dS = +0.0759m
%   * BL14888 dS = +0.1524m
%   * QM14891 dS = 0.0472m
%   * IM14895 -> IM14890 to match controls IM14890 dS = -0.1095m
%   * Q14901 dS = -0.2048m
%   * Add in PR14892 (non-functional, mystery PROF device)
%   * Add VV14940 SLC fast valve
% ------------------------------------------------------------------------------
% 01-DEC-2022, G. White
%  * re-matched QM14891 - Q15601 to keep S15 strengths on quad strings more equal
% 28-FEB-2020, G. White - changes after visual inspection of beamline
%  * Moved YC14820 upstream of BCX14796 and changed unit number to 780
%  * Moved IM14890 toroid to between QM14891 & Q14901 and changed unit # to 895
%  * Moved VV14890 to just after BCX14883 and changed unit # to 885
%  * Changed CE14802 unit number to 815
% 23-AUG-2018, M. Woodley
%  * quadrupole K1 values from FACET2e_baseline.mat
%  * CQ's tweaked for Yuri-style match
%  * add VV14711 and VV14890
% ------------------------------------------------------------------------------
% 31-JAN-2017, M. Woodley
%  * from FACET2e.xsif
% ------------------------------------------------------------------------------
% ==============================================================================
% SBEN
% ------------------------------------------------------------------------------
% BC14 (electron side)
% - use series approximation for sinc(x)=sin(x)/x to allow AB14=0
% GB14  : chicane bend gap height (m)
% ZB14  : chicane bend "Z" length (m)
% AB14  : chicane bend angle (rad)
% LB14  : chicane bend path length (m)
% AB14s : "short" half chicane bend angle (rad)
% LB14s : "short" half chicane bend path length (m)
% AB14l : "long" half chicane bend angle (rad)
% LB14l : "long" half chicane bend path length (m)
GB14 =  0.03335;
ZB14 =  0.549;
AB14 =  0.04192;
AB14_2 =  AB14*AB14;
AB14_4 =  AB14_2*AB14_2;
AB14_6 =  AB14_4*AB14_2;
SINCAB14 =  1-AB14_2/6+AB14_4/120-AB14_6/5040 ;%~sinc(AB14)=sin(AB14)/AB14
LB14 =  ZB14/SINCAB14;
AB14S =  asin(sin(AB14)/2);
AB14S_2 =  AB14S*AB14S;
AB14S_4 =  AB14S_2*AB14S_2;
AB14S_6 =  AB14S_4*AB14S_2;
SINCAB14S =  1-AB14S_2/6+AB14S_4/120-AB14S_6/5040 ;%~sinc(AB14s)=sin(AB14s)/AB14s
LB14S =  (ZB14/2)/SINCAB14S;
AB14L =  AB14-AB14S;
LB14L =  LB14-LB14S;
EBC14E =  0.199858196713E-2;
BCX14720A={'be' 'BCX14720' LB14S [-AB14S GB14/2 0 0 0.633 0 0]}';
BCX14720B={'be' 'BCX14720' LB14L [-AB14L GB14/2 0 EBC14E 0 0.633 0]}';
BCX14796A={'be' 'BCX14796' LB14L [+AB14L GB14/2 +AB14 0 0.633 0 0]}';
BCX14796B={'be' 'BCX14796' LB14S [+AB14S GB14/2 0 0 0 0.633 0]}';
BCX14808A={'be' 'BCX14808' LB14S [+AB14S GB14/2 0 0 0.633 0 0]}';
BCX14808B={'be' 'BCX14808' LB14L [+AB14L GB14/2 0 +AB14 0 0.633 0]}';
BCX14883A={'be' 'BCX14883' LB14L [-AB14L GB14/2 EBC14E 0 0.633 0 0]}';
BCX14883B={'be' 'BCX14883' LB14S [-AB14S GB14/2 0 0 0 0.633 0]}';
% define unsplit SBENs for BMAD ... not used by MAD
BCX14720={'be' 'BCX14720' LB14 [-AB14 GB14/2 0 EBC14E 0.633 0.633 0]}';
BCX14796={'be' 'BCX14796' LB14 [+AB14 GB14/2 +AB14 0 0.633 0.633 0]}';
BCX14808={'be' 'BCX14808' LB14 [+AB14 GB14/2 0 +AB14 0.633 0.633 0]}';
BCX14883={'be' 'BCX14883' LB14 [-AB14 GB14/2 EBC14E 0 0.633 0.633 0]}';
% BC14 (positron side)
EBC14P =  -EBC14E;
BCX141720A={'be' 'BCX141720' LB14S [+AB14S GB14/2 0 0 0.633 0 0]}';
BCX141720B={'be' 'BCX141720' LB14L [+AB14L GB14/2 0 EBC14P 0 0.633 0]}';
BCX141796A={'be' 'BCX141796' LB14L [-AB14L GB14/2 -AB14 0 0.633 0 0]}';
BCX141796B={'be' 'BCX141796' LB14S [-AB14S GB14/2 0 0 0 0.633 0]}';
BCX141808A={'be' 'BCX141808' LB14S [-AB14S GB14/2 0 0 0.633 0 0]}';
BCX141808B={'be' 'BCX141808' LB14L [-AB14L GB14/2 0 -AB14 0 0.633 0]}';
BCX141883A={'be' 'BCX141883' LB14L [+AB14L GB14/2 EBC14P 0 0.633 0 0]}';
BCX141883B={'be' 'BCX141883' LB14S [+AB14S GB14/2 0 0 0 0.633 0]}';
% define unsplit SBENs for BMAD ... not used by MAD
BCX141720={'be' 'BCX141720' LB14 [+AB14 GB14/2 0 EBC14P 0.633 0.633 0]}';
BCX141796={'be' 'BCX141796' LB14 [-AB14 GB14/2 -AB14 0 0.633 0.633 0]}';
BCX141808={'be' 'BCX141808' LB14 [-AB14 GB14/2 0 -AB14 0.633 0.633 0]}';
BCX141883={'be' 'BCX141883' LB14 [+AB14 GB14/2 EBC14P 0 0.633 0.633 0]}';
% ==============================================================================
% QUAD
% ------------------------------------------------------------------------------
% common
KQ14701 =  -1.665165076061;
KQM14715 =   1.569247133195;
KQM14891 =   1.54788;
KQ14901 =  -1.954;
Q14701={'qu' 'Q14701' LQE/2 [QSIGN*(KQ14701) 0]}';
QM14715={'qu' 'QM14715' LQE/2 [QSIGN*(KQM14715) 0]}';
QM14891={'qu' 'QM14891' LQE/2 [QSIGN*(KQM14891) 0]}';
Q14901={'qu' 'Q14901' LQE/2 [QSIGN*(KQ14901) 0]}';
% electron
KCQ14 =   0;
KCQ14738 =  -0.228768471644E-2;
KCQ14866 =  -0.228768471644E-2;
CQ14738={'qu' 'CQ14738' LQC2/2 [KCQ14738 0]}';
CQ14866={'qu' 'CQ14866' LQC2/2 [KCQ14866 0]}';
% positron
KCQ141738 =  -0.228768471644E-2;
KCQ141866 =  -0.228768471644E-2;
CQ141738={'qu' 'CQ141738' LQC2/2 [KCQ141738 0]}';
CQ141866={'qu' 'CQ141866' LQC2/2 [KCQ141866 0]}';
% ==============================================================================
% DRIF
% ------------------------------------------------------------------------------
DM20={'dr' '' 0.0342 []}';
DM21={'dr' '' 1.9339006 []}';%0.6096+0.193+0.6340002+0.4973004
DM21A={'dr' '' 1.4366002 []}';
DM21B={'dr' '' DM21{3}-DM21A{3} []}';
DM22={'dr' '' 0.509664791 []}';%0.316596+0.193068791
ZDBQ2A =  1.924996725395;
ZD210A1 =  6.908547661855;
ZD210A2 =  0.852976041214;
ZD210B1 =  0.962912148546;
ZD210B2 =  6.863711554523;
ZDBQ2B =  1.859896725395;
DBQ2A={'dr' '' ZDBQ2A/cos(AB14) []}';
D210A1={'dr' '' ZD210A1/cos(AB14) []}';
D210A2={'dr' '' ZD210A2/cos(AB14) []}';
DDG21={'dr' '' 0.360994711244 []}';
DDG22={'dr' '' 0.144666 []}';
DDG23={'dr' '' 0.58687672 []}';
D210B1={'dr' '' ZD210B1/cos(AB14) []}';
D210B2={'dr' '' ZD210B2/cos(AB14) []}';
DBQ2B={'dr' '' ZDBQ2B/cos(AB14) []}';
DM23A={'dr' '' 0.589065836686 []}';%0.589065382243
DM23B={'dr' '' 0.190735 []}';%0.152434617756996
DM23C={'dr' '' 0.122265 []}';%0.160565382243004
DM24={'dr' '' 0.5552 []}';
DM24A={'dr' '' 0.1498 []}';
DM24B={'dr' '' 0.1534 []}';
DM25={'dr' '' 0.0446 []}';
DM25A={'dr' '' 0.254 []}';
DM25B={'dr' '' 0.2797 []}';
% ==============================================================================
% XCORs and YCORs
% ------------------------------------------------------------------------------
% electron
YC14780={'mo' 'YC14780' 0 []}';
% positron
YC141780={'mo' 'YC141780' 0 []}';
% ==============================================================================
% diagnostics
% ------------------------------------------------------------------------------
% BPMs (LCLS-II type designations)
% common
BPM14701={'mo' 'BPM14701' 0 []}';
BPM14715={'mo' 'BPM14715' 0 []}';
BPM14891={'mo' 'BPM14891' 0 []}';
BPM14901={'mo' 'BPM14901' 0 []}';
% electron
BPM14801={'mo' 'BPM14801' 0 []}';
% positron
BPM141801={'mo' 'BPM141801' 0 []}';
% misc
% common
BL14888={'mo' 'BL14888' 0 []}';
IM14890={'mo' 'IM14890' 0 []}';
% electron
PR14803={'mo' 'PR14803' 0 []}';
PR14892={'mo' 'PR14892' 0 []}';% NON-FUNCTIONAL
% positron
PR141803={'mo' 'PR141803' 0 []}';
% ==============================================================================
% collimators
% ------------------------------------------------------------------------------
% electron
CE14815={'dr' 'CE14815' 0 []}';%energy collimator
% positron
CE141815={'dr' 'CE141815' 0 []}';%energy collimator
% ==============================================================================
% vacuum valves
% ------------------------------------------------------------------------------
VV14711={'mo' 'VV14711' 0 []}';%BC14 upstream vacuum valve
VV14887={'mo' 'VV14887' 0 []}';%BC14 downstream vacuum valve
VV14940={'mo' 'VV14940' 0 []}';%SLC fast-valve
% ==============================================================================
% MARK
% ------------------------------------------------------------------------------
BEGBC14_1={'mo' 'BEGBC14_1' 0 []}';
ENDBC14_1={'mo' 'ENDBC14_1' 0 []}';
BEGBC14E={'mo' 'BEGBC14E' 0 []}';
ENDBC14E={'mo' 'ENDBC14E' 0 []}';
BEGBC14P={'mo' 'BEGBC14P' 0 []}';
ENDBC14P={'mo' 'ENDBC14P' 0 []}';
BEGBC14_2={'mo' 'BEGBC14_2' 0 []}';
CNT2B={'mo' 'CNT2B' 0 []}';
ENDBC14_2={'mo' 'ENDBC14_2' 0 []}';
% ==============================================================================
% BEAMLINEs
% ------------------------------------------------------------------------------
Q14701_FULL=[Q14701,BPM14701,Q14701];
QM14715_FULL=[QM14715,BPM14715,QM14715];
BC14_1=[BEGBC14_1,DM20,Q14701_FULL,DM21A,VV14711,DM21B,QM14715_FULL,DM22,ENDBC14_1];
BCX14720_FULL=[BCX14720A,BCX14720B];
CQ14738_FULL=[CQ14738,CQ14738];
BCX14796_FULL=[BCX14796A,BCX14796B];
BCX14808_FULL=[BCX14808A,BCX14808B];
CQ14866_FULL=[CQ14866,CQ14866];
BCX14883_FULL=[BCX14883A,BCX14883B];
BC14E=[BEGBC14E,BCX14720_FULL,DBQ2A,CQ14738_FULL,D210A1,YC14780,D210A2,BCX14796_FULL,DDG21,BPM14801,DDG22,PR14803,DDG23,BCX14808_FULL,D210B1,CE14815,D210B2,CQ14866_FULL,DBQ2B,BCX14883_FULL,ENDBC14E];
BCX141720_FULL=[BCX141720A,BCX141720B];
CQ141738_FULL=[CQ141738,CQ141738];
BCX141796_FULL=[BCX141796A,BCX141796B];
BCX141808_FULL=[BCX141808A,BCX141808B];
CQ141866_FULL=[CQ141866,CQ141866];
BCX141883_FULL=[BCX141883A,BCX141883B];
BC14P=[BEGBC14P,BCX141720_FULL,DBQ2A,CQ141738_FULL,D210A1,YC141780,D210A2,BCX141796_FULL,DDG21,BPM141801,DDG22,PR141803,DDG23,BCX141808_FULL,D210B1,CE141815,D210B2,CQ141866_FULL,DBQ2B,BCX141883_FULL,ENDBC14P];
QM14891_FULL=[QM14891,BPM14891,QM14891];
Q14901_FULL=[Q14901,BPM14901,Q14901];
BC14_2=[BEGBC14_2,CNT2B,DM23A,VV14887,DM23B,BL14888,DM23C,QM14891_FULL,DM24A,IM14890,DM24B,Q14901_FULL,DM25,PR14892,DM25A,VV14940,DM25B,ENDBC14_2];
% ==============================================================================

% *** OPTICS=FACET2-28JAN2026 ***
% ==============================================================================
% Modification History
% ------------------------------------------------------------------------------
% 15-FEB-2024, M. Woodley
%  * restore 19-8a 10' structure (unpowered), along with XC19802 and YC19803
% ------------------------------------------------------------------------------
% 26-APR-2023, G. White
%  * moved elements to reflect measurements from Georg:
%   * XC14900 dz +0.0993m
%   * YC14900 dz +0.1642m
% 02-MAY-2022, G. White
%  * moved XC19900 & YC19900 to reflect in-tunnel measurements by L. Alsburg
% 04-NOV-2021, M. Woodley
%  * remove unused K15_2d definition
% ------------------------------------------------------------------------------
% 16-APR-2021, G. White
%  * 19-8a removed, XC19900 & YC19900 moved to match as-installed
% 04-MAR-2019, M. Woodley
%  * restore 19-8a accelerator section
% 21-FEB-2019, M. Woodley
%  * rename some LI19 quads: Q1979x->Q19801, Q19801->Q19851, Q19901->Q19871
%  * move Q19801 back to it's nominal (SLC) location ... ~46 cm d/s
%  * restore original FACET correctors and toroid on 19-8 and 19-9
%  * rename LI18 collimators per L. Alsberg: CX1896->CX18960, CY1896->CY18960
% ------------------------------------------------------------------------------
% 23-AUG-2018, M. Woodley
%  * gradL3, phiL3, gradFB3, and phiFB3 values from FACET2e_baseline.mat;
%    set KlossL3=0
%  * quadrupole K1 values from FACET2e_baseline.mat
%  * remove TCAV isolation valves in LI15
% ------------------------------------------------------------------------------
% 22-DEC-2017, M. Woodley
%  * revert to original FACET (v35) layout and optics in LI18/LI19
%  * rematch into BC20 with Q18601, Q18701, Q18801, and Q18901
%  * reinstate TCAV3 (TCY15280) as an LCAV; add frequency (Sband) and
%    TYPE ("LOLA")
%  * shorten 18-1d, 18-2d, 18-3d, and 18-4d to 9.4' ... add wire scanners
% ------------------------------------------------------------------------------
% 31-JAN-2017, M. Woodley
%  * from FACET2e.xsif
% ------------------------------------------------------------------------------
% ==============================================================================
% accelerating structures
% ------------------------------------------------------------------------------
% the L3 S-band linac consists of: 146 x 10'   structure @ 25% power
%                                    1 x 10'   structure @ 50% power
%                                    3 x  9.4' structure @ 25% power
%                                    1 x  7'   structure @ 25% power
% ------------------------------------------------------------------------------
FL3 =   0.981867111676;
GRADL3 =  12.734495002141*FL3 ;%MeV/m
PHIL3 =   0.0/360            ;%rad/2pi
GRADFB3 =  18.48619126048*FL3  ;%MeV/m
PHIFB3 =  60.0/360            ;%rad/2pi
KLOSSL3 =   0                  ;%V/C/m
G15_1 =  GRADL3  ;
G15_2 =  GRADL3  ;
G15_3 =  GRADL3  ;
G15_4 =  GRADL3;
G15_5 =  GRADL3  ;
G15_6 =  GRADL3  ;
G15_7 =  GRADL3  ;
G15_8 =  GRADL3;
G16_1 =  GRADL3  ;
G16_2 =  GRADL3  ;
G16_3 =  GRADL3  ;
G16_4 =  GRADL3;
G16_5 =  GRADL3  ;
G16_6 =  GRADL3  ;
G16_7 =  GRADL3  ;
G16_8 =  GRADL3;
G17_1 =  GRADL3  ;
G17_2 =  GRADL3  ;
G17_3 =  GRADL3  ;
G17_4 =  GRADL3;
G17_5 =  GRADL3  ;
G17_6 =  GRADL3  ;
G17_7 =  GRADL3  ;
G17_8 =  GRADL3;
G18_1 =  GRADL3  ;
G18_2 =  GRADL3  ;
G18_3 =  GRADL3  ;
G18_4 =  GRADL3;
G18_5 =  GRADL3  ;
G18_6 =  GRADL3  ;
G18_7 =  GRADL3  ;
G18_8 =  GRADL3;
G19_1 =  GRADFB3 ;
G19_2 =  GRADFB3 ;
G19_3 =  GRADFB3 ;
G19_4 =  GRADFB3;
G19_5 =  GRADFB3 ;
G19_6 =  GRADFB3;
K15_1A1={'lc' 'K15_1A' 0.4293 [SBANDF P25*G15_1*0.4293 PHIL3*TWOPI]}';
K15_1A2={'lc' 'K15_1A' 0.3149 [SBANDF P25*G15_1*0.3149 PHIL3*TWOPI]}';
K15_1A3={'lc' 'K15_1A' 2.2999 [SBANDF P25*G15_1*2.2999 PHIL3*TWOPI]}';
K15_1B={'lc' 'K15_1B' DLWL10 [SBANDF P25*G15_1*DLWL10 PHIL3*TWOPI]}';
K15_1C={'lc' 'K15_1C' DLWL10 [SBANDF P25*G15_1*DLWL10 PHIL3*TWOPI]}';
K15_1D={'lc' 'K15_1D' DLWL10 [SBANDF P25*G15_1*DLWL10 PHIL3*TWOPI]}';
K15_2A1={'lc' 'K15_2A' 0.3256 [SBANDF P25*G15_2*0.3256 PHIL3*TWOPI]}';
K15_2A2={'lc' 'K15_2A' 0.25 [SBANDF P25*G15_2*0.25 PHIL3*TWOPI]}';
K15_2A3={'lc' 'K15_2A' 2.4685 [SBANDF P25*G15_2*2.4685 PHIL3*TWOPI]}';
K15_2B={'lc' 'K15_2B' DLWL10 [SBANDF P25*G15_2*DLWL10 PHIL3*TWOPI]}';
K15_2C={'lc' 'K15_2C' DLWL10 [SBANDF P50*G15_2*DLWL10 PHIL3*TWOPI]}';
K15_3A1={'lc' 'K15_3A' 0.3312 [SBANDF P25*G15_3*0.3312 PHIL3*TWOPI]}';
K15_3A2={'lc' 'K15_3A' 0.4044 [SBANDF P25*G15_3*0.4044 PHIL3*TWOPI]}';
K15_3A3={'lc' 'K15_3A' 2.3085 [SBANDF P25*G15_3*2.3085 PHIL3*TWOPI]}';
K15_3B={'lc' 'K15_3B' DLWL10 [SBANDF P25*G15_3*DLWL10 PHIL3*TWOPI]}';
K15_3C={'lc' 'K15_3C' DLWL10 [SBANDF P25*G15_3*DLWL10 PHIL3*TWOPI]}';
K15_3D={'lc' 'K15_3D' DLWL10 [SBANDF P25*G15_3*DLWL10 PHIL3*TWOPI]}';
K15_4A1={'lc' 'K15_4A' 0.3268 [SBANDF P25*G15_4*0.3268 PHIL3*TWOPI]}';
K15_4A2={'lc' 'K15_4A' 0.25 [SBANDF P25*G15_4*0.25 PHIL3*TWOPI]}';
K15_4A3={'lc' 'K15_4A' 2.4673 [SBANDF P25*G15_4*2.4673 PHIL3*TWOPI]}';
K15_4B={'lc' 'K15_4B' DLWL10 [SBANDF P25*G15_4*DLWL10 PHIL3*TWOPI]}';
K15_4C={'lc' 'K15_4C' DLWL10 [SBANDF P25*G15_4*DLWL10 PHIL3*TWOPI]}';
K15_4D={'lc' 'K15_4D' DLWL10 [SBANDF P25*G15_4*DLWL10 PHIL3*TWOPI]}';
K15_5A1={'lc' 'K15_5A' 0.3324 [SBANDF P25*G15_5*0.3324 PHIL3*TWOPI]}';
K15_5A2={'lc' 'K15_5A' 0.3937 [SBANDF P25*G15_5*0.3937 PHIL3*TWOPI]}';
K15_5A3={'lc' 'K15_5A' 2.318 [SBANDF P25*G15_5*2.318 PHIL3*TWOPI]}';
K15_5B={'lc' 'K15_5B' DLWL10 [SBANDF P25*G15_5*DLWL10 PHIL3*TWOPI]}';
K15_5C={'lc' 'K15_5C' DLWL10 [SBANDF P25*G15_5*DLWL10 PHIL3*TWOPI]}';
K15_5D={'lc' 'K15_5D' DLWL10 [SBANDF P25*G15_5*DLWL10 PHIL3*TWOPI]}';
K15_6A1={'lc' 'K15_6A' 0.328 [SBANDF P25*G15_6*0.328 PHIL3*TWOPI]}';
K15_6A2={'lc' 'K15_6A' 0.25 [SBANDF P25*G15_6*0.25 PHIL3*TWOPI]}';
K15_6A3={'lc' 'K15_6A' 2.4661 [SBANDF P25*G15_6*2.4661 PHIL3*TWOPI]}';
K15_6B={'lc' 'K15_6B' DLWL10 [SBANDF P25*G15_6*DLWL10 PHIL3*TWOPI]}';
K15_6C={'lc' 'K15_6C' DLWL10 [SBANDF P25*G15_6*DLWL10 PHIL3*TWOPI]}';
K15_6D={'lc' 'K15_6D' DLWL10 [SBANDF P25*G15_6*DLWL10 PHIL3*TWOPI]}';
K15_7A1={'lc' 'K15_7A' 0.3336 [SBANDF P25*G15_7*0.3336 PHIL3*TWOPI]}';
K15_7A2={'lc' 'K15_7A' 0.25 [SBANDF P25*G15_7*0.25 PHIL3*TWOPI]}';
K15_7A3={'lc' 'K15_7A' 2.4605 [SBANDF P25*G15_7*2.4605 PHIL3*TWOPI]}';
K15_7B={'lc' 'K15_7B' DLWL10 [SBANDF P25*G15_7*DLWL10 PHIL3*TWOPI]}';
K15_7C={'lc' 'K15_7C' DLWL10 [SBANDF P25*G15_7*DLWL10 PHIL3*TWOPI]}';
K15_7D={'lc' 'K15_7D' DLWL10 [SBANDF P25*G15_7*DLWL10 PHIL3*TWOPI]}';
K15_8A1={'lc' 'K15_8A' 0.3292 [SBANDF P25*G15_8*0.3292 PHIL3*TWOPI]}';
K15_8A2={'lc' 'K15_8A' 0.4413 [SBANDF P25*G15_8*0.4413 PHIL3*TWOPI]}';
K15_8A3={'lc' 'K15_8A' 2.2736 [SBANDF P25*G15_8*2.2736 PHIL3*TWOPI]}';
K15_8B={'lc' 'K15_8B' DLWL10 [SBANDF P25*G15_8*DLWL10 PHIL3*TWOPI]}';
K15_8C={'lc' 'K15_8C' DLWL10 [SBANDF P25*G15_8*DLWL10 PHIL3*TWOPI]}';
K15_8D1={'lc' 'K15_8D' 2.3869 [SBANDF P25*G15_8*2.3869 PHIL3*TWOPI]}';
K15_8D2={'lc' 'K15_8D' 0.25 [SBANDF P25*G15_8*0.25 PHIL3*TWOPI]}';
K15_8D3={'lc' 'K15_8D' 0.4072 [SBANDF P25*G15_8*0.4072 PHIL3*TWOPI]}';
K16_1A={'lc' 'K16_1A' DLWL10 [SBANDF P25*G16_1*DLWL10 PHIL3*TWOPI]}';
K16_1B={'lc' 'K16_1B' DLWL10 [SBANDF P25*G16_1*DLWL10 PHIL3*TWOPI]}';
K16_1C={'lc' 'K16_1C' DLWL10 [SBANDF P25*G16_1*DLWL10 PHIL3*TWOPI]}';
K16_1D={'lc' 'K16_1D' DLWL10 [SBANDF P25*G16_1*DLWL10 PHIL3*TWOPI]}';
K16_2A1={'lc' 'K16_2A' 0.3256 [SBANDF P25*G16_2*0.3256 PHIL3*TWOPI]}';
K16_2A2={'lc' 'K16_2A' 0.3655 [SBANDF P25*G16_2*0.3655 PHIL3*TWOPI]}';
K16_2A3={'lc' 'K16_2A' 2.353 [SBANDF P25*G16_2*2.353 PHIL3*TWOPI]}';
K16_2B={'lc' 'K16_2B' DLWL10 [SBANDF P25*G16_2*DLWL10 PHIL3*TWOPI]}';
K16_2C={'lc' 'K16_2C' DLWL10 [SBANDF P25*G16_2*DLWL10 PHIL3*TWOPI]}';
K16_2D={'lc' 'K16_2D' DLWL10 [SBANDF P25*G16_2*DLWL10 PHIL3*TWOPI]}';
K16_3A1={'lc' 'K16_3A' 0.3312 [SBANDF P25*G16_3*0.3312 PHIL3*TWOPI]}';
K16_3A2={'lc' 'K16_3A' 0.25 [SBANDF P25*G16_3*0.25 PHIL3*TWOPI]}';
K16_3A3={'lc' 'K16_3A' 2.4629 [SBANDF P25*G16_3*2.4629 PHIL3*TWOPI]}';
K16_3B={'lc' 'K16_3B' DLWL10 [SBANDF P25*G16_3*DLWL10 PHIL3*TWOPI]}';
K16_3C={'lc' 'K16_3C' DLWL10 [SBANDF P25*G16_3*DLWL10 PHIL3*TWOPI]}';
K16_3D={'lc' 'K16_3D' DLWL10 [SBANDF P25*G16_3*DLWL10 PHIL3*TWOPI]}';
K16_4A1={'lc' 'K16_4A' 0.3268 [SBANDF P25*G16_4*0.3268 PHIL3*TWOPI]}';
K16_4A2={'lc' 'K16_4A' 0.3675 [SBANDF P25*G16_4*0.3675 PHIL3*TWOPI]}';
K16_4A3={'lc' 'K16_4A' 2.3498 [SBANDF P25*G16_4*2.3498 PHIL3*TWOPI]}';
K16_4B={'lc' 'K16_4B' DLWL10 [SBANDF P25*G16_4*DLWL10 PHIL3*TWOPI]}';
K16_4C={'lc' 'K16_4C' DLWL10 [SBANDF P25*G16_4*DLWL10 PHIL3*TWOPI]}';
K16_4D={'lc' 'K16_4D' DLWL10 [SBANDF P25*G16_4*DLWL10 PHIL3*TWOPI]}';
K16_5A1={'lc' 'K16_5A' 0.3324 [SBANDF P25*G16_5*0.3324 PHIL3*TWOPI]}';
K16_5A2={'lc' 'K16_5A' 0.3651 [SBANDF P25*G16_5*0.3651 PHIL3*TWOPI]}';
K16_5A3={'lc' 'K16_5A' 2.3466 [SBANDF P25*G16_5*2.3466 PHIL3*TWOPI]}';
K16_5B={'lc' 'K16_5B' DLWL10 [SBANDF P25*G16_5*DLWL10 PHIL3*TWOPI]}';
K16_5C={'lc' 'K16_5C' DLWL10 [SBANDF P25*G16_5*DLWL10 PHIL3*TWOPI]}';
K16_5D={'lc' 'K16_5D' DLWL10 [SBANDF P25*G16_5*DLWL10 PHIL3*TWOPI]}';
K16_6A1={'lc' 'K16_6A' 0.328 [SBANDF P25*G16_6*0.328 PHIL3*TWOPI]}';
K16_6A2={'lc' 'K16_6A' 0.4139 [SBANDF P25*G16_6*0.4139 PHIL3*TWOPI]}';
K16_6A3={'lc' 'K16_6A' 2.3022 [SBANDF P25*G16_6*2.3022 PHIL3*TWOPI]}';
K16_6B={'lc' 'K16_6B' DLWL10 [SBANDF P25*G16_6*DLWL10 PHIL3*TWOPI]}';
K16_6C={'lc' 'K16_6C' DLWL10 [SBANDF P25*G16_6*DLWL10 PHIL3*TWOPI]}';
K16_6D={'lc' 'K16_6D' DLWL10 [SBANDF P25*G16_6*DLWL10 PHIL3*TWOPI]}';
K16_7A1={'lc' 'K16_7A' 0.3336 [SBANDF P25*G16_7*0.3336 PHIL3*TWOPI]}';
K16_7A2={'lc' 'K16_7A' 0.3893 [SBANDF P25*G16_7*0.3893 PHIL3*TWOPI]}';
K16_7A3={'lc' 'K16_7A' 2.3212 [SBANDF P25*G16_7*2.3212 PHIL3*TWOPI]}';
K16_7B={'lc' 'K16_7B' DLWL10 [SBANDF P25*G16_7*DLWL10 PHIL3*TWOPI]}';
K16_7C={'lc' 'K16_7C' DLWL10 [SBANDF P25*G16_7*DLWL10 PHIL3*TWOPI]}';
K16_7D={'lc' 'K16_7D' DLWL10 [SBANDF P25*G16_7*DLWL10 PHIL3*TWOPI]}';
K16_8A1={'lc' 'K16_8A' 0.3292 [SBANDF P25*G16_8*0.3292 PHIL3*TWOPI]}';
K16_8A2={'lc' 'K16_8A' 0.25 [SBANDF P25*G16_8*0.25 PHIL3*TWOPI]}';
K16_8A3={'lc' 'K16_8A' 2.4649 [SBANDF P25*G16_8*2.4649 PHIL3*TWOPI]}';
K16_8B={'lc' 'K16_8B' DLWL10 [SBANDF P25*G16_8*DLWL10 PHIL3*TWOPI]}';
K16_8C={'lc' 'K16_8C' DLWL10 [SBANDF P25*G16_8*DLWL10 PHIL3*TWOPI]}';
K16_8D1={'lc' 'K16_8D' 2.3869 [SBANDF P25*G16_8*2.3869 PHIL3*TWOPI]}';
K16_8D2={'lc' 'K16_8D' 0.25 [SBANDF P25*G16_8*0.25 PHIL3*TWOPI]}';
K16_8D3={'lc' 'K16_8D' 0.4072 [SBANDF P25*G16_8*0.4072 PHIL3*TWOPI]}';
K17_1A={'lc' 'K17_1A' DLWL10 [SBANDF P25*G17_1*DLWL10 PHIL3*TWOPI]}';
K17_1B={'lc' 'K17_1B' DLWL10 [SBANDF P25*G17_1*DLWL10 PHIL3*TWOPI]}';
K17_1C={'lc' 'K17_1C' DLWL10 [SBANDF P25*G17_1*DLWL10 PHIL3*TWOPI]}';
K17_1D={'lc' 'K17_1D' DLWL10 [SBANDF P25*G17_1*DLWL10 PHIL3*TWOPI]}';
K17_2A1={'lc' 'K17_2A' 0.3256 [SBANDF P25*G17_2*0.3256 PHIL3*TWOPI]}';
K17_2A2={'lc' 'K17_2A' 0.3814 [SBANDF P25*G17_2*0.3814 PHIL3*TWOPI]}';
K17_2A3={'lc' 'K17_2A' 2.3371 [SBANDF P25*G17_2*2.3371 PHIL3*TWOPI]}';
K17_2B={'lc' 'K17_2B' DLWL10 [SBANDF P25*G17_2*DLWL10 PHIL3*TWOPI]}';
K17_2C={'lc' 'K17_2C' DLWL10 [SBANDF P25*G17_2*DLWL10 PHIL3*TWOPI]}';
K17_2D={'lc' 'K17_2D' DLWL10 [SBANDF P25*G17_2*DLWL10 PHIL3*TWOPI]}';
K17_3A1={'lc' 'K17_3A' 0.3312 [SBANDF P25*G17_3*0.3312 PHIL3*TWOPI]}';
K17_3A2={'lc' 'K17_3A' 0.25 [SBANDF P25*G17_3*0.25 PHIL3*TWOPI]}';
K17_3A3={'lc' 'K17_3A' 2.4629 [SBANDF P25*G17_3*2.4629 PHIL3*TWOPI]}';
K17_3B={'lc' 'K17_3B' DLWL10 [SBANDF P25*G17_3*DLWL10 PHIL3*TWOPI]}';
K17_3C={'lc' 'K17_3C' DLWL10 [SBANDF P25*G17_3*DLWL10 PHIL3*TWOPI]}';
K17_3D={'lc' 'K17_3D' DLWL10 [SBANDF P25*G17_3*DLWL10 PHIL3*TWOPI]}';
K17_4A1={'lc' 'K17_4A' 0.3268 [SBANDF P25*G17_4*0.3268 PHIL3*TWOPI]}';
K17_4A2={'lc' 'K17_4A' 0.25 [SBANDF P25*G17_4*0.25 PHIL3*TWOPI]}';
K17_4A3={'lc' 'K17_4A' 2.4673 [SBANDF P25*G17_4*2.4673 PHIL3*TWOPI]}';
K17_4B={'lc' 'K17_4B' DLWL10 [SBANDF P25*G17_4*DLWL10 PHIL3*TWOPI]}';
K17_4C={'lc' 'K17_4C' DLWL10 [SBANDF P25*G17_4*DLWL10 PHIL3*TWOPI]}';
K17_4D={'lc' 'K17_4D' DLWL10 [SBANDF P25*G17_4*DLWL10 PHIL3*TWOPI]}';
K17_5A1={'lc' 'K17_5A' 0.3324 [SBANDF P25*G17_5*0.3324 PHIL3*TWOPI]}';
K17_5A2={'lc' 'K17_5A' 0.3841 [SBANDF P25*G17_5*0.3841 PHIL3*TWOPI]}';
K17_5A3={'lc' 'K17_5A' 2.3276 [SBANDF P25*G17_5*2.3276 PHIL3*TWOPI]}';
K17_5B={'lc' 'K17_5B' DLWL10 [SBANDF P25*G17_5*DLWL10 PHIL3*TWOPI]}';
K17_5C={'lc' 'K17_5C' DLWL10 [SBANDF P25*G17_5*DLWL10 PHIL3*TWOPI]}';
K17_5D={'lc' 'K17_5D' DLWL10 [SBANDF P25*G17_5*DLWL10 PHIL3*TWOPI]}';
K17_6A1={'lc' 'K17_6A' 0.328 [SBANDF P25*G17_6*0.328 PHIL3*TWOPI]}';
K17_6A2={'lc' 'K17_6A' 0.3949 [SBANDF P25*G17_6*0.3949 PHIL3*TWOPI]}';
K17_6A3={'lc' 'K17_6A' 2.3212 [SBANDF P25*G17_6*2.3212 PHIL3*TWOPI]}';
K17_6B={'lc' 'K17_6B' DLWL10 [SBANDF P25*G17_6*DLWL10 PHIL3*TWOPI]}';
K17_6C={'lc' 'K17_6C' DLWL10 [SBANDF P25*G17_6*DLWL10 PHIL3*TWOPI]}';
K17_6D={'lc' 'K17_6D' DLWL10 [SBANDF P25*G17_6*DLWL10 PHIL3*TWOPI]}';
K17_7A1={'lc' 'K17_7A' 0.3336 [SBANDF P25*G17_7*0.3336 PHIL3*TWOPI]}';
K17_7A2={'lc' 'K17_7A' 0.25 [SBANDF P25*G17_7*0.25 PHIL3*TWOPI]}';
K17_7A3={'lc' 'K17_7A' 2.4605 [SBANDF P25*G17_7*2.4605 PHIL3*TWOPI]}';
K17_7B={'lc' 'K17_7B' DLWL10 [SBANDF P25*G17_7*DLWL10 PHIL3*TWOPI]}';
K17_7C={'lc' 'K17_7C' DLWL10 [SBANDF P25*G17_7*DLWL10 PHIL3*TWOPI]}';
K17_7D={'lc' 'K17_7D' DLWL10 [SBANDF P25*G17_7*DLWL10 PHIL3*TWOPI]}';
K17_8A1={'lc' 'K17_8A' 0.3292 [SBANDF P25*G17_8*0.3292 PHIL3*TWOPI]}';
K17_8A2={'lc' 'K17_8A' 0.25 [SBANDF P25*G17_8*0.25 PHIL3*TWOPI]}';
K17_8A3={'lc' 'K17_8A' 2.4649 [SBANDF P25*G17_8*2.4649 PHIL3*TWOPI]}';
K17_8B={'lc' 'K17_8B' DLWL10 [SBANDF P25*G17_8*DLWL10 PHIL3*TWOPI]}';
K17_8C={'lc' 'K17_8C' DLWL10 [SBANDF P25*G17_8*DLWL10 PHIL3*TWOPI]}';
K17_8D1={'lc' 'K17_8D' 2.2856 [SBANDF P25*G17_8*2.2856 PHIL3*TWOPI]}';
K17_8D2={'lc' 'K17_8D' 0.3513 [SBANDF P25*G17_8*0.3513 PHIL3*TWOPI]}';
K17_8D3={'lc' 'K17_8D' 0.4072 [SBANDF P25*G17_8*0.4072 PHIL3*TWOPI]}';
K18_1A={'lc' 'K18_1A' DLWL10 [SBANDF P25*G18_1*DLWL10 PHIL3*TWOPI]}';
K18_1B={'lc' 'K18_1B' DLWL10 [SBANDF P25*G18_1*DLWL10 PHIL3*TWOPI]}';
K18_1C={'lc' 'K18_1C' DLWL10 [SBANDF P25*G18_1*DLWL10 PHIL3*TWOPI]}';
K18_1D={'lc' 'K18_1D' DLWL10 [SBANDF P25*G18_1*DLWL10 PHIL3*TWOPI]}';
K18_2A1={'lc' 'K18_2A' 0.3256 [SBANDF P25*G18_2*0.3256 PHIL3*TWOPI]}';
K18_2A2={'lc' 'K18_2A' 0.3878 [SBANDF P25*G18_2*0.3878 PHIL3*TWOPI]}';
K18_2A3={'lc' 'K18_2A' 2.3307 [SBANDF P25*G18_2*2.3307 PHIL3*TWOPI]}';
K18_2B={'lc' 'K18_2B' DLWL10 [SBANDF P25*G18_2*DLWL10 PHIL3*TWOPI]}';
K18_2C={'lc' 'K18_2C' DLWL10 [SBANDF P25*G18_2*DLWL10 PHIL3*TWOPI]}';
K18_2D={'lc' 'K18_2D' DLWL10 [SBANDF P25*G18_2*DLWL10 PHIL3*TWOPI]}';
K18_3A1={'lc' 'K18_3A' 0.3312 [SBANDF P25*G18_3*0.3312 PHIL3*TWOPI]}';
K18_3A2={'lc' 'K18_3A' 0.25 [SBANDF P25*G18_3*0.25 PHIL3*TWOPI]}';
K18_3A3={'lc' 'K18_3A' 2.4629 [SBANDF P25*G18_3*2.4629 PHIL3*TWOPI]}';
K18_3B={'lc' 'K18_3B' DLWL10 [SBANDF P25*G18_3*DLWL10 PHIL3*TWOPI]}';
K18_3C={'lc' 'K18_3C' DLWL10 [SBANDF P25*G18_3*DLWL10 PHIL3*TWOPI]}';
K18_3D={'lc' 'K18_3D' DLWL10 [SBANDF P25*G18_3*DLWL10 PHIL3*TWOPI]}';
K18_4A1={'lc' 'K18_4A' 0.3268 [SBANDF P25*G18_4*0.3268 PHIL3*TWOPI]}';
K18_4A2={'lc' 'K18_4A' 0.25 [SBANDF P25*G18_4*0.25 PHIL3*TWOPI]}';
K18_4A3={'lc' 'K18_4A' 2.4673 [SBANDF P25*G18_4*2.4673 PHIL3*TWOPI]}';
K18_4B={'lc' 'K18_4B' DLWL10 [SBANDF P25*G18_4*DLWL10 PHIL3*TWOPI]}';
K18_4C={'lc' 'K18_4C' DLWL10 [SBANDF P25*G18_4*DLWL10 PHIL3*TWOPI]}';
K18_4D={'lc' 'K18_4D' DLWL10 [SBANDF P25*G18_4*DLWL10 PHIL3*TWOPI]}';
K18_5A1={'lc' 'K18_5A' 0.3324 [SBANDF P25*G18_5*0.3324 PHIL3*TWOPI]}';
K18_5A2={'lc' 'K18_5A' 0.4095 [SBANDF P25*G18_5*0.4095 PHIL3*TWOPI]}';
K18_5A3={'lc' 'K18_5A' 2.3022 [SBANDF P25*G18_5*2.3022 PHIL3*TWOPI]}';
K18_5B={'lc' 'K18_5B' DLWL10 [SBANDF P25*G18_5*DLWL10 PHIL3*TWOPI]}';
K18_5C={'lc' 'K18_5C' DLWL10 [SBANDF P25*G18_5*DLWL10 PHIL3*TWOPI]}';
K18_5D={'lc' 'K18_5D' DLWL10 [SBANDF P25*G18_5*DLWL10 PHIL3*TWOPI]}';
K18_6A1={'lc' 'K18_6A' 0.328 [SBANDF P25*G18_6*0.328 PHIL3*TWOPI]}';
K18_6A2={'lc' 'K18_6A' 0.3885 [SBANDF P25*G18_6*0.3885 PHIL3*TWOPI]}';
K18_6A3={'lc' 'K18_6A' 2.3276 [SBANDF P25*G18_6*2.3276 PHIL3*TWOPI]}';
K18_6B={'lc' 'K18_6B' DLWL10 [SBANDF P25*G18_6*DLWL10 PHIL3*TWOPI]}';
K18_6C={'lc' 'K18_6C' DLWL10 [SBANDF P25*G18_6*DLWL10 PHIL3*TWOPI]}';
K18_6D={'lc' 'K18_6D' DLWL10 [SBANDF P25*G18_6*DLWL10 PHIL3*TWOPI]}';
K18_7A1={'lc' 'K18_7A' 0.3336 [SBANDF P25*G18_7*0.3336 PHIL3*TWOPI]}';
K18_7A2={'lc' 'K18_7A' 0.3798 [SBANDF P25*G18_7*0.3798 PHIL3*TWOPI]}';
K18_7A3={'lc' 'K18_7A' 2.3307 [SBANDF P25*G18_7*2.3307 PHIL3*TWOPI]}';
K18_7B={'lc' 'K18_7B' DLWL10 [SBANDF P25*G18_7*DLWL10 PHIL3*TWOPI]}';
K18_7C={'lc' 'K18_7C' DLWL10 [SBANDF P25*G18_7*DLWL10 PHIL3*TWOPI]}';
K18_7D={'lc' 'K18_7D' DLWL10 [SBANDF P25*G18_7*DLWL10 PHIL3*TWOPI]}';
K18_8A1={'lc' 'K18_8A' 0.3292 [SBANDF P25*G18_8*0.3292 PHIL3*TWOPI]}';
K18_8A2={'lc' 'K18_8A' 0.4127 [SBANDF P25*G18_8*0.4127 PHIL3*TWOPI]}';
K18_8A3={'lc' 'K18_8A' 2.3022 [SBANDF P25*G18_8*2.3022 PHIL3*TWOPI]}';
K18_8B={'lc' 'K18_8B' DLWL10 [SBANDF P25*G18_8*DLWL10 PHIL3*TWOPI]}';
K18_8C={'lc' 'K18_8C' DLWL10 [SBANDF P25*G18_8*DLWL10 PHIL3*TWOPI]}';
K18_8D1={'lc' 'K18_8D' 2.2729 [SBANDF P25*G18_8*2.2729 PHIL3*TWOPI]}';
K18_8D2={'lc' 'K18_8D' 0.364 [SBANDF P25*G18_8*0.364 PHIL3*TWOPI]}';
K18_8D3={'lc' 'K18_8D' 0.4072 [SBANDF P25*G18_8*0.4072 PHIL3*TWOPI]}';
K19_1A={'lc' 'K19_1A' DLWL10 [SBANDF P25*G19_1*DLWL10 -PHIFB3*TWOPI]}';
K19_1B={'lc' 'K19_1B' DLWL10 [SBANDF P25*G19_1*DLWL10 -PHIFB3*TWOPI]}';
K19_1C={'lc' 'K19_1C' DLWL10 [SBANDF P25*G19_1*DLWL10 -PHIFB3*TWOPI]}';
K19_1D={'lc' 'K19_1D' DLWL9 [SBANDF P25*G19_1*DLWL9 -PHIFB3*TWOPI]}';
K19_2A1={'lc' 'K19_2A' 0.3256 [SBANDF P25*G19_2*0.3256 -PHIFB3*TWOPI]}';
K19_2A2={'lc' 'K19_2A' 0.3592 [SBANDF P25*G19_2*0.3592 -PHIFB3*TWOPI]}';
K19_2A3={'lc' 'K19_2A' 2.3593 [SBANDF P25*G19_2*2.3593 -PHIFB3*TWOPI]}';
K19_2B={'lc' 'K19_2B' DLWL10 [SBANDF P25*G19_2*DLWL10 -PHIFB3*TWOPI]}';
K19_2C={'lc' 'K19_2C' DLWL10 [SBANDF P25*G19_2*DLWL10 -PHIFB3*TWOPI]}';
K19_2D={'lc' 'K19_2D' DLWL9 [SBANDF P25*G19_2*DLWL9 -PHIFB3*TWOPI]}';
K19_3A1={'lc' 'K19_3A' 0.3312 [SBANDF P25*G19_3*0.3312 -PHIFB3*TWOPI]}';
K19_3A2={'lc' 'K19_3A' 0.3695 [SBANDF P25*G19_3*0.3695 -PHIFB3*TWOPI]}';
K19_3A3={'lc' 'K19_3A' 2.3434 [SBANDF P25*G19_3*2.3434 -PHIFB3*TWOPI]}';
K19_3B={'lc' 'K19_3B' DLWL10 [SBANDF P25*G19_3*DLWL10 -PHIFB3*TWOPI]}';
K19_3C={'lc' 'K19_3C' DLWL10 [SBANDF P25*G19_3*DLWL10 -PHIFB3*TWOPI]}';
K19_3D={'lc' 'K19_3D' DLWL9 [SBANDF P25*G19_3*DLWL9 -PHIFB3*TWOPI]}';
K19_4A1={'lc' 'K19_4A' 0.3268 [SBANDF P25*G19_4*0.3268 +PHIFB3*TWOPI]}';
K19_4A2={'lc' 'K19_4A' 0.3897 [SBANDF P25*G19_4*0.3897 +PHIFB3*TWOPI]}';
K19_4A3={'lc' 'K19_4A' 2.3276 [SBANDF P25*G19_4*2.3276 +PHIFB3*TWOPI]}';
K19_4B={'lc' 'K19_4B' DLWL10 [SBANDF P25*G19_4*DLWL10 +PHIFB3*TWOPI]}';
K19_4C={'lc' 'K19_4C' DLWL10 [SBANDF P25*G19_4*DLWL10 +PHIFB3*TWOPI]}';
K19_4D1={'lc' 'K19_4D' 1.595 [SBANDF P25*G19_4*1.595 +PHIFB3*TWOPI]}';
K19_4D2={'lc' 'K19_4D' 0.3207 [SBANDF P25*G19_4*0.3207 +PHIFB3*TWOPI]}';
K19_4D3={'lc' 'K19_4D' 0.2537 [SBANDF P25*G19_4*0.2537 +PHIFB3*TWOPI]}';
K19_5A1={'lc' 'K19_5A' 0.3324 [SBANDF P25*G19_5*0.3324 +PHIFB3*TWOPI]}';
K19_5A2={'lc' 'K19_5A' 0.4032 [SBANDF P25*G19_5*0.4032 +PHIFB3*TWOPI]}';
K19_5A3={'lc' 'K19_5A' 2.3085 [SBANDF P25*G19_5*2.3085 +PHIFB3*TWOPI]}';
K19_5B={'lc' 'K19_5B' DLWL10 [SBANDF P25*G19_5*DLWL10 +PHIFB3*TWOPI]}';
K19_5C={'lc' 'K19_5C' DLWL10 [SBANDF P25*G19_5*DLWL10 +PHIFB3*TWOPI]}';
K19_5D={'lc' 'K19_5D' DLWL10 [SBANDF P25*G19_5*DLWL10 +PHIFB3*TWOPI]}';
K19_6A1={'lc' 'K19_6A' 0.328 [SBANDF P25*G19_6*0.328 +PHIFB3*TWOPI]}';
K19_6A2={'lc' 'K19_6A' 0.3854 [SBANDF P25*G19_6*0.3854 +PHIFB3*TWOPI]}';
K19_6A3={'lc' 'K19_6A' 2.3307 [SBANDF P25*G19_6*2.3307 +PHIFB3*TWOPI]}';
K19_6B={'lc' 'K19_6B' DLWL10 [SBANDF P25*G19_6*DLWL10 +PHIFB3*TWOPI]}';
K19_6C={'lc' 'K19_6C' DLWL10 [SBANDF P25*G19_6*DLWL10 +PHIFB3*TWOPI]}';
K19_6D1={'lc' 'K19_6D' 0.8425 [SBANDF P25*G19_6*0.8425 +PHIFB3*TWOPI]}';
K19_6D2={'lc' 'K19_6D' 1.4632 [SBANDF P25*G19_6*1.4632 +PHIFB3*TWOPI]}';
K19_6D3={'lc' 'K19_6D' 0.7384 [SBANDF P25*G19_6*0.7384 +PHIFB3*TWOPI]}';
% 19-8b,c,d removed ... 19-8a remains (unpowered)
K19_8A1={'lc' 'K19_8A' 0.3895 [SBANDF 0 0*TWOPI]}';
K19_8A2={'lc' 'K19_8A' 0.2985 [SBANDF 0 0*TWOPI]}';
K19_8A3={'lc' 'K19_8A' 2.3561 [SBANDF 0 0*TWOPI]}';
% define unsplit LCAVs for BMAD ... not used by MAD
K15_1A={'lc' 'K15_1A' DLWL10 [SBANDF P25*G15_1*DLWL10 PHIL3*TWOPI]}';
K15_2A={'lc' 'K15_2A' DLWL10 [SBANDF P25*G15_2*DLWL10 PHIL3*TWOPI]}';
K15_3A={'lc' 'K15_3A' DLWL10 [SBANDF P25*G15_3*DLWL10 PHIL3*TWOPI]}';
K15_4A={'lc' 'K15_4A' DLWL10 [SBANDF P25*G15_4*DLWL10 PHIL3*TWOPI]}';
K15_5A={'lc' 'K15_5A' DLWL10 [SBANDF P25*G15_5*DLWL10 PHIL3*TWOPI]}';
K15_6A={'lc' 'K15_6A' DLWL10 [SBANDF P25*G15_6*DLWL10 PHIL3*TWOPI]}';
K15_7A={'lc' 'K15_7A' DLWL10 [SBANDF P25*G15_7*DLWL10 PHIL3*TWOPI]}';
K15_8A={'lc' 'K15_8A' DLWL10 [SBANDF P25*G15_8*DLWL10 PHIL3*TWOPI]}';
K15_8D={'lc' 'K15_8D' DLWL10 [SBANDF P25*G15_8*DLWL10 PHIL3*TWOPI]}';
K16_2A={'lc' 'K16_2A' DLWL10 [SBANDF P25*G16_2*DLWL10 PHIL3*TWOPI]}';
K16_3A={'lc' 'K16_3A' DLWL10 [SBANDF P25*G16_3*DLWL10 PHIL3*TWOPI]}';
K16_4A={'lc' 'K16_4A' DLWL10 [SBANDF P25*G16_4*DLWL10 PHIL3*TWOPI]}';
K16_5A={'lc' 'K16_5A' DLWL10 [SBANDF P25*G16_5*DLWL10 PHIL3*TWOPI]}';
K16_6A={'lc' 'K16_6A' DLWL10 [SBANDF P25*G16_6*DLWL10 PHIL3*TWOPI]}';
K16_7A={'lc' 'K16_7A' DLWL10 [SBANDF P25*G16_7*DLWL10 PHIL3*TWOPI]}';
K16_8A={'lc' 'K16_8A' DLWL10 [SBANDF P25*G16_8*DLWL10 PHIL3*TWOPI]}';
K16_8D={'lc' 'K16_8D' DLWL10 [SBANDF P25*G16_8*DLWL10 PHIL3*TWOPI]}';
K17_2A={'lc' 'K17_2A' DLWL10 [SBANDF P25*G17_2*DLWL10 PHIL3*TWOPI]}';
K17_3A={'lc' 'K17_3A' DLWL10 [SBANDF P25*G17_3*DLWL10 PHIL3*TWOPI]}';
K17_4A={'lc' 'K17_4A' DLWL10 [SBANDF P25*G17_4*DLWL10 PHIL3*TWOPI]}';
K17_5A={'lc' 'K17_5A' DLWL10 [SBANDF P25*G17_5*DLWL10 PHIL3*TWOPI]}';
K17_6A={'lc' 'K17_6A' DLWL10 [SBANDF P25*G17_6*DLWL10 PHIL3*TWOPI]}';
K17_7A={'lc' 'K17_7A' DLWL10 [SBANDF P25*G17_7*DLWL10 PHIL3*TWOPI]}';
K17_8A={'lc' 'K17_8A' DLWL10 [SBANDF P25*G17_8*DLWL10 PHIL3*TWOPI]}';
K17_8D={'lc' 'K17_8D' DLWL10 [SBANDF P25*G17_8*DLWL10 PHIL3*TWOPI]}';
K18_2A={'lc' 'K18_2A' DLWL10 [SBANDF P25*G18_2*DLWL10 PHIL3*TWOPI]}';
K18_3A={'lc' 'K18_3A' DLWL10 [SBANDF P25*G18_3*DLWL10 PHIL3*TWOPI]}';
K18_4A={'lc' 'K18_4A' DLWL10 [SBANDF P25*G18_4*DLWL10 PHIL3*TWOPI]}';
K18_5A={'lc' 'K18_5A' DLWL10 [SBANDF P25*G18_5*DLWL10 PHIL3*TWOPI]}';
K18_6A={'lc' 'K18_6A' DLWL10 [SBANDF P25*G18_6*DLWL10 PHIL3*TWOPI]}';
K18_7A={'lc' 'K18_7A' DLWL10 [SBANDF P25*G18_7*DLWL10 PHIL3*TWOPI]}';
K18_8A={'lc' 'K18_8A' DLWL10 [SBANDF P25*G18_8*DLWL10 PHIL3*TWOPI]}';
K18_8D={'lc' 'K18_8D' DLWL10 [SBANDF P25*G18_8*DLWL10 PHIL3*TWOPI]}';
K19_2A={'lc' 'K19_2A' DLWL10 [SBANDF P25*G19_2*DLWL10 -PHIFB3*TWOPI]}';
K19_3A={'lc' 'K19_3A' DLWL10 [SBANDF P25*G19_3*DLWL10 -PHIFB3*TWOPI]}';
K19_4A={'lc' 'K19_4A' DLWL10 [SBANDF P25*G19_4*DLWL10 +PHIFB3*TWOPI]}';
K19_4D={'lc' 'K19_4D' DLWL7 [SBANDF P25*G19_4*DLWL7 +PHIFB3*TWOPI]}';
K19_5A={'lc' 'K19_5A' DLWL10 [SBANDF P25*G19_5*DLWL10 +PHIFB3*TWOPI]}';
K19_6A={'lc' 'K19_6A' DLWL10 [SBANDF P25*G19_6*DLWL10 +PHIFB3*TWOPI]}';
K19_6D={'lc' 'K19_6D' DLWL10 [SBANDF P25*G19_6*DLWL10 +PHIFB3*TWOPI]}';
% 19-8b,c,d removed ... 19-8a remains (unpowered)
K19_8A={'lc' 'K19_8A' DLWL10 [SBANDF 0 0*TWOPI]}';
% transverse deflecting cavity (TCAV)
TCY15280={'tc' 'TCY15280' 2.438/2 [SBANDF 0 0*TWOPI]}';
%TCY15280 : DRIF, L=2.438/2
% ==============================================================================
% QUADs
% ------------------------------------------------------------------------------
KQL3 =  0.81745071251 ;%65 degree cells
QFL3={'qu' 'QFL3' LQE/2 [+KQL3 0]}';
QDL3={'qu' 'QDL3' LQE/2 [-KQL3 0]}';
KQ15201 =   0.854425;
KQ15301 =  -0.833455;
KQ15401 =   0.635191;
KQ15501 =  -0.443457;
KQ15601 =   0.547422;
KQ15701 =  -0.623667834891;
KQ15801 =   0.791763004861;
KQ15901 =  -0.801101447545;
KQ16201 =   0.821463787359;
KQ16301 =  -0.850401234411;
KQ16401 =   0.788733005474;
KQ16501 =  -0.769929455281;
KQ16601 =   0.770004947138;
KQ16701 =  -0.795410224465;
KQ16801 =   0.850112860734;
KQ16901 =  -0.817211247334;
KQ17201 =   0.823948541708;
KQ17301 =  -0.848144311116;
KQ17401 =   0.786314789976;
KQ17501 =  -0.768808328938;
KQ17601 =   0.769910102282;
KQ17701 =  -0.796241135217;
KQ17801 =   0.851384850293;
KQ17901 =  -0.818176454295;
KQ18201 =   0.823948541708;
KQ18301 =  -0.848144311116;
KQ18401 =   0.786314789976;
KQ18501 =  -0.768808328938;
KQ18601 =   0.769910102282;
KQ18701 =  -0.796241135217;
KQ18801 =   0.851384850293;
KQ18901 =  -0.818176454295;
% LI19 bulk (201-801) I = 47.018 A, max boost < 15 A
KQ19201 =   0.824748618649;
KQ19301 =  -0.927074994776;
KQ19401 =   0.84178984371;
KQ19501 =  -1.024647998816 ;%-0.77122157 ... scavenger extraction
KQ19601 =   0.773964264856 ;% 0.77143041 ... scavenger extraction
KQ19701 =  -0.749016000184 ;%-0.7786437  ... scavenger extraction
KQ19801 =   1.002139999604;
KQ19851 =  -1.799181572264;
KQ19871 =   1.261375590325;
Q15201={'qu' 'Q15201' LQE/2 [QSIGN*(KQ15201) 0]}';
Q15301={'qu' 'Q15301' LQE/2 [QSIGN*(KQ15301) 0]}';
Q15401={'qu' 'Q15401' LQE/2 [QSIGN*(KQ15401) 0]}';
Q15501={'qu' 'Q15501' LQE/2 [QSIGN*(KQ15501) 0]}';
Q15601={'qu' 'Q15601' LQE/2 [QSIGN*(KQ15601) 0]}';
Q15701={'qu' 'Q15701' LQE/2 [QSIGN*(KQ15701) 0]}';
Q15801={'qu' 'Q15801' LQE/2 [QSIGN*(KQ15801) 0]}';
Q15901={'qu' 'Q15901' LQE/2 [QSIGN*(KQ15901) 0]}';
Q16201={'qu' 'Q16201' LQE/2 [QSIGN*(KQ16201) 0]}';
Q16301={'qu' 'Q16301' LQE/2 [QSIGN*(KQ16301) 0]}';
Q16401={'qu' 'Q16401' LQE/2 [QSIGN*(KQ16401) 0]}';
Q16501={'qu' 'Q16501' LQE/2 [QSIGN*(KQ16501) 0]}';
Q16601={'qu' 'Q16601' LQE/2 [QSIGN*(KQ16601) 0]}';
Q16701={'qu' 'Q16701' LQE/2 [QSIGN*(KQ16701) 0]}';
Q16801={'qu' 'Q16801' LQE/2 [QSIGN*(KQ16801) 0]}';
Q16901={'qu' 'Q16901' LQE/2 [QSIGN*(KQ16901) 0]}';
Q17201={'qu' 'Q17201' LQE/2 [QSIGN*(KQ17201) 0]}';
Q17301={'qu' 'Q17301' LQE/2 [QSIGN*(KQ17301) 0]}';
Q17401={'qu' 'Q17401' LQE/2 [QSIGN*(KQ17401) 0]}';
Q17501={'qu' 'Q17501' LQE/2 [QSIGN*(KQ17501) 0]}';
Q17601={'qu' 'Q17601' LQE/2 [QSIGN*(KQ17601) 0]}';
Q17701={'qu' 'Q17701' LQE/2 [QSIGN*(KQ17701) 0]}';
Q17801={'qu' 'Q17801' LQE/2 [QSIGN*(KQ17801) 0]}';
Q17901={'qu' 'Q17901' LQE/2 [QSIGN*(KQ17901) 0]}';
Q18201={'qu' 'Q18201' LQE/2 [QSIGN*(KQ18201) 0]}';
Q18301={'qu' 'Q18301' LQE/2 [QSIGN*(KQ18301) 0]}';
Q18401={'qu' 'Q18401' LQE/2 [QSIGN*(KQ18401) 0]}';
Q18501={'qu' 'Q18501' LQE/2 [QSIGN*(KQ18501) 0]}';
Q18601={'qu' 'Q18601' LQE/2 [QSIGN*(KQ18601) 0]}';
Q18701={'qu' 'Q18701' LQE/2 [QSIGN*(KQ18701) 0]}';
Q18801={'qu' 'Q18801' LQE/2 [QSIGN*(KQ18801) 0]}';
Q18901={'qu' 'Q18901' LQE/2 [QSIGN*(KQ18901) 0]}';
Q19201={'qu' 'Q19201' LQE/2 [QSIGN*(KQ19201) 0]}';
Q19301={'qu' 'Q19301' LQE/2 [QSIGN*(KQ19301) 0]}';
Q19401={'qu' 'Q19401' LQE/2 [QSIGN*(KQ19401) 0]}';
Q19501={'qu' 'Q19501' LQE/2 [QSIGN*(KQ19501) 0]}';
Q19601={'qu' 'Q19601' LQE/2 [QSIGN*(KQ19601) 0]}';
Q19701={'qu' 'Q19701' LQE/2 [QSIGN*(KQ19701) 0]}';
Q19801={'qu' 'Q19801' LQE/2 [QSIGN*(KQ19801) 0]}';
Q19851={'qu' 'Q19851' LQE/2 [QSIGN*(KQ19851) 0]}';
Q19871={'qu' 'Q19871' LQE/2 [QSIGN*(KQ19871) 0]}';
% ==============================================================================
% drifts
% ------------------------------------------------------------------------------
DL3FODO={'dr' '' 12.2376 []}';
D152DA={'dr' '' 0.1806+0.15 []}';
D152DB={'dr' '' 0.1255+0.15 []}';
DAQ4D={'dr' '' 0.5 []}';
DAQ4E={'dr' '' DAQ4{3}-DAQ4D{3} []}';
DAQ4K={'dr' '' 0.222952 []}';
DAQ4L={'dr' '' 0.443738 []}';
DAQ4M={'dr' '' 0.31891 []}';
DAQ4N={'dr' '' LDAQ4-(DAQ4K{3}+DAQ4L{3}+DAQ4M{3}) []}';
DAQ5={'dr' '' DLWL10-DLWL9 []}';
DAQ5A={'dr' '' 0.1 []}';
DAQ5B={'dr' '' DAQ5{3}-DAQ5A{3} []}';
DBKY170={'dr' '' 0.5/2 []}';
DAA7M={'dr' '' 0.585-(ZKS+LQE)/2 []}';%LDAA7-DAA7l[L]
DAA7L={'dr' '' LDAA7-ZKS-DAA7M{3} []}';%0.4378
D10={'dr' '' DLWL10 []}';
D10A={'dr' '' 0.9196 []}';
D10B={'dr' '' D10{3}-D10A{3} []}';
LMDW =  25.801374683414 ;%Q19701 exit to ENDL3F (from FACET v35)
D199A={'dr' '' 3.2062 []}';%6.2773
D199B={'dr' '' 2.0867 []}';
D199B1={'dr' '' 0.7277 []}';
D199B2={'dr' '' 0.343 []}';
D199B3={'dr' '' D199B{3}-D199B1{3}-D199B2{3} []}';%1.016
D199C={'dr' '' 4.570884 []}';
%D199d  : DRIF, L=Lmdw-(DAQ2[L]+4*DLWL10+DAQ1[L]+LQE+DAQ2[L]+&
%                 DLWL10+3.2062+LQE+2.0867+LQE+D199c[L])
D199D={'dr' '' LMDW-(DAQ2{3}+4*DLWL10+DAQ1{3}+LQE+DAQ2{3}+DLWL10+D199A{3}+LQE+D199B{3}+LQE+D199C{3}) []}';
% ==============================================================================
% XCORs and YCORs
% ------------------------------------------------------------------------------
XC14900={'mo' 'XC14900' 0 []}';
XC15202={'mo' 'XC15202' 0 []}';
XC15302={'mo' 'XC15302' 0 []}';
XC15402={'mo' 'XC15402' 0 []}';
XC15502={'mo' 'XC15502' 0 []}';
XC15602={'mo' 'XC15602' 0 []}';
XC15702={'mo' 'XC15702' 0 []}';
XC15802={'mo' 'XC15802' 0 []}';
XC15900={'mo' 'XC15900' 0 []}';
XC16202={'mo' 'XC16202' 0 []}';
XC16302={'mo' 'XC16302' 0 []}';
XC16402={'mo' 'XC16402' 0 []}';
XC16502={'mo' 'XC16502' 0 []}';
XC16602={'mo' 'XC16602' 0 []}';
XC16702={'mo' 'XC16702' 0 []}';
XC16802={'mo' 'XC16802' 0 []}';
XC16900={'mo' 'XC16900' 0 []}';
XC17202={'mo' 'XC17202' 0 []}';
XC17302={'mo' 'XC17302' 0 []}';
XC17402={'mo' 'XC17402' 0 []}';
XC17502={'mo' 'XC17502' 0 []}';
XC17602={'mo' 'XC17602' 0 []}';
XC17702={'mo' 'XC17702' 0 []}';
XC17802={'mo' 'XC17802' 0 []}';
XC17900={'mo' 'XC17900' 0 []}';
XC18202={'mo' 'XC18202' 0 []}';
XC18302={'mo' 'XC18302' 0 []}';
XC18402={'mo' 'XC18402' 0 []}';
XC18502={'mo' 'XC18502' 0 []}';
XC18602={'mo' 'XC18602' 0 []}';
XC18702={'mo' 'XC18702' 0 []}';
XC18802={'mo' 'XC18802' 0 []}';
XC18900={'mo' 'XC18900' 0 []}';
XC19202={'mo' 'XC19202' 0 []}';
XC19302={'mo' 'XC19302' 0 []}';
XC19402={'mo' 'XC19402' 0 []}';
XC19502={'mo' 'XC19502' 0 []}';
XC19602={'mo' 'XC19602' 0 []}';
XC19700={'mo' 'XC19700' 0 []}';
XC19802={'mo' 'XC19802' 0 []}';
XC19900={'mo' 'XC19900' 0 []}';
YC14900={'mo' 'YC14900' 0 []}';
YC15203={'mo' 'YC15203' 0 []}';
YC15303={'mo' 'YC15303' 0 []}';
YC15403={'mo' 'YC15403' 0 []}';
YC15503={'mo' 'YC15503' 0 []}';
YC15603={'mo' 'YC15603' 0 []}';
YC15703={'mo' 'YC15703' 0 []}';
YC15803={'mo' 'YC15803' 0 []}';
YC15900={'mo' 'YC15900' 0 []}';
YC16203={'mo' 'YC16203' 0 []}';
YC16303={'mo' 'YC16303' 0 []}';
YC16403={'mo' 'YC16403' 0 []}';
YC16503={'mo' 'YC16503' 0 []}';
YC16603={'mo' 'YC16603' 0 []}';
YC16703={'mo' 'YC16703' 0 []}';
YC16803={'mo' 'YC16803' 0 []}';
YC16900={'mo' 'YC16900' 0 []}';
YC17203={'mo' 'YC17203' 0 []}';
YC17303={'mo' 'YC17303' 0 []}';
YC17403={'mo' 'YC17403' 0 []}';
YC17503={'mo' 'YC17503' 0 []}';
YC17603={'mo' 'YC17603' 0 []}';
YC17703={'mo' 'YC17703' 0 []}';
YC17803={'mo' 'YC17803' 0 []}';
YC17900={'mo' 'YC17900' 0 []}';
YC18203={'mo' 'YC18203' 0 []}';
YC18303={'mo' 'YC18303' 0 []}';
YC18403={'mo' 'YC18403' 0 []}';
YC18503={'mo' 'YC18503' 0 []}';
YC18603={'mo' 'YC18603' 0 []}';
YC18703={'mo' 'YC18703' 0 []}';
YC18803={'mo' 'YC18803' 0 []}';
YC18900={'mo' 'YC18900' 0 []}';
YC19203={'mo' 'YC19203' 0 []}';
YC19303={'mo' 'YC19303' 0 []}';
YC19403={'mo' 'YC19403' 0 []}';
YC57145={'mo' 'YC57145' 0 []}';%scavenger e- DC extraction
YC57146={'mo' 'YC57146' 0 []}';%scavenger e- DC extraction
YC19503={'mo' 'YC19503' 0 []}';
YC19603={'mo' 'YC19603' 0 []}';
YC19700={'mo' 'YC19700' 0 []}';
YC19803={'mo' 'YC19803' 0 []}';
YC19900={'mo' 'YC19900' 0 []}';
% ==============================================================================
% diagnostics
% ------------------------------------------------------------------------------
% BPMs (LCLS-II type designations)
BPM15201={'mo' 'BPM15201' 0 []}';
BPM15301={'mo' 'BPM15301' 0 []}';
BPM15401={'mo' 'BPM15401' 0 []}';
BPM15501={'mo' 'BPM15501' 0 []}';
BPM15601={'mo' 'BPM15601' 0 []}';
BPM15701={'mo' 'BPM15701' 0 []}';
BPM15801={'mo' 'BPM15801' 0 []}';
BPM15901={'mo' 'BPM15901' 0 []}';
BPM16201={'mo' 'BPM16201' 0 []}';
BPM16301={'mo' 'BPM16301' 0 []}';
BPM16401={'mo' 'BPM16401' 0 []}';
BPM16501={'mo' 'BPM16501' 0 []}';
BPM16601={'mo' 'BPM16601' 0 []}';
BPM16701={'mo' 'BPM16701' 0 []}';
BPM16801={'mo' 'BPM16801' 0 []}';
BPM16901={'mo' 'BPM16901' 0 []}';
BPM17201={'mo' 'BPM17201' 0 []}';
BPM17301={'mo' 'BPM17301' 0 []}';
BPM17401={'mo' 'BPM17401' 0 []}';
BPM17501={'mo' 'BPM17501' 0 []}';
BPM17601={'mo' 'BPM17601' 0 []}';
BPM17701={'mo' 'BPM17701' 0 []}';
BPM17801={'mo' 'BPM17801' 0 []}';
BPM17901={'mo' 'BPM17901' 0 []}';
BPM18201={'mo' 'BPM18201' 0 []}';
BPM18301={'mo' 'BPM18301' 0 []}';
BPM18401={'mo' 'BPM18401' 0 []}';
BPM18501={'mo' 'BPM18501' 0 []}';
BPM18601={'mo' 'BPM18601' 0 []}';
BPM18701={'mo' 'BPM18701' 0 []}';
BPM18801={'mo' 'BPM18801' 0 []}';
BPM18901={'mo' 'BPM18901' 0 []}';
BPM19201={'mo' 'BPM19201' 0 []}';
BPM19301={'mo' 'BPM19301' 0 []}';
BPM19401={'mo' 'BPM19401' 0 []}';
BPM19501={'mo' 'BPM19501' 0 []}';
BPM19601={'mo' 'BPM19601' 0 []}';
BPM19701={'mo' 'BPM19701' 0 []}';
BPM19801={'mo' 'BPM19801' 0 []}';
BPM19851={'mo' 'BPM19851' 0 []}';
BPM19871={'mo' 'BPM19871' 0 []}';
% misc
PR15944={'mo' 'PR15944' 0 []}';%post-BC14 longitudinal diagnostics
%NOTE: dMUY from TCY15280 is ~135 deg (should be 90)
BL18900={'mo' 'BL18900' 0 []}';%bunch length monitor / OTR at 18-9
WS18944={'mo' 'WS18944' 0 []}';
WS19144={'mo' 'WS19144' 0 []}';
WS19244={'mo' 'WS19244' 0 []}';
WS19344={'mo' 'WS19344' 0 []}';
IM1988={'mo' 'IM1988' 0 []}';
% ==============================================================================
% collimators
% ------------------------------------------------------------------------------
CX18960={'dr' 'CX18960' 0 []}';
CY18960={'dr' 'CY18960' 0 []}';
% ==============================================================================
% MARK
% ------------------------------------------------------------------------------
BEGL3F_1={'mo' 'BEGL3F_1' 0 []}';
LI15BEG={'mo' 'LI15BEG' 0 []}';
LI15END={'mo' 'LI15END' 0 []}';
LI16BEG={'mo' 'LI16BEG' 0 []}';
LI16END={'mo' 'LI16END' 0 []}';
LI17BEG={'mo' 'LI17BEG' 0 []}';
LI17END={'mo' 'LI17END' 0 []}';
LI18BEG={'mo' 'LI18BEG' 0 []}';
LI18END={'mo' 'LI18END' 0 []}';
LI19BEG={'mo' 'LI19BEG' 0 []}';
MSCAVEXT={'mo' 'MSCAVEXT' 0 []}';%match to SCAVX19.TRANS here
ENDL3F_1={'mo' 'ENDL3F_1' 0 []}';
BEGL3F_2={'mo' 'BEGL3F_2' 0 []}';
MSCAV={'mo' 'MSCAV' 0 []}';%match to SCAV20A.TRANS here
HLAM={'mo' 'HLAM' 0 []}';%entrance face of HLAM172
LI19TERM={'mo' 'LI19TERM' 0 []}';
ENDL3F_2={'mo' 'ENDL3F_2' 0 []}';
% ==============================================================================
% BEAMLINEs
% ------------------------------------------------------------------------------
L3C=[QFL3,DL3FODO,QDL3,QDL3,DL3FODO,QFL3];
% ------------------------------------------------------------------------------
K15_1A_FULL=[K15_1A1,XC14900,K15_1A2,YC14900,K15_1A3];
K15_1B_FULL=[K15_1B];
K15_1C_FULL=[K15_1C];
K15_1D_FULL=[K15_1D];
K15_2A_FULL=[K15_2A1,XC15202,K15_2A2,YC15203,K15_2A3];
K15_2B_FULL=[K15_2B];
K15_2C_FULL=[K15_2C];
K15_3A_FULL=[K15_3A1,XC15302,K15_3A2,YC15303,K15_3A3];
K15_3B_FULL=[K15_3B];
K15_3C_FULL=[K15_3C];
K15_3D_FULL=[K15_3D];
K15_4A_FULL=[K15_4A1,XC15402,K15_4A2,YC15403,K15_4A3];
K15_4B_FULL=[K15_4B];
K15_4C_FULL=[K15_4C];
K15_4D_FULL=[K15_4D];
K15_5A_FULL=[K15_5A1,XC15502,K15_5A2,YC15503,K15_5A3];
K15_5B_FULL=[K15_5B];
K15_5C_FULL=[K15_5C];
K15_5D_FULL=[K15_5D];
K15_6A_FULL=[K15_6A1,XC15602,K15_6A2,YC15603,K15_6A3];
K15_6B_FULL=[K15_6B];
K15_6C_FULL=[K15_6C];
K15_6D_FULL=[K15_6D];
K15_7A_FULL=[K15_7A1,XC15702,K15_7A2,YC15703,K15_7A3];
K15_7B_FULL=[K15_7B];
K15_7C_FULL=[K15_7C];
K15_7D_FULL=[K15_7D];
K15_8A_FULL=[K15_8A1,XC15802,K15_8A2,YC15803,K15_8A3];
K15_8B_FULL=[K15_8B];
K15_8C_FULL=[K15_8C];
K15_8D_FULL=[K15_8D1,XC15900,K15_8D2,YC15900,K15_8D3];
TCY15280_FULL=[TCY15280,TCY15280];
K15_1=[K15_1A_FULL,K15_1B_FULL,K15_1C_FULL,K15_1D_FULL];
K15_2=[K15_2A_FULL,K15_2B_FULL,K15_2C_FULL            ];
K15_3=[K15_3A_FULL,K15_3B_FULL,K15_3C_FULL,K15_3D_FULL];
K15_4=[K15_4A_FULL,K15_4B_FULL,K15_4C_FULL,K15_4D_FULL];
K15_5=[K15_5A_FULL,K15_5B_FULL,K15_5C_FULL,K15_5D_FULL];
K15_6=[K15_6A_FULL,K15_6B_FULL,K15_6C_FULL,K15_6D_FULL];
K15_7=[K15_7A_FULL,K15_7B_FULL,K15_7C_FULL,K15_7D_FULL];
K15_8=[K15_8A_FULL,K15_8B_FULL,K15_8C_FULL,K15_8D_FULL];
Q15201_FULL=[Q15201,BPM15201,Q15201];
Q15301_FULL=[Q15301,BPM15301,Q15301];
Q15401_FULL=[Q15401,BPM15401,Q15401];
Q15501_FULL=[Q15501,BPM15501,Q15501];
Q15601_FULL=[Q15601,BPM15601,Q15601];
Q15701_FULL=[Q15701,BPM15701,Q15701];
Q15801_FULL=[Q15801,BPM15801,Q15801];
Q15901_FULL=[Q15901,BPM15901,Q15901];
LI15=[LI15BEG,K15_1,DAQ1,Q15201_FULL,DAQ2,K15_2,D152DA,TCY15280_FULL,D152DB,DAQ1,Q15301_FULL,DAQ2,K15_3,DAQ1,Q15401_FULL,DAQ2,K15_4,DAQ1,Q15501_FULL,DAQ2,K15_5,DAQ1,Q15601_FULL,DAQ2,K15_6,DAQ1,Q15701_FULL,DAQ2,K15_7,DAQ1,Q15801_FULL,DAQ2,K15_8,DAQ3,Q15901_FULL,DAQ4D,PR15944,DAQ4E,LI15END];
% ------------------------------------------------------------------------------
K16_1A_FULL=[K16_1A];
K16_1B_FULL=[K16_1B];
K16_1C_FULL=[K16_1C];
K16_1D_FULL=[K16_1D];
K16_2A_FULL=[K16_2A1,XC16202,K16_2A2,YC16203,K16_2A3];
K16_2B_FULL=[K16_2B];
K16_2C_FULL=[K16_2C];
K16_2D_FULL=[K16_2D];
K16_3A_FULL=[K16_3A1,XC16302,K16_3A2,YC16303,K16_3A3];
K16_3B_FULL=[K16_3B];
K16_3C_FULL=[K16_3C];
K16_3D_FULL=[K16_3D];
K16_4A_FULL=[K16_4A1,XC16402,K16_4A2,YC16403,K16_4A3];
K16_4B_FULL=[K16_4B];
K16_4C_FULL=[K16_4C];
K16_4D_FULL=[K16_4D];
K16_5A_FULL=[K16_5A1,XC16502,K16_5A2,YC16503,K16_5A3];
K16_5B_FULL=[K16_5B];
K16_5C_FULL=[K16_5C];
K16_5D_FULL=[K16_5D];
K16_6A_FULL=[K16_6A1,XC16602,K16_6A2,YC16603,K16_6A3];
K16_6B_FULL=[K16_6B];
K16_6C_FULL=[K16_6C];
K16_6D_FULL=[K16_6D];
K16_7A_FULL=[K16_7A1,XC16702,K16_7A2,YC16703,K16_7A3];
K16_7B_FULL=[K16_7B];
K16_7C_FULL=[K16_7C];
K16_7D_FULL=[K16_7D];
K16_8A_FULL=[K16_8A1,XC16802,K16_8A2,YC16803,K16_8A3];
K16_8B_FULL=[K16_8B];
K16_8C_FULL=[K16_8C];
K16_8D_FULL=[K16_8D1,XC16900,K16_8D2,YC16900,K16_8D3];
K16_1=[K16_1A_FULL,K16_1B_FULL,K16_1C_FULL,K16_1D_FULL];
K16_2=[K16_2A_FULL,K16_2B_FULL,K16_2C_FULL,K16_2D_FULL];
K16_3=[K16_3A_FULL,K16_3B_FULL,K16_3C_FULL,K16_3D_FULL];
K16_4=[K16_4A_FULL,K16_4B_FULL,K16_4C_FULL,K16_4D_FULL];
K16_5=[K16_5A_FULL,K16_5B_FULL,K16_5C_FULL,K16_5D_FULL];
K16_6=[K16_6A_FULL,K16_6B_FULL,K16_6C_FULL,K16_6D_FULL];
K16_7=[K16_7A_FULL,K16_7B_FULL,K16_7C_FULL,K16_7D_FULL];
K16_8=[K16_8A_FULL,K16_8B_FULL,K16_8C_FULL,K16_8D_FULL];
Q16201_FULL=[Q16201,BPM16201,Q16201];
Q16301_FULL=[Q16301,BPM16301,Q16301];
Q16401_FULL=[Q16401,BPM16401,Q16401];
Q16501_FULL=[Q16501,BPM16501,Q16501];
Q16601_FULL=[Q16601,BPM16601,Q16601];
Q16701_FULL=[Q16701,BPM16701,Q16701];
Q16801_FULL=[Q16801,BPM16801,Q16801];
Q16901_FULL=[Q16901,BPM16901,Q16901];
LI16=[LI16BEG,K16_1,DAQ1,Q16201_FULL,DAQ2,K16_2,DAQ1,Q16301_FULL,DAQ2,K16_3,DAQ1,Q16401_FULL,DAQ2,K16_4,DAQ1,Q16501_FULL,DAQ2,K16_5,DAQ1,Q16601_FULL,DAQ2,K16_6,DAQ1,Q16701_FULL,DAQ2,K16_7,DAQ1,Q16801_FULL,DAQ2,K16_8,DAQ3,Q16901_FULL,DAQ4,LI16END];
% ------------------------------------------------------------------------------
K17_1A_FULL=[K17_1A];
K17_1B_FULL=[K17_1B];
K17_1C_FULL=[K17_1C];
K17_1D_FULL=[K17_1D];
K17_2A_FULL=[K17_2A1,XC17202,K17_2A2,YC17203,K17_2A3];
K17_2B_FULL=[K17_2B];
K17_2C_FULL=[K17_2C];
K17_2D_FULL=[K17_2D];
K17_3A_FULL=[K17_3A1,XC17302,K17_3A2,YC17303,K17_3A3];
K17_3B_FULL=[K17_3B];
K17_3C_FULL=[K17_3C];
K17_3D_FULL=[K17_3D];
K17_4A_FULL=[K17_4A1,XC17402,K17_4A2,YC17403,K17_4A3];
K17_4B_FULL=[K17_4B];
K17_4C_FULL=[K17_4C];
K17_4D_FULL=[K17_4D];
K17_5A_FULL=[K17_5A1,XC17502,K17_5A2,YC17503,K17_5A3];
K17_5B_FULL=[K17_5B];
K17_5C_FULL=[K17_5C];
K17_5D_FULL=[K17_5D];
K17_6A_FULL=[K17_6A1,XC17602,K17_6A2,YC17603,K17_6A3];
K17_6B_FULL=[K17_6B];
K17_6C_FULL=[K17_6C];
K17_6D_FULL=[K17_6D];
K17_7A_FULL=[K17_7A1,XC17702,K17_7A2,YC17703,K17_7A3];
K17_7B_FULL=[K17_7B];
K17_7C_FULL=[K17_7C];
K17_7D_FULL=[K17_7D];
K17_8A_FULL=[K17_8A1,XC17802,K17_8A2,YC17803,K17_8A3];
K17_8B_FULL=[K17_8B];
K17_8C_FULL=[K17_8C];
K17_8D_FULL=[K17_8D1,XC17900,K17_8D2,YC17900,K17_8D3];
K17_1=[K17_1A_FULL,K17_1B_FULL,K17_1C_FULL,K17_1D_FULL];
K17_2=[K17_2A_FULL,K17_2B_FULL,K17_2C_FULL,K17_2D_FULL];
K17_3=[K17_3A_FULL,K17_3B_FULL,K17_3C_FULL,K17_3D_FULL];
K17_4=[K17_4A_FULL,K17_4B_FULL,K17_4C_FULL,K17_4D_FULL];
K17_5=[K17_5A_FULL,K17_5B_FULL,K17_5C_FULL,K17_5D_FULL];
K17_6=[K17_6A_FULL,K17_6B_FULL,K17_6C_FULL,K17_6D_FULL];
K17_7=[K17_7A_FULL,K17_7B_FULL,K17_7C_FULL,K17_7D_FULL];
K17_8=[K17_8A_FULL,K17_8B_FULL,K17_8C_FULL,K17_8D_FULL];
Q17201_FULL=[Q17201,BPM17201,Q17201];
Q17301_FULL=[Q17301,BPM17301,Q17301];
Q17401_FULL=[Q17401,BPM17401,Q17401];
Q17501_FULL=[Q17501,BPM17501,Q17501];
Q17601_FULL=[Q17601,BPM17601,Q17601];
Q17701_FULL=[Q17701,BPM17701,Q17701];
Q17801_FULL=[Q17801,BPM17801,Q17801];
Q17901_FULL=[Q17901,BPM17901,Q17901];
LI17=[LI17BEG,K17_1,DAQ1,Q17201_FULL,DAQ2,K17_2,DAQ1,Q17301_FULL,DAQ2,K17_3,DAQ1,Q17401_FULL,DAQ2,K17_4,DAQ1,Q17501_FULL,DAQ2,K17_5,DAQ1,Q17601_FULL,DAQ2,K17_6,DAQ1,Q17701_FULL,DAQ2,K17_7,DAQ1,Q17801_FULL,DAQ2,K17_8,DAQ3,Q17901_FULL,DAQ4,LI17END];
% ------------------------------------------------------------------------------
K18_1A_FULL=[K18_1A];
K18_1B_FULL=[K18_1B];
K18_1C_FULL=[K18_1C];
K18_1D_FULL=[K18_1D];
K18_2A_FULL=[K18_2A1,XC18202,K18_2A2,YC18203,K18_2A3];
K18_2B_FULL=[K18_2B];
K18_2C_FULL=[K18_2C];
K18_2D_FULL=[K18_2D];
K18_3A_FULL=[K18_3A1,XC18302,K18_3A2,YC18303,K18_3A3];
K18_3B_FULL=[K18_3B];
K18_3C_FULL=[K18_3C];
K18_3D_FULL=[K18_3D];
K18_4A_FULL=[K18_4A1,XC18402,K18_4A2,YC18403,K18_4A3];
K18_4B_FULL=[K18_4B];
K18_4C_FULL=[K18_4C];
K18_4D_FULL=[K18_4D];
K18_5A_FULL=[K18_5A1,XC18502,K18_5A2,YC18503,K18_5A3];
K18_5B_FULL=[K18_5B];
K18_5C_FULL=[K18_5C];
K18_5D_FULL=[K18_5D];
K18_6A_FULL=[K18_6A1,XC18602,K18_6A2,YC18603,K18_6A3];
K18_6B_FULL=[K18_6B];
K18_6C_FULL=[K18_6C];
K18_6D_FULL=[K18_6D];
K18_7A_FULL=[K18_7A1,XC18702,K18_7A2,YC18703,K18_7A3];
K18_7B_FULL=[K18_7B];
K18_7C_FULL=[K18_7C];
K18_7D_FULL=[K18_7D];
K18_8A_FULL=[K18_8A1,XC18802,K18_8A2,YC18803,K18_8A3];
K18_8B_FULL=[K18_8B];
K18_8C_FULL=[K18_8C];
K18_8D_FULL=[K18_8D1,XC18900,K18_8D2,YC18900,K18_8D3];
K18_1=[K18_1A_FULL,K18_1B_FULL,K18_1C_FULL,K18_1D_FULL];
K18_2=[K18_2A_FULL,K18_2B_FULL,K18_2C_FULL,K18_2D_FULL];
K18_3=[K18_3A_FULL,K18_3B_FULL,K18_3C_FULL,K18_3D_FULL];
K18_4=[K18_4A_FULL,K18_4B_FULL,K18_4C_FULL,K18_4D_FULL];
K18_5=[K18_5A_FULL,K18_5B_FULL,K18_5C_FULL,K18_5D_FULL];
K18_6=[K18_6A_FULL,K18_6B_FULL,K18_6C_FULL,K18_6D_FULL];
K18_7=[K18_7A_FULL,K18_7B_FULL,K18_7C_FULL,K18_7D_FULL];
K18_8=[K18_8A_FULL,K18_8B_FULL,K18_8C_FULL,K18_8D_FULL];
Q18201_FULL=[Q18201,BPM18201,Q18201];
Q18301_FULL=[Q18301,BPM18301,Q18301];
Q18401_FULL=[Q18401,BPM18401,Q18401];
Q18501_FULL=[Q18501,BPM18501,Q18501];
Q18601_FULL=[Q18601,BPM18601,Q18601];
Q18701_FULL=[Q18701,BPM18701,Q18701];
Q18801_FULL=[Q18801,BPM18801,Q18801];
Q18901_FULL=[Q18901,BPM18901,Q18901];
LI18=[LI18BEG,K18_1,DAQ1,Q18201_FULL,DAQ2,K18_2,DAQ1,Q18301_FULL,DAQ2,K18_3,DAQ1,Q18401_FULL,DAQ2,K18_4,DAQ1,Q18501_FULL,DAQ2,K18_5,DAQ1,Q18601_FULL,DAQ2,K18_6,DAQ1,Q18701_FULL,DAQ2,K18_7,DAQ1,Q18801_FULL,DAQ2,K18_8,DAQ3,Q18901_FULL,DAQ4K,BL18900,DAQ4L,WS18944,DAQ4M,CX18960,CY18960,DAQ4N,LI18END];
% ------------------------------------------------------------------------------
K19_1A_FULL=[K19_1A];
K19_1B_FULL=[K19_1B];
K19_1C_FULL=[K19_1C];
K19_1D_FULL=[K19_1D];
K19_2A_FULL=[K19_2A1,XC19202,K19_2A2,YC19203,K19_2A3];
K19_2B_FULL=[K19_2B];
K19_2C_FULL=[K19_2C];
K19_2D_FULL=[K19_2D];
K19_3A_FULL=[K19_3A1,XC19302,K19_3A2,YC19303,K19_3A3];
K19_3B_FULL=[K19_3B];
K19_3C_FULL=[K19_3C];
K19_3D_FULL=[K19_3D];
K19_4A_FULL=[K19_4A1,XC19402,K19_4A2,YC19403,K19_4A3];
K19_4B_FULL=[K19_4B];
K19_4C_FULL=[K19_4C];
K19_4D_FULL=[K19_4D1,YC57145,K19_4D2,YC57146,K19_4D3];
K19_5A_FULL=[K19_5A1,XC19502,K19_5A2,YC19503,K19_5A3];
K19_5B_FULL=[K19_5B];
K19_5C_FULL=[K19_5C];
K19_5D_FULL=[K19_5D];
K19_6A_FULL=[K19_6A1,XC19602,K19_6A2,YC19603,K19_6A3];
K19_6B_FULL=[K19_6B];
K19_6C_FULL=[K19_6C];
K19_6D_FULL=[K19_6D1,YC19700,K19_6D2,XC19700,K19_6D3];
K19_8A_FULL=[K19_8A1,XC19802,K19_8A2,YC19803,K19_8A3];
K19_1=[K19_1A_FULL,K19_1B_FULL,K19_1C_FULL,K19_1D_FULL];
K19_2=[K19_2A_FULL,K19_2B_FULL,K19_2C_FULL,K19_2D_FULL];
K19_3=[K19_3A_FULL,K19_3B_FULL,K19_3C_FULL,K19_3D_FULL];
K19_4=[K19_4A_FULL,K19_4B_FULL,K19_4C_FULL,K19_4D_FULL];
K19_5=[K19_5A_FULL,K19_5B_FULL,K19_5C_FULL,K19_5D_FULL];
K19_6=[K19_6A_FULL,K19_6B_FULL,K19_6C_FULL,K19_6D_FULL];
D19_7=[D10A,HLAM,D10B,D10,D10,D10];
K19_8=[K19_8A_FULL];
Q19201_FULL=[Q19201,BPM19201,Q19201];
Q19301_FULL=[Q19301,BPM19301,Q19301];
Q19401_FULL=[Q19401,BPM19401,Q19401];
Q19501_FULL=[Q19501,BPM19501,Q19501];
Q19601_FULL=[Q19601,BPM19601,Q19601];
Q19701_FULL=[Q19701,BPM19701,Q19701];
Q19801_FULL=[Q19801,BPM19801,Q19801];
Q19851_FULL=[Q19851,BPM19851,Q19851];
Q19871_FULL=[Q19871,BPM19871,Q19871];
LI19A=[LI19BEG,K19_1,DAQ5A,WS19144,DAQ5B,DAQ1,Q19201_FULL,DAQ2,K19_2,DAQ5A,WS19244,DAQ5B,DAQ1,Q19301_FULL,DAQ2,K19_3,DAQ5A,WS19344,DAQ5B,DAQ1,Q19401_FULL,DAQ2,K19_4,DAA7L,MSCAVEXT];
LI19B=[DBKY170,DBKY170,DAA7M,DAQ1,Q19501_FULL,DAQ2,K19_5,DAQ1,Q19601_FULL,DAQ2,K19_6,DAQ1,Q19701_FULL];
LI19C=[DAQ2,D19_7,DAQ1,Q19801_FULL,DAQ2,K19_8,D199A,Q19851_FULL,D199B1,XC19900,D199B2,YC19900,D199B3,Q19871_FULL,D199C,IM1988,D199D,LI19TERM];
LI19=[LI19A,LI19B,LI19C];
% ------------------------------------------------------------------------------
L3F_1=[BEGL3F_1,LI15,LI16,LI17,LI18,LI19A,ENDL3F_1];
L3F_2=[BEGL3F_2,LI19B,LI19C,ENDL3F_2];
L3F=[L3F_1,L3F_2];
% ==============================================================================

%CALL, FILENAME="BC20W.xsif"  FACET "W" chicane
%CALL, FILENAME="FF20W.xsif"  FACET FF/EXPT/SPECT
%CALL, FILENAME="BC20E.xsif"  FACET2 Sector 20 upgrade
%CALL, FILENAME="FF20E.xsif"  FACET2 Sector 20 upgrade
%CALL, FILENAME="FF20H.xsif"  FACET2 "hybrid"
% *** OPTICS=FACET2-28JAN2026 ***
% ==============================================================================
% FACET e- optics
% ------------------------------------------------------------------------------
% 07-JAN-2026, M. Woodley
%  * add PAX chicane per C. Emma
% ------------------------------------------------------------------------------
% 17-DEC-2021, M. Woodley
%  * change WIGE "_full" definitions for Bmad translation per C. Mayes
% 04-NOV-2021, M. Woodley
%  * concatenation of BC20W.xsif and FF20H.xsif for fall 2021 operatons
% ------------------------------------------------------------------------------
% 25-AUG-2021, M. Woodley
%  * from LI20.xsif (18MAY21 release): extract W-chicane definitions
%  * XTCAV moved into Final Focus
% ------------------------------------------------------------------------------
% ==============================================================================
% BC20W
% ==============================================================================
% ------------------------------------------------------------------------------
% global parameters
% ------------------------------------------------------------------------------
WSCL =  1 ;%wiggler field scale (=1 for wiggler ON, =0 for wiggler OFF)
% Initial parameters at BEGBC20
TW20=struct('BETX',12.2509,'ALFX',0.6685,'BETY',22.3869,'ALFY',1.1657);

% ==============================================================================
% BENDs
% ------------------------------------------------------------------------------
% chicane bends
AB1 =  0.02258935           ;%full bend angle
AB1H =  AB1/2                ;%half bend angle
ZB1 =  1.063                ;%full effective Z-length
LB1 =  ZB1*(AB1H/sin(AB1H)) ;%full path length
GB1 =  0.023                ;%full gap height
B1LE1={'be' 'B1LE' LB1/2 [AB1H GB1/2 AB1H 0 0.5 0 0]}';
B1LE2={'be' 'B1LE' LB1/2 [AB1H GB1/2 0 AB1H 0 0.5 0]}';
B1RE1={'be' 'B1RE' LB1/2 [AB1H GB1/2 AB1H 0 0.5 0 0]}';
B1RE2={'be' 'B1RE' LB1/2 [AB1H GB1/2 0 AB1H 0 0.5 0]}';
AB2 =  -AB1*1.45245;
LB2 =  1.8249;
GB2 =  0.0127;
AB3 =  -(AB1+AB2);
LB3 =  0.5287;
GB3 =  0.02065;
B2LE1={'be' 'B2LE' LB2/2 [AB2/2 GB2/2 AB2/2 0 0.5 0 0]}';
B2LE2={'be' 'B2LE' LB2/2 [AB2/2 GB2/2 0 AB2/2 0 0.5 0]}';
B3LE1={'be' 'B3LE' LB3/2 [AB3/2 GB3/2 AB3/2 0 0.5 0 0]}';
B3LE2={'be' 'B3LE' LB3/2 [AB3/2 GB3/2 0 AB3/2 0 0.5 0]}';
B3RE1={'be' 'B3RE' LB3/2 [AB3/2 GB3/2 AB3/2 0 0.5 0 0]}';
B3RE2={'be' 'B3RE' LB3/2 [AB3/2 GB3/2 0 AB3/2 0 0.5 0]}';
B2RE1={'be' 'B2RE' LB2/2 [AB2/2 GB2/2 AB2/2 0 0.5 0 0]}';
B2RE2={'be' 'B2RE' LB2/2 [AB2/2 GB2/2 0 AB2/2 0 0.5 0]}';
% define unsplit SBENs for BMAD ... not used by MAD
B1LE={'be' 'B1LE' LB1 [AB1 GB1/2 AB1H AB1H 0.5 0.5 0]}';
B2LE={'be' 'B2LE' LB2 [AB2 GB2/2 AB2/2 AB2/2 0.5 0.5 0]}';
B3LE={'be' 'B3LE' LB3 [AB3 GB3/2 AB3/2 AB3/2 0.5 0.5 0]}';
B3RE={'be' 'B3RE' LB3 [AB3 GB3/2 AB3/2 AB3/2 0.5 0.5 0]}';
B2RE={'be' 'B2RE' LB2 [AB2 GB2/2 AB2/2 AB2/2 0.5 0. 0]}';
B1RE={'be' 'B1RE' LB1 [AB1 GB1/2 AB1H AB1H 0.5 0.5 0]}';
% SLC BSY vertical wiggler (described in SLAC-PUB-3945)
% (use series approximation for sinc(x)=sin(x)/x, in case Wscl=0)
AWE =  -0.0025*WSCL ;%bend angle
ZWE =  0.244        ;%half-pole Z length
GWE =  0.02032      ;%gap height
ZDWE =  0.126525     ;%pole-to-pole Z spacing (per G. Gassner)
ZWIG =  4*ZWE+2*ZDWE ;%total wiggler Z length
AWE_2 =  AWE*AWE;
AWE_4 =  AWE_2*AWE_2;
AWE_6 =  AWE_4*AWE_2;
SINCA =  1-AWE_2/6+AWE_4/120-AWE_6/5040;
LWE =  ZWE/SINCA    ;%half-pole path length
AWES =  asin(sin(AWE)/2) ;%"short half" half-pole bend angle
AWES_2 =  AWES*AWES;
AWES_4 =  AWES_2*AWES_2;
AWES_6 =  AWES_4*AWES_2;
SINCAS =  1-AWES_2/6+AWES_4/120-AWES_6/5040;
LWES =  (ZWE/2)/SINCAS   ;%"short half" half-pole path length
AWEL =  AWE-AWES         ;%"long half" half-pole bend angle
LWEL =  LWE-LWES         ;%"long half" half-pole path length
WIGE11={'be' 'WIGE1' LWES [AWES GWE/2 0 0 0.5 0 pi/2]}';
WIGE12={'be' 'WIGE1' LWEL [AWEL GWE/2 0 AWE 0 0.5 pi/2]}';
WIGE21={'be' 'WIGE2' LWE [-AWE GWE/2 -AWE 0 0.5 0 pi/2]}';
WIGE22={'be' 'WIGE2' LWE [-AWE GWE/2 0 -AWE 0 0.5 pi/2]}';
WIGE31={'be' 'WIGE3' LWEL [AWEL GWE/2 AWE 0 0.5 0 pi/2]}';
WIGE32={'be' 'WIGE3' LWES [AWES GWE/2 0 0 0 0.5 pi/2]}';
% define unsplit SBENs for BMAD ... not used by MAD
WIGE1={'be' 'WIGE' LWE [AWE GWE/2 0 AWE 0.5 0.5 pi/2]}';
WIGE2={'be' 'WIGE' 2*LWE [-2*AWE GWE/2 -AWE -AWE 0.5 0.5 pi/2]}';
WIGE3={'be' 'WIGE3' LWE [AWE GWE/2 AWE 0 0.5 0.5 pi/2]}';
LDWE =  ZDWE/cos(AWE);
% ==============================================================================
% QUADs
% ------------------------------------------------------------------------------
% chicane quadrupoles
% NOTE: quad apertures are bore radius minus 2 mm
LQ1 =  0.5962;
AQ1 =  8.325E-3;
LQ2 =  1.0;
AQ2 =  25.0E-3;
LQ3 =  0.7142;
AQ3 =  18.6375E-3;
LQ4 =  0.7142;
AQ4 =  18.6375E-3;
LQ5 =  0.4284;
AQ5 =  18.6375E-3;
LQ6 =  0.31;
AQ6 =  8.0E-3;
% low-beta, R56 = -7 mm (symmetric)
KQ1EL =   0.811267045415;
KQ2EL =  -0.462440331546;
KQ3EL =   0.462720996981 ;%2 magnets
KQ4EL =   0.55529780432  ;%3 magnets
KQ5EL =  -0.163592075808;
KQ6E =  -1.376084837022;
KQ5ER =  -0.163597934794 ;%KQ5EL
KQ4ER =   KQ4EL          ;%3 magnets
KQ3ER =   KQ3EL          ;%2 magnets
KQ2ER =  -0.46243862549  ;%KQ2EL
KQ1ER =   0.811228495466 ;%KQ1EL
Q1EL={'qu' 'Q1EL' LQ1/2 [KQ1EL 0]}';
Q2EL={'qu' 'Q2EL' LQ2/2 [KQ2EL 0]}';
Q3EL_1={'qu' 'Q3EL_1' LQ3/2 [KQ3EL 0]}';
Q3EL_2={'qu' 'Q3EL_2' LQ3/2 [KQ3EL 0]}';
Q4EL_1={'qu' 'Q4EL_1' LQ4/2 [KQ4EL 0]}';
Q4EL_2={'qu' 'Q4EL_2' LQ4/2 [KQ4EL 0]}';
Q4EL_3={'qu' 'Q4EL_3' LQ4/2 [KQ4EL 0]}';
Q5EL={'qu' 'Q5EL' LQ5/2 [KQ5EL 0]}';
Q6E={'qu' 'Q6E' LQ6/2 [KQ6E 0]}';
Q5ER={'qu' 'Q5ER' LQ5/2 [KQ5ER 0]}';
Q4ER_1={'qu' 'Q4ER_1' LQ4/2 [KQ4ER 0]}';
Q4ER_2={'qu' 'Q4ER_2' LQ4/2 [KQ4ER 0]}';
Q4ER_3={'qu' 'Q4ER_3' LQ4/2 [KQ4ER 0]}';
Q3ER_1={'qu' 'Q3ER_1' LQ3/2 [KQ3ER 0]}';
Q3ER_2={'qu' 'Q3ER_2' LQ3/2 [KQ3ER 0]}';
Q2ER={'qu' 'Q2ER' LQ2/2 [KQ2ER 0]}';
Q1ER={'qu' 'Q1ER' LQ1/2 [KQ1ER 0]}';
% skew quadrupole
KSQ1 =  0;
SQ1={'qu' 'SQ1' LSQ/2 [KSQ1 pi/4]}';
% ==============================================================================
% SEXTs
% ------------------------------------------------------------------------------
% chicane sextupoles
% NOTE: sext apertures are bore radius minus 2 mm
LS1 =  0.250;
LS2 =  0.762;
LS3 =  0.250 ;%x2 magnets
AS =  18.6375E-3;
% from G. White's optics2013.xlsx
KS1E =    9.651773976528;
KS2E =   -8.063275568009;
KS3E =  -12.223414013976;
S1EL={'dr' 'S1EL' LS1/2 []}';
S2EL={'dr' 'S2EL' LS2/2 []}';
S3EL_1={'dr' 'S3EL_1' LS3/2 []}';
S3EL_2={'dr' 'S3EL_2' LS3/2 []}';
S3ER_1={'dr' 'S3ER_1' LS3/2 []}';
S3ER_2={'dr' 'S3ER_2' LS3/2 []}';
S2ER={'dr' 'S2ER' LS2/2 []}';
S1ER={'dr' 'S1ER' LS1/2 []}';
% ==============================================================================
% DRIFTs
% ------------------------------------------------------------------------------
% NOTE: there is a slight asymmetry in the locations of the S1E and S2E
%       left and right side sextupoles with respect to chicane center,
%       introduced in v27 at the request of the engineers
DSB1O =   0.158696902041E-7 ;%adjust z-position for rotated B1
DSB1I =  -0.433155152299E-8;
DLB1 =  (LB1-1.04)/2 ;%length reduction for drifts adjacent to B1
LD1E =  4.98-DLB1-0.0051-DSB1I;
LDMQ =  0.1;
LD1EM =  LD1E-LDMQ+0.0051;
LDTM =  1.0;
LD1ET1 =  LD1EM-LDTM-0.002515;
LDTPM =  0.169+0.002515;
LDM1Q =  LDMQ+0.038489;
LD2EA =  0.917799 ;%1.144-0.0051
LDCNHC =  0.417207 ;%0.25
LDHCYC1 =  1.554794 ;%0.35-0.160
LDCC =  0.15;
LDYC1XC1 =  LDCC-0.036;
LDCB =  0.15;
LDXC1B =  LDCB+0.11755;
LD3E =  0.25-0.05640;
LD4EA =  0.375-0.01345;
LD4EA1 =  LD4EA-0.167894;
LD4EB =  0.25-0.0103;
LDM3Q =  LDMQ+0.085052;
LD4EB1 =  LD4EB-LDM3Q+0.167894;
LDQ3E =  0.25-0.0206;
LD5EA =  0.25-0.2663;
LD5EA1 =  LD5EA+0.257175;
LD5EB =  0.15-0.256;
LD5EB1 =  LD5EB+0.324886;
LC22 =  0.224;
LD5EC =  2.46665-LD5EB-LC22-0.256;
LD5EC1 =  LD5EC-0.618061;
LDS3E =  0.25-0.0103;
LDQ4E =  0.25-0.0206;
LD6EA =  0.25-0.0103;
LDM4Q =  LD6EA;
LD6EA1 =  LD6EA-LDM4Q;
LTP =  49.113411+2*DLB1-4E-6+DSB1O;
DLTE =  -52.7E-3;
LTE =  LTP+DLTE;
AB2P =  2*asin(sin(AB1)/2);
LB2P =  LB1*AB2P/AB1;
DLB2P =  (LB2P-1.04)/2;
LD1P =  LD1E+LQ1+0.21-DLB2P-0.0051;
LD2P =  0.25-DLB2P;
LQ1P =  0.8394;
DLB3P =  DLB2P;
LD3P =  0.25-DLB3P;
LB3P =  LB2P;
LD4PA =  LB2+0.21-0.0365-DLB3P-0.00645;
LD2E =  LD1P-LD1E-LQ1+LB2P+LD2P+LQ1P+LD3P+LB3P+LD4PA-LB2;
LD6EC =  0.13317775-0.01035;
LD7E =  0.25-0.02070;
LD8E =  0.25715-0.015352;
LD6EB =  LTE/2-LB1-LB2-LB3 -LQ1-LQ2-2*LQ3-3*LQ4-LQ5-LQ6/2 -LS1-LS2-2*LS3 -LD1E-LD2E-LD3E-LD4EA-LD4EB-LD5EA-LD5EB-LD5EC -LD6EA-LD6EC-LD7E-LD8E-LDQ3E-2*LDQ4E-LDS3E-2*LC22;
LD6EB1 =  LD6EB+0.0544;
LC136 =  0.136 ;%3D4.0 corrector effective length
LD6EC1 =  LD6EC+0.0336;
LDM5Q =  0;
LD7E1 =  LD7E-LDM5Q;
LDM6Q =  LDMQ-0.005043;
LD8EM =  LD8E-LDM6Q;
LDM7Q =  0;
LD7E2 =  LD7E-LDM7Q;
LD6EC2 =  LD6EC+0.063493;
LD6EB2 =  LD6EB+0.024507;
LDM8Q =  LD6EA;
LD6EA2 =  LD6EA-LDM8Q;
LD5EC2 =  LD5EC-0.6172;
LD5EB2 =  LD5EB+0.3272;
LD5EA2 =  LD5EA+0.254;
LDM9Q =  LDMQ+0.085895;
LD4EB2 =  LD4EB-LDM9Q+0.158750;
LD4EA2 =  LD4EA-LDMQ-0.158750;
LB2XC4 =  0.264;
LXC4TCA =  0.029862;
LXC4TCB =  0.1915;
LTCWIG =  0.303373944002;
LDYAGW =  0.168162610455 ;%0.25
LDYAGQ =  0.326737389545 ;%0.25-0.0051
LDM11Q =  LDMQ+0.055967;
LDTM11 =  LDTM-0.061067;
LD1ET2 =  LD1EM-LDTM;
D1ET1A={'dr' '' 2.821013 []}';
D1ET1B={'dr' '' LD1ET1-D1ET1A{3} []}';
DTPM={'dr' '' LDTPM+0.433425 []}';
DPMM1={'dr' '' 0.339192 []}';
DM1Q={'dr' '' LDM1Q+0.014794 []}';
D2EA={'dr' '' LD2EA []}';
DCNHC={'dr' '' LDCNHC []}';
DHCYC1={'dr' '' LDHCYC1 []}';
DHCYC1A={'dr' '' 0.615+2.0*IN2M []}';
DHCYC1B={'dr' '' DHCYC1{3}-LSQ-DHCYC1A{3} []}';
DYC1E={'dr' '' LC260/2 []}';
DYC1XC1={'dr' '' LDYC1XC1 []}';
DXC1E={'dr' '' LC260/2 []}';
DXC1B={'dr' '' LDXC1B []}';
D3E={'dr' '' LD3E []}';
D4EA1={'dr' '' LD4EA1 []}';
D4EB1={'dr' '' LD4EB1+0.02787 []}';
DM3Q={'dr' '' LDM3Q-0.02787 []}';
DQ3E={'dr' '' LDQ3E []}';
D5EA1={'dr' '' LD5EA1 []}';
D5EB1={'dr' '' LD5EB1 []}';
DXC2E={'dr' '' LC260/2 []}';
D5EC1A={'dr' '' 0.745056 []}';
D5EC1B={'dr' '' LD5EC1-D5EC1A{3} []}';
DS3E={'dr' '' LDS3E []}';
DQ4E={'dr' '' LDQ4E []}';
DM4Q={'dr' '' LDM4Q []}';
D6EA1={'dr' '' LD6EA1 []}';
D6EB1={'dr' '' LD6EB1 []}';
DYC2E={'dr' '' LC136/2 []}';
D6EC1={'dr' '' LD6EC1 []}';
DM5Q={'dr' '' LDM5Q []}';
D7E1={'dr' '' LD7E1 []}';
D8EM={'dr' '' LD8EM+0.058314 []}';
DM6Q={'dr' '' LDM6Q-0.058314 []}';
D8E={'dr' '' LD8E []}';
D7E2={'dr' '' LD7E2 []}';
DM7Q={'dr' '' LDM7Q []}';
D6EC2={'dr' '' LD6EC2 []}';
DYC3E={'dr' '' LC136/2 []}';
D6EB2={'dr' '' LD6EB2 []}';
D6EA2={'dr' '' LD6EA2 []}';
DM8Q={'dr' '' LDM8Q []}';
D5EC2A={'dr' '' 0.890355 []}';
D5EC2B={'dr' '' LD5EC2-D5EC2A{3} []}';
DXC3E={'dr' '' LC260/2 []}';
D5EB2={'dr' '' LD5EB2 []}';
D5EA2={'dr' '' LD5EA2 []}';
DM9Q={'dr' '' LDM9Q []}';
D4EB2={'dr' '' LD4EB2 []}';
D4EA2={'dr' '' LD4EA2 []}';
DMQ={'dr' '' LDMQ []}';
DB2XC4={'dr' '' LB2XC4 []}';
DXC4E={'dr' '' LC260/2 []}';
DXC4TCA={'dr' '' LXC4TCA []}';
DXC4TCB={'dr' '' LXC4TCB []}';
DXTCAVF={'dr' '' 40.687*IN2M []}';
DTCWIG={'dr' '' LTCWIG []}';
DWE={'dr' '' LDWE []}';
DYAGW={'dr' '' LDYAGW+0.098819 []}';
DYAGQ={'dr' '' LDYAGQ-0.113604744002 []}';
DM11Q={'dr' '' LDM11Q+0.022716 []}';
DTM11={'dr' '' LDTM11-0.022716 []}';
D1ET2A={'dr' '' 2.452948 []}';
D1ET2B={'dr' '' LD1ET2-D1ET2A{3} []}';
% ==============================================================================
% XCORs
% ------------------------------------------------------------------------------
XC1996={'mo' 'XC1996' 0 []}';%2.031" gap
XC1E={'mo' 'XC1E' 0 []}';%0.815" gap
XCB2LE={'mo' 'XCB2LE' 0 []}';
XC2E={'mo' 'XC2E' 0 []}';%0.815" gap
XCB3LE={'mo' 'XCB3LE' 0 []}';
XCB3RE={'mo' 'XCB3RE' 0 []}';
XC3E={'mo' 'XC3E' 0 []}';%0.815" gap
XCB2RE={'mo' 'XCB2RE' 0 []}';
XC4E={'mo' 'XC4E' 0 []}';%0.815" gap
XC2460={'mo' 'XC2460' 0 []}';%2.031" gap
% ==============================================================================
% YCORs
% ------------------------------------------------------------------------------
YC1E={'mo' 'YC1E' 0 []}';%0.815" gap
YC2181={'mo' 'YC2181' 0 []}';%2.031" gap
YC2E={'mo' 'YC2E' 0 []}';%1.181" gap
YC3E={'mo' 'YC3E' 0 []}';%1.181" gap
YC2321={'mo' 'YC2321' 0 []}';%2.031" gap
YCWIGE={'mo' 'YCWIGE' 0 []}';
% ==============================================================================
% BPMs
% ------------------------------------------------------------------------------
M1E={'mo' 'M1E' 0 []}';
M3E={'mo' 'M3E' 0 []}';
MS2EL={'mo' 'MS2EL' 0 []}';
M4E={'mo' 'M4E' 0 []}';
M5E={'mo' 'M5E' 0 []}';
M6E={'mo' 'M6E' 0 []}';
M7E={'mo' 'M7E' 0 []}';
M8E={'mo' 'M8E' 0 []}';
MS2ER={'mo' 'MS2ER' 0 []}';
M9E={'mo' 'M9E' 0 []}';
M11E={'mo' 'M11E' 0 []}';
% ==============================================================================
% diagnostics, collimators, MARKERs, etc.
% ------------------------------------------------------------------------------
% profile monitors
PMON={'mo' 'PMON' 0 []}';%P202042T
SYAG={'mo' 'SYAG' 0 []}';%P202432T
% toroids
IM2040={'mo' 'IM2040' 0 []}';%T202040T
IM2452={'mo' 'IM2452' 0 []}';%T202452T
% collimators
CN2069={'dr' 'CN2069' 0 []}';%notch collimator
CX2085={'dr' 'CX2085' 0 []}';%horizontal jaw collimator
% sextupole movers
AS1EL={'mo' 'AS1EL' 0 []}';
AS2EL={'mo' 'AS2EL' 0 []}';
AS2ER={'mo' 'AS2ER' 0 []}';
AS1ER={'mo' 'AS1ER' 0 []}';
% other points of interest (INSTs go into Oracle database)
BEGBC20={'mo' 'BEGBC20' 0 []}';
CB1LE={'mo' 'CB1LE' 0 []}';
CB2LE={'mo' 'CB2LE' 0 []}';
MSEP1E={'mo' 'MSEP1E' 0 []}';
CB3LE={'mo' 'CB3LE' 0 []}';
MCE={'mo' 'MCE' 0 []}';
CB3RE={'mo' 'CB3RE' 0 []}';
MSEP2E={'mo' 'MSEP2E' 0 []}';
CB2RE={'mo' 'CB2RE' 0 []}';
CB1RE={'mo' 'CB1RE' 0 []}';
ENDBC20={'mo' 'ENDBC20' 0 []}';
% ==============================================================================
% BEAMLINEs
% ------------------------------------------------------------------------------
B1LE_FULL=[B1LE1,B1LE2];
B2LE_FULL=[B2LE1,XCB2LE,B2LE2];
B3LE_FULL=[B3LE1,XCB3LE,B3LE2];
B3RE_FULL=[B3RE1,XCB3RE,B3RE2];
B2RE_FULL=[B2RE1,XCB2RE,B2RE2];
B1RE_FULL=[B1RE1,B1RE2];
WIGE1_FULL=[WIGE11,WIGE12];
WIGE2_FULL=[WIGE21,YCWIGE,WIGE22];
WIGE3_FULL=[WIGE31,WIGE32];
WIGE_LINE=[WIGE1_FULL,DWE,WIGE2_FULL,DWE,WIGE3_FULL];
Q1EL_FULL=[Q1EL,Q1EL];
SQ1_FULL=[SQ1,SQ1];
Q2EL_FULL=[Q2EL,Q2EL];
Q3EL_1_FULL=[Q3EL_1,Q3EL_1];
Q3EL_2_FULL=[Q3EL_2,Q3EL_2];
Q4EL_1_FULL=[Q4EL_1,Q4EL_1];
Q4EL_2_FULL=[Q4EL_2,Q4EL_2];
Q4EL_3_FULL=[Q4EL_3,Q4EL_3];
Q5EL_FULL=[Q5EL,Q5EL];
Q6E_FULL=[Q6E,MCE,Q6E];
Q5ER_FULL=[Q5ER,Q5ER];
Q4ER_1_FULL=[Q4ER_1,Q4ER_1];
Q4ER_2_FULL=[Q4ER_2,Q4ER_2];
Q4ER_3_FULL=[Q4ER_3,Q4ER_3];
Q3ER_1_FULL=[Q3ER_1,Q3ER_1];
Q3ER_2_FULL=[Q3ER_2,Q3ER_2];
Q2ER_FULL=[Q2ER,Q2ER];
Q1ER_FULL=[Q1ER,Q1ER];
S1EL_FULL=[S1EL,AS1EL,S1EL];
S2EL_FULL=[S2EL,AS2EL,S2EL];
S3EL_1_FULL=[S3EL_1,S3EL_1];
S3EL_2_FULL=[S3EL_2,S3EL_2];
S3ER_1_FULL=[S3ER_1,S3ER_1];
S3ER_2_FULL=[S3ER_2,S3ER_2];
S2ER_FULL=[S2ER,AS2ER,S2ER];
S1ER_FULL=[S1ER,AS1ER,S1ER];
CHICANE1=[B1LE_FULL,CB1LE,D1ET1A,XC1996,D1ET1B,IM2040,DTPM,PMON,DPMM1,M1E,DM1Q,Q1EL_FULL,D2EA,CN2069,DCNHC,CX2085,DHCYC1A,SQ1_FULL,DHCYC1B,DYC1E,YC1E,DYC1E,DYC1XC1,DXC1E,XC1E,DXC1E,DXC1B,B2LE_FULL,CB2LE,MSEP1E,D3E,Q2EL_FULL,D4EA1,S1EL_FULL,D4EB1,M3E,DM3Q,Q3EL_1_FULL,DQ3E,Q3EL_2_FULL,D5EA1,MS2EL,S2EL_FULL,D5EB1,DXC2E,XC2E,DXC2E,D5EC1A,YC2181,D5EC1B,S3EL_1_FULL,DS3E,Q4EL_1_FULL,DQ4E,Q4EL_2_FULL,DQ4E,Q4EL_3_FULL,DM4Q,M4E,D6EA1,S3EL_2_FULL,D6EB1,DYC2E,YC2E,DYC2E,D6EC1,Q5EL_FULL,DM5Q,M5E,D7E1,B3LE_FULL,CB3LE,D8EM,M6E,DM6Q,Q6E];
CHICANE2=[Q6E,D8E,B3RE_FULL,CB3RE,D7E2,M7E,DM7Q,Q5ER_FULL,D6EC2,DYC3E,YC3E,DYC3E,D6EB2,S3ER_1_FULL,D6EA2,M8E,DM8Q,Q4ER_1_FULL,DQ4E,Q4ER_2_FULL,DQ4E,Q4ER_3_FULL,DS3E,S3ER_2_FULL,D5EC2A,YC2321,D5EC2B,DXC3E,XC3E,DXC3E,D5EB2,MS2ER,S2ER_FULL,D5EA2,Q3ER_1_FULL,DQ3E,Q3ER_2_FULL,DM9Q,M9E,D4EB2,S1ER_FULL,D4EA2,DMQ,Q2ER_FULL,D3E,MSEP2E,B2RE_FULL,CB2RE,DB2XC4,DXC4E,XC4E,DXC4E,DXC4TCA,DXC4TCB,DXTCAVF,DTCWIG,WIGE_LINE,DYAGW,SYAG,DYAGQ,Q1ER_FULL,DM11Q,M11E,DTM11,IM2452,D1ET2A,XC2460,D1ET2B,B1RE_FULL,CB1RE];
BC20W=[BEGBC20,B1LE_FULL,CB1LE,D1ET1A,XC1996,D1ET1B,IM2040,DTPM,PMON,DPMM1,M1E,DM1Q,Q1EL_FULL,D2EA,CN2069,DCNHC,CX2085,DHCYC1A,SQ1_FULL,DHCYC1B,DYC1E,YC1E,DYC1E,DYC1XC1,DXC1E,XC1E,DXC1E,DXC1B,B2LE_FULL,CB2LE,MSEP1E,D3E,Q2EL_FULL,D4EA1,S1EL_FULL,D4EB1,M3E,DM3Q,Q3EL_1_FULL,DQ3E,Q3EL_2_FULL,D5EA1,MS2EL,S2EL_FULL,D5EB1,DXC2E,XC2E,DXC2E,D5EC1A,YC2181,D5EC1B,S3EL_1_FULL,DS3E,Q4EL_1_FULL,DQ4E,Q4EL_2_FULL,DQ4E,Q4EL_3_FULL,DM4Q,M4E,D6EA1,S3EL_2_FULL,D6EB1,DYC2E,YC2E,DYC2E,D6EC1,Q5EL_FULL,DM5Q,M5E,D7E1,B3LE_FULL,CB3LE,D8EM,M6E,DM6Q,Q6E_FULL,D8E,B3RE_FULL,CB3RE,D7E2,M7E,DM7Q,Q5ER_FULL,D6EC2,DYC3E,YC3E,DYC3E,D6EB2,S3ER_1_FULL,D6EA2,M8E,DM8Q,Q4ER_1_FULL,DQ4E,Q4ER_2_FULL,DQ4E,Q4ER_3_FULL,DS3E,S3ER_2_FULL,D5EC2A,YC2321,D5EC2B,DXC3E,XC3E,DXC3E,D5EB2,MS2ER,S2ER_FULL,D5EA2,Q3ER_1_FULL,DQ3E,Q3ER_2_FULL,DM9Q,M9E,D4EB2,S1ER_FULL,D4EA2,DMQ,Q2ER_FULL,D3E,MSEP2E,B2RE_FULL,CB2RE,DB2XC4,DXC4E,XC4E,DXC4E,DXC4TCA,DXC4TCB,DXTCAVF,DTCWIG,WIGE_LINE,DYAGW,SYAG,DYAGQ,Q1ER_FULL,DM11Q,M11E,DTM11,IM2452,D1ET2A,XC2460,D1ET2B,B1RE_FULL,CB1RE,ENDBC20];
% ==============================================================================
% FF20H
% ==============================================================================
% ==============================================================================
% FACET e- optics
% ------------------------------------------------------------------------------
% 01-FEB-2022, G. White
%  * add XC1FF & YC1FF correctors (X203026 & Y203017), positions according to
%    measurements by C. Clarke
%  * MQ4FF -> M2FF for consistent naming, MONI list match spreadsheet by D. Storey
% ------------------------------------------------------------------------------
% 26-AUG-2021, M. Woodley
%  * Hybrid Final Focus for fall 2021 operation
% 23-AUG-2021, G. White
%  * Update to expected as-built condition for full Sector 20 Upgrade
%  * Removed SQ2, XC1FF, YC1FF, MS2EL, MS2ER
%  * Added Q5-Q3FF as independent quads
%  * Added second IMOVN in FFS (IM20FF)
%  * Merged with changes to master deck to describe existing experimental areas
% 23-JUN-2021, G. White
%  * Changed QFF1 & QFF2 -> Q5FF, Q4FF to match SCP
% 17-MAY-2021, G. White
%  * Removed QFF4, inserted Q0FF, Q1FF & Q2FF
%    - z locations measured by metrology: see FACET elog 05/05/2021
% 17-APR-2021, G. White
%  * Fixed as-installed Q0D, Q1D, Q2D locations: measurements by Georg
%  * Added all current expt table and spectrometer table devices according to
%    walk-through by Christine: https://docs.google.com/spreadsheets/d/
%    1Qw85KBUfSJ6Jt8tArqjcGTVlpZCtUz2hDWX8EuCRcb4/edit?usp=sharing
% 17-JUN-2020, G. White
%  * Added new beamline components to IP area as per engineering drawing from
%    D. Storey
% ------------------------------------------------------------------------------
% ==============================================================================
% transverse deflecting structure
% ------------------------------------------------------------------------------
LXTCAV =  40.687*IN2M;
XTCAVF={'tc' 'XTCAVF' LXTCAV/2 [11424 0 0*TWOPI]}';
% ==============================================================================
% BENDs
% ------------------------------------------------------------------------------
% PAX chicane
% - coils/yoke to the left (wall side)
% - offsets beam toward +X (to the left)
% - use series approximation for sinc(x)=sin(x)/x to allow ABP=0
% GBP  : chicane bend gap height (m)
% ZBP  : chicane bend "Z" length (m)
% FBP  : chicane bend field integral (measured)
% BBP  : nominal bend field (kG)
% RBP  : nominal bend radius (m) ... brho/B
% LBP  : chicane bend path length (m)
% ABP  : chicane bend angle (rad)
% ABPs : "short" half chicane bend angle (rad)
% LBPs : "short" half chicane bend path length (m)
% ABPl : "long" half chicane bend angle (rad)
% LBPl : "long" half chicane bend path length (m)
GBP =  0.008;
ZBP =  13.08*IN2M ;%0.332
FBP =  1.1112;
BBP =  11.4;
RBP =  (CB*E20)/BBP;
LBP =  RBP*asin(ZBP/RBP);
ABP =  LBP/RBP;
ABP_2 =  ABP*ABP;
ABP_4 =  ABP_2*ABP_2;
ABP_6 =  ABP_4*ABP_2;
SINCABP =  1-ABP_2/6+ABP_4/120-ABP_6/5040 ;%~sinc(ABP)=sin(ABP)/ABP
ABPS =  asin(sin(ABP)/2);
ABPS_2 =  ABPS*ABPS;
ABPS_4 =  ABPS_2*ABPS_2;
ABPS_6 =  ABPS_4*ABPS_2;
SINCABPS =  1-ABPS_2/6+ABPS_4/120-ABPS_6/5040 ;%~sinc(ABPs)=sin(ABPs)/ABPs
LBPS =  ZBP/2/SINCABPS;
ABPL =  ABP-ABPS;
LBPL =  LBP-LBPS;
%VALUE, Cb*E20*ABP BL (kG-m)
BCX203280A={'be' 'BCX203280' LBPS [-ABPS GBP/2 0 0 FBP 0 0]}';
BCX203280B={'be' 'BCX203280' LBPL [-ABPL GBP/2 0 -ABP 0 FBP 0]}';
BCX203282A={'be' 'BCX203282' LBPL [+ABPL GBP/2 +ABP 0 FBP 0 0]}';
BCX203282B={'be' 'BCX203282' LBPS [+ABPS GBP/2 0 0 0 FBP 0]}';
BCX203284A={'be' 'BCX203284' LBPS [+ABPS GBP/2 0 0 FBP 0 0]}';
BCX203284B={'be' 'BCX203284' LBPL [+ABPL GBP/2 0 +ABP 0 FBP 0]}';
BCX203286A={'be' 'BCX203286' LBPL [-ABPL GBP/2 -ABP 0 FBP 0 0]}';
BCX203286B={'be' 'BCX203286' LBPS [-ABPS GBP/2 0 0 0 FBP 0]}';
% define unsplit SBENs for BMAD ... not used by MAD
BCX203280={'be' 'BCX203280' LBP [-ABP GBP/2 0 -ABP FBP FBP 0]}';
BCX203282={'be' 'BCX20328' LBP [+ABP GBP/2 +ABP 0 FBP FBP 0]}';
BCX203284={'be' 'BCX203284' LBP [+ABP GBP/2 0 +ABP FBP FBP 0]}';
BCX203286={'be' 'BCX203286' LBP [-ABP GBP/2 -ABP 0 FBP FBP 0]}';
% spectrometer / dump bend
ABD =  0.006;
LBD =  0.9779;
GBD =  0.0635;
B5D361={'be' 'B5D36' LBD/2 [ABD/2 GBD/2 ABD/2 0 0.5 0 pi/2]}';
B5D362={'be' 'B5D36' LBD/2 [ABD/2 GBD/2 0 ABD/2 0 0.5 pi/2]}';
% define unsplit SBENs for BMAD ... not used by MAD
B5D36={'be' 'B5D36' LBD [ABD GBD/2 ABD/2 ABD/2 0.5 0.5 pi/2]}';
% ==============================================================================
% QUADs
% ------------------------------------------------------------------------------
% Final Focus
KQ5FF =  -0.467268800263;
KQ4FF =  -0.34106223812;
KQ3FF =   0.416509781058;
KQ2FF =   0.53037956706;
KQ1FF =  -0.987354552247;
KQ0FF =   KQ2FF;
Q5FF={'qu' 'Q5FF' 0.23045 [KQ5FF 0]}';
Q4FF={'qu' 'Q4FF' 0.3571 [KQ4FF 0]}';
Q3FF={'qu' 'Q3FF' 0.3571 [KQ3FF 0]}';
Q2FF={'qu' 'Q2FF' 0.3571 [KQ2FF 0]}';
Q1FF={'qu' 'Q1FF' 0.3571 [KQ1FF 0]}';
Q0FF={'qu' 'Q0FF' 0.3571 [KQ0FF 0]}';
% spectrometer
AQS =  27.0E-3;
KQ0D =  -0.328856645709;
KQ1D =   0.542213589911;
KQ2D =   KQ0D;
Q0D={'qu' 'Q0D' 0.5 [KQ0D 0]}';
Q1D={'qu' 'Q1D' 0.5 [KQ1D 0]}';
Q2D={'qu' 'Q2D' 0.5 [KQ2D 0]}';
% ==============================================================================
% DRIFTs
% ------------------------------------------------------------------------------
% hybrid Final Focus (from git/master/MAD/LI20.xsif)
D11FF={'dr' '' 0.383821 []}';
D10FF={'dr' '' 0.137748 []}';
D9FF={'dr' '' 0.304419 []}';
D8FF={'dr' '' 0.2694 []}';
D8AFF={'dr' '' 0.3877 []}';
D8BFF={'dr' '' 0.237877 []}';
D7FF={'dr' '' 0.2057 []}';
D6FF={'dr' '' 0.1957 []}';
DKRK1={'dr' '' 0.691936 []}';
DKRK2={'dr' '' 0.231287 []}';
D4FF={'dr' '' 0.484 []}';
D3AFF={'dr' '' 0.817 []}';
D3FF={'dr' '' 0.49065 []}';
DXTC1={'dr' '' 0.750525 []}';
DXTC2={'dr' '' 1.182251 []}';
D2AFF={'dr' '' 0.034973 []}';
D2FF={'dr' '' 0.829867 []}';
D1AFF={'dr' '' 0.025933 []}';
D1FF={'dr' '' 0.837967 []}';
D0FF={'dr' '' 0.027833 []}';
DMQ0FF={'dr' '' 0.17701988 []}';%0.177019
D2FFA={'dr' '' 0.591951477022 []}';%0.829867
D2FFB={'dr' '' D2FF{3}-D2FFA{3} []}';
D1FFA={'dr' '' 0.601951477023 []}';%0.837967
D1FFB={'dr' '' D1FF{3}-D1FFA{3} []}';
% experiment table
DEX20_1={'dr' '' 0.072982138228781 []}';
DEX20_2={'dr' '' 0.229949458794 []}';%0.26
DEX20_3={'dr' '' 0.06 []}';%0.1
DEX20_4={'dr' '' 0.130050541206 []}';%0.06
DEX20_5={'dr' '' 0.12 []}';
DEX20_6={'dr' '' 0.72 []}';
DEX20_7={'dr' '' 0.4017 []}';
DEX20_8={'dr' '' 0.052 []}';
DEX20_9={'dr' '' 0.4663 []}';
DEX20_10={'dr' '' 0.09 []}';
DEX20_10A={'dr' '' 0.04 []}';
DEX20_11={'dr' '' 0.04 []}';
DEX20_12={'dr' '' 1.13 []}';
DEX20_12A={'dr' '' 0.05 []}';
DEX20_13={'dr' '' 0.9 []}';
DEX20_14={'dr' '' 0.11 []}';
DEX20_15={'dr' '' 0.05 []}';
DEX20_16={'dr' '' 0.13 []}';
% spectrometer
DMQ0D={'dr' '' 0.202439 []}';
DM1QEX={'dr' '' 0.286595 []}';
D1D={'dr' '' 0.580966 []}';
DMQ1D={'dr' '' 0.356564 []}';
D2D={'dr' '' 0.754657 []}';
DMQ2D={'dr' '' 0.183182 []}';
D3D={'dr' '' 0.056305 []}';
D4D={'dr' '' 3.177152733425373 []}';
DM3BEX={'dr' '' 0.357498/cos(ABD) []}';
D5D={'dr' '' 1.22355893395943 []}';
D6D={'dr' '' 1.36 []}';
D7D={'dr' '' 1.13+1.0E-06 []}';
D8D={'dr' '' 0.5 []}';
D9D={'dr' '' 4.26-1.0E-06 []}';
D10D={'dr' '' 0.05 []}';
D11D={'dr' '' 0.23 []}';
D12D={'dr' '' 0.09 []}';
D13D={'dr' '' 0.590168 []}';
D14D={'dr' '' 1.31 []}';
LAIRG={'dr' '' 1.675292943463 []}';
% PAX chicane, and chicane bypass
ZB2BO =  0.188468 ;%0.189 outer bend-to-bend drift
ZB2BI =  ZB2BO/2         ;%half inner bend-to-bend drift (to chicane center)
LPAX =  4*ZBP+3*ZB2BO   ;%total PAX chicane Z-length
DBPO={'dr' '' ZB2BO/cos(ABP) []}';
DBPI={'dr' '' ZB2BI []}';
D3DP={'dr' '' 0.369872158794 []}';
D4DP={'dr' '' D3D{3}+D4D{3}-LPAX-D3DP{3} []}';
DPBYP1={'dr' '' ZBP+ZB2BO+ZBP+ZB2BO+ZBP/2 []}';%same Z as 3rd PAX dipole
DPBYP2={'dr' '' LPAX-DPBYP1{3} []}';
% ==============================================================================
% XCORs
% ------------------------------------------------------------------------------
XC1FF={'mo' 'XC1FF' 0 []}';%X203026
XC3FF={'mo' 'XC3FF' 0 []}';%X203086
XC4FF={'mo' 'XC4FF' 0 []}';%X203116
XC1EX={'mo' 'XC1EX' 0 []}';%X203276
% ==============================================================================
% YCORs
% ------------------------------------------------------------------------------
YC1FF={'mo' 'YC1FF' 0 []}';%Y203017
YC2FF={'mo' 'YC2FF' 0 []}';%Y203057
YC4FF={'mo' 'YC4FF' 0 []}';%Y203147
% ==============================================================================
% BPMs
% ------------------------------------------------------------------------------
M1FF={'mo' 'M1FF' 0 []}';
M2FF={'mo' 'M2FF' 0 []}';
M3FF={'mo' 'M3FF' 0 []}';
M4FF={'mo' 'M4FF' 0 []}';
M5FF={'mo' 'M5FF' 0 []}';
M0EX={'mo' 'M0EX' 0 []}';
M1EX={'mo' 'M1EX' 0 []}';
M2EX={'mo' 'M2EX' 0 []}';
M3EX={'mo' 'M3EX' 0 []}';
% ==============================================================================
% diagnostics, collimators, MARKERs, etc.
% ------------------------------------------------------------------------------
% profile monitors
USTHZ={'mo' 'USTHZ' 0 []}';%OTRS:LI20:3070 CAMR:LI20:106
DSTHZ={'mo' 'DSTHZ' 0 []}';%PD203075
IPOTR1P={'mo' 'IPOTR1P' 0 []}';%OTRS:LI20:3175 (in plasma oven)
IPOTR1={'mo' 'IPOTR1' 0 []}';%OTRS:LI20:3180 (in bypass line) CAMR:LI20:102
IPWS1={'mo' 'IPWS1' 0 []}';%WIRE:LI20:3179 (in bypass line)
IPOTR2={'mo' 'IPOTR2' 0 []}';%OTRS:LI20:3202 (in bypass line)
IPWS3={'mo' 'IPWS3' 0 []}';%WIRE:LI20:3229
DSOTR={'mo' 'DSOTR' 0 []}';%OTRS:LI20:3206 CAMR:LI20:103
WDSOTR={'mo' 'WDSOTR' 0 []}';%OTRS:LI20:3239 CAMR:LI20:104
DTOTR={'mo' 'DTOTR' 0 []}';
PGAM1={'mo' 'PGAM1' 0 []}';%PROF:LI20:3500 CMOS:LI20:3490 (Gamma 1 screen)
CNEAR={'mo' 'CNEAR' 0 []}';%CMOS:LI20:3490
PDUMP={'mo' 'PDUMP' 0 []}';%P203475T
% toroids
IQMON20={'mo' 'IQMON20' 0 []}';%TORO:LI20:3163 "Resonant charge monitor" in S20 experimental region
IM3255={'mo' 'IM3255' 0 []}';%TORO:LI20:3255
% bunch length monitors
BL20_4={'mo' 'BL20_4' 0 []}';%FF bunch length monitor
% other points of interest (INSTs go into Oracle database)
BEGFF20={'mo' 'BEGFF20' 0 []}';
MFFF={'mo' 'MFFF' 0 []}';%Q5FF entrance
KRK={'mo' 'KRK' 0 []}';%Kraken chamber focus point
DBMARK67={'mo' 'DBMARK67' 0 []}';%USTHz
ENDFF20={'mo' 'ENDFF20' 0 []}';%Q0FF exit
BEGEXPT20={'mo' 'BEGEXPT20' 0 []}';
EXTHOLE1={'mo' 'EXTHOLE1' 0 []}';%Extension table (first hole)
BEWIN1={'mo' 'BEWIN1' 0 []}';%1st Beryllium window
LCUBE={'mo' 'LCUBE' 0 []}';%Laser Injection Cube
PIC_CENT={'mo' 'PIC_CENT' 0 []}';%Center location of "picnic basket"
FILG={'mo' 'FILG' 0 []}';%Filamentation experiment gas target
FILS={'mo' 'FILS' 0 []}';%Filamentation experiment solid target
PENT={'mo' 'PENT' 0 []}';%E300 plasma entrance (Gate valve A)
MIP={'mo' 'MIP' 0 []}';%Default IP location (for optics reference)
PEXT={'mo' 'PEXT' 0 []}';%Plasma oven exit (Gate valve B)
BEWIN2={'mo' 'BEWIN2' 0 []}';%2nd Beryllium window
ENDEXPT20={'mo' 'ENDEXPT20' 0 []}';
BEGSPECT20_1={'mo' 'BEGSPECT20_1' 0 []}';
ENDSPECT20_1={'mo' 'ENDSPECT20_1' 0 []}';
BEGPAXBYP={'mo' 'BEGPAXBYP' 0 []}';
ENDPAXBYP={'mo' 'ENDPAXBYP' 0 []}';
BEGSPECT20_2={'mo' 'BEGSPECT20_2' 0 []}';
PDCBEG={'mo' 'PDCBEG' 0 []}';%Upstream end of PDC chamber
PDCEND={'mo' 'PDCEND' 0 []}';%Downstream end of PDC chamber
EDCBEG={'mo' 'EDCBEG' 0 []}';%Upstream end of EDC chamber
EDCEND={'mo' 'EDCEND' 0 []}';%Downstream end of EDC chamber
BFLYMID={'mo' 'BFLYMID' 0 []}';%Middle of Butterfly chamber
EXTWIN={'mo' 'EXTWIN' 0 []}';%Exit window (5mm thick Al)
MAINDUMP={'mo' 'MAINDUMP' 0 []}';%dump face
DBMARK30={'mo' 'DBMARK30' 0 []}';
ENDSPECT20_2={'mo' 'ENDSPECT20_2' 0 []}';
BEGPAX={'mo' 'BEGPAX' 0 []}';
PAXMID={'mo' 'PAXMID' 0 []}';
ENDPAX={'mo' 'ENDPAX' 0 []}';
% ==============================================================================
% BEAMLINEs
% ------------------------------------------------------------------------------
XTCAVF_FULL=[XTCAVF,XTCAVF];
BCX203280_FULL=[BCX203280A,BCX203280B];
BCX203282_FULL=[BCX203282A,BCX203282B];
BCX203284_FULL=[BCX203284A,BCX203284B];
BCX203286_FULL=[BCX203286A,BCX203286B];
B5D36_FULL=[B5D361,B5D362];
Q5FF_FULL=[Q5FF,Q5FF];
Q4FF_FULL=[Q4FF,Q4FF];
Q3FF_FULL=[Q3FF,Q3FF];
Q2FF_FULL=[Q2FF,Q2FF];
Q1FF_FULL=[Q1FF,Q1FF];
Q0FF_FULL=[Q0FF,Q0FF];
Q0D_FULL=[Q0D,Q0D];
Q1D_FULL=[Q1D,Q1D];
Q2D_FULL=[Q2D,Q2D];
FF=[BEGFF20,D11FF,MFFF,Q5FF_FULL,D10FF,M1FF,D9FF,BL20_4,D8FF,XC1FF,D8AFF,YC1FF,D8BFF,M2FF,D7FF,Q4FF_FULL,D6FF,Q3FF_FULL,DKRK1,KRK,DKRK2,YC2FF,D4FF,USTHZ,DBMARK67,D3AFF,DSTHZ,D3FF,XC3FF,DXTC1,XTCAVF_FULL,DXTC2,Q2FF_FULL,D2AFF,M3FF,D2FFA,YC4FF,D2FFB,Q1FF_FULL,D1AFF,M4FF,D1FFA,XC4FF,D1FFB,Q0FF_FULL,D0FF,M5FF,DMQ0FF,ENDFF20];
EXPT=[BEGEXPT20,DEX20_1,EXTHOLE1,DEX20_2,BEWIN1,DEX20_3,IQMON20,DEX20_4,DEX20_5,LCUBE,DEX20_6,PIC_CENT,DEX20_7,FILG,DEX20_8,FILS,DEX20_9,IPOTR1P,DEX20_10,IPOTR1,DEX20_10A,PENT,MIP,DEX20_11,IPWS1,DEX20_12,PEXT,DEX20_12A,IPOTR2,DEX20_13,IM3255,DEX20_14,BEWIN2,DEX20_15,IPWS3,DEX20_16,DSOTR,ENDEXPT20];
SPECT_1=[BEGSPECT20_1,DMQ0D,Q0D_FULL,DM1QEX,M0EX,D1D,WDSOTR,DMQ1D,Q1D_FULL,DM1QEX,M1EX,D2D,DMQ2D,Q2D_FULL,DM1QEX,M2EX,D3DP,ENDSPECT20_1];
PAXBYPASS=[BEGPAXBYP,DPBYP1,XC1EX,DPBYP2,ENDPAXBYP];
SPECT_2=[BEGSPECT20_2,D4DP,B5D36_FULL,DM3BEX,M3EX,D5D,PDCBEG,D6D,PDCEND,D7D,EDCBEG,D8D,EDCEND,D9D,DTOTR,D10D,BFLYMID,D11D,EXTWIN,D12D,PGAM1,D13D,CNEAR,D14D,PDUMP,LAIRG,MAINDUMP,DBMARK30,ENDSPECT20_2];
FF20H=[FF,EXPT,SPECT_1,PAXBYPASS,SPECT_2];
PAXCHICANE=[BEGPAX,BCX203280_FULL,DBPO,BCX203282_FULL,DBPI,PAXMID,DBPI,BCX203284_FULL,DBPO,BCX203286_FULL,ENDPAX];
% ==============================================================================
SECTOR20=[BC20W,FF20H];
SECTOR20_PAX=[BC20W,FF,EXPT,SPECT_1,PAXCHICANE,SPECT_2];
% ------------------------------------------------------------------------------

% *** OPTICS=FACET2-28JAN2026 ***
% ==============================================================================
% Modification History
% ------------------------------------------------------------------------------
% 31-DEC-2018, M. Woodley
%  * use original SCAVX19 kicker strength
%  * off-axis quadrupole strengths from FACET2e_patch.mad8 ... geometry not
%    rematched (assume we can steer onto target)
% ------------------------------------------------------------------------------
% 31-JAN-2017, M. Woodley
%  * from facet_v28E_posi_target.mad8
%  * change "skeleton deck" element names to LCLS-II/FACET-II format
%  * include extraction from LI19
% ------------------------------------------------------------------------------
% 16-FEB-2011, Y. Nosochkov
%    Rematch quad strengths for FACET v.28
% 24-AUG-2010, Y. Nosochkov
%    Rematch quad strengths for updated FACET v.27
% 03-MAY-2010, Y. Nosochkov
%    Rematch quad strengths for FACET v.27.
% 08-NOV-1989, M. Woodley [MDW]
%    From SCAV20A TRANS on IBM/VM disk RHH 193; surveyed positions for
%    BPMS EP01 170, PROF EP01 171, XCOR EP01 175, YCOR EP01 184,
%    XCOR EP01 185, and YCOR EP01 187; turn off ROB, RO8, -RO8, and
%    0ROL rolls
% 08-MAY-1990, MDW
%    Remove YCOR EP01 187, XCOR EP01 380, and YCOR EP01 380; add
%    YCOR EP01 183 and BPMS EP01 185; move XCOR EP01 175,
%    YCOR EP01 184, and XCOR EP01 185 ... locations of new and moved
%    diagnostic and correction devices per A. Kulikov; turn rolls back
%    on; set roll values and field for BNDS EP01 275 as per
%    current version of SCAV20A TRANS on IBM/VM-XA disk RHH 193
% 18-JAN-1991, MDW
%    New lattice for increased energy acceptance from RHH deck
%    EXTR02 TRANS; locations of devices between VLAM-174 and QF-186,
%    and COLL-186, STOP-193, and STOP-194 taken from drawings
%    ID-234-108-15 and SA-234-101-35; TORO-199 and pre-target
%    diagnostic device locations from visual survey
% 05-FEB-1991, MDW
%    Add PROF EP02 390 at location specified by A. Kulikov
% ------------------------------------------------------------------------------
% ==============================================================================
% accelerating structures
% ------------------------------------------------------------------------------
LK19_5X =  4*DLWL10;
LK19_6X =  4*DLWL10+1.0E-6 ;%kicked beam
K19_5X={'lc' 'K19_5X' LK19_5X [SBANDF P25*G19_5*LK19_5X +PHIFB3*TWOPI]}';
K19_6X={'lc' 'K19_6X' LK19_6X [SBANDF P25*G19_6*LK19_6X +PHIFB3*TWOPI]}';
% ==============================================================================
% SBEN
% ------------------------------------------------------------------------------
% LI19 scavenger kicker
GKS =  0.0254;
AKS =  0.204084364491E-3 ;%0.1969E-3
AKS_2 =  AKS*AKS;
AKS_4 =  AKS_2*AKS_2;
AKS_6 =  AKS_4*AKS_2;
SINCAKS =  1-AKS_2/6+AKS_4/120-AKS_6/5040 ;%~sinc(AKs)=sin(AKs)/AKs
LKS =  ZKS/SINCAKS;
AKSS =  asin(sin(AKS)/2);
AKSS_2 =  AKSS*AKSS;
AKSS_4 =  AKSS_2*AKSS_2;
AKSS_6 =  AKSS_4*AKSS_2;
SINCAKSS =  1-AKSS_2/6+AKSS_4/120-AKSS_6/5040 ;%~sinc(AKsS)=sin(AKsS)/AKsS
LKSS =  ZKS/(2*SINCAKSS);
AKSL =  AKS-AKSS;
LKSL =  LKS-LKSS;
BKY170A={'be' 'BKY170' LKSS [AKSS GKS/2 0 0 0.5 0 pi/2]}';
BKY170B={'be' 'BKY170' LKSL [AKSL GKS/2 0 AKS 0 0.5 pi/2]}';
% define unsplit SBENs for BMAD ... not used by MAD
BKY170={'be' 'BKY170' LKS [AKS GKS/2 0 AKS 0.5 0.5 pi/2]}';
% EP01 (bitid 57) bends
BLX57172A={'be' 'BLX57172' 1.0 [-0.91770196E-02 0.00334 0.0 0 0 0 0]}';
BLX57172B={'be' 'BLX57172' 1.0 [-0.91770196E-02 0.00334 0 -0.18360864E-01 0 0 0]}';
BLY57174A={'be' 'BLY57174' 1.0 [0.37402743E-02 0.00334 0.0 0 0 0 pi/2]}';
BLY57174B={'be' 'BLY57174' 1.0 [0.37402743E-02 0.00334 0 0.74874625E-02 0 0 pi/2]}';
BY57202A={'be' 'BY57202' 0.7798 [0.11953120E-01 0.0129 0.11955505E-01 0 0 0 pi/2]}';
BY57202B={'be' 'BY57202' 0.7798 [0.11953120E-01 0.0129 0 0.11955505E-01 0 0 pi/2]}';
BX57205A={'be' 'BX57205' 1.6383 [-0.25112588E-01 0.0129 -0.25115288E-01 0 0 0 0]}';
BX57205B={'be' 'BX57205' 1.6383 [-0.25112588E-01 0.0129 0 -0.25115288E-01 0 0 0]}';
BX57215A={'be' 'BX57215' 1.6383 [-0.25112588E-01 0.0129 -0.25115288E-01 0 0 0 0]}';
BX57215B={'be' 'BX57215' 1.6383 [-0.25112588E-01 0.0129 0 -0.25115288E-01 0 0 0]}';
BX57225A={'be' 'BX57225' 1.6383 [-0.25112588E-01 0.0129 -0.25115288E-01 0 0 0 0]}';
BX57225B={'be' 'BX57225' 1.6383 [-0.25112588E-01 0.0129 0 -0.25115288E-01 0 0 0]}';
BX57235A={'be' 'BX57235' 1.6383 [-0.25112588E-01 0.0129 -0.25115288E-01 0 0 0 0]}';
BX57235B={'be' 'BX57235' 1.6383 [-0.25112588E-01 0.0129 0 -0.25115288E-01 0 0 0]}';
BX57245A={'be' 'BX57245' 1.6383 [-0.25112588E-01 0.0129 -0.25115288E-01 0 0 0 0]}';
BX57245B={'be' 'BX57245' 1.6383 [-0.25112588E-01 0.0129 0 -0.25115288E-01 0 0 0]}';
BX57255A={'be' 'BX57255' 1.6383 [-0.25112588E-01 0.0129 -0.25115288E-01 0 0 0 0]}';
BX57255B={'be' 'BX57255' 1.6383 [-0.25112588E-01 0.0129 0 -0.25115288E-01 0 0 0]}';
BX57265A={'be' 'BX57265' 1.6383 [-0.25112588E-01 0.0129 -0.25115288E-01 0 0 0 0]}';
BX57265B={'be' 'BX57265' 1.6383 [-0.25112588E-01 0.0129 0 -0.25115288E-01 0 0 0]}';
BX57275A={'be' 'BX57275' 1.6383 [-0.25125566E-01 0.0129 -0.25115288E-01 0 0 0 0]}';
BX57275B={'be' 'BX57275' 1.6383 [-0.25125566E-01 0.0129 0 -0.25115288E-01 0 0 0]}';
% define unsplit SBENs for BMAD ... not used by MAD
LBLX57 =  2*BLX57172A{3} ;
ABLX57 =  2*BLX57172A{4};
LBLY57 =  2*BLY57174A{3} ;
ABLY57 =  2*BLY57174A{4};
LBY57 =  2*BY57202A{3}  ;
ABY57 =  2*BY57202A{4};
LBX57 =  2*BX57205A{3}  ;
ABX57 =  2*BX57205A{4};
ABX57E =  2*BX57275A{4};
BLX57172={'be' 'BLX5717' LBLX57 [ABLX57 0.00334 0.0 -0.18360864E-01 0 0 0]}';
BLY57174={'be' 'BLY57174' LBLY57 [ABLY57 0.00334 0.0 +0.74874625E-02 0 0 pi/2]}';
BY57202={'be' 'BY5720' LBY57 [ABY57 0.0129 +0.11955505E-01 +0.11955505E-01 0 0 pi/2]}';
BX57205={'be' 'BX57205' LBLX57 [ABX57 0.0129 -0.25115288E-01 -0.25115288E-01 0 0 0]}';
BX57215={'be' 'BX57215' LBLX57 [ABX57 0.0129 -0.25115288E-01 -0.25115288E-01 0 0 0]}';
BX57225={'be' 'BX57225' LBLX57 [ABX57 0.0129 -0.25115288E-01 -0.25115288E-01 0 0 0]}';
BX57235={'be' 'BX57235' LBLX57 [ABX57 0.0129 -0.25115288E-01 -0.25115288E-01 0 0 0]}';
BX57245={'be' 'BX57245' LBLX57 [ABX57 0.0129 -0.25115288E-01 -0.25115288E-01 0 0 0]}';
BX57255={'be' 'BX57255' LBLX57 [ABX57 0.0129 -0.25115288E-01 -0.25115288E-01 0 0 0]}';
BX57265={'be' 'BX57265' LBLX57 [ABX57 0.0129 -0.25115288E-01 -0.25115288E-01 0 0 0]}';
BX57275={'be' 'BX57275' LBLX57 [ABX57E 0.0129 -0.25115288E-01 -0.25115288E-01 0 0 0]}';
% ==============================================================================
% QUAD (off-axis)
% ------------------------------------------------------------------------------
% run MGEO1 match; then use offAxisQuad.m to compute initial angles;
% then use MGEO2 match to fine-tune the kicker angle and Q19701 bend angle
AKQ501 =  -0.697627545E-5;
AKQ601 =  +0.1010490572E-3;
AKQ701 =  -0.339144696197E-3 ;%-0.3388633894E-3
%comment
Q19501X={'bg' 'Q19501X' Q19501{3} [AKQ501 0 0 0 0 0 pi/2 -Q19501{4}(1)]}';
Q19601X={'bg' 'Q19601X' Q19601{3} [AKQ601 0 0 0 0 0 pi/2 -Q19601{4}(1)]}';
Q19701X={'bg' 'Q19701X' Q19701{3} [AKQ701 0 0 0 0 0 pi/2 -Q19701{4}(1)]}';
%endcomment
% 
% Q19501X={'bg' 'Q19501X' LQE/2 [AKQ501 0 0 0 0 0 pi/2 -QSIGN*(KQ19501)]}';
% Q19601X={'bg' 'Q19601X' LQE/2 [AKQ601 0 0 0 0 0 pi/2 -QSIGN*(KQ19601)]}';
% Q19701X={'bg' 'Q19701X' LQE/2 [AKQ701 0 0 0 0 0 pi/2 -QSIGN*(KQ19701)]}';

% ==============================================================================
% QUAD
% ------------------------------------------------------------------------------
KQ57186 =   0.434861284035             ;
KQ57190 =  -0.573614963712              ;
KQ57204 =   0.453180832142               ;
KQ57210 =  -0.508482870713               ;
KQ57220 =   0.738426550257              ;
KQ57280 =   0.686348881151               ;
Q57186={'qu' 'Q57186' 0.2486 [KQ57186 0]}';
Q57190={'qu' 'Q57190' 0.2799 [KQ57190 0]}';
Q57204={'qu' 'Q57204' 0.3077 [KQ57204 0]}';
Q57210={'qu' 'Q57210' 0.3077 [KQ57210 0]}';
Q57220={'qu' 'Q57220' 0.3077 [KQ57220 0]}';
Q57230={'qu' 'Q57230' 0.3077 [KQ57210 0]}';
Q57240={'qu' 'Q57240' 0.3077 [KQ57220 0]}';
Q57250={'qu' 'Q57250' 0.3077 [KQ57210 0]}';
Q57260={'qu' 'Q57260' 0.3077 [KQ57220 0]}';
Q57270={'qu' 'Q57270' 0.3077 [KQ57210 0]}';
Q57280={'qu' 'Q57280' 0.3077 [KQ57280 0]}';
% ==============================================================================
% DRIF
% ------------------------------------------------------------------------------
DRI34001={'dr' '' 0.2506 []}';
DRI34002={'dr' '' 0.2413 []}';
DRI34003={'dr' '' 0.4547 []}';
DRI34004={'dr' '' 0.3556 []}';
DRI34005={'dr' '' 0.99 []}';
DRI34006={'dr' '' 1.475 []}';
DRI34007={'dr' '' 0.1651 []}';
DRI34008={'dr' '' 0.2093 []}';
DRI34009={'dr' '' 1.7501 []}';
DRI34010={'dr' '' 3.6591 []}';
DRI34011={'dr' '' 0.4008 []}';
DRI34012={'dr' '' 0.3073 []}';
DRI34013={'dr' '' 0.3907 []}';
DRI34014={'dr' '' 0.3747 []}';
DRI34015={'dr' '' 0.2248 []}';
DRI34016={'dr' '' 0.1476 []}';
DRI34017={'dr' '' 0.1518 []}';
DRI34018={'dr' '' 0.0875 []}';
DRI34019={'dr' '' 0.3868 []}';
DRI34020={'dr' '' 2.8567 []}';
DRI34021={'dr' '' 1.8288 []}';
DRI34022={'dr' '' 1.2081 []}';
DRI34023={'dr' '' 0.2838 []}';
DRI34024={'dr' '' 0.2765 []}';
DRI34025={'dr' '' 0.5051 []}';
DRI34026={'dr' '' 0.2921 []}';
DRI34027={'dr' '' 0.2159 []}';
DRI34028={'dr' '' 1.2637 []}';
DRI34029={'dr' '' 0.1969 []}';
DRI34030={'dr' '' 0.1524 []}';
DRI34031={'dr' '' 0.7366 []}';
DRI34032={'dr' '' 0.254 []}';
DRI34033={'dr' '' 0.6486 []}';
DRI34034={'dr' '' 0.6 []}';
% ==============================================================================
% SROT
% ------------------------------------------------------------------------------
ROB={'ro' 'ROB' 0 [-(0.77600829E-01)]}';
RO8={'ro' 'RO8' 0 [-(0.38221140E-01)]}';
NEG_RO8={'ro' 'NEG_RO8' 0 [-(-0.38221140E-01)]}';
ZERO_ROL={'ro' 'ZERO_ROL' 0 [-(-0.85532653E-01)]}';
% ==============================================================================
% XCORs and YCORs
% ------------------------------------------------------------------------------
XC57175={'mo' 'XC57175' 0 []}';
XC57185={'mo' 'XC57185' 0 []}';
XC57205={'mo' 'XC57205' 0 []}';
XC57215={'mo' 'XC57215' 0 []}';
XC57225={'mo' 'XC57225' 0 []}';
XC57235={'mo' 'XC57235' 0 []}';
XC57245={'mo' 'XC57245' 0 []}';
XC57255={'mo' 'XC57255' 0 []}';
XC57265={'mo' 'XC57265' 0 []}';
XC57275={'mo' 'XC57275' 0 []}';
XC57282={'mo' 'XC57282' 0 []}';
YC57183={'mo' 'YC57183' 0 []}';
YC57184={'mo' 'YC57184' 0 []}';
YC57210={'mo' 'YC57210' 0 []}';
YC57220={'mo' 'YC57220' 0 []}';
YC57230={'mo' 'YC57230' 0 []}';
YC57240={'mo' 'YC57240' 0 []}';
YC57250={'mo' 'YC57250' 0 []}';
YC57260={'mo' 'YC57260' 0 []}';
YC57270={'mo' 'YC57270' 0 []}';
YC57280={'mo' 'YC57280' 0 []}';
YC57282={'mo' 'YC57282' 0 []}';
% ==============================================================================
% diagnostics
% ------------------------------------------------------------------------------
% BPMs (LCLS-II type designations)
BPM57170={'mo' 'BPM57170' 0 []}';
BPM57175={'mo' 'BPM57175' 0 []}';
BPM57185={'mo' 'BPM57185' 0 []}';
BPM57186={'mo' 'BPM57186' 0 []}';
BPM57190={'mo' 'BPM57190' 0 []}';
BPM57204={'mo' 'BPM57204' 0 []}';
BPM57210={'mo' 'BPM57210' 0 []}';
BPM57220={'mo' 'BPM57220' 0 []}';
BPM57230={'mo' 'BPM57230' 0 []}';
BPM57240={'mo' 'BPM57240' 0 []}';
BPM57250={'mo' 'BPM57250' 0 []}';
BPM57260={'mo' 'BPM57260' 0 []}';
BPM57270={'mo' 'BPM57270' 0 []}';
BPM57280={'mo' 'BPM57280' 0 []}';
BPM57383={'mo' 'BPM57383' 0 []}';
BPM57400={'mo' 'BPM57400' 0 []}';
% misc
PR57171={'mo' 'PR57171' 0 []}';
PC57175={'dr' 'PC57175' 0 []}';
IM57175={'mo' 'IM57175' 0 []}';
PC57178={'dr' 'PC57178' 0.9525 []}';
WS57184={'mo' 'WS57184' 0 []}';
PR57185={'mo' 'PR57185' 0 []}';
PC57186A={'dr' 'PC57186A' 0.9563 []}';
PC57186C={'dr' 'PC57186C' 0.4176 []}';
ST57193={'mo' 'ST57193' 0 []}';
ST57194={'mo' 'ST57194' 0 []}';
IM57199={'mo' 'IM57199' 0 []}';
SP57281={'mo' 'SP57281' 0 []}';
BZ57372={'mo' 'BZ57372' 0 []}';
IM57375={'mo' 'IM57375' 0 []}';
IM57376={'mo' 'IM57376' 0 []}';
PR57385={'mo' 'PR57385' 0 []}';
PC57388={'dr' 'PC57388' 0 []}';
PR58390={'mo' 'PR58390' 0 []}';
% ==============================================================================
% MARK
% ------------------------------------------------------------------------------
BEGSCAV={'mo' 'BEGSCAV' 0 []}';
DBMARK32={'mo' 'DBMARK32' 0 []}';%e+ production target
ENDSCAV={'mo' 'ENDSCAV' 0 []}';
% ==============================================================================
% BEAMLINEs
% ------------------------------------------------------------------------------
BKY170_FULL=[BKY170A,BKY170B];
Q19501X_FULL=[Q19501X,BPM19501,Q19501X];
Q19601X_FULL=[Q19601X,BPM19601,Q19601X];
Q19701X_FULL=[Q19701X,BPM19701,Q19701X];
SCAV19X=[BEGSCAV,BKY170_FULL,DAA7M,DAQ1,Q19501X_FULL,DAQ2,K19_5X,DAQ1,Q19601X_FULL,DAQ2,K19_6X,DAQ1,Q19701X_FULL,MSCAV];
BLX57172_FULL=[BLX57172A,BLX57172B];
BLY57174_FULL=[BLY57174A,BLY57174B];
BY57202_FULL=[BY57202A,BY57202B];
BX57205_FULL=[BX57205A,XC57205,BX57205B];
BX57215_FULL=[BX57215A,XC57215,BX57215B];
BX57225_FULL=[BX57225A,XC57225,BX57225B];
BX57235_FULL=[BX57235A,XC57235,BX57235B];
BX57245_FULL=[BX57245A,XC57245,BX57245B];
BX57255_FULL=[BX57255A,XC57255,BX57255B];
BX57265_FULL=[BX57265A,XC57265,BX57265B];
BX57275_FULL=[BX57275A,XC57275,BX57275B];
Q57186_FULL=[Q57186,Q57186];
Q57190_FULL=[Q57190,Q57190];
Q57204_FULL=[Q57204,Q57204];
Q57210_FULL=[Q57210,YC57210,Q57210];
Q57220_FULL=[Q57220,YC57220,Q57220];
Q57230_FULL=[Q57230,YC57230,Q57230];
Q57240_FULL=[Q57240,YC57240,Q57240];
Q57250_FULL=[Q57250,YC57250,Q57250];
Q57260_FULL=[Q57260,YC57260,Q57260];
Q57270_FULL=[Q57270,YC57270,Q57270];
Q57280_FULL=[Q57280,YC57280,Q57280];
SCAV20A=[DRI34001,BPM57170,DRI34002,PR57171,DRI34003,BLX57172_FULL,DRI34004,BLY57174_FULL,DRI34005,XC57175,DRI34006,PC57175,DRI34007,BPM57175,DRI34008,IM57175,DRI34009,PC57178,DRI34010,YC57183,DRI34011,YC57184,DRI34012,WS57184,DRI34013,XC57185,DRI34014,PR57185,DRI34015,BPM57185,DRI34016,BPM57186,Q57186_FULL,DRI34017,PC57186A,DRI34018,PC57186C,DRI34019,BPM57190,Q57190_FULL,DRI34020,ST57193,DRI34021,ST57194,DRI34022,IM57199,DRI34023,BY57202_FULL,DRI34024,ROB,BPM57204,Q57204_FULL,DRI34024,BX57205_FULL,DRI34024,BPM57210,Q57210_FULL,DRI34024,BX57215_FULL,DRI34024,BPM57220,Q57220_FULL,DRI34024,BX57225_FULL,DRI34024,BPM57230,Q57230_FULL,DRI34024,BX57235_FULL,DRI34024,BPM57240,Q57240_FULL,DRI34024,BX57245_FULL,DRI34024,BPM57250,Q57250_FULL,DRI34024,BX57255_FULL,DRI34024,BPM57260,Q57260_FULL,DRI34024,BX57265_FULL,DRI34024,BPM57270,Q57270_FULL,DRI34024,RO8,BX57275_FULL,NEG_RO8,DRI34024,BPM57280,Q57280_FULL,ZERO_ROL,DRI34025,SP57281,DRI34026,YC57282,DRI34027,XC57282,DRI34028,BZ57372,DRI34029,IM57375,DRI34030,IM57376 ,DRI34031,BPM57383,DRI34030,PR57385,DRI34032,PC57388,DRI34033,PR58390,DRI34034,BPM57400,DBMARK32,ENDSCAV];
SCAV=[SCAV19X,SCAV20A];
% ------------------------------------------------------------------------------
% extraction
% ------------------------------------------------------------------------------
% use MGEO1 to compute
KKY170={'mo' 'KKY170' ZKS []}';%-0.19169E-3
KQ701={'mo' 'KQ701' 0 []}';% 0
LI19X=[BEGSCAV,KKY170,DAA7M,DAQ1,Q19501,Q19501,DAQ2,K19_5,DAQ1,Q19601,Q19601,DAQ2,K19_6,DAQ1,Q19701,KQ701,Q19701,MSCAV];
% ==============================================================================

% *** OPTICS=FACET2-28JAN2026 ***
% F2_ELEC     : BEGINJ /ENDSPECT20
% F2_SCAV     : BEGSCAV/ENDSCAV
% F2_ELEC_PAX : PAXBEG /PAXEND












































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































BC14=[BC14_1,BC14E,BC14_2];%electron side
FACET2E=[DL10,L1F,BC11,L2F,BC14,L3F,SECTOR20];
FACET2E_PAX=[DL10,L1F,BC11,L2F,BC14,L3F,SECTOR20_PAX];
FACET2S=[DL10,L1F,BC11,L2F,BC14,L3F_1,SCAV];
% beam path definitions
%F2_ELEC     : e- gun to LI20 dump (PAX chicane bypassed)
%F2_ELEC_PAX : e- gun to LI20 dump
%F2_SCAV     : e- gun to e+ production target
%F2_S10AIP   : e- gun line for AIP injector tests/comissioning
F2_ELEC=[INJ,FACET2E];
F2_ELEC_PAX=[INJ,FACET2E_PAX];
F2_SCAV=[INJ,FACET2S];
F2_S10AIP=[INJS10AIP];
% ------------------------------------------------------------------------------
% SURVEY coordinates
% ------------------------------------------------------------------------------
% at CATHODEF (Gun moved 50.47 cm closer to gun c.f. LCLS-I)
LLL =  7.51*0.3048-1.42  ;%loadlock length [m]
XLL =    10.693567344496 ;%X at loadlock start [m]
ZLL =  1001.562110341    ;%Z at loadlock start [m]
XC =  XLL+LLL*sin(ADL1) ;%X at cathode [m]  10.12329735 (LCLS-I = 10.448934873335)
YC =  0                 ;%Y at cathode [m]
ZC =  ZLL+LLL*cos(ADL1) ;%Z at cathode [m]  1002.376541 (LCLS-I = 1001.911433068)
THETAC =  ADL1                                ;%-35*RADDEG
PHIC =  0;
PSIC =  0;
% at BEGDL10
LINJ =  7.955897298 ;%was 8.398441604
XI =  XC+LINJ*sin(ADL1);
YI =  YC;
ZI =  ZC+LINJ*cos(ADL1);
THETAI =  THETAC;
PHII =  PHIC;
PSII =  PSIC;
% at MSCAVEXT
Z19 =  1877.228;
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% BETA0 block definitions
% ------------------------------------------------------------------------------
TW0=struct('BETX',BX0,'ALFX',AX0,'BETY',BY0,'ALFY',AY0);
TWI=struct('BETX',BXI,'ALFX',AXI,'BETY',BYI,'ALFY',AYI);
% misc


% ------------------------------------------------------------------------------
% BEAM definition (from FACET2e_baseline.mat)
% ------------------------------------------------------------------------------


% ------------------------------------------------------------------------------
% commands
% ------------------------------------------------------------------------------



%CALL, FILENAME="FACET2e_match.mad8"
%CALL, FILENAME="FACET2e_makeSymbols.mad8"
%CALL, FILENAME="FACET2e_makeElegant.mad8"
% ------------------------------------------------------------------------------
%STOP
% ------------------------------------------------------------------------------
% %for testing the Matlab model
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

% ------------------------------------------------------------------------------
%COMMENT standard output







% 
% 
% 
% 
% 
% 
% 
% 


% 
% 
% 
% 
% 
% 
% 
% 




















%ENDCOMMENT
% ------------------------------------------------------------------------------

