% FACET2E_LUCRETIA Parse FACET2 MAD electron deck into Lucretia
% Run from MAD folder

global BEAMLINE WF
% Clear existing beamline database if exists
if ~isempty(BEAMLINE)
  BEAMLINE={};
end
close all

% Read in deck from mad8 definitions
DT=DeckTool('XSIF');
DT.ReadDeck('FACET2e.mad8','F2_ELEC','TW0','BEAM0');

% Beam defnition and inital match from GPT tracking simulation
% load ../Lucretia/beams/FACET2e_op.mat Beam Qmatch Initial % 30-deg schottky phase, current operational gun params
load ../Lucretia/beams/FACET2e_2bunch_op.mat Beam Qmatch Initial % 30-deg schottky phase, current operational gun params for double-pulse on cathode
for iquad=1:length(Qmatch.Quads)
  ele=findcells(BEAMLINE,'Name',char(Qmatch.Quads(iquad)));
  for iele=ele
    BEAMLINE{iele}.B=Qmatch.B(iquad)/2;
  end
end

% Set survey coordinates
L0a=findcells(BEAMLINE,'Name','L0AF__1');
dL=1.474748004-BEAMLINE{L0a}.S;
ANGINJ=-deg2rad(35);
X0=10.448934873335+dL*sin(ANGINJ);
Z0=1001.911433068+dL*cos(ANGINJ);
SetFloorCoordinates(1,length(BEAMLINE),[X0 0 Z0 ANGINJ 0 0]);

% Match optics to get correct Twiss parameters in main Linac
M=Match;
M.beam=MakeBeam6DGauss(Initial,1e3,3,1);
M.iInitial=1;
M.initStruc=Initial;
M.verbose=false; % see optimizer output or not
M.optim='fminsearch';
M.optimDisplay='off';
M.useParallel=false;
M.addVariable('INITIAL',1,'betax',0.1,150);
M.addVariable('INITIAL',1,'betay',0.1,150);
M.addVariable('INITIAL',1,'alphax',-50,50);
M.addVariable('INITIAL',1,'alphay',-50,50);
ENDBC11=findcells(BEAMLINE,'Name','BC11CEND');
M.addMatch(ENDBC11,'beta_x',3,1e-4);
M.addMatch(ENDBC11,'beta_y',3,1e-4);
M.addMatch(ENDBC11,'alpha_x',0,1e-5);
M.addMatch(ENDBC11,'alpha_y',0,1e-5);
M.doMatch;
display(M)
Initial=M.initStruc;
Initial.SigPUncorrel=0.125*0.1e-2;

% Set default Sector 20 optics
match_S20(Initial,"pwfa_50cm",1,1);
T=TwissPlot(1,length(BEAMLINE),Initial,[1 1 0]);

% S-band structure apertures
for iele=findcells(BEAMLINE,'Class','LCAV')
  BEAMLINE{iele}.aper=0.00955;
end

% Database initialization
SetSPositions(1,length(BEAMLINE),0);
% SetElementSlices(1,length(BEAMLINE));
% SetElementBlocks(1,length(BEAMLINE));

% Check end co-ordinates
dX_end = BEAMLINE{end}.Coordf-[0 -0.0796021612227596 2019.20524109153];
dANG_end = BEAMLINE{end}.Anglef-[0 -0.006 0];
disp('END BEAMLINE POSITION ERROR:')
disp(dX_end)
disp('END BEAMLINE ANGLE ERROR:')
disp(dANG_end)

% Save lattices
% - BMAD
% DT=DeckTool('BMAD');
% DT.WriteDeck(Initial,'FACET2e.bmad','FACET2e',true) ;
% - Luctretia
% convert 0 angle bends into drifts- breaks Lucretia tracking
for iele=1:length(BEAMLINE)
  if strcmp(BEAMLINE{iele}.Class,'SBEN') && BEAMLINE{iele}.Angle==0
    BLe=BEAMLINE{iele};
    BEAMLINE{iele}=DrifStruc(BEAMLINE{iele}.L,BEAMLINE{iele}.Name);
    cf={'P' 'Coordi' 'Anglei' 'Coordf' 'Anglef' 'S'};
    for icf=1:length(cf); BEAMLINE{iele}.(cf{icf}) = BLe.(cf{icf}) ; end
  end
end
save FACET2e.mat BEAMLINE WF Initial Beam

