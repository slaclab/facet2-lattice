function match_S20(Initial,ConfigName,Nbunch,SextupoleMatch)
%MATCHS20 Configure Sector 20 magnets by matching in-memory Lucretia file
% First load Lucretia file from facet2-lattice/Lucretia/models repository directory.
% Available configurations depicted in presentation files on:
% https://www.slac.stanford.edu/~whitegr/F2S2E/Sector20_optics.pdf
% https://www.slac.stanford.edu/~whitegr/F2_S2E/F2_S2E.pdf
%
% Use DeckTool to convert to new XSIF or BMAD lattices
%
% (Matlab optimization toolbox required)
%
%match_s20(Initial,ConfigurationName,Nbunch)
%
% Initial : Lucretia Initial structire from loaded repository model
%
% ConfigurationName :
% "PWFA_5cm" - PWFA oven IP (PENT), beta*_x,y = 5cm (default in stored lattice, for 1 bunch mode)
% "PWFA_10cm" - PWFA oven IP (PENT), beta*_x,y = 10cm
% "PWFA_50cm" - PWFA oven IP (PENT), beta*_x,y = 50cm
% "PWFA_100cm" - PWFA oven IP (PENT), beta*_x,y = 100cm
% "TCAV" - Optics optimized for TCAV measurements
% "SFQED" - Optics optimized for SFQED experiment, beta* = 10m, E=13 GeV
% "Kraken_RoundBeam" - Beam waist at Kraken chamber, round focus (29 x 29 cm)
% "Kraken_FlatBeam" - Beam waist at Kraken chamber, flat beam at focus (200 x 20 cm)
% "Filamentation_Solid" - Optimized for Filamentation solid target experiment
% "Filamentation_Gas" - Optimized for Filamentation gas target experiment
%
% BC20 chicane matching:
%  "..._R56=N" - Append R56 request to re-match chicane (-10:10 (mm) is matchable range) - request R56 is in mm units
%
% Nbunch = 1 or 2
%  If 2, then optimize for second ("witness") bunch in 2-bunch configuration (with energy determined from tracking studies)
%
%match_s20(Initial,ConfigurationName,Nbunch,true)
%
% Also match sextupoles in Sectro 20 chicane to minimize chromatic abberations at the IP
global BEAMLINE PS

if nargin<3
  error('Check input arguments');
end

% Assumed PS not assigned in base lattice
if ~isempty(PS)
  error('PS not empty')
end

if ~exist('SextupoleMatch','var')
  SextupoleMatch=false;
end

odisp='iter';
% odisp='off';

ipname='PENT';
if ~exist('E0','var')
  E0=[];
end
if Nbunch==2
  E0=9.921; % Witness bunch energy
end
de=0.8e-2; % rms energy spread in S20
psno=1:5;
ConfigName = lower(string(ConfigName)) ;

% Pull off chicane match request first
t=regexp(ConfigName,'_r56=(.+)','tokens','once');
if isempty(t)
  dochicane=false;
else
  dochicane=true;
  r56=str2double(t{1});
  ConfigName=regexprep(ConfigName,'_r56=.+$','');
end

switch ConfigName
  case "pwfa_5cm"
    ipbeta=[0.05 0.05];
  case "pwfa_10cm"
    ipbeta=[0.1 0.1];
  case "pwfa_50cm"
    ipbeta=[0.5 0.5];
  case "pwfa_100cm"
    ipbeta=[1 1];
  case "tcav"
    ipbeta=[0.05 0.05];
    ipname='DSOTR';
  case "sfqed"
    psno=2:5;
    ipbeta=[10 10];
    E0=13;
    de=0.1e-2;
  case "kraken_roundbeam"
    ipname='KRK';
    ipbeta=[0.29 0.29];
    psno=3:5;
  case "kraken_flatbeam"
    ipname='KRK';
    ipbeta=[2 0.2];
    psno=3:5;
  case "filamentation_solid"
    ipbeta=[0.05 0.05];
    ipname='FILS';
  case "filamentation_gas"
    ipbeta=[0.7 0.7];
    ipname='FILG';
  otherwise
    error('Unknown configuration');
end

% Required beamline indices
pent=findcells(BEAMLINE,'Name','PENT');
iscr=findcells(BEAMLINE,'Name','PDUMP');
i1=findcells(BEAMLINE,'Name','BEGBC20');
l3_1=findcells(BEAMLINE,'Name','BEGL3F_1');
% l3_2=findcells(BEAMLINE,'Name','ENDL3F_1');
ipele=findcells(BEAMLINE,'Name',ipname);

% Set specific design energy if requested
if ~isempty(E0)
  SetDesignMomentumProfile(l3_1,length(BEAMLINE),Initial.Q,BEAMLINE{l3_1}.P,E0);
end

% Get Initial structure at start of S20
[~,T]=GetTwiss(1,i1,Initial.x.Twiss,Initial.y.Twiss);
I=TwissToInitial(T,i1,Initial);
I.Q=2e-9;
I.x.NEmit=3.0e-6; I.y.NEmit=3e-6;
I.SigPUncorrel=I.Momentum.*de;

if dochicane
  M=Match;
  chquad_name={{'Q1EL' 'Q1ER'} {'Q2EL' 'Q2ER'} ...
    {'Q3EL_1' 'Q3EL_2' 'Q3ER_1' 'Q3ER_2'} ...
    {'Q4EL_1' 'Q4EL_2' 'Q4ER_2' 'Q4ER_3'} ...
    {'Q5EL' 'Q5ER'} {'Q6E'}};
  pslist=[];
  for ich=1:length(chquad_name)
    bl=[];
    for iq=1:length(chquad_name{ich})
      bl=[bl findcells(BEAMLINE,'Name',chquad_name{ich}{iq})]; %#ok<*AGROW> 
    end
    AssignToPS(  bl, length(PS)+1 );
    pslist(end+1)=length(PS);
    MovePhysicsVarsToPS(length(PS));
  end
  BMAX=[400 -400 330 320 -120 -300].*[2 2 4 4 2 1]; % from db in kG (LGPS)
  BMIN=[0 0 0 0 0 0]; % from db in kG (LGPS)
  B0=[388.633 -310.417 422.236 669.562 -51.9686 -142.578];
  for ips=1:length(pslist)
    lim1=BMIN(ips)/10;
    lim2=BMAX(ips)/10;
    if lim2>=lim1
      M.addVariable('PS', pslist(ips),'Ampl',lim1,lim2);
    else
      M.addVariable('PS', pslist(ips),'Ampl',lim2,lim1);
    end
  end
  M.beam=MakeBeam6DGauss(I,1000,3,1);
  M.iInitial=i1;
  M.initStruc=I;
  M.verbose=false; % see optimizer output or not
  M.optim='lsqnonlin';
  M.optimDisplay='iter';
  mce=findcells(BEAMLINE,'Name','MCE');
  cend=findcells(BEAMLINE,'Name','MFFF');
  bx=findcells(BEAMLINE,'Name','DQ3E'); bx=bx(1);
  by=findcells(BEAMLINE,'Name','Q2EL'); by=by(1);
  M.addMatch(bx,'beta_x',0,200); %#ok<*UNRCH>
  M.addMatch(by,'beta_y',0,200);
  M.addMatch(mce,'alpha_x',0,0.001);
  M.addMatch(mce,'alpha_y',0,0.001);
  M.addMatch(mce,'beta_y',0,200);
  M.addMatch(pent,'R',r56*1e-3,1e-5,'56');
  M.addMatch(cend,'eta_x',0,1e-5);
  M.addMatch(cend,'etap_x',0,1e-5);
  M.addMatch(mce,'etap_x',0,1e-5);
  M.doMatch;
  disp(M);
  for ips=1:length(PS)
    for iele=PS(ips).Element
      if GetTrueStrength(iele)==0
        BEAMLINE{iele}.B=0;
      else
        RenormalizePS(ips);
      end
      BEAMLINE{iele}.PS=0;
    end
  end
  PS=[];
end

% Form power supplies for matching magnet strengths
% PS(1) = [Q0 Q2]; PS(2:5) = [Q1, Q3, Q4, Q5]
qm={'Q0FF' 'Q1FF' 'Q2FF' 'Q3FF' 'Q4FF' 'Q5FF'};
iele=[findcells(BEAMLINE,'Name','Q0FF') findcells(BEAMLINE,'Name','Q2FF')];
AssignToPS( iele, 1 ) ;
for iquad=[2 4 5 6]
  iele=findcells(BEAMLINE,'Name',qm{iquad});
  AssignToPS( iele, length(PS)+1 ) ;
end
MovePhysicsVarsToPS(1:length(PS));

% Match IP
% - If SFQED, then operate final triplet as a doublet
if ConfigName == "sfqed"
  PS(1).Ampl=0; PS(2).Ampl=0;
end
if any(ipbeta>0.5) || startsWith(ConfigName,'kraken')
  lval=[-44 -44 -44 0 0]; uval=[0 44 0 44 20.3];
  optim='fminsearch';
else
  lval=[-44 -44 -44 -44 -20.3]; uval=[44 44 44 44 0];
  optim='lsqnonlin';
end
M=Match;
for ips=psno
  M.addVariable('PS',ips,'Ampl',lval(ips),uval(ips)); 
end
M.beam=MakeBeam6DGauss(I,1e3,3,1);
M.iInitial=i1;
M.initStruc=I;
M.verbose=false; % see optimizer output or not
M.optim=optim;
M.optimDisplay=odisp;
M.addMatch(ipele,'alpha_x',0,0.0001);
M.addMatch(ipele,'alpha_y',0,0.0001);
M.addMatch(ipele,'beta_x',ipbeta(1),0.0001);
M.addMatch(ipele,'beta_y',ipbeta(2),0.0001);
M.doMatch();
if strcmp(odisp,'iter')
  disp(M);
end

% If Kraken, then re-match waist at dump
if startsWith(ConfigName,"kraken")
  % Switch off spectrometer quads
  qd=[findcells(BEAMLINE,'Name','Q0D') findcells(BEAMLINE,'Name','Q1D') findcells(BEAMLINE,'Name','Q2D')];
  for iq=qd
    BEAMLINE{iq}.B=0;
  end
  M=Match;
  M.beam=MakeBeam6DGauss(I,1e3,3,1);
  M.iInitial=i1;
  M.initStruc=I;
  M.verbose=false; % see optimizer output or not
  M.optim='fminsearch';
  M.optimDisplay=odisp;
  M.addMatch(iscr,'beta_x',0,100);
  M.addMatch(iscr,'beta_y',0,100);
  M.addVariable('PS',1,'Ampl',lval(1),uval(1)); 
  M.addVariable('PS',2,'Ampl',lval(2),uval(2)); 
  M.doMatch();
  if strcmp(odisp,'iter')
    disp(M);
  end
end

% Match dump optics from IP
if ~startsWith(ConfigName,"kraken")
  if ConfigName == "sfqed"
    qd=findcells(BEAMLINE,'Name','Q2D');
    qf=findcells(BEAMLINE,'Name','Q1D');
    for iele=findcells(BEAMLINE,'Name','Q0D')
      BEAMLINE{iele}.B=0;
    end
  else
    qd=[findcells(BEAMLINE,'Name','Q0D') findcells(BEAMLINE,'Name','Q2D')];
    qf=findcells(BEAMLINE,'Name','Q1D');
  end
  for iele=qf;BEAMLINE{iele}.B=1; end
  for iele=qd;BEAMLINE{iele}.B=-1; end
  AssignToPS(qf,length(PS)+1); psqf=length(PS);
  AssignToPS(qd,length(PS)+1); psqd=length(PS);
  MovePhysicsVarsToPS([psqf psqd]);
  qmax=100;
  M=Match;
  M.beam=MakeBeam6DGauss(I,1e3,3,1);
  M.iInitial=i1;
  M.initStruc=I;
  M.verbose=false; % see optimizer output or not
  M.optim='lsqnonlin';
  M.optimDisplay=odisp;
  if ConfigName == "sfqed"
%     M.addMatch(iscr,'R',0,1e-14,'11');
%     M.addMatch(iscr,'R',0,1e-14,'33');
    M.addMatch(iscr,'beta_x',0,ipbeta(1));
    M.addMatch(iscr,'beta_y',0,ipbeta(2));
  else
    M.addMatch(iscr,'R',0,1e-9,'12');
    M.addMatch(iscr,'R',0,1e-9,'34');
  end
  M.addVariable('PS', psqf,'Ampl',0,qmax*2);
  M.addVariable('PS', psqd,'Ampl',-qmax*2,0);
  M.doMatch();
  if strcmp(odisp,'iter')
    disp(M);
  end
end

% Match sextupoles
if SextupoleMatch
  disp('Optimizing Sextupole Strengths...');
  ss = findcells(BEAMLINE,'Class','SEXT') ;
  AssignToPS(ss([1:2 15:16]),length(PS)+1); sps(1)=length(PS); MovePhysicsVarsToPS(sps(1));
  AssignToPS(ss([3:4 13:14]),length(PS)+1); sps(2)=length(PS); MovePhysicsVarsToPS(sps(2));
  AssignToPS(ss(5:12),length(PS)+1); sps(3)=length(PS); MovePhysicsVarsToPS(sps(2));
  for ips=sps
    PS(ips).Ampl=0; PS(ips).SetPt=0;
  end
  SetTrackFlags('ZMotion',1,1,length(BEAMLINE));
  if ConfigName == "sfqed"
    tind=iscr;
  else
    tind=ipele;
  end
  M=Match;
  M.beam=MakeBeam6DGauss(I,1e4,3,1);
  M.iInitial=i1;
  M.initStruc=I;
  M.verbose=false; % see optimizer output or not
  M.optim='fminsearch';
  M.optimDisplay='iter';
  M.addVariable('PS', sps(1),'Ampl',-2170.6,2170.6); % S1's
  M.addVariable('PS', sps(2),'Ampl',-776.8,776.8); % S2's
  M.addVariable('PS', sps(3),'Ampl',-776.8,776.8); % S3's
  M.addMatch(tind,'NEmit_x',0,I.x.NEmit);
  M.addMatch(tind,'NEmit_y',0,I.y.NEmit);
  M.doMatch();
  disp(M);
end

% Display matched magnet values and restore database to BEAMLINE
qm=[qm {'Q0D' 'Q1D' 'Q2D'}];
for ips=1:length(PS)
  for iele=PS(ips).Element
    if GetTrueStrength(iele)==0
      BEAMLINE{iele}.B=0;
    else
      RenormalizePS(ips);
    end
    BEAMLINE{iele}.PS=0;
  end
end
PS=[];
TwissPlot(i1,length(BEAMLINE),I,[1 1 1]);
clight=2.99792458e8; % speed of light (m/sec)
Cb=1e9/clight;       % rigidity conversion (T-m/GeV)
for ips=1:length(qm)
  iele=findcells(BEAMLINE,'Name',qm{ips});
  fprintf('K%s := %g\n',BEAMLINE{iele(1)}.Name,BEAMLINE{iele(1)}.B/(BEAMLINE{iele(1)}.L*Cb*BEAMLINE{iele(1)}.P));
end
bmax=[44 44 44 44 44 20.3 44 44 44].*10;
for ips=1:length(qm)
  iele=findcells(BEAMLINE,'Name',qm{ips});
  fprintf('BDES %s := %.1f (BMAX = %.1f) \n',BEAMLINE{iele(1)}.Name,10*sum(arrayfun(@(x) BEAMLINE{x}.B,iele)),bmax(ips));
end
if SextupoleMatch
  ss = findcells(BEAMLINE,'Class','SEXT') ;
  snames=unique(arrayfun(@(x) BEAMLINE{x}.Name,ss,'UniformOutput',false));
  for isext=1:length(snames)
    iele=findcells(BEAMLINE,'Name',snames{isext});
    fprintf('K%s := %g\n',BEAMLINE{iele(1)}.Name,BEAMLINE{iele(1)}.B/(BEAMLINE{iele(1)}.L*Cb*BEAMLINE{iele(1)}.P));
  end
end
% Show FFS magnets in format for import into FFS_magnets.xlsx
iv=IVB; % current lookup object
magnames={'Q5FF' 'Q4FF' 'Q3FF' 'Q2FF' 'Q1FF' 'Q0FF' 'Q0D' 'Q1D' 'Q2D'};
disp(magnames)
for ips=1:length(magnames)
  iele=findcells(BEAMLINE,'Name',magnames{ips});
  BDES = sum(arrayfun(@(x) BEAMLINE{x}.B,iele)) ;
  IDES=GetI(iv,BEAMLINE{iele(1)}.Type,BDES);
  fprintf('%g %g ',BDES*10,IDES);
end
fprintf('\n');
