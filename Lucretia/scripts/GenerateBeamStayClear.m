function StayClearTable = GenerateBeamStayClear(Initial,tablename)
% Generate beam stay clear definitions for FACET2 lattice
% Based on : FACET-II Technical Note #007: "Aperture Restrictions in Sector 20"

% Script operates on Twiss lattice using Lucretia deck loaded from facet2-lattice/Lucretia/models/
% cell array "BEAMLINE" must be in global memory and Lucretia "Initial" structure variable supplied

% Beam stay clear x, and y half-widths and radii are calculated corresponding to the start of the Beamline element in question and quoted in meters

global BEAMLINE

[~,T]=GetTwiss(1,length(BEAMLINE),Initial.x.Twiss,Initial.y.Twiss);

% Stay-clear defined as Naper X beta size + NaperE * rms dispersive size
Nbeta=10; % # sigma apertures to define (Naper * rms beam size determined from Twiss parameters)
minclear=1/sqrt(2)*2.25*1e-3; % Minimum clearance to use (in x and y)
nemit=5e-6; % Max expected operating normalized transverse emittance
% Include TCAV streaked images for stay-clear in relevant sections
tcavname={'TCY10490' 'TCY15280' 'XTCAVF'};
tcavendname={'ENDBC11_2' 'ENDL3F_1' 'MAINDUMP'}; % Considered end section for each streak diagnostic
tcavsz=[2.5e-3 2e-3 200e-6]; % Largest longitudinal offset to consider for streaked beam (from 2-bunch compression profile)
tcavV=[3 10 20]; % max tcav operating voltage / MV
tcavF=[2.8560e+09 2.8560e+09 11.424e9]; % TCAV frequencies / Hz
tcavdim=[2 1 1]; % 1=horizontal streaking, 2=vertical streaking
tcavind=zeros(1,length(tcavname)); tcavendind=tcavind;
for itcav=1:length(tcavname)
  ind = findcells(BEAMLINE,'Name',tcavname{itcav});
  tcavind(itcav)=ind(end) ;
  tcavendind(itcav)=findcells(BEAMLINE,'Name',tcavendname{itcav});
end
c=299792458; % speed of light in vacuum [m/s]

% Define max energy offset based on 2-bunch profile configuration
de=[0.0014    0.0084    0.2700    0.3000 0.3].*1.1; % Max energy offset (GeV) defined at exit of [injector,L1,L2,L3] - use linear interpolation for intermediate locations
dename={'BEGL1F' 'ENDL1F' 'ENDL2F' 'ENDL3F_2' 'MAINDUMP'};
deind=zeros(1,length(dename));
for iele=1:length(dename)
  deind(iele) = findcells(BEAMLINE,'Name',dename{iele});
end

% Calculate the stay-clear data
i1=findcells(BEAMLINE,'Name','ENDINJ'); % Only fill in beam stay clear beyong this point
S=arrayfun(@(x) BEAMLINE{x}.S,1:length(BEAMLINE));
for iele=1:length(BEAMLINE)
  StayClear_x(iele).beta=0; StayClear_y(iele).beta=0; StayClear_r(iele).beta=0 ;
  StayClear_x(iele).eta=0; StayClear_y(iele).eta=0; StayClear_r(iele).eta=0 ;
  StayClear_x(iele).tcav=0; StayClear_y(iele).tcav=0; StayClear_r(iele).tcav=0 ;
  StayClear_x(iele).all=0; StayClear_y(iele).all=0; StayClear_r(iele).all=0 ;
end
for iele=1:length(BEAMLINE)
  if iele<=i1
    continue
  end
  if iele<=deind(1)
    derel = de(1) / BEAMLINE{iele}.P ;
  else
    for iele2=2:length(deind)
      if iele>deind(iele2-1) && iele<=deind(iele2)
        de1 = interp1([BEAMLINE{deind(iele2-1)}.S BEAMLINE{deind(iele2)}.S] , de(iele2-1:iele2),BEAMLINE{iele}.S) ;
        derel = de1 / BEAMLINE{iele}.P ;
      end
    end
  end
  gamma=BEAMLINE{iele}.P/0.511e-3;
  emit=nemit/gamma;
  StayClear_x(iele).beta = sqrt(emit*T.betax(iele) ) * Nbeta ; %#ok<*SAGROW>
  StayClear_y(iele).beta = sqrt(emit*T.betay(iele) ) * Nbeta ;
  StayClear_r(iele).beta = sqrt( StayClear_x(iele).beta^2 + StayClear_y(iele).beta^2 ) ;
  
  StayClear_x(iele).eta = abs(derel*T.etax(iele)) ;
  StayClear_y(iele).eta = abs(derel*T.etay(iele)) ;
  StayClear_r(iele).eta = sqrt( StayClear_x(iele).eta^2 + StayClear_y(iele).eta^2 ) ;
  
  StayClear_x(iele).tcav = 0 ;
  StayClear_y(iele).tcav = 0 ;
  StayClear_r(iele).tcav = 0 ;
  
  for itcav=1:length(tcavname)
    if iele>tcavind(itcav) && iele<=tcavendind(itcav)
      [~,R]=RmatAtoB(tcavind(itcav),iele-1);
      E0=BEAMLINE{iele}.P*1e3; % Beam energy / MeV
      kp=(2*pi*tcavF(itcav))/c/E0;
      Sx=tcavV(itcav)*kp;
      if tcavdim(itcav)==1
        StayClear_x(iele).tcav = tcavsz(itcav)*Sx*abs(R(1,2)) ;
        StayClear_r(iele).tcav = StayClear_x(iele).tcav ;
      else
        StayClear_y(iele).tcav = tcavsz(itcav)*Sx*abs(R(3,4)) ;
        StayClear_r(iele).tcav = StayClear_y(iele).tcav ;
      end
    end
  end
  
  StayClear_x(iele).all = max( [StayClear_x(iele).beta StayClear_x(iele).eta StayClear_x(iele).tcav ] ) ;
  StayClear_y(iele).all = max( [StayClear_y(iele).beta StayClear_y(iele).eta StayClear_y(iele).tcav ] ) ;
  StayClear_x(iele).all = max([StayClear_x(iele).all minclear]) ;
  StayClear_y(iele).all = max([StayClear_y(iele).all minclear]) ;
  StayClear_r(iele).all = sqrt( StayClear_x(iele).all.^2 + StayClear_y(iele).all^2 ) ;
  
end
close all
plot(S,[StayClear_r(:).beta].*1e3,S,[StayClear_r(:).eta].*1e3,S,[StayClear_r(:).tcav].*1e3,S,[StayClear_r(:).all].*1e3,'--');
legend({'Betatron stay-clear' 'Dispersive stay-clear' 'TCAV streak stay-clear' 'Overall stay-clear'});
xlabel('S [m]'); ylabel('Stay-Clear Radius [mm]');
% -- add apertures
q2ff = findcells(BEAMLINE,'Name','Q2FF');
q1ff = findcells(BEAMLINE,'Name','Q1FF');
q0ff = findcells(BEAMLINE,'Name','Q0FF');
q0d = findcells(BEAMLINE,'Name','Q0D');
q1d = findcells(BEAMLINE,'Name','Q1D');
q2d = findcells(BEAMLINE,'Name','Q2D');
usStrawBeg = findcells(BEAMLINE,'Name','LCUBE');
dsStrawEnd = findcells(BEAMLINE,'Name','IM3255');
S_Ap = [BEAMLINE{q2ff(1)}.S+[-.1 .6] 0 BEAMLINE{q1ff(1)}.S+[-.1 .6]  0 BEAMLINE{q0ff(1)}.S+[-.1 .6] 0 ...
        BEAMLINE{usStrawBeg}.S+[0 .1] 0 BEAMLINE{dsStrawEnd}.S+[-.1 0] 0 ...
        BEAMLINE{q0d(1)}.S+[0 1] 0 BEAMLINE{q1d(1)}.S+[0 1] 0 BEAMLINE{q2d(1)}.S+[0 1]];
 
Ap = [9 9 NaN 9 9 NaN 9 9 NaN ...
      2.5 2.5 NaN 2.5 2.5 NaN ...
      24 24 NaN 24 24 NaN 24 24]; % mm
hold on
plot(S_Ap,Ap, 'k-','LineWidth',5);
hold off
AddMagnetPlot(1,length(BEAMLINE));
vnames={'BeamlineIndex' 'S' 'x_all' 'x_beta' 'x_eta' 'x_tcav' 'y_all' 'y_beta' 'y_eta' 'y_tcav' 'r_all' 'r_beta' 'r_eta' 'r_tcav'};
StayClearTable=table((1:length(BEAMLINE))',S(:),...
  [StayClear_x(:).all]',[StayClear_x(:).beta]',[StayClear_x(:).eta]',[StayClear_x(:).tcav]',...
  [StayClear_y(:).all]',[StayClear_y(:).beta]',[StayClear_y(:).eta]',[StayClear_y(:).tcav]',...
  [StayClear_r(:).all]',[StayClear_r(:).beta]',[StayClear_r(:).eta]',[StayClear_r(:).tcav]','VariableNames',vnames);
StayClearTable.Properties.Description=char(tablename);
writetable(StayClearTable,sprintf('BeamStayClear_%s.xlsx',tablename));
save(sprintf('BeamStayClear_%s',tablename),'StayClearTable');
