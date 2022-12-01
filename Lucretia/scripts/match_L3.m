function T=match_L3(Initial)
% Match through BC14 into L3
global BEAMLINE PS

% Keep Q15501 & Q16601 closer to other magnets in their strings
qm={'QM14891' 'Q14901' 'Q15201' 'Q15301' 'Q15401' 'Q15501'  'Q15601'};
qlim={[1 3]    [-4 -2] [0.8 1.5]   [-1.5 -0.8]  [0.8 1.5]  [-1.5 -0.8]  [0.8 1.5]} ;
i1=findcells(BEAMLINE,'Name',qm{1}); i1=i1(1);
i2=findcells(BEAMLINE,'Name',qm{end}); i2=i2(2);
[~,T]=GetTwiss(1,i2,Initial.x.Twiss,Initial.y.Twiss); I=TwissToInitial(T,i1,Initial);
TwissPlot(i1,i2,I,[1 1 0],0.01);

for iquad=1:length(qm)
  iele=findcells(BEAMLINE,'Name',qm{iquad});
  AssignToPS( iele, length(PS)+1 ) ;
end
MovePhysicsVarsToPS(1:length(PS));
M=Match;
M.beam=MakeBeam6DGauss(I,1e3,3,1);
M.iInitial=i1;
M.initStruc=I;
M.verbose=false; % see optimizer output or not
M.optim='lsqnonlin';
M.optimDisplay='iter';
M.addMatch(i2,'beta_x',T.betax(end),1e-4);
M.addMatch(i2,'beta_y',T.betay(end),1e-4);
M.addMatch(i2,'alpha_x',T.alphax(end),1e-4);
M.addMatch(i2,'alpha_y',T.alphay(end),1e-4);
% M.addMatch(PS(11).Element(1),'beta_x',0,100);
% M.addMatch(PS(9).Element(1),'beta_y',0,100);
% M.addMatch(PS(5).Element(1),'beta_y',0,100);
% M.addMatch(PS(6).Element(1),'beta_y',0,100);
qmatch=1:length(qm);
for ips=qmatch
  M.addVariable('PS',ips,'Ampl',qlim{ips}(1),qlim{ips}(2));
end
M.doMatch();
disp(M);
TwissPlot(i1,i2,I,[1 1 0],0.01);

for ips=1:length(qm)
  iele=findcells(BEAMLINE,'Name',qm{ips});
  BDES = sum(arrayfun(@(x) BEAMLINE{x}.B*PS(ips).Ampl,iele)) ;
  fprintf('%s: BDES= %g\n',qm{ips},BDES*10);
end

% restore DB
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
T=TwissPlot(1,length(BEAMLINE),Initial,[1 1 0],0.01);

clight=299792458; Brho=1e9/clight;
for ips=1:length(qm)
  iele=findcells(BEAMLINE,'Name',qm{ips}); iele=iele(1);
  KDES = BEAMLINE{iele}.B/(Brho*BEAMLINE{iele}.P)/BEAMLINE{iele}.L ;
  fprintf('%s := %g\n',"K"+string(qm{ips}),KDES);
end