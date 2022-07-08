function T=match_L2(Initial)
% Match through BC11 into L2
global BEAMLINE PS

% Alternative scheme with Q11401 OFF
% qm={'QA11132' 'Q11201' 'QA11265' 'Q11301' 'QM11312' 'QM11358'  'QM11362' 'QM11393' 'Q11401' 'Q11501' 'Q11601' 'Q11701' 'Q11801' 'Q11901' 'Q12201'};
qm={'QA11132' 'Q11201' 'QA11265' 'Q11301' 'QM11312' 'QM11358'  'QM11362' 'QM11393' 'Q11401' 'Q11501' 'Q11601' 'Q11701' 'Q11801'};
i1=findcells(BEAMLINE,'Name',qm{1}); i1=i1(1);
i2=findcells(BEAMLINE,'Name',qm{end}); i2=i2(2);
[~,T]=GetTwiss(1,i2,Initial.x.Twiss,Initial.y.Twiss); I=TwissToInitial(T,i1,Initial);
TwissPlot(i1,i2,I,[1 1 0],0.01);

for iquad=1:length(qm)
  iele=findcells(BEAMLINE,'Name',qm{iquad});
  AssignToPS( iele, length(PS)+1 ) ;
end
MovePhysicsVarsToPS(1:length(PS));
qOFF=[4 5 6 9];
% qOFF=[2 4 5 6 8];
for iq=qOFF
  PS(iq).Ampl=0; PS(iq).SetPt=0;
end
qFLIP=[3 8];
for iq=qFLIP
  PS(iq).Ampl=-PS(iq).Ampl; PS(iq).SetPt=-PS(iq).SetPt;
end
PS(3).Ampl=-0.25;
PS(4).Ampl=0.25;
PS(7).Ampl=-0.15;
PS(8).Ampl=0.2;
M=Match;
M.beam=MakeBeam6DGauss(I,1e3,3,1);
M.iInitial=i1;
M.initStruc=I;
M.verbose=false; % see optimizer output or not
M.optim='fminsearch';
% M.optim='fmincon';
M.optimDisplay='iter';
M.addMatch(i2,'beta_x',T.betax(end),1e-4);
M.addMatch(i2,'beta_y',T.betay(end),1e-4);
M.addMatch(i2,'alpha_x',T.alphax(end),1e-4);
M.addMatch(i2,'alpha_y',T.alphay(end),1e-4);
% M.addMatch(PS(11).Element(1),'beta_x',0,100);
% M.addMatch(PS(9).Element(1),'beta_y',0,100);
% M.addMatch(PS(5).Element(1),'beta_y',0,100);
% M.addMatch(PS(6).Element(1),'beta_y',0,100);
qmatch=[10:length(qm)];
for ips=qmatch
  if PS(ips).Ampl>0
    M.addVariable('PS',ips,'Ampl',0,200);
  else
    M.addVariable('PS',ips,'Ampl',-200,0);
  end
%   M.addVariable('PS',ips,'Ampl',-200,200);
end
M.doMatch();
disp(M);
% TwissPlot(i1,PS(10).Element(end),I,[1 1 0],0.01);
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