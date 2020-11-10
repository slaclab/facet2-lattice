% Generate master Beam Stay clear file
%  Use max stay clear definition from all BeamStayClear_* files in local directory

datfile=dir("BeamStayClear_*.mat");

tab={};
for idat=1:length(datfile)
  id=string(regexp(regexprep(datfile(idat).name,"(\.mat)",""),"BeamStayClear_(\S+)","tokens","once"));
  if id~="master"
    dat=load(datfile(idat).name,'StayClearTable');
    tab{end+1}=dat; %#ok<*SAGROW>
  end
end
mastertab=tab{1}.StayClearTable;
if length(tab)<=1
  return
end
for itab=2:length(tab)
  masterdat=mastertab{:,3:end};
  newdat=tab{itab}.StayClearTable{:,3:end};
  seldat = newdat>masterdat ;
  if any(seldat(:))
    masterdat(seldat) = newdat(seldat) ;
    mastertab{:,3:end} = masterdat ;
  end
end

StayClearTable=mastertab;
StayClearTable.Properties.Description='MasterBeamClearDefinition';
writetable(StayClearTable,'BeamStayClear_master.xlsx');
save('BeamStayClear_master','StayClearTable');

figure
plot(mastertab{:,'S'},mastertab{:,'r_beta'}.*1e3,mastertab{:,'S'},mastertab{:,'r_eta'}.*1e3,...
  mastertab{:,'S'},mastertab{:,'r_tcav'}.*1e3,mastertab{:,'S'},mastertab{:,'r_all'}.*1e3,'--');
legend({'Betatron stay-clear' 'Dispersive stay-clear' 'TCAV streak stay-clear' 'Overall stay-clear'});
xlabel('S [m]'); ylabel('Stay-Clear Radius [mm]');
AddMagnetPlot(1,length(BEAMLINE));