% Concatinate multiple injector beam definitions into single bunch

dirname='data/inj/errFiles';
fname='BeamErr';

files=dir(fullfile(dirname,sprintf('%s*',fname)));

for ifile=1:length(files)
  d=load(fullfile(files(ifile).folder,files(ifile).name),'ImpactData');
  B=d.ImpactData.ldata.Beam;
  coef=polyfit(B.Bunch.x(6,~B.Bunch.stop),B.Bunch.x(2,~B.Bunch.stop),1);
  coef(2)=0;
  B.Bunch.x(2,~B.Bunch.stop)=B.Bunch.x(2,~B.Bunch.stop)-polyval(coef,B.Bunch.x(6,~B.Bunch.stop));
  for ii=1:5
    if ii<5
      coef=polyfit(B.Bunch.x(6,~B.Bunch.stop),B.Bunch.x(ii,~B.Bunch.stop),1);
      coef(2)=0;
      B.Bunch.x(ii,~B.Bunch.stop)=B.Bunch.x(ii,~B.Bunch.stop)-polyval(coef,B.Bunch.x(6,~B.Bunch.stop));
    end
    B.Bunch.x(ii,~B.Bunch.stop)=B.Bunch.x(ii,~B.Bunch.stop)-mean(B.Bunch.x(ii,~B.Bunch.stop));
  end
  % Put on-energy
  B.Bunch.x(6,:)=(B.Bunch.x(6,:)-mean(B.Bunch.x(6,:)))+0.135;
  if ifile>1
    Beam.Bunch.x=[Beam.Bunch.x B.Bunch.x];
    Beam.Bunch.stop=[Beam.Bunch.stop B.Bunch.stop];
    Beam.Bunch.Q=[Beam.Bunch.Q B.Bunch.Q];
  else
    Beam=B;
  end
  fprintf('Processed %d / %d files.\n',ifile,length(files))
end
Beam.Bunch.Q=Beam.Bunch.Q./length(files);
fprintf('Final Beam size= %.2g macro-particles\n',length(Beam.Bunch.Q))
beamImage(Beam)