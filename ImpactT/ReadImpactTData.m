function data=ReadImpactTData(dirname)
% Read in data from an executed ImpactT run directory

% Open required data files for reading
datfiles=openFiles(dirname);

try
  
  % Exctract data from master deck file
  Q=parseInputDeck(datfiles.run);
  
  % Get energy profile
  edat=parseEData(datfiles.edat);
  
  % Get X/Y/Z data
  pdat=parsePosData(datfiles.xdat,datfiles.ydat,datfiles.zdat);
  
  % Get initial and final bunch data
  try
    bdat=parseBunchData(datfiles.ibunch,datfiles.fbunch,Q);
  catch
    bdat=[];
  end
  
catch ME % Close open files before error
  
  fn=fieldnames(datfiles);
  for ifn=1:length(fn)
    fclose(datfiles.(fn{ifn}));
  end
  rethrow(ME);
  
end

% Package data
data.Q=Q; data.E=edat; data.POS=pdat; data.Beam=bdat;

% Close files
fn=fieldnames(datfiles);
for ifn=1:length(fn)
  fclose(datfiles.(fn{ifn}));
end

function dat=parseBunchData(fidI,fidF,Q)
di=textscan(fidI,'%f'); df=textscan(fidF,'%f');
Bi=CreateBlankBeam(1,length(di{1})/6,1,0);
Bi.Bunch.Q=ones(size(Bi.Bunch.Q)).*(Q/length(Bi.Bunch.Q));
Bi.Bunch.x(1,:)=di{1}(1:6:end);
Bi.Bunch.x(2,:)=atan(di{1}(2:6:end)./di{1}(6:6:end));
Bi.Bunch.x(3,:)=di{1}(3:6:end);
Bi.Bunch.x(4,:)=atan(di{1}(4:6:end)./di{1}(6:6:end));
Bi.Bunch.x(5,:)=di{1}(5:6:end);
Bi.Bunch.x(6,:)=sqrt(0.511e-3^2*(1+(di{1}(2:6:end).^2+di{1}(4:6:end).^2+di{1}(6:6:end).^2)));
Bf=CreateBlankBeam(1,length(df{1})/6,1,0);
Bf.Bunch.Q=ones(size(Bf.Bunch.Q)).*(Q/length(Bf.Bunch.Q));
Bf.Bunch.x(1,:)=df{1}(1:6:end);
Bf.Bunch.x(2,:)=atan(df{1}(2:6:end)./df{1}(6:6:end));
Bf.Bunch.x(3,:)=df{1}(3:6:end);
Bf.Bunch.x(4,:)=atan(df{1}(4:6:end)./df{1}(6:6:end));
Bf.Bunch.x(5,:)=df{1}(5:6:end); Bf.Bunch.x(5,:)=-(Bf.Bunch.x(5,:)-median(Bf.Bunch.x(5,:)));
Bf.Bunch.x(6,:)=sqrt(0.511e-3^2*(1+(df{1}(2:6:end).^2+df{1}(4:6:end).^2+df{1}(6:6:end).^2)));
dat.Binit=Bi; dat.Bfin=Bf;

function dat=parsePosData(fidX,fidY,fidZ)
dX=textscan(fidX,'%f'); dY=textscan(fidY,'%f'); dZ=textscan(fidZ,'%f');
dat.Z=dZ{1}(2:7:end);
dat.X=dX{1}(3:8:end);
dat.Y=dY{1}(3:8:end);
dat.rmsX=dX{1}(4:8:end);
dat.rmsY=dY{1}(4:8:end);
dat.rmsZ=dZ{1}(3:7:end);
dat.emitX=dX{1}(8:8:end);
dat.emitY=dY{1}(8:8:end);
dat.emitZ=dZ{1}(7:7:end);

function dat=parseEData(fid)
d=textscan(fid,'%f');
dat.S=d{1}(2:7:end);
dat.E=d{1}(4:7:end);
dat.dE=d{1}(7:7:end);

function Q=parseInputDeck(fid)
lc=0;
while 1
  tline=fgetl(fid);
  if ~ischar(tline), break, end
  if ~strcmp(tline(1),'!')
    lc=lc+1;
  end
  if lc==9
    % charge = KE * current / scaleFreq
    % line 9 format: current, kinetic energy, particle rest energy, particle charge, scale frequency, initial  cavity phase !
    vals=textscan(tline,'%f'); vals=vals{1};
    Q=(vals(1)*vals(2))/vals(5); % Charge / nC
    break
  end
end

function [datfiles]=openFiles(dirname)
% Open required data files for reading
% - File containing run details
runfile=fopen(fullfile(dirname,'ImpactT.in'));
if runfile<0; error('Cannot open ImpactT.in file'); end;
% - Time/distance/E etc
datfile=fopen(fullfile(dirname,'fort.18'));
if datfile<0
  fclose(runfile);
  error('Cannot open fort.18 file')
end
% - X rms information
xdatfile=fopen(fullfile(dirname,'fort.24'));
if xdatfile<0
  fclose(runfile); fclose(datfile);
  error('Cannot open fort.24 file')
end
% - Y rms information
ydatfile=fopen(fullfile(dirname,'fort.25'));
if ydatfile<0
  fclose(runfile); fclose(datfile); fclose(xdatfile);
  error('Cannot open fort.25 file')
end
% - Z rms information
zdatfile=fopen(fullfile(dirname,'fort.26'));
if zdatfile<0
  fclose(runfile); fclose(datfile); fclose(xdatfile); fclose(ydatfile);
  error('Cannot open fort.26 file')
end
% - Initial bunch
ibunchfile=fopen(fullfile(dirname,'fort.40'));
if ibunchfile<0
  fclose(runfile); fclose(datfile); fclose(xdatfile); fclose(ydatfile); fclose(zdatfile);
  error('Cannot open fort.40 file')
end
% - Final bunch
fbunchfile=fopen(fullfile(dirname,'fort.50'));
if fbunchfile<0
  fclose(runfile); fclose(datfile); fclose(xdatfile); fclose(ydatfile); fclose(zdatfile); fclose(ibunchfile);
  error('Cannot open fort.50 file')
end
datfiles.run=runfile; datfiles.edat=datfile; datfiles.xdat=xdatfile; datfiles.ydat=ydatfile;
datfiles.zdat=zdatfile; datfiles.ibunch=ibunchfile; datfiles.fbunch=fbunchfile;