function runData=RunImpactT(npart,GunPhase,Sol1,Sol2,dGunDrift,laserR,laserZ,laserModFreq)
% RunImpactT(npart,GunPhase,Sol1,Sol2,dGunDrift,laserR,laserZ)
% Setup ImpactT run file, assuming existence of 'ImpactT_master.in'. Then
% run the program.
% Input arguments optional, defaults taken from ImpactT_master.in file
%   Set unwanted arguments to []
% npart: number of macro-particles. Default is 20k.
% GunPhase: phase of RF gun / degrees relative to charge 0-crossing
%   (Setting this variable causes L0a phase to be set to maximize energy
%   post-L0a)
% Sol1: Scale factor to apply to Solenoid 1 field
% Sol2: Scale factor to apply to Solenoid 2 field
% dGunDrift: Change to Gun exit->L0a drift distance (m)
%   (Setting this variable causes L0a phase to be set to maximize energy
%   post-L0a) 
% LaserR: spot size radius of laser on cathode / m
% LaserZ: length of flat-top of laser pulse / s
% Initialize arguments
if ~exist('npart','var')
  npart=2e4;
end
if ~exist('GunPhase','var') || isempty(GunPhase)
  GunPhase=30;
end
if ~exist('Sol1','var')
  Sol1=[];
end
if ~exist('Sol2','var')
  Sol2=[];
end
if ~exist('dGunDrift','var')
  dGunDrift=[];
end
if ~exist('laserR','var') || isempty(laserR)
  laserR=2.68e-3;
end
if ~exist('laserZ','var') || isempty(laserZ)
  laserZ=7e-12;
end

% Random number seed
rng('shuffle');

% Scale laser radius to compensate for change in gun phase, assuming charge
% scales linearly as 1/30 e9 / degree with 1e9 charge @ 30 degrees
% laserR=laserR*sqrt(1/(GunPhase/30));
% GunPhase=GunPhase+270;
GunPhase=(GunPhase-5.6)+270;


% Specify which lines elements of interest are to be found
% Commented lines are ignored (deck starts on line 10)
mline.Gun=10;
mline.Sol1=11;
mline.Sol2=13;
mline.GunDrift=12;
mline.L0a=15:18;
mline.laserR=6:7;
mline.laserZ=8;% if ~isempty(laserR)
%   for iline=1:2
%     vals=textscan(Iline{mline.laserR(iline)},'%s'); vals=vals{1};
%     vals{1}=num2str(laserR,6); vals{2}=num2str(laserR,6);
%     Iline{mline.laserR(iline)}=writeImpactLine(vals);
%   end
% end

mline.npart=3;

% - Read in master ImpactT file
fid=fopen('ImpactT_master.in');
if fid<0
  error('ImpactT_master.in file not found!');
end
lc=0;
Iline={};
while 1
  tline=fgetl(fid);
  if ~ischar(tline), break, end
  if ~strcmp(tline(1),'!')
    lc=lc+1;
    Iline{lc}=tline;
  end
end
fclose(fid);

% Set Gun Phase
if ~isempty(GunPhase)
  vals=textscan(Iline{mline.Gun},'%s'); vals=vals{1};
  vals{8}=num2str(GunPhase);
  Iline{mline.Gun}=writeImpactLine(vals);
end

% Set Sol1 strength
if ~isempty(Sol1)
  vals=textscan(Iline{mline.Sol1},'%s'); vals=vals{1};
  vals{16}=num2str(Sol1);
  Iline{mline.Sol1}=writeImpactLine(vals);
end

% Set Sol2 strength
if ~isempty(Sol2)
  vals=textscan(Iline{mline.Sol2},'%s'); vals=vals{1};
  vals{16}=num2str(Sol2);
  Iline{mline.Sol2}=writeImpactLine(vals);
end

% Set laser pulse radius
% if ~isempty(laserR)
%   for iline=1:2
%     vals=textscan(Iline{mline.laserR(iline)},'%s'); vals=vals{1};
%     vals{1}=num2str(laserR,6); vals{2}=num2str(laserR,6);
%     Iline{mline.laserR(iline)}=writeImpactLine(vals);
%   end
% end

% Set laser flat top length
% if ~isempty(laserZ)
%   vals=textscan(Iline{mline.laserZ},'%s'); vals=vals{1};
%   vals{5}=num2str(laserZ,6);
%   Iline{mline.laserZ}=writeImpactLine(vals);
% end

% Set Gun -> L0a drift distance
if ~isempty(dGunDrift)
  vals=textscan(Iline{mline.GunDrift},'%s'); vals=vals{1};
  z1=str2double(vals{5})+str2double(vals{1}); % initial d/s location of drift
  vals{1}=num2str(str2double(vals{1})+dGunDrift);
  % - move z location of everything that starts after gun drift
  for ival=10:length(Iline)
    vals2=textscan(Iline{ival},'%s'); vals2=vals2{1};
    if str2double(vals2{5})>=z1
      vals2{5}=num2str(str2double(vals2{5})+dGunDrift);
      Iline{ival}=writeImpactLine(vals2);
    end
  end
  Iline{mline.GunDrift}=writeImpactLine(vals);
end

% If setting new Gun Phase or gun drift length, scan L0a phase to maximize energy gain at L0a
% exit using single particle tracking in ImpactT and then set L0a phase in
% generated input file
if ~isempty(GunPhase) || ~isempty(dGunDrift)
  % - find L0a phase which maximizes energy by tracking in 1 bunch mode
  minphase=fminbnd(@(x) ScanL0aFun(x,Iline,mline),-180,180,optimset('TolX',0.1,'TolFun',0.01,'Display','iter'));
  % - Set new L0a phase
  for irow=mline.L0a
    vals=textscan(Iline{irow},'%s'); vals=vals{1};
    vals{8}=num2str(str2double(vals{8})+minphase);
    Iline{irow}=writeImpactLine(vals);
  end
end

% Set # macro particles
if exist('laserModFreq','var') && ~isempty(laserModFreq)
  MakeImpactBeam(npart,laserR,laserZ,1,laserModFreq);
else
  MakeImpactBeam(npart,laserR,laserZ,1);
end
vals=textscan(Iline{mline.npart},'%s'); vals=vals{1};
vals{2}=num2str(npart);
Iline{mline.npart}=writeImpactLine(vals);

% Write ImpactT run file
fid=fopen('ImpactT.in','w');
if fid<0
  error('Cannot open ImpactT.in file for writing!')
end
for il=1:length(Iline)
  fprintf(fid,'%s\n',Iline{il});
end
fclose(fid);

% Run Program and return data
stxt=evalc('sret=system(''./impact'')');if sret;error(stxt);end;
idata=ReadImpactTData('.');

% Use Lucretia to configure L0b and match into DL1 and to track beam to DL1
clear global BEAMLINE KLYSTRON PS
ldata=TrackInj(idata,0);

% Package output data
runData.idata=idata;
runData.ldata=ldata;
if ~isempty(minphase)
  runData.L0aPhase=minphase;
end

function SetImpactTSingleBunchFile(Iline1,mline,dL0aPhase)

% -- set single macro particle
MakeImpactBeam('single');
vals=textscan(Iline1{3},'%s'); vals=vals{1};
vals{2}='1';
Iline1{3}=writeImpactLine(vals);
% -- Set appropriate dist type
% vals=textscan(Iline1{5},'%s'); vals=vals{1};
% vals{1}='17';
% Iline1{5}=writeImpactLine(vals);
% -- Set zero beam size
% for irow=6:8
%   vals=textscan(Iline1{irow},'%s'); vals=vals{1};
%   for icol=1:3
%     vals{icol}='0.0';
%   end
%   Iline1{irow}=writeImpactLine(vals);
% end
% -- Set zero current
vals=textscan(Iline1{9},'%s'); vals=vals{1};
vals{1}='0';
Iline1{9}=writeImpactLine(vals);
% -- Set L0a phase
for irow=mline.L0a
  vals=textscan(Iline1{irow},'%s'); vals=vals{1};
  vals{8}=num2str(str2double(vals{8})+dL0aPhase);
  Iline1{irow}=writeImpactLine(vals);
end
% Write 1-bunch impact file
fid=fopen('ImpactT.in','w');
if fid<0
  error('Cannot open ImpactT.in file for writing!')
end
for il=1:length(Iline1)
  fprintf(fid,'%s\n',Iline1{il});
end
fclose(fid);

function optval=ScanL0aFun(x,Iline,mline)
% _ generate ImpactT input file for given L0a phase change x
SetImpactTSingleBunchFile(Iline,mline,x);
% - run ImpactT
stxt=evalc('sret=system(''./impact'')');if sret;error(stxt);end;
% - gather output
data=ReadImpactTData('.');
% Min function is final energy
optval=-data.E.E(end);

function tline=writeImpactLine(txtcell)
% format cell array of strings for ImpactT.in line output
tline=txtcell{1};
if length(txtcell)>1
  for iline=2:length(txtcell)
    tline=[tline sprintf(' %s',txtcell{iline})];
  end
end
