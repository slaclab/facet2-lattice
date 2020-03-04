function Bout=GPTtrackExample()
% Example function to take pre-generated Lucretia beam and 3D field map
% and track through GPT and recover the rays at discrete locations through
% the map
%
% Output is cell array containing particle coordinates at screen
% positions set in GPTtrack sub-function
%
% Required files are:
%  LucretiaTestBeam.mat (Matlab file with example Lucretia Beam)
%  BfieldExample.gdf (example 3D field map in GPT input format)
%  save_struct_to_gdf_file.m (to convert matlab data to gdf format)
%  load_gdf.m (to convert from gdf format to Matlab data)

% - Temp directory for running GPT - this is fine for linux or mac, change
% to something more windows-like if Windows OS required
% Code generates a temporary subdirectory here for running GPT and then
% deletes afterwards- this is useful if you want to take this example and
% use in a parallized routine
tmpdir='/tmp';

% Load a pre-defined Lucretia beam
% Beam.Bunch.x is a [6xN] array of ray coordinates: [x,xp,y,yp,-ct,P]
% (note P is absolute momentum in GeV NOT dp/p), other coordinates are in
% [m, rad] units
% This beam has 1,000 macro-particles defined
ld = load('LucretiaTestBeam.mat','Beam') ;
Beam=ld.Beam;

% Track this beam through BfieldExample.gdf field (change locations of
% screens in GPTtrack function to alter where data is read out)
Bout = GPTtrack(Beam,tmpdir) ;

% As an example, plot rms horizontal beam size
% and position at each screen location
for iscreen=1:length(Bout)
  sigma_x(iscreen) = std(Bout{iscreen}.Bunch.x(1,:));
  pos_x(iscreen) = mean(Bout{iscreen}.Bunch.x(1,:));
end
figure
subplot(2,1,1),plot(1:length(Bout),sigma_x);
xlabel('Screen Number'); ylabel('\sigma_x [m]');
subplot(2,1,2),plot(1:length(Bout),pos_x);
xlabel('Screen Number'); ylabel('x [m]');

end

function Bout=GPTtrack(BeamIn,tmpdir)
% Track provided Lucretia bunch through 3D field map. Output at pre-stored s
% locations by placement of screens in GPT code.
% Output argument Bout is cell array of Lucretia bunches. Each element
% of cell array corresponds to pre-defined s locations.
% ===================================================
% - some constants
clight=299792458; % speed of light in vacuum / m/s
qe=1.60217653e-19; % electron charge / C
me=0.511e-3; % Rest Mass of electron / GeV
% - simulation parameters
BfieldFileName='BfieldExample.gdf'; % example B field map
nscreen=10; % Number of screen elements to place in field map
xscr=zeros(1,10); % x position of screens
xpscr=zeros(1,10); % horizontal angle of screens
yscr=zeros(1,10); % y position of screens
ypscr=zeros(1,10); % vertical angle of screens
zscr=linspace(0,0.5,10); % z position of screens
t0=0; % GPT sim start time
t1=1.4e-9; % GOT tracking end time
tprec=1e-10; % GPT tracking time interval
% - temp folder for running GPT
gptdir=fullfile(tmpdir,sprintf('gpt_%s',num2str(sum(clock).*1e6,12)));
mkdir(gptdir);
Bout={}; % initialize output bunch cell array
% - Write Luctia bunch particles out in GPT input file
xGPT=BeamIn.Bunch.x;
gamma=BeamIn.Bunch.x(6,:)./0.511e-3;
beta=sqrt(1-gamma.^-2);
beta_x=beta.*sin(xGPT(2,:)); beta_y=beta.*sin(xGPT(4,:)); beta_z=sqrt(1-(beta_x./beta).^2-(beta_y./beta).^2).*beta;
data=struct; data.d.x=xGPT(1,:); data.d.y=xGPT(3,:); % GPT data structure for pos coordinates
% Convert z distribution into particle release times
data.d.z=xGPT(5,:);
v=beta*clight;
data.d.t=xGPT(5,:)./v;
data.d.t=data.d.t;
data.d.Bx=beta_x'; data.d.By=beta_y'; data.d.Bz=beta_z'; data.d.nmacro=BeamIn.Bunch.Q./qe';
data.d.G=1./sqrt(1-(beta_x.^2+beta_y.^2+beta_z.^2));
save_struct_to_gdf_file(fullfile(gptdir,'bunchData.gdf'), data); % save Matlab data file to GPT .gdf file
% - write GPT run file and execute GPT tracking
fid=fopen(fullfile(gptdir,'GPT.in'),'w'); % open txt file for generating GPT run instructions
fprintf(fid,'accuracy(6);\n'); % GPT accuracy parameter
fprintf(fid,'m=me ;\n'); % electron mass
fprintf(fid,'q=qe ;\n'); % electron charge
fprintf(fid,'setfile("beam","%s") ;\n',fullfile(gptdir,'bunchData.gdf')) ; % the bunch data file we generated above
x0=0; xp0=0;
xx=cos(xp0); xz=sin(xp0); % x position/angle offset
zoff=0; % z offset
fprintf(fid,'settransform("wcs",%g,0,%g,%g,0,%g,0,1,0,"beam");\n',x0,zoff,xx,xz) ; % transform the input beam coordinates if desired
fprintf(fid,'map3D_B("wcs",0,0,0,1,0,0,0,1,0,"%s","x","y","z","bx","by","bz",1) ;\n',fullfile(pwd,BfieldFileName)); % load pre-generated B field map
% Place screen elements where we want the beam coordinates to be read out
for iscreen=1:nscreen
  ox=xscr(iscreen); oy=yscr(iscreen); oz=zscr(iscreen);
  xp=xpscr(iscreen); yp=ypscr(iscreen);
  xx=cos(xp); xz=sin(xp); yy=cos(yp); yz=sin(yp);
  fprintf(fid,'screen("wcs",%g,%g,%g,%g,0,%g,0,%g,%g,0) ;\n',ox,oy,oz,xx,xz,yy,yz);
end
fprintf(fid,'snapshot(%g,%g,%g) ;',t0,t1,tprec);
fclose(fid);
% Execute GPT tracking - results go into "result.gdf" file
sid=system(sprintf('gpt -o %s %s',fullfile(gptdir,'result.gdf'),fullfile(gptdir,'GPT.in')));
if sid; error('GPT run error'); end;
% - read in GPT tracking results at screen locations and store in
% cell array of Lucretia beams
g=load_gdf(fullfile(gptdir,'result.gdf'));
scr=flip(find(arrayfun(@(x) isfield(g(x).p,'position'),1:length(g))));
for iscr=1:nscreen
  id=(g(scr(iscr)).d.ID); % particle order in GPT file entries
  % - flag any missing missing particles and set these as NaN's in output
  badid=[];
  if length(id)~=length(BeamIn.Bunch.Q)
    badid=find(~ismember(1:length(BeamIn.Bunch.Q),id));
  end
  % - make cell array of output beams in Lucretia format
  Bout{iscr}=BeamIn;
  xang(id)=atan(g(scr(iscr)).d.Bx./g(scr(iscr)).d.Bz);
  yang(id)=atan(g(scr(iscr)).d.By./g(scr(iscr)).d.Bz);
  % If required, can set the z coordinate compared with some reference
  % trajectory, here just set to zero
%   dt=g(scr(iscr)).d.t-obj.refTraj.t(iscr);
  dt=0;
  beta=sqrt(g(scr(iscr)).d.Bx.^2+g(scr(iscr)).d.By.^2+g(scr(iscr)).d.Bz.^2);
  z(id)=dt.*clight.*beta;
  gamma=1./sqrt(1-beta.^2);
  E(id)=gamma.*me;
  xv(id)=g(scr(iscr)).d.x;
  yv(id)=g(scr(iscr)).d.y;
  if ~isempty(badid)
    xv(badid)=NaN; xang(badid)=NaN; yv(badid)=NaN; yang(badid)=NaN; z(badid)=NaN; E(badid)=NaN;
  end
  Bout{iscr}.Bunch.x=[xv(:)'; xang(:)'; yv(:)'; yang(:)'; z(:)'; E(:)'];
end
% - Delete temp directory
rmdir(gptdir,'s');
% - if required, can uncomment below to get output at each tprec timestep
%       xv=arrayfun(@(x) g(x).d.x,1:length(g)-obj.nslice-1,'UniformOutput',false);
%       xpv=arrayfun(@(x) g(x).d.Bx,1:length(g)-obj.nslice-1,'UniformOutput',false);
%       zv=arrayfun(@(x) g(x).d.z,1:length(g)-obj.nslice-1,'UniformOutput',false);
%       zpv=arrayfun(@(x) g(x).d.Bz,1:length(g)-obj.nslice,'UniformOutput',false);
%       xang=arrayfun(@(x) mean(atan(xpv{x}./zpv{x})),1:length(xv));
%       mx=cellfun(@(x) mean(x),xv);
%       mz=cellfun(@(x) mean(x),zv);
end