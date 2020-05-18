function MakeImpactBeam(npart,sigT,Len,KE,modFreq)
% sigT in m (rms), Len in s (FWHM), KE in eV

sigcut=3;
zoff=0;
pzoff=0.001978;
temitmult=1; % 1 for LCLS gun (0.9 um-rad /mm)
% modFreq: 2e11:5e12
if exist('modFreq','var')
  modAmp=0.1; % Modulation amplitude
else
  modAmp=0;
  modFreq=1;
end

fid=fopen('partcl.data','w');
if isequal(npart,'single')
  fprintf(fid,'1\n');
  fprintf(fid,'0.0 0.0 0.0 0.0 %g %g\n',zoff,pzoff);
  fclose(fid);
  return
end

% Set Gaussian width fraction of cut radius
transcut=sigT;
sigT=sigT*2;

xvals=randn(npart*10,1).*sigT;
yvals=randn(npart*10,1).*sigT;
R=sqrt(xvals.^2+yvals.^2);
xvals(R>transcut)=[];
xvals=xvals(1:npart);
yvals(R>transcut)=[];
yvals=yvals(1:npart);
gamma=1+KE/0.511e6;
beta=sqrt(1-gamma^-2);
v=299792458*beta;
Len=Len/(2*sqrt(2*log(2)));
t=linspace(-Len*sigcut,Len*sigcut,1e7);
modsig=sin(2.*pi.*modFreq.*t);
gsig=gauss(t,0,Len);
zsig=gsig+modsig.*gsig.*modAmp;
cdf=cumsum(zsig); cdf=cdf./max(cdf);
zvals=interp1(cdf,t,rand(1,npart)).*v+zoff;
% zvals=randn(npart,1).*Len.*v+zoff;
px=randn(npart,1).*0.0009*temitmult;
py=randn(npart,1).*0.0009*temitmult;
pz=abs(randn(npart,1).*2.9e-6)+pzoff;
% histogram(zvals)
fprintf(fid,'%d\n',npart);
% ================================================
% MC RANDOM CATHODE JITTER
% ================================================
xvals=xvals+randnt(3)*std(xvals)*0.03;
yvals=yvals+randnt(3)*std(yvals)*0.03;
zvals=zvals+10e-15*v*randnt(3);
for ip=1:npart
  fprintf(fid,'%g %g %g %g %g %g\n',xvals(ip),px(ip),yvals(ip),py(ip),zvals(ip),pz(ip));
end
fclose(fid);
