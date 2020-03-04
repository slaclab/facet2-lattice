classdef GPT_FacetInjector < handle & physConsts & batchLSF
  % GPT_FACETINJECTOR
  % run on batch submission machine with >> FI=GPT_FacetInjector(nseed,'bsetup');FI.run
  
  properties
    tmpDirectory='/scratch';
    jobName='F2Injector';
    dataDirectory='/home/glen/OneDrive/Work/FACET2/data/F2Injector';
    runDirectory='/home/glen/OneDrive/Work/FACET2/data/F2Injector';
    batchQueueName='short'; % nmacro<=1e5 (short) 1e5<nmacro<1e7 (medium)
    mcrRoot='/nfs/slac/g/nlcrd/u01/whitegr/mcr/v96';
    srcdir='/home/glen/OneDrive/Work/FACET2/GPT';
    gptloc='/usr/local/bin/';
    gunEz=120; % on-crest gun accelerating gradient (MV/m)
    gunPhase=280; % gun RF phase offset (degrees)
    gunSolB=0.25; % Peak gun solenoid field strength (T)
    L0a_phase=133; % L0a phase (degrees) [peak accelerating phase] %133
    L0_phaseOffset=0;
    L0a_V=20; % L0a accelerating gradient (MV/m)
    nmacro=1e6; % # of macro-particles
    nsigcut=3; % Cut threhold for gaussian random number generator
    dojitter=false; % apply jitter to initial laser distribution or not
    laserjitter=[0.03,0.03]; % n std jitter of laser on cathode
    qjitter=0.01; % Relative charge jitter
    tjitter=200e-15;  % timing jitter
%     Q0=[0.5e-9 1.5e-9]; % Initial bunch charge
    Q0=2e-9; % Initial bunch charge
    laserR=1.5e-3; % cut radius, distribution is rms with 2XlaserR sigma (2.68)
    dt=11.62e-12; % no LH=11e-12 LH,500=11.62e-12
%     sigt=[2e-12 5.25e-12]; % Laser pulse duration per bunch (FWHM in s)
    sigt=7e-12; % Laser pulse duration per bunch (FWHM in s)
    BeamOut % Storage for beam at exit of simulation
    SimData
    StoreVec=false; % Store particle vectors
    nscan=30; % number of scan points (nscan^2 for 2d scans)
    scanpoints
    verbose=0;
    batchmethod='linscan'; % linscan | mltrain | rgen
    speckleWeight=0; % if >0, speckle the source distribution with relative fluctations at this level (uniform random aplitude +/-)
  end
  properties
    Beam
    BeamSingle
  end
  properties(SetAccess=private)
    nnet
  end
  properties(Access=private)
    tjitterval
  end
  properties(Constant)
    datafiles={'rfgundata.gdf' 'soldata.gdf' 'bunchData.gdf' 'bunchDataSingle.gdf' 'slwake.gdf' 'INJL0.saveline'};
  end
  
  methods
    function obj = GPT_FacetInjector(varargin)
      if nargin==0
        return
      end
      obj.runArgs='rhel60';
      obj.storeJobOutput=true; % keep screen output of job?
      obj.storeErrOutput=true; % keep error output stream?
      obj.nscan=varargin{1};
      obj.submitWaitTime=0.1;
      obj.maxfails=5;
      obj.wallTime=60; % nmacro=1e5 : 60 1e6 : 600
      obj.soLibs={sprintf('/nfs/slac/g/nlcrd/u01/whitegr/LucretiaLibs/%s',mexext)};
      switch obj.batchmethod
        case 'linscan'
          obj.nseed=obj.nscan^2;
          gp=linspace(270,300,obj.nscan);
          sb=linspace(0.2,0.27,obj.nscan);
          [GP,SB]=ndgrid(gp,sb);
          obj.scanpoints.A=GP;
          obj.scanpoints.B=SB;
        case 'mltrain'
          obj.nseed=varargin{1};
          obj.scanpoints.gunEz=[90 130];
          obj.scanpoints.gunPhase=[270 300];
          obj.scanpoints.gunSolB=[0.2 0.27];
          obj.scanpoints.laserR=[2e-3 3.5e-3];
          obj.scanpoints.sigt=[3e-12 10e-12];
        case 'rgen'
          % Nothing to do, just generating bunches with different random seeds
          obj.nseed=varargin{1};
        otherwise
          error('Unknown batch method');
      end
      if nargin>1
        for iarg=2:2:nargin
          if isprop(obj,varargin{iarg})
            obj.(varargin{iarg})=varargin{iarg+1};
          elseif strcmp(varargin{iarg},'bsetup')
            obj.dataDirectory='/nfs/slac/g/nlcrd/u01/whitegr/F2Injector';
            obj.runDirectory='/nfs/slac/g/nlcrd/u01/whitegr/F2Injector';
            obj.srcdir='/media/Data/OneDrive/FACET2/GPT';
            obj.gptloc='/var/GPT/gpt337_31_jan_2018_rhel_6_9_gcc_4_8_3/bin';
          else
            error('No such property: %s',varargin{iarg})
          end
        end
      end
      obj.MakeCathodeBeam;
      obj.writeBeamFile(obj.Beam);
      obj.writeBeamFile(obj.BeamSingle,'single');
      for ifile=1:length(obj.datafiles)
        copyfile(fullfile(obj.srcdir,obj.datafiles{ifile}),obj.runDirectory);
      end
      copyfile(fullfile(obj.gptloc,'gpt'),obj.runDirectory);
    end
    function B=run1(obj)
      % single bunch ops
      iseed=1;
      % Set gun phase & solenoid strength
      obj.gunPhase=288.6; % 8ps: 290
      obj.gunSolB=0.2507; % 8ps: 0.249
      % Set L0a phase for max acceleration
      xmin=fminsearch(@(x) L0aOptFun(obj,x,iseed),obj.L0a_phase,optimset('Display','off'));
      obj.L0a_phase=xmin+obj.L0_phaseOffset;
      % Write run file and run GPT (beam to end of L0a @ 66 MeV)
      td=pwd;
      cd(obj.runDirectory)
      obj.writeRunFile(iseed);
      obj.runGPT(iseed);
      % Track beam to end of injector @135 MeV with Lucretia
      B=obj.readData(iseed); B=B{1};
      tdat=obj.TrackInj(B);
      B=tdat.Beam;
      cd(td);
      [nx,ny]=GetNEmit90FromBeam(B,1);
      fprintf('Emittance: %g / %g Energy: %g\n',nx,ny,mean(B.Bunch.x(6,~B.Bunch.stop)))
      fprintf('sigz: %g mm sige: %g %%\n',std(B.Bunch.x(5,~B.Bunch.stop))*1e3,(std(B.Bunch.x(6,~B.Bunch.stop))/mean(B.Bunch.x(6,~B.Bunch.stop)))*100)
    end
    function data=run2(obj,iseed)
      % 2-bunch ops
      if ~exist('iseed','var')
        iseed=1;
      end
      % Set gun phase & solenoid strength
      obj.gunPhase=285;
      obj.gunSolB=0.244;
      obj.L0_phaseOffset=-9;
      % Set L0a phase for max acceleration
      xmin=fminsearch(@(x) L0aOptFun(obj,x,iseed),obj.L0a_phase,optimset('Display','off'));
      obj.L0a_phase=xmin+obj.L0_phaseOffset;
      % Write run file and run GPT (beam to end of L0a @ 66 MeV)
      td=pwd;
      cd(obj.runDirectory)
      obj.writeRunFile(iseed);
      obj.runGPT(iseed);
      % Track beam to end of injector @135 MeV with Lucretia
      B=obj.readData(iseed); B=B{1};
      tdat=obj.TrackInj(B);
      B=tdat.Beam;
      [nx1,ny1,nx2,ny2,nz1,nz2,dz,sz1,sz2]=obj.Get2bunchData(B);
      data.nx=[nx1 nx2]; data.ny=[ny1 ny2]; data.nz=[nz1 nz2]; data.dz=dz; data.sz=[sz1 sz2];
      cd(td);
      data.Beam=B;
      data.wID=1:obj.nmacro; data.dID=obj.nmacro+1:obj.nmacro*2;
    end
    function data=run1_env(obj)
      % single bunch ops - Get Envelope data
      iseed=1;
      obj.StoreVec=true;
      % Set gun phase & solenoid strength
      obj.gunPhase=280; % 8ps: 290
      obj.gunSolB=0.25; % 8ps: 0.249
      % Set L0a phase for max acceleration
      xmin=fminsearch(@(x) L0aOptFun(obj,x,iseed),obj.L0a_phase,optimset('Display','off'));
      obj.L0a_phase=xmin+obj.L0_phaseOffset;
      % Write run file and run GPT (beam to end of L0a @ 66 MeV)
      td=pwd;
      cd(obj.runDirectory)
      obj.writeRunFile(iseed);
      obj.runGPT(iseed);
      % Track beam to end of injector @135 MeV with Lucretia
      [Bout,xv,yv,zv,ev]=obj.readData(iseed);
      data.B=Bout; data.xv=xv; data.yv=yv; data.zv=zv; data.ev=ev;
    end
    function data=run2_env(obj,iseed)
      % Get envelope data
      if ~exist('iseed','var')
        iseed=1;
      end
      obj.StoreVec=true;
      % Set gun phase & solenoid strength
      obj.gunPhase=285;
      obj.gunSolB=0.244;
      obj.L0_phaseOffset=-9;
      % Set L0a phase for max acceleration
      xmin=fminsearch(@(x) L0aOptFun(obj,x,iseed),obj.L0a_phase,optimset('Display','off'));
      obj.L0a_phase=xmin+obj.L0_phaseOffset;
      % Write run file and run GPT (beam to end of L0a @ 66 MeV)
      td=pwd;
      cd(obj.runDirectory)
      obj.writeRunFile(iseed);
      obj.runGPT(iseed);
      % Track beam to end of injector @135 MeV with Lucretia
      [Bout,xv,yv,zv,ev]=obj.readData(iseed);
      data.B=Bout; data.xv=xv; data.yv=yv; data.zv=zv; data.ev=ev;
    end
    function data=submitFunc(obj,iseed)
      %#function applyLSC applyCSR
      rng(iseed);
      switch obj.batchmethod
        case 'rgen'
          obj.MakeCathodeBeam;
          obj.writeBeamFile(obj.Beam);
          obj.writeBeamFile(obj.BeamSingle,'single');
        case 'linscan'
          % Set gun phase & solenoid strength
          obj.gunPhase=obj.scanpoints.A(iseed);
          obj.gunSolB=obj.scanpoints.B(iseed);
        case 'mltrain'
          % Randomly set control variables requested by fieldnames of scanpoints property
          fn=fieldnames(obj.scanpoints);
          for ifn=1:length(fn)
            vals=obj.scanpoints.(fn{ifn});
            obj.(fn{ifn}) = vals(1) + rand*range(vals) ;
          end
      end
      % Set L0a phase for max acceleration
      xmin=fminsearch(@(x) L0aOptFun(obj,x,iseed),obj.L0a_phase,optimset('Display','off'));
      obj.L0a_phase=xmin+obj.L0_phaseOffset;
      % Write run file and run GPT (beam to end of L0a @ 66 MeV)
      if ~isdeployed
        td=pwd;
        cd(obj.runDirectory)
      end
      try
        obj.writeRunFile(iseed);
        obj.runGPT(iseed);
      catch ME
        delete(sprintf('GPT_%d.in',iseed));
        if ~isdeployed; cd(td); end
        rethrow(ME)
      end
      delete(sprintf('GPT_%d.in',iseed));
      % Track beam to end of injector @135 MeV with Lucretia
      B=obj.readData(iseed); B=B{1}; B_L0a=B;
      tdat=obj.TrackInj(B);
      B=tdat.Beam;
      if length(obj.Q0)==2
        [nx1,ny1,nx2,ny2,nz1,nz2,dz,sz1,sz2]=obj.Get2bunchData(B);
        data.nx=[nx1 nx2]; data.ny=[ny1 ny2]; data.nz=[nz1 nz2]; data.dz=dz; data.sz=[sz1 sz2];
      else
        [nx,ny,nz]=GetNEmitFromBeam(B,1);
        [nx90,ny90,nz90]=GetNEmit90FromBeam(B,1);
        [nxS,nyS,Q,~,deS]=plotSliceEmit(B,10e-6,false); [~,mI]=max(Q); nxS=nxS(mI); nyS=nyS(mI);
        sigz=std(B.Bunch.x(5,~B.Bunch.stop));
        sige=(std(B.Bunch.x(6,~B.Bunch.stop))/mean(B.Bunch.x(6,~B.Bunch.stop)))*100;
        data.nx=nx; data.ny=ny; data.sigz=sigz; data.sige=sige; data.E=mean(B.Bunch.x(6,~B.Bunch.stop));
        data.nz=nz; data.nz90=nz90;
        data.sigeS=deS(mI);
        data.nx90=nx90; data.ny90=ny90;
        data.nxS=nxS; data.nyS=nyS;
        switch obj.batchmethod
          case 'rgen'
            data.Beam=B;
%             data.Beam0=obj.Beam;
%             data.BeamL0a=B_L0a;
          case 'mltrain'
            fn=fieldnames(obj.scanpoints);
            for ifn=1:length(fn)
              data.trainvar.(fn{ifn}) = obj.(fn{ifn}) ;
            end
        end
      end
      if ~isdeployed; cd(td); end
    end
    function [nx,nx90,nxS,nz] = analyze_rgen(obj)
      nx=nan(1,obj.nseed); nx90=nx; nxS=nx; nz=nx;
      for iseed=1:obj.nseed
        if ~isempty(obj.retdata{iseed})
          nx(iseed)=obj.retdata{iseed}.nx;
          nx90(iseed)=obj.retdata{iseed}.nx90;
          nxS(iseed)=obj.retdata{iseed}.nxS;
          nz(iseed)=obj.retdata{iseed}.nz;
        end
      end
      figure
      subplot(2,2,1), histogram(nx);xlabel('\eta_x');
      subplot(2,2,2), histogram(nx90);xlabel('\eta_x (90%)');
      subplot(2,2,3), histogram(nxS);xlabel('\eta_x (Slice)');
      subplot(2,2,4), histogram(nx);xlabel('\eta_z');
    end
    function analyze2(obj)
      gp=obj.scanpoints.A;
      sb=obj.scanpoints.B;
      emit1=nan(size(gp)); emit2=emit1; dz=emit1; sz1=emit1; sz2=emit1;
      for iseed=1:obj.nseed
        if ~isempty(obj.retdata{iseed})
          emit1(iseed)=max([obj.retdata{iseed}.nx(1) obj.retdata{iseed}.ny(1)]);
          emit2(iseed)=max([obj.retdata{iseed}.nx(2) obj.retdata{iseed}.ny(2)]);
          dz(iseed)=obj.retdata{iseed}.dz;
          sz1(iseed)=obj.retdata{iseed}.sz(1);
          sz2(iseed)=obj.retdata{iseed}.sz(2);
        end
      end
      figure
      subplot(3,1,1)
      mesh(gp,sb,emit1)
      subplot(3,1,2)
      mesh(gp,sb,emit2)
      subplot(3,1,3)
      mesh(gp,sb,dz)
      figure
      subplot(3,1,1)
      mesh(gp,sb,sz1)
      subplot(3,1,2)
      mesh(gp,sb,sz2)
    end
    function analyze(obj)
      gp=obj.scanpoints.A;
      sb=obj.scanpoints.B;
      emit=nan(size(gp)); emitS=emit; sz=emit;
      for iseed=1:obj.nseed
        if ~isempty(obj.retdata{iseed})
          emit(iseed)=max([obj.retdata{iseed}.nx90 obj.retdata{iseed}.ny90]);
          emitS(iseed)=max([obj.retdata{iseed}.nxS obj.retdata{iseed}.nyS]);
          sz(iseed)=obj.retdata{iseed}.sigz;
        end
      end
      figure
      subplot(3,1,1)
      mesh(gp,sb,emit)
      subplot(3,1,2)
      mesh(gp,sb,emitS)
      subplot(3,1,3)
      mesh(gp,sb,sz)
    end
    function [nx1,ny1,nx2,ny2,nz1,nz2,dz,sz1,sz2]=Get2bunchData(obj,Beam)
      ib1=1:obj.nmacro; ib2=obj.nmacro+1:obj.nmacro*2;
      B1=Beam; B2=B1;
      B1.Bunch.x=B1.Bunch.x(:,ib1); B1.Bunch.Q=B1.Bunch.Q(ib1); B1.Bunch.stop=B1.Bunch.stop(ib1);
      B2.Bunch.x=B2.Bunch.x(:,ib2); B2.Bunch.Q=B2.Bunch.Q(ib2); B2.Bunch.stop=B2.Bunch.stop(ib2);
      [nx1,ny1,nz1]=GetNEmit90FromBeam(B1,1);
      [nx2,ny2,nz2]=GetNEmit90FromBeam(B2,1);
      dz=mean(B2.Bunch.x(5,~B2.Bunch.stop))-mean(B1.Bunch.x(5,~B1.Bunch.stop));
      sz1=std(B1.Bunch.x(5,~B1.Bunch.stop)); sz2=std(B2.Bunch.x(5,~B2.Bunch.stop));
      if obj.verbose>0
        fprintf('Emit B1: %g / %g um-rad B2: %g / %g um-rad [x/y]\n',nx1*1e6,ny1*1e6,nx2*1e6,ny2*1e6)
        fprintf('Bunch spacing dz= %g mm\n',dz*1e3)
        fprintf('Bunch length B1: %g B2: %g [mm]\n',sz1*1e3,sz2*1e3)
      end
    end
    function val=L0aOptFun(obj,x,iseed)
      obj.L0a_phase=x;
      obj.writeRunFile(iseed,'single');
      try
        obj.runGPT(iseed);
        Bout=obj.readData(iseed,'single');
        val=1/mean(Bout{1}.Bunch.x(6,:));
      catch
        val=10000000000;
      end
    end
    function MakeCathodeBeam(obj)
      % sigT in m (rms), Len in s (FWHM)
      for ibunch=1:length(obj.sigt)
        sigT=obj.laserR;
        Len=obj.sigt(ibunch);
        sigcut=obj.nsigcut;
        temitmult=1; % 1 for LCLS gun (0.9 um-rad /mm)
        npart=obj.nmacro;
        KE=1; % eV
        if ibunch==1
          obj.BeamSingle=CreateBlankBeam(1, length(obj.sigt), KE*1e-9, 1);
          obj.BeamSingle.Bunch.Q=obj.Q0;
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
        gamma=1+(KE/1e9)/obj.emass;
        beta=sqrt(1-gamma^-2);
        v=obj.clight*beta;
        Len=Len/(2*sqrt(2*log(2)));
        zvals=randn(npart,1).*Len.*v;
        px=randnt(sigcut,npart,1).*0.0009*temitmult;
        py=randnt(sigcut,npart,1).*0.0009*temitmult;
        pz=abs(randnt(sigcut,npart,1).*2.9e-6)+0.001978;
%         xvals=xvals+randnt(1)*std(xvals)*obj.laserjitter(1); % do this with offsets in GPT.in file instead
%         yvals=yvals+randnt(1)*std(yvals)*obj.laserjitter(2);
        zoff=v*((ibunch-1)*obj.dt);
        if ibunch==1
          B=CreateBlankBeam(1, npart, KE*1e-9, 1);
          X=B.Bunch.x;
          X(1:5,:)=[xvals'; atan(px./pz)'; yvals'; atan(py./pz)'; zvals'];
          B.Bunch.Q=ones(1,npart).*(obj.Q0(ibunch)/npart);
        else
          X=[xvals'; atan(px./pz)'; yvals'; atan(py./pz)'; zoff+zvals'; B.Bunch.x(6,:)];
          B.Bunch.Q=[B.Bunch.Q ones(1,npart).*(obj.Q0(ibunch)/npart)];
          B.Bunch.stop=zeros(size(B.Bunch.Q));
        end
        % Speckle the cathode distribution? (model QE non-uniformity)
        if obj.speckleWeight>0
          X([1 3],:)=obj.generateSpeckledGaussian(obj.laserR,ceil(sqrt(npart*10)),2,0.1,obj.speckleWeight,20,npart,0) ;
        end
        if ibunch==1
          B.Bunch.x=X;
        else
          B.Bunch.x=[B.Bunch.x X];
        end
        obj.Beam=B;
      end
    end
    function runGPT(obj,iseed)
      % Execute GPT program
      sid=system(sprintf('%s -o %s %s > /dev/null 2>&1',fullfile(obj.runDirectory,'gpt'),sprintf('result_%d.gdf',iseed),sprintf('GPT_%d.in',iseed)));
      if sid
        delete(sprintf('GPT_%d.in',iseed));
        error('GPT run error');
      end
      delete(sprintf('GPT_%d.in',iseed));
    end
    function [Bout,xv,yv,zv,ev]=readData(obj,iseed,cmd)
      % - read in GPT tracking results at screen location and store in
      % Lucretia Beam format
      nscreen=1;
      g=load_gdf(sprintf('result_%d.gdf',iseed));
      scr=flip(find(arrayfun(@(x) isfield(g(x).p,'position'),1:length(g))));
      if exist('cmd','var') && isequal(cmd,'single')
        BeamIn=obj.BeamSingle;
      else
        BeamIn=obj.Beam;
      end
      if ~isempty(scr)
        Bout=cell(1,nscreen);
        for iscr=1:nscreen
          id=(g(scr(iscr)).d.ID); % particle order in GPT file entries
          % - flag any missing missing particles and set these as stopped
          % particles in Lucretia beam definition data structure
          goodray=ismember(1:length(BeamIn.Bunch.Q),id) ;
          % - make cell array of output beams in Lucretia format
          if exist('cmd','var') && isequal(cmd,'single')
            Bout{iscr}=obj.BeamSingle;
          else
            Bout{iscr}=obj.Beam;
          end
          xang=atan(g(scr(iscr)).d.Bx./g(scr(iscr)).d.Bz);
          yang=atan(g(scr(iscr)).d.By./g(scr(iscr)).d.Bz);
          gamma=g(scr(iscr)).d.G;
          beta=sqrt(1-(1./gamma.^2));
          v=beta.*obj.clight;
          z=g(scr(iscr)).d.t.*v-g(scr(iscr)).d.z;
          E=obj.emass.*sqrt(gamma.^2-1);
          xv=g(scr(iscr)).d.x;
          yv=g(scr(iscr)).d.y;
          Bout{iscr}.Bunch.x(:,goodray)=[xv'; xang'; yv'; yang'; z'; E'];
          Bout{iscr}.Bunch.stop(~goodray)=1;
          Bout{iscr}.Bunch.x(5,goodray)=Bout{iscr}.Bunch.x(5,goodray)-mean(Bout{iscr}.Bunch.x(5,goodray));
          Bout{iscr}.Bunch.x(:,~goodray)=0;
        end
      else
        Bout=[];
      end
      % Read snapshot locations into vectors
      if nargout>1
        xv=nan(length(g),length(obj.Beam.Bunch.Q)); yv=xv; zv=xv; ev=xv;
        for ig=1:length(g)
          if any(ig==scr)
            continue;
          end
          id=(g(ig).d.ID); % particle order in GPT file entries
          xv(ig,id)=g(ig).d.x;
          yv(ig,id)=g(ig).d.y;
          zv(ig,id)=g(ig).d.z;
          beta=sqrt(g(ig).d.Bx.^2+g(ig).d.By.^2+g(ig).d.Bz.^2);
          gamma=1./sqrt(1-beta.^2);
          ev(ig,id)=gamma.*obj.emass;
        end
      end
%       delete(sprintf('result_%d.gdf',iseed));
    end
    function writeBeamFile(obj,BeamIn,cmd)
      % - Write Luctia bunch particles out in GPT input file
      xGPT=BeamIn.Bunch.x;
      gamma=1+BeamIn.Bunch.x(6,:)./0.511e-3;
      beta=sqrt(1-gamma.^-2);
      beta_x=beta.*sin(xGPT(2,:)); beta_y=beta.*sin(xGPT(4,:)); beta_z=sqrt(1-(beta_x./beta).^2-(beta_y./beta).^2).*beta;
      data=struct; data.d.x=xGPT(1,:); data.d.y=xGPT(3,:); % GPT data structure for pos coordinates
      % Convert z distribution into particle release times
      data.d.z=randn(1,length(BeamIn.Bunch.Q)).*0;
      v=beta*obj.clight;
      data.d.t=xGPT(5,:)./v;
      data.d.t=data.d.t-min(data.d.t);
      data.d.Bx=beta_x'; data.d.By=beta_y'; data.d.Bz=beta_z'; data.d.nmacro=BeamIn.Bunch.Q./obj.eQ';
      data.d.G=1./sqrt(1-(beta_x.^2+beta_y.^2+beta_z.^2));
      if exist('cmd','var') && isequal(cmd,'single')
        save_struct_to_gdf_file(fullfile(obj.runDirectory,'bunchDataSingle.gdf'), data); % save Matlab data file to GPT .gdf file
      else
        save_struct_to_gdf_file(fullfile(obj.runDirectory,'bunchData.gdf'), data); % save Matlab data file to GPT .gdf file
      end
    end
    function writeRunFile(obj,iseed,cmd)
      % Generate GPT lattice and run file
      % - require beam to have been setup
      if isempty(obj.Beam)
        error('No beam defined');
      end
      fid=fopen(sprintf('GPT_%d.in',iseed),'w'); % open txt file for generating GPT run instructions
      %       fprintf(fid,'accuracy(6);\n'); % GPT accuracy parameter
      fprintf(fid,'m=me ;\n'); % electron mass
      fprintf(fid,'q=qe ;\n'); % electron charge
      if exist('cmd','var') && isequal(cmd,'single')
        fprintf(fid,'setfile("beam","%s") ;\n',fullfile(obj.runDirectory,'bunchDataSingle.gdf')) ; % the bunch data file (single particle)
      else
        fprintf(fid,'setfile("beam","%s") ;\n',fullfile(obj.runDirectory,'bunchData.gdf')) ; % the bunch data file
      end
      if obj.dojitter
        fprintf(fid,'settotalcharge("beam",%g) ;\n',-sum(obj.Beam.Bunch.Q)*(1+randnt(1)*obj.qjitter));
      end
      if obj.dojitter
        boff=randnt(1,1,2).*[std(obj.Beam.Bunch.x(1,:)).*obj.laserjitter(1) std(obj.Beam.Bunch.x(3,:)).*obj.laserjitter(2)];
        fprintf(fid,'settransform("wcs",%.8g,%.8g,0,1,0,0,0,1,0,"beam");\n',boff) ; % offset beam
      end
      fprintf(fid,'Spacecharge3Dmesh("Cathode");\n'); % Cathode
%       fprintf(fid,'Spacecharge3Dmesh("Cathode","MeshNfac",1.5,"MeshAdapt",0.25,"SolverAcc",0.005,"MeshBoxSize",5);\n'); % Cathode
      fprintf(fid,'map1D_TM("wcs","z",0,"%s","z","Ez",%g,%g,1.79447772373049e10) ;\n',fullfile(obj.runDirectory,'rfgundata.gdf'),obj.gunEz*1e6/2.4523,deg2rad(obj.gunPhase)); % RF Gun
      fprintf(fid,'map1D_B("wcs","z",0,"%s","z","Bz",%g) ;\n',fullfile(obj.runDirectory,'soldata.gdf'),obj.gunSolB); % Gun solenoid
      gamma=0.004734/0.511e-3; % Design for 6 MeV
      rfL=3.0429; %3.0429;
      a0=0.1;
      P=(obj.L0a_V*1e6)^2/(2*50e6*a0);
      gunDrift=0.9696;
      % Convert time jitter to degrees of s-band phase
      if obj.dojitter
        obj.tjitterval=randnt(1)*obj.tjitter;
        dpha=(obj.tjitterval)/(1/2.8560e+09/360);
      else
        dpha=0;
      end
      fprintf(fid,'trwlinac("wcs","z",%g,%g,50e6,%g,%g,%g,0,%g,1.79447772373049e10,%g) ;\n',gunDrift+rfL/2,a0,P,P,gamma,deg2rad(obj.L0a_phase+dpha),rfL) ; % L0a
      if ~exist('cmd','var') || ~isequal(cmd,'single')
        fprintf(fid,'wakefield("wcs","z",%g,%g,%g,"%s","z","","","Wz") ;\n',gunDrift+rfL/2,rfL,50,fullfile(obj.runDirectory,'slwake.gdf')); % Wakefield
      end
      fprintf(fid,'screen("wcs","I",%g);\n',gunDrift+rfL); % Beam output location
      if obj.StoreVec
        fprintf(fid,'snapshot(%g,%g,%g) ;\n',0,14e-9,0.1e-10);
      end
      fclose(fid);
    end
    function datDL1BEG=TrackInj(obj,B)
      global BEAMLINE KLYSTRON PS WF GIRDER
      
      % Physics switches
      LSCon=1;
      LSCLEN=0.01;
      
      % Apply timing jitter to beam
      if obj.dojitter
        clight=299792458;
        B.Bunch.x(5,:)=B.Bunch.x(5,:)+clight*obj.tjitterval;
      end
      
      % Load injector beamline
      BL_INIT=BEAMLINE; PS_INIT=PS; WF_INIT=WF; KLY_INIT=KLYSTRON; GIR_INIT=GIRDER;
      if ~isempty(BEAMLINE)
        BEAMLINE=[];
        PS=[]; WF=[]; KLYSTRON=[]; GIRDER={}; %#ok<NASGU>
      end
      XSIFToLucretia('INJL0.saveline','INJL0');
      SetSPositions(1, length(BEAMLINE), 0);
      
      % assign element slices and blocks
      SetElementSlices(1,length(BEAMLINE));
      SetElementBlocks(1,length(BEAMLINE));
      
      % Identify elements
      inds.i1=findcells(BEAMLINE,'Name','L0BBEG');
      inds.i2=length(BEAMLINE);
      inds.L0b=findcells(BEAMLINE,'Name','L0B*');
      inds.L0b=inds.L0b(ismember(inds.L0b,findcells(BEAMLINE,'Class','LCAV')));
      inds.qm=[findcells(BEAMLINE,'Name','QA01') findcells(BEAMLINE,'Name','QA02') ...
        findcells(BEAMLINE,'Name','QE01B') findcells(BEAMLINE,'Name','QE02B')];
      
      % Assign PS's
      inds.qmps=[];
      for iqm=inds.qm(1:2:end)
        AssignToPS(BEAMLINE{iqm}.Slices,length(PS)+1);
        inds.qmps(end+1)=length(PS);
      end
      
      % define wakefields
      decm=100;
      bw=0.01;
      load srwf_long_sband.mat wf wfz
      WF.ZSR(1).z=wfz(1:decm:end); WF.ZSR(1).K=wf(1:decm:end); WF.ZSR(1).BinWidth=bw;
      load srwf_trans_sband.mat wf wfz
      WF.TSR(1).z=wfz(1:decm:end); WF.TSR(1).K=wf(1:decm:end); WF.TSR(1).BinWidth=bw;
      
      % assign wakefields
      id=findcells(BEAMLINE,'Class','LCAV',1,length(BEAMLINE));
      for n=1:length(id)
        BEAMLINE{id(n)}.Wakes=[1,1];
      end
      
      % assign klystrons
      AssignToKlystron(inds.L0b,1);
      
      % Set L0b amplitude and phase to deliver unchirped beam at design
      % energy at L0b exit
      goodray=~B.Bunch.stop;
      for idim=1:5
        B.Bunch.x(idim,goodray)=B.Bunch.x(idim,goodray)-mean(B.Bunch.x(idim,goodray));
      end
      T=Track(B);
      T.startInd=inds.i1;
      T.finishInd=inds.i2;
      SetDesignMomentumProfile( inds.i1, inds.i2, sum(B.Bunch.Q(goodray)), mean(B.Bunch.x(6,goodray))) ;
      if obj.verbose>1
        dtxt='iter';
      else
        dtxt='off';
      end
      xmin=fminbnd(@(x) obj.eopt(x,T,0.135),-180,180,optimset('Display',dtxt));
      KLYSTRON(1).Phase=xmin+obj.L0_phaseOffset;
      L0bPhaseOpt=xmin+obj.L0_phaseOffset;
      SetDesignMomentumProfile( inds.i1, inds.i2, sum(B.Bunch.Q(goodray)), mean(B.Bunch.x(6,goodray)), 0.135) ;
      
      % Match beam to required twiss at exit of QE02B
      % bx=6.9889, ax=-1.4, by=9.3226, ay=-1.7619
      MovePhysicsVarsToPS( inds.qmps );
      xmin=lsqnonlin(@(x) obj.topt(x,inds.qmps,T,6.9889,-1.4,9.3226,-1.7619),...
        arrayfun(@(x) PS(x).Ampl,inds.qmps),-2.*ones(1,4),2.*ones(1,4),optimset('Display',dtxt));
      for ips=1:length(inds.qmps)
        PS(inds.qmps(ips)).Ampl=xmin(ips);
        PS(inds.qmps(ips)).SetPt=xmin(ips);
        RenormalizePS( inds.qmps(ips) );
      end
      
      % Matching quad settings
      ind=0;
      clight=2.99792458e8; % speed of light (m/sec)
      Cb=1e9/clight;       % rigidity conversion (T-m/GeV)
      if obj.verbose>0; disp('Matching QUAD settings:'); end
      for ips=inds.qmps
        ind=ind+1;
        K(ind)=PS(ips).Ampl*sum(arrayfun(@(x) BEAMLINE{x}.B(1),PS(ips).Element))/...
          ((Cb*BEAMLINE{PS(ips).Element(1)}.P)*sum(arrayfun(@(x) BEAMLINE{x}.L,PS(ips).Element)));
        BDES(ind)=PS(ips).Ampl*sum(arrayfun(@(x) BEAMLINE{x}.B(1),PS(ips).Element))*10;
        if obj.verbose>0
          fprintf('%s: BDES= %g K1 = %g\n',BEAMLINE{PS(ips).Element(1)}.Name,...
            BDES(ind),K(ind))
        end
      end
      
      % -- Space charge
      if LSCon
        for iele=findcells(BEAMLINE,'TrackFlag')
          BEAMLINE{iele}.TrackFlag.LSC=1;
          BEAMLINE{iele}.TrackFlag.LSC_storeData=0;
          BEAMLINE{iele}.TrackFlag.LSC_npow2Bins=11;
          BEAMLINE{iele}.TrackFlag.LSC_smoothFactor=3;
          % Set NBPM on LCAV elements or Splits elsewhere
          % to ensure 0.3mm drift sections for application of LSC
          if strcmp(BEAMLINE{iele}.Class,'LCAV')
            BEAMLINE{iele}.NBPM=ceil(BEAMLINE{iele}.L/LSCLEN);
            BEAMLINE{iele}.GetSBPMData=1;
            BEAMLINE{iele}.GetInstData=1;
          else
            BEAMLINE{iele}.TrackFlag.Split=ceil(BEAMLINE{iele}.L/LSCLEN);
          end
        end
      end
      
      % Beam at DL1BEG
      T.finishInd=findcells(BEAMLINE,'Name','DL1BEGB');
      T.trackThru();
      B_dl1=T.beamOut();
      goodray=~B_dl1.Bunch.stop;
      [nx,ny,nz]=GetNEmitFromBeam(B_dl1,1);
      [nx90,ny90,nz90,Tx,Ty]=GetNEmit90FromBeam(B_dl1);
      E=mean(B_dl1.Bunch.x(6,goodray));
      dE=100.*std(B_dl1.Bunch.x(6,goodray))/E;
      sigz=std(B_dl1.Bunch.x(5,goodray));
      obj.BeamOut=B_dl1;
      % de-trended dE/E
      [P,~,MU] = polyfit(B_dl1.Bunch.x(5,goodray),B_dl1.Bunch.x(6,goodray),9);
      dE_min=100*std(B_dl1.Bunch.x(6,goodray) - polyval(P,B_dl1.Bunch.x(5,goodray),[],MU))/E;
      
      if obj.verbose>0
        disp('Tracked Bunch @ DL1BEG:')
        fprintf('Emit (um): x= %g y= %g z= %g\n',nx*1e6,ny*1e6,nz*1e6)
        fprintf('Emit90 (um): x= %g y= %g z= %g\n',nx90*1e6,ny90*1e6,nz90*1e6)
        fprintf('E= %g dE/E= %g dE/E_min= %g sigz= %g\n',E,dE,dE_min,sigz)
        fprintf('betax= %g alphax= %g betay= %g alphay= %g\n',Tx.beta,Tx.alpha,Ty.beta,Ty.alpha)
        fprintf('L0b Phase = %g deg\n',BEAMLINE{inds.L0b(1)}.Phase)
      end
      if obj.verbose>1
        beamImage(B_dl1,0);
      end
      
      % Organize output data
      datDL1BEG.Beam=B_dl1;
      datDL1BEG.Tx=Tx; datDL1BEG.Ty=Ty;
      datDL1BEG.emit=[nx ny nz];
      datDL1BEG.emit90=[nx90 ny90 nz90];
      datDL1BEG.dE=dE;
      datDL1BEG.dE_min=dE_min;
      datDL1BEG.QM_K1=K;
      datDL1BEG.QM_BDES=BDES;
      datDL1BEG.L0bPhase=L0bPhaseOpt;
      
      % Put back original beamline if there was one
      BEAMLINE=BL_INIT; PS=PS_INIT; WF=WF_INIT; KLYSTRON=KLY_INIT; GIRDER=GIR_INIT;
    end
    function [xopt,fval]=SurrogateOpt(obj,controlvar,optvar)
      %SURROGATEOPT Optimize injector parameters based on Surrogate model
      % controlvar: vector of control variable indices to use (index is order of scanpoints property fields)
      % optvar: Vector of ovservables to read out (index is into retdata property trainvar fields)
      
      % Set upper and lower bounds for optimizer based on scanpoints property fields
      fn=fieldnames(obj.scanpoints); fn=fn(controlvar);
      lb=zeros(1,length(fn)); ub=lb;
      for ifn=1:length(fn)
        lb(ifn)=obj.scanpoints.(fn{ifn})(1);
        ub(ifn)=obj.scanpoints.(fn{ifn})(2);
      end
      opts = optimoptions('gamultiobj','PlotFcn',@gaplotpareto,'Display','iter');
      [xopt,fval] = gamultiobj(@(x) obj.SurrogateRun(x,controlvar,optvar),length(fn),-ones(length(lb),length(controlvar)),zeros(length(lb),1),[],[],lb,ub,[],opts);
      
    end
    function out = SurrogateRun(obj,inX,Xvar,optvar)
      %SURROGATERUN Run surrogate model with desired input parameters
      % inX: Input parameters corresponding to scanpoints property fields and ranges
      % Xvar: Vector of control variables used (index is order of scanpoints property fields), others are set according to default object values
      % optvar: Vector of ovservables to read out (index is into retdata property trainvar fields)
      % Check input
      
      if numel(inX) ~= numel(Xvar)
        error('Incorrect input variables')
      end
      
      % Set non provided scan variables to object default values
      fn=fieldnames(obj.scanpoints);
      xval=nan(length(fn),1);
      xval(Xvar)=inX;
      for ifn=1:length(fn)
        if isnan(xval(ifn))
          xval(ifn)=obj.(fn{ifn});
        end
      end
      
      % Run neural net
      out=obj.nnet(xval);
      out=abs(out(optvar));
      
    end
    function [data_in,data_out,names_in,names_out,performance,net]=SurrogateTrain(obj)
      % DOTRAIN Train NN surrogate model based on scanpoints data
      
      % Unpack data into format for NN regression learner
      fn=fieldnames(obj.retdata{1});
      names_out=fn(~ismember(fn,'trainvar'));
      names_in=fieldnames(obj.retdata{1}.trainvar);
      data_in=nan(obj.nseed,length(names_in));
      data_out=nan(obj.nseed,length(names_out));
      for idata=1:obj.nseed
        for iname=1:length(names_in)
          if ~isempty(obj.retdata{idata}) && isfield(obj.retdata{idata}.trainvar,names_in{iname}) && ~isempty(obj.retdata{idata}.trainvar.(names_in{iname}))
            data_in(idata,iname)=obj.retdata{idata}.trainvar.(names_in{iname});
          end
        end
        for iname=1:length(names_out)
          if ~isempty(obj.retdata{idata}) && isfield(obj.retdata{idata},names_out{iname}) && ~isempty(obj.retdata{idata}.(names_out{iname}))
            data_out(idata,iname)=obj.retdata{idata}.(names_out{iname});
          end
        end
      end
      data_in(any(isnan(data_in)'),:)=[];
      data_out(any(isnan(data_out)'),:)=[]; data_out=abs(data_out);
      
      % Train the neural network
      x = data_in';
      t = data_out';
      
      % Choose a Training Function
      % For a list of all training functions type: help nntrain
      % 'trainlm' is usually fastest.
      % 'trainbr' takes longer but may be better for challenging problems.
      % 'trainscg' uses less memory. Suitable in low memory situations.
      trainFcn = 'trainbr';  % Bayesian Regularization backpropagation.
      
      % Create a Fitting Network
      hiddenLayerSize = 30;
      net = fitnet(hiddenLayerSize,trainFcn);
      
      % Setup Division of Data for Training, Validation, Testing
      net.divideParam.trainRatio = 70/100;
      net.divideParam.valRatio = 15/100;
      net.divideParam.testRatio = 15/100;
      
      % Train the Network
      net.trainParam.max_fail=10;
      [net,~] = train(net,x,t);
      
      % Test the Network
      y = net(x);
      e = gsubtract(t,y);
      performance = perform(net,t,y);
      
      % View the Network
      view(net)
      
      % Plots
      figure, plotregression(t,y)
      
      % Store network
      obj.nnet=net;
    end
  end
  methods(Static,Hidden)
    function xyrays = generateSpeckledGaussian(R,nbin,nsigma,speckleSize,speckleWeight,nspeckle,nray,doPlot)
      % Generates a Gaussian Distribution with speckles
      % Speckle amplitude is implemented as a random gaussian spot (size speckleSize relative
      %  to underlying Gaussian), placed randomly within generated Gaussian area (nspeckle times)
      % nray = number of macro particles to generate
      % R= physical radius (m)
      
      % Generate grid
      xran=linspace(-R,R,nbin); dx=diff(xran(1:2));
      [X,Y] = ndgrid(xran, xran); XY = [X(:) Y(:)]; [~,Ri]=cart2pol(X,Y);
      
      % Generate background Gaussian
      gaussDist = mvnpdf(XY,[0 0],[R*nsigma R*nsigma].^2); gaussDist=gaussDist./max(gaussDist(:)); gaussDist = reshape(gaussDist,size(X));
      gaussDist(Ri>R) = 0;
      speckleDist = gaussDist;
      
      % Generate speckles, width, amplitude and position relative to background Gaussian
      for n=1:nspeckle
        x0 = (rand*2-1)*R ;
        y0 = (rand*2-1)*R ;
        xwidth = speckleSize*R;
        wt=interp2(xran,xran,gaussDist,x0,y0)*speckleWeight*(1-2*rand); % weight is relative to Gaussian function at point of placement
        speckle = mvnpdf(XY,[x0 y0],ones(1,2).*xwidth^2); speckle = wt .* speckle./max(speckle(:)) ; speckle = reshape(speckle,size(X));
        speckleDist = speckleDist + speckle;
      end
      
      % randomly generate nray samples from speckleDist
      speckleDist(speckleDist<0) = 0 ;
      speckleDist = abs(speckleDist)./sum(abs(speckleDist(:))) ;
      snum=randsample(1:numel(speckleDist),nray,true,speckleDist(:));
      xyrays=[X(snum); Y(snum)];
      % Randomly smear each ray inside the bin it was generated to prevent putting multiple rays on top of each other
      xyrays=xyrays+(2.*rand(size(xyrays))-1).*(dx/2);
      
      if doPlot
        figure;
        subplot(1,3,1)
        mesh(X,Y,gaussDist);
        minval=min(gaussDist(gaussDist>0));
        set(gca,'Zlim',[minval 1]);
        caxis([minval 1]);
        title('Ideal Gaussian')
        colorbar
        subplot(1,3,2)
        mesh(X,Y,speckleDist);
        speckleDist(Ri>R)=0;
        minval=min(speckleDist(speckleDist>0));
        set(gca,'Zlim',[minval max(speckleDist(:))]);
        caxis([minval max(speckleDist(:))]);
        title('Gaussian with Speckles')
        colorbar
        subplot(1,3,3)
        histogram2(X(snum),Y(snum)); xlabel('X'); ylabel('Y');
      end
    end
    function ret=topt(x,qmps,T,bx,ax,by,ay)
      global PS
      for ips=1:length(qmps)
        PS(qmps(ips)).Ampl=x(ips);
      end
      T.trackThru(); bo=T.beamOut;
      try
        [~,~,~,Tx,Ty]=GetNEmit90FromBeam(bo);
      catch
        ret=[100 100 100 100];
        return
      end
      ret = [Tx.beta-bx Tx.alpha-ax Ty.beta-by Ty.alpha-ay];
    end
    function ret=eopt(x,T,E0)
      global KLYSTRON
      KLYSTRON(1).Phase=x;
      T.trackThru(); bo=T.beamOut;
      if T.trackStatus{1}~=1
        ret=1e5;
      else
        ret = std(bo.Bunch.x(6,~bo.Bunch.stop))/E0;
      end
    end
  end
end

