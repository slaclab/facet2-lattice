function [nx,ny,nz,Tx,Ty]=GetNEmit90FromBeam(B, bid)

if ~exist('bid','var')
  bid=1;
end
rejfrac=0.1;
goodray=~B.Bunch.stop;
if ~any(goodray)
  error('All particles stopped')
end
nrej=floor(sum(goodray)*rejfrac);
for idim=[1 3 5]
  Bx=B;
  Bx.Bunch(bid).x=Bx.Bunch(bid).x(:,goodray); Bx.Bunch(bid).Q=Bx.Bunch.Q(goodray); Bx.Bunch(bid).stop=Bx.Bunch(bid).stop(goodray);
  [~,Ix]=sort(abs(Bx.Bunch(bid).x(idim,:)-mean(Bx.Bunch(bid).x(idim,:))),'descend');
  [~,Ixp]=sort(abs(Bx.Bunch(bid).x(idim+1,:)-mean(Bx.Bunch(bid).x(idim+1,:))),'descend');
  I=unique([Ix(nrej+1:end) Ixp(nrej+1:end)]);
  Bx.Bunch(bid).x=Bx.Bunch(bid).x(:,I); Bx.Bunch(bid).Q=Bx.Bunch(bid).Q(I); Bx.Bunch(bid).stop=Bx.Bunch(bid).stop(I);
  if idim==1
    Tx = GetUncoupledTwissFromBeamPars(Bx,bid);
    nx=GetNEmitFromBeam(Bx,bid);
  elseif idim==3
    [~,Ty] = GetUncoupledTwissFromBeamPars(Bx,bid);
    [~,ny]=GetNEmitFromBeam(Bx,bid);
  else
    [~,~,nz]=GetNEmitFromBeam(Bx,bid);
  end
end