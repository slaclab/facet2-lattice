global BEAMLINE
% Associate bmad wakefield definition files with all s-band structures
for iele=findcells(BEAMLINE,'Freq',2856)
  BEAMLINE{iele}.SRWFfile='sband.wake';
end
% Write out bmad lattice, un-splitting all magnetic elements
DT=DeckTool('BMAD',-1);
DT.WriteDeck(Initial,'FACET2e.bmad','facet2e',true);