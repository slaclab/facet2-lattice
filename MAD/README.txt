-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
07/21/2017, G. White: Re-formatted decks by MDW and re-configured injector for FACET-II requirements
05/02/2017, G. White: Updated to v1.4.4 DR optics (added steering correctors and BPMs to lattice)
04/28/2017, G. White: Added BC20E optics & updated DRTBC11 optics for v1.43 DR optics
04/27/2017, G. White: Added BC20P optics
07/13/2017, G. White: Updated structure for main Linac decks compiled by M. Woodley
04/26/2017, G. White: Update PRLTDR optics to be compatible with v1.43 DR optics
03/23/2016, G. White: Update to v1.43 of DR optics (sliced magnet models)
11/15/2016, G. White: Add v1.4 of DR optics (Harmonic #=51, with 2 or 4 RF cavities)
04/13/2016, G.White: TDR version
02/08/2015, G.White: Initial CVS release.
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

This archive file contains lattice files in MAD8 (or XSIF) and Lucretia format for
the FACET2 electron and positron beam lines.

See doc directory for some documentation.

===============
CONTENTS
===============
* doc/ : Documentation directory
* FACET2e.mad8 : Master deck file for parsing electron beamlines
* FACET2e.mat : Lucretia file for electron beamlines
* FACET2p.mad8 : Master deck file for parsing positron beamlines
* INJ / DL10 / BC11 / BC14 / L1 / L2 / L3 / LI20 .xsif : Subsystem decks
* FACET2p_PRLTDR.xsif : XSIF lattice for positron return line to DR injection.
* FACET2p_DRTBC11.xsif : Positron system from DR extraction, through vertical extraction
                        dogleg, through first e+ bunch compressor chicane, beam diagnostic
                        waist and horizontal dogleg injection to main beamline in last
                        BC11 bend magnet.
* FACET2p_DR.xsif : Positron damping ring v1.4.4
* FACET2p_DR.bmad : Positron damping ring v1.4.3 (BMAD format)
* FACET2p_DR.mat : AT and Lucretia lattice files for positron damping ring
* QDDSQ*.[xsif|bmad] : Sliced quad-sextupole arc models for DR
* BA/BD.[xsif|bmad] : Sliced models for bend-quad-sextupole arc magnets for DR
* sband_l.dat : Longitudinal short-range wakefield definition file for s-band structures.
* sband_t.dat : Transverse short-range wakefield definition file for s-band structures.
* xband_l.dat : Longitudinal short-range wakefield definition file for x-band structures.
* xband_t.dat : Transverse short-range wakefield definition file for x-band structures.
-----> BC20 optics below are upgrade optics for use in phase 3 with simultaneouse e-/e+ running (original W-chicane is contained in above optics)
* BC20P.xsif : BC20 chicane optics: positron arm
* BC20E.xsif : BC20 chicane optics: electron arm

