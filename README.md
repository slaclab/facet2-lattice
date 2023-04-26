# facet2-lattice
Lattice files for the FACET-II accelerator

## Setup 

Many of the paths will use the `$FACET2_LATTICE` environmental variable, which should point to this repository, as in:

`export FACET2_LATTICE=/path/to/this/repo`




## Changelog 
* Apr-26-2023, G. White: BC14 elements moved per measurements by Georg (see changes in L2.xsif, BC14.xsif & L3.xsif)
* Oct-17-2022, G. White: Implemented BC11 modifications per Brendan & survey measurements by Georg (See MAD/BC11.xsif for changelog)
* Oct-28-2021, G. White: Removed XTCAV from chicane, tweaked FFS & SPECT elements to match as-built reported by D. Storey 9-SEPT-2021 (LI20.xsif modified)
* Sept-3-2021, G. White: As-built layout for FFS changes in LI20.xsif (all singleton Q5,Q4,Q3,Q2,Q1,Q0FF quads)
* 05/05/2021, G. White: Added Q0FF, Q1FF, Q2FF, removed QFF4, locations according to metrology mesurements
* 04/17/2021, G. White: Fixed LI19 corrector locations, removed 19-8a and updated S20 experiment object locations @ Q0-2D quad locations
* 02/25/2021, G. White & M. Woodley: changes to reflect in-tunnel hardware for INJ, BC11, BC14 & S20 (see MAD xsif files for details)
* 01/22/2021, G. White: Move BC11 magnets and devices per B. O'Shea for edge radiation equipmemt installation
* 11/05/2020, G. White: Added model-independent data directory for storing e.g. Beam stay clear data
* 05/18/2020, G. White: Imported (older) ImpactT files for FACET-II injector (warning, doesn't exactly match GPT files)
* 03/03/2020, G. White: Files imported into slac lab / facet2-lattice github repo
* 07/21/2017, G. White: Re-formatted decks by MDW and re-configured injector for FACET-II requirements
* 05/02/2017, G. White: Updated to v1.4.4 DR optics (added steering correctors and BPMs to lattice)
* 04/28/2017, G. White: Added BC20E optics & updated DRTBC11 optics for v1.43 DR optics
* 04/27/2017, G. White: Added BC20P optics
* 07/13/2017, G. White: Updated structure for main Linac decks compiled by M. Woodley
* 04/26/2017, G. White: Update PRLTDR optics to be compatible with v1.43 DR optics
* 03/23/2016, G. White: Update to v1.43 of DR optics (sliced magnet models)
* 11/15/2016, G. White: Add v1.4 of DR optics (Harmonic #=51, with 2 or 4 RF cavities)
* 04/13/2016, G. White: TDR version
* 02/08/2015, G. White: Initial CVS release.

The "source" description of all lattices should be considered to be the MAD8 files, which are mirrored in the SLAC CVS repository.
This archive file contains lattice files in MAD8 (or XSIF), Lucretia and BMAD format for the FACET2 electron and positron beam lines.
The positron damping ring is also in AT format.
The injector, up to the exit of L0a is modelled in GPT, for which a GPT input and suporting files exists.

CONTENTS:
* Data/ : Model independent data : Beam stay clear files
* MAD/doc/ : Documentation directory
* MAD/FACET2e.mad8 : Master deck file for parsing electron beamlines
* Lucretia/models/FACET2e/FACET2e.mat : Lucretia file for electron beamlines
* MAD/FACET2p.mad8 : Master deck file for parsing positron beamlines
*   INJ / DL10 / BC11 / BC14 / L1 / L2 / L3 / LI20 .xsif : Subsystem decks
* MAD/FACET2p_PRLTDR.xsif : XSIF lattice for positron return line to DR injection.
* MAD/FACET2p_DRTBC11.xsif : Positron system from DR extraction, through vertical extraction
                        dogleg, through first e+ bunch compressor chicane, beam diagnostic
                        waist and horizontal dogleg injection to main beamline in last
                        BC11 bend magnet.
* MAD/FACET2p_DR.xsif : Positron damping ring v1.4.4
* MAD/FACET2p_DR.bmad : Positron damping ring v1.4.3 (BMAD format)
* Lucretia/models/FACET2p_DR.mat : Lucretia lattice files for positron damping ring
* AT/FACET2p_DR.mat : AT positron damping ring model
* QDDSQ*.[xsif|bmad] : Sliced quad-sextupole arc models for DR
* BA/BD.[xsif|bmad] : Sliced models for bend-quad-sextupole arc magnets for DR
* GPT: Matlab scripts and supporting ImpactT data files to generate FACET-II electron injector lattice up to exit of L0a
