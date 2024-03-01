1. Note: a local bmad instance is necessary for this conversion.  The conda bmad / tao does not have all the required files.
  - e.g. cd ~/forks/bmad-ecosystem/
  - e.g. . util/dist_source_me
  - Check existance of: $ACC_ROOT_DIR/util_programs/mad_to_bmad/mad8_to_bmad.py
2. Obtain updated release from CVS, either by mounting a CVS repo, or checking out CVS locally.
  - `cvs update` from local CVS/optics directory
3. Create a local branch of lcls-lattice to stage the update.
  - git checkout master
  - git checkout -b DDMMYYYY_conversion
4. Copy files from CVS into facet2-lattice/mad (omit mad/CVS, bmad/CVS directory from copy)
  - cp ~/CVS/optics/etc/lattice/facet2/* facet2-lattice/MAD
    - Omit directories.
  - copy CVS/optics/script/elementdevices.dat facet2-lattice/MAD
5. Obtain updated lcls_elements.csv
  - Go to [https://oraweb.slac.stanford.edu/apex/slacprod/f?p=116:600](https://oraweb.slac.stanford.edu/apex/slacprod/f?p=116:600)
  - Click on Actions and Download
  - Save file to lcls-lattice/bmad/conversion/from_oracle
6. Open a command prompt for jupyter notebook
  - export FACET2_LATTICE=YOUR LOCAL FACET2 BRANCH
  - check that $LCLS_LATTICE is set to the location of the DDMMYY_conversion branch of lcls-lattice.
    - Yes, LCLS_LATTICE should be set to the location of the lcls-lattice repo, not the location of the facet2-lattice repo.
  - Start a jupyter notebook session in facet2-lattice directory
7. Within jupyter notebook cd to bmad/conversion and open facet2_slac_to_bmad.ipynb
  - Check hard-coded paths in slac_to_bmad.ipynb match facet2-lattice repo conversion branch.
  - Run all cells in facet2_slac_to_bmad.ipynb
8. Within jupyter nodebook cd to bmad/conversion/device_mapping and open device_mapping.ipynb
  - Check LCLS_LATTICE environment variable points to conversion branch of lcls-lattice repo.
  - Run all cells.
  - This generates the lcls-lattice/bmad/master/*_devicenames.bmad files
7. Within jupyter nodebook cd to bmad/conversion/device_mapping and open device_mapping.ipynb
  - Check LCLS_LATTICE environment variable points to conversion branch of lcls-lattice repo.
  - Run all cells.
  - This generates the lcls-lattice/bmad/master/*_devicenames.bmad files
  - If tao fails to start, it may be necessary to comment out the call to FACET2e_devicenames from f2_elec.lat.bmad
8. Compare optics between Bmad and Lucretia or mad8s.
  - Common sources of discrepancy
    - extraneous quad settings in the top level bmad file (e.g. f2_elec.lat.bmad)
    - lcavity tracking at low energy.
    - Suggestion:  Obtain the beginning optics by matching at MRK0F or another location > 100 MeV.
