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