Setup
=====
Before launching jupyter notebook run "source set_env.sh" being in the directory containing this file. "set_env.sh" is tested on mac machine and should work on linux too on a bash shell
The code requires both python 2.7 and python 3.0+ to be present in the system and might also require ruby to be present

Main code and data files
========================
1) data_preprocessing.ipynb - Does data preprocessing and stores them as pickle files Requires python 2.7
-- This generates preprocessed data in "processed_data" folder in current directory
2) PSN.ipynb - The proposed Phrase Structure Network working on the pickle data files in "processed_data" folder (Requires python 3.0+)

Original Dataset and TEES code
==============================
GENIA dataset is present in the folder "TEES/data/corpora" in the form of xml files

TEES code base(Subsection of download from https://github.com/jbjorne/TEES) is located in TEES folder in current directory

Other info
==========
PS: Since the preprocessed pickle files are also saved, PSN.ipynb can be directly run. If you would like to regenerate the preprocessed pickle file you can run data_preprocessing.ipynb which reads original GENIA dataset and generates preprocessed pickle file.
