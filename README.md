
Main code and data files
1) data_preprocessing.ipynb - Does data preprocessing and stores them as pickle files Requires python 2.7
-- This generates preprocessed data in "processed_data" folder in current directory
2) PSN.ipynb - The proposed Phrase Structure Network working on the pickle data files in "processed_data" folder (Requires python 3.0+)

Original Dataset:
GENIA dataset is present in the folder "TEES/data/corpora" in the form of xml files

TEES code base:
It is located in TEES folder in current directory


PS: Since the preprocessed pickle files are also saved, PSN.ipynb can be directly run. If you would like to regenerate the preprocessed pickle file you can run data_preprocessing.ipynb which reads original GENIA dataset and generates preprocessed pickle file.
