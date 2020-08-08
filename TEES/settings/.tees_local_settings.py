
# Edit these settings to configure TEES. A variable must have a value 
# other than None for it to be usable. This file is interpreted as
# a Python module, so Python code can be used.

DATAPATH = os.getenv('DATAPATH'); # Main directory for datafiles

# Tools
SVM_MULTICLASS_DIR = DATAPATH + '/tools/SVMMultiClass'
BANNER_DIR = DATAPATH + '/tools/BANNER'
GENIA_SENTENCE_SPLITTER_DIR = DATAPATH + '/tools/geniass'
RUBY_PATH = "ruby" # Command to run Ruby (used only by the GENIA Sentence Splitter)
BLLIP_PARSER_DIR = DATAPATH + '/tools/BLLIP/dmcc-bllip-parser-558adf6'
MCCLOSKY_BIOPARSINGMODEL_DIR = DATAPATH + '/tools/BLLIP/biomodel'
STANFORD_PARSER_DIR = DATAPATH + '/tools/stanford-parser-2012-03-09'

# Data
CORPUS_DIR = DATAPATH + '/corpora'
MODEL_DIR = DATAPATH + '/models'
BIONLP_EVALUATOR_DIR = DATAPATH + '/tools/evaluators'
BIONLP_EVALUATOR_GOLD_DIR = DATAPATH + '/tools/evaluators/gold'
