from argparse import Namespace

constants = Namespace()

# Tokenizer Constants
constants.MOD_BEGIN = "("
constants.MOD_END = ")"
constants.GLYCAN_BEGIN = "{"
constants.GLYCAN_END = "}"
constants.N_TERM_DEFAULT = "H"
constants.C_TERM_DEFAULT = "OH"

# N-glycan Sequon detection
constants.ALLOW_MODIFIED_ASPARAGINE = False
# Sequence Fragment Constants
constants.FRAG_OFFSET = 1
constants.PARTIAL_HEXNAC_LOSS = True
constants.EXCLUDE_B1 = True
