try:
    from glycopeptidepy._c.structure.constants import Configuration
    constants = Configuration()
except ImportError:
    from argparse import Namespace
    constants = Namespace()

# Tokenizer Constants
constants.N_TERM_DEFAULT = "H"
constants.C_TERM_DEFAULT = "OH"

# N-glycan Sequon detection
constants.ALLOW_MODIFIED_ASPARAGINE = False

# Sequence Fragment Constants
constants.FRAG_OFFSET = 1
constants.PARTIAL_HEXNAC_LOSS = True
constants.EXCLUDE_B1 = True
