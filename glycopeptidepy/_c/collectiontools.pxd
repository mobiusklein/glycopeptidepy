from glycopeptidepy._c.count_table cimport CountTable


ctypedef fused mapping_t:
    object
    CountTable
    # dict

cpdef list descending_combination_counter(mapping_t counter)
