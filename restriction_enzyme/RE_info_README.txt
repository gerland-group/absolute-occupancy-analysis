middle_shift =  middle of the recognition sequence
site_shift = how far left the cut site is from the middle for the upper strand

example HindIII
5' - A|A G C T T - 3'  (plus strand)
3' - T T C G A|A - 5'  (minus strand)

If the 3' end is shorter than the 5' end after cutting, the 3' end is elongated to match the 5' end (by polymerase).
If the 3' end is longer  that the 5' end after cutting, the 3' end is digested  to match the 5' end.

Thus site_shift is also how many bases to add from the middle of the sequence when looking for reads that were cut.

cut_sites = recognition_site start on the plus strand + middle_shift

cut_fragment_position = cut_sites - site_shift (for reads starting at the cut site using the genomic position)
cut_fragment_position = cut_sites + site_shift (for reads ending at the cut site using the genomic position)

