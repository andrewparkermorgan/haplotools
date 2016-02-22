# haptools

Some rudimentary Python classes and scripts for dealing with collections of haplotype segments, as for indviduals from multi-founder experimental populations like the [Collaborative Cross](http://csbio.unc.edu/CCstatus) and [Diversity Outbred](http://do.jax.org).

Dependencies:

* Python >= 2.7
* `bx-python` library (https://pypi.python.org/pypi/bx-python/0.7.3)

# Usage

```
import genome

# manually construct genome as mosaic of haplotype segments
mouse = genome.Genome("CC001/Unc", sex = genome.SEX_MALE)
mouse.add_segment("chr1",0,100000,0, genome.PHASE_PATERNAL, "A")
mouse.add_segment("chr1",100000,200000,0, genome.PHASE_PATERNAL, "C")
mouse.add_segment("chr1",0,200000, genome.PHASE_MATERNAL, "B")

# iterate over haplotype segments
for block in mouse:
	print(block)

# query haplotypes at a locus
mouse.get_genotype_at("chr1", 500)

```