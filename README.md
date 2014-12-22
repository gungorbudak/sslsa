Structural Superimposition of Local Sequence Alignment

CENG465 HOMEWORK, 2014, Middle East Technical University

GUNGOR BUDAK, 1720184, Bioinformatics MSc

A program which finds out whether a local sequence
alignment of two protein sequences also implies structural
similarity of the aligned parts

The script takes 4 parameters:
1. First PDB ID – no default value, must be given
2. Second PDB ID – no default value, must be given
3. Gap open penalty (a nonpositive number) – default: -10
4. Gap extension penalty (a nonpositive number) – default: -5

Requirements: Python 2.7.x, NumPy 1.9.1, BioPython 1.64

How to run (on Windows CMD) and example output:

D:\projects\sslsa>python sslsa.py 1xyx 1xu0
Downloading PDB structure '1xyx'...
Downloading PDB structure '1xu0'...
LGGYMLGSAMSRPMIHFGNDWEDRYYRENMYRYPNQVYYRPV---DQYSNQNNFVHDCVNITIKQHTVTTTTKGEN--FT
ETDVKMMERVVEQMCVTQYQKES
IGGYMLGNAVGRMSYQFNNPMESRYYNDYYNQMPNRVY-RPMYRGEEYVSEDRFVRDCYNMSVTEYIIKPAEGKNNSELN
QLDTTVKSQIIREMCITEYRRGS
Alignment length: 103
RMSD: 2.87702543237

or

D:\projects\sslsa>python sslsa.py 1xyx 1xu0 -5 -5
Downloading PDB structure '1xyx'...
Downloading PDB structure '1xu0'...
LGGYMLGSAMSRPMIHFGNDWEDRYYRENMY-RYPNQVYYRPV---DQYSNQNNFVHDCVNITIKQHTVTTTTKGENFTE
-T--DVKMMERVVEQMCVTQYQKES
IGGYMLGNAVGRMSYQFNNPMESRYYND-YYNQMPNRVY-RPMYRGEEYVSEDRFVRDCYNMSVTEYIIKPA-EGKNNSE
LNQLDTTVKSQIIREMCITEYRRGS
Alignment length: 105
RMSD: 2.96977891523

Interpretation:

Similar sequences yield similar lower RMSD because they are actually similar
but since several amino acids can be grouped according to their structural
properties, even if the sequences might not be very similar we might get lower
RMSD due to structural biochemical similarity. This shows that structure is more
conserved during the course of evolution, whereas sequences are subject to changes
due to the mechanisms of evolution such as mutation, genetic drift, gene flow.