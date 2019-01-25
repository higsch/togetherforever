# Collection of dependencies

* gffread
* threeFrameTranslator.py
* awk
* OpenMS (in PATH!)
* seqtk
* piDeepNet (R with libraries)

# Run command
nextflow search.nf --gtfs "data/*.gtf" --genome resources/Homo_sapiens.GRCh37.dna.primary_assembly.fa --canonical resources/Homo_sapiens.GRCh38.ENS90.pep.all.fa -resume
