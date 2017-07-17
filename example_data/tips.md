In order to extract all loci from a GFF3, a command like the following can work. Tweak the *awk* parts as needed if your GFF3 is setup differently. This step is required when providing a **subset_list** for the extract_sequences.cwl step.


```
grep '\tgene\t' PlasmoDB-24_Pfalciparum3D7.gff | rev | cut -f1 | rev | awk -F';' '{print $1}' | awk -F':' '{print $1}' | awk -F'=' '{print $2}'
```