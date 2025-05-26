Guide to files:

sars_spike_protein_assembled.fna: 
The final assembled spike protein sequence

raw_input/sars_spike_protein_raw_reads.fastq:
The simulated raw shotgun sequencing data (contains sequencing errors, in case you want to do your own data cleaning)

trimmed_input/sars_spike_protein_raw_reads.fastq:
The simulated data, after quality control, with most of the sequencing errors and low-quality reads removed. (contains the quality scores, in case you want to use them in your algorithm)

trimmed_input/sars_spike_protein_raw_reads.fna:
The simulated data, after quality control, just the clean DNA sequences. (Does not contain the quality scores – simplest to work with).

