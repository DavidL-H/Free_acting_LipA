# Free_acting_LipA
Ph.D. Chapter III - Efforts towards developing a free-acting Lipoyl Synthase
For easy over-view, the files and scripts will be described here in the order in which they appear as figures in the Chapter

## Figure 4: Heatmap of the lipoic acid biosynthetic genes present in bacterial genomes.
This data originates from Microbesonline (Dehal, P. S., Joachimiak, M. P., Price, M. N., Bates, J. T., Baumohl, J. K., Chivian, D., … Arkin, A. P. (2010). MicrobesOnline: an integrated portal for comparative and functional genomics. Nucleic Acids Research, 38(Database), D396–D400. https://doi.org/10.1093/nar/gkp919)

### LipA_B_lplA_LipM_LipL_gcvH orthologs.xlsx
Main data file, containing the exported orthologs and homologs of lipoic acid biosynthetic genes in the genes available in Microbesonline
### MO_bacteria_metadata.csv
Metadata file, giving information on each of the genomes in LipA_B_lplA_LipM_LipL_gcvH orthologs.xlsx, including the phylum.
### Heatmap of lipoic acid biosynthetic genes.R
R script for handling, cleaning and visualizing the data, resulting in Figure 4 heatmap

## Figure 5: Unrooted guide tree with neighborhood-joining, based on multiple sequence alignment of the Pfam BPL_LplA_LipB (PF03099) family
This is a simple MSA tree, generated as describes in the Metods. The files used for generating it are as follows:
### uniprot_pf03099_reviewed_bacteria.fa
Fasta file containing the sequences used for the guide-tree.

## Figure 8: difference in residue conservation score for complementing vs. non-complementing LipAs & FigureS2: mtLipA 3d strucutre with more conserved residues highlighted
### hetLipAseq_ConservationAnalysis.R
R-script for handling and combining alignment and ConSurf data files, visualizing the data for the figure.
### alignment_all_LipAs.clustal_num
Clustal-Omega MSA of the LipA sequences used for the analysis
### Allalign_to_mtLipA_alignscores.tab
Output from ConSurf, giving a conservation score and confidence interval of each residue (rooted in the mtLipA residues)
### CompAlign_to_mtLipA_alignscores.tab
Output from ConSurf, giving a conservation score and confidence interval of each residue (rooted in the mtLipA residues), only complementing LipAs
### Conservation_CompvsNoncomp_mtLipA.csv
Final Data-file indicating the differences in conservations of LipA residues, used for plots in Figure 8.

### mtLipA5EXK_functionalMoreConsResidues_highlighted.pse
PyMOL session file showing the structure and highlighted residues used for figure S2

## Figure 9: epPCR mutability of ecLipA study
A deep-sequence analysis of missense mutations of ecLipA allowed for in vivo functionality. Data originates from NGS data supplied by Eurofins Genomics as describes in the Methods
### epPCR_NGS.R
R-script for the handling, cleaning and analysis of the epPCR data
### Selected epPCR ecoLipA deep sequencing.csv
Output of the CLC genomics Workbench basic variant analysis tool for the selected epPCR LipA library
### Unselected epPCR ecoLipA deep sequencing.csv
Output of the CLC genomics Workbench basic variant analysis tool for the un-selected epPCR LipA library

## Figure 10: Combining conservation and epPCR model scores
Generating a logistical regression scoring function for predicting the chance of a heterologous LipA complementing 

#
