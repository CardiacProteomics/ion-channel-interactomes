# ECG plotter adaption
This is an adaptation of Niek Verweij's ECG plotter developed for the paper **Outlining cardiac ion channel protein interactors and their signature in the human electrocardiogram** by Maurya et al., published in Nature Cardiovascular Research in 2023: https://doi.org/10.1038/s44161-023-00294-y 
  
The original version of the tool can be found here: https://github.com/niekverw/ecgenetics  
and it's original publication here:  https://doi.org/10.1016/j.cels.2020.08.005  

The purpose of this adaptation to use ECG plotter and automated fashion in a
pipeline to analyze multiple genes at once. It extracts the most significant
SNP and the corresponding p-value for each gene in the input for a region of
interest in the ECG. This way one can quickly check many genes and see which
ones have a strong influence on the ECG. The script additionaly determines the
LD block of each gene looks for SNPs within LD instead of a specified size window
arround the gene as in the original ECG plotter version.  
  
  
# Detailed description
The script takes a csv file with a column named Gene
with human gene names as input, for each gene it determines the rsid of the
first and the last snp within the gene, then uses those to determine the LD
block upstream of the first SNP and downstream of the last SNP using the
LDlink API. An LD block is calculated as the region in which MAF > 0.05 and 
D prime > 0.5 (default). Next, all the SNPs in that region are extracted and the SNP
with the most significant absolute p value in the area of interest of the ECG
is identified.  
  
The output consits of four columns:  
**Gene**: input gene list  
**rsid**: rsID of the most significant SNP in the LD block of the gene and the area of interest of
the ECG plot  
**minP**: minimum p-value of the most significant SNP in the LD block
of the gene and the area of interest of the ECG plot  
**Distance**: distance
(either positive or negative) from the ends of the gene (or from the snps
closest to the ends) to the reported SNP.  
  
The script also creates a pdf file with all either all ECG-plots for each
of the SNPs of each input gene LD block (can be very slow, not reccomended) or only for the most
significant SNP per gene (fast(er)).  

# Instructions
To run our adapted version you need to download the original version and the accompanying data files (note: only runs on linux systems). If the original verison is running, just swap out the helpers.R file with the one provided here and adapt the file paths in standalone_final_v_ion-channel-interactomes.R script to your system and input files.
