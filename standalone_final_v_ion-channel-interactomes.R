# ECG plotter standalone script for ion channel interactomes
#
# Date: 2021_03_19 Authors: Mads Kurtzhals & Konstantin Kahnert
#
#
# This is an adaptation of Niek Verweij's ECG plotter
# (https://github.com/niekverw/ecgenetics)
#
# The purpose of this adaptation to use ECG plotter and automated fashion in a
# pipeline to analyze multiple genes at once. It extracts the most significant
# SNP and the corresponding p-value for each gene in the input for a region of
# interest in the ECG. This way one can quickly check many genes and see which
# ones have a strong influence on the ECG. The script additionaly determines the
# LD block of each gene looks for SNPs within LD instead of a specified size window
# arround the gene as in the original ECG plotter version.
#
#
# Detailed description: The script takes a csv file with a column named Gene
# with human gene names as input, for each gene it determines the rsid of the
# first and the last snp within the gene, then uses those to determine the LD
# block upstream of the first SNP and downstream of the last SNP using the
# LDlink API. An LD block is calculated as the region in which MAF > 0.05 and D'
# > 0.5 (default). Next, all the SNPs in that region are extracted and the SNP
# with the most significant absolute p value in the area of interest of the ECG
# is identified.
#
# The output consits of four columns: Gene: input gene list rsid: rsID of the
# most significant SNP in the LD block of the gene and the area of interest of
# the ECG plot minP: minimum p-value of the most significant SNP in the LD block
# of the gene and the area of interest of the ECG plot Distance: distance
# (either positive or negative) from the ends of the gene (or from the snps
# closest to the ends) to the reported SNP
#
# The script also creates a pdf file with all either all ECG-plots for each
# of the SNPs of each input gene LD block (slow!) or only for the most
# significant SNP per gene (fast(er)).
#


# set parameter for running ECG plotter
file_path_input = "/media/sf_E_DRIVE/Data/Cardiac_ionchannels/IP dat_human_orthologs_Gja1.tsv"

search_dataset = "GW" # can be "tophits" or "GW"
phenotype = "stretch" # can be "unadjusted" or "stretch"
ecg_region = 1:500 # can be any slice between 1 and 500 (for choosing a region of interest of the ECG)
dprime_thresh = 0.5 # default is 0.5, the higher the cutoff the more stringent i.e. the smaller the LD block
print_all_ecg_plots = "one" # can be "all", "one", or "none"
max_query = 1000
pops="EUR" # Set population of the group the data should be retrieved from (LD link parameter)
include_LD = T

# set working directory
setwd("/media/sf_E_DRIVE/Software/ECG_plotter/ecgenetics/")

# load required files
source("global.R")
source("helpers.R")
source("helpers.heatmap.R")
source("helpers.regionalplot.R")
source("get_nearest_gene.r")

# load libraries
library(ggpubr)

# define required functions
flipeffects <- function(row) {
  if (is.na(sum(row)) == TRUE){
    return(row)
  } else if (abs(min(row)) > max(row)){
    #print("flip")
    return(-1*row)
  } else{
    #print("don't flip")
    return(1*row)
  }
}

rescale_amplitudes<-function(col){
  colrescale<-col
  colrescale[colrescale>0]=(col[col>0] - 0)/(max(col) - 0)
  colrescale[colrescale<0]=-1*(abs(col[col<0]) - 0)/(max(col) - 0)
  return(colrescale)
}

rescale_amplitudes<-function(col){
  colrescale<-col
  colrescale[colrescale>-log10(5E-08)]=-log10(5E-08)
  colrescale[colrescale<log10(5E-08)]=log10(5E-08)
  return(colrescale)
}


make_heatmap_plot <- function(df_snp_p = data$df_snp_p,
                              vct_snp_info = data$df_snp_info,
                              file) {
  
  #colorscale = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)
  colorscale = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  mx <- data.matrix(df_snp_p)
  rownames(mx)<-1:nrow(mx)
  
  #mx <- t(apply(mx, 1, flipeffects))
  mx <- t(apply(mx, 1, rescale_amplitudes))
  
  rownames <-paste0(vct_snp_info$SNP," (",input_file$`Human gene names`,")")
  #rownames <- stringr::str_wrap(rownames,width=30)
  #mx <- t(apply(mx, 1, abs))
  
  hm <- heatmaply(mx,labRow=rownames,
                  showticklabels =c(FALSE,TRUE),
                  color=colorscale, 
                  #limits = c(-30, 30),
                  Colv=FALSE,
                  dendrogram="none",
                  fontsize_row=7,
                  margins = c(60,400,40,20)) %>% ggplotly(height=max(500,25*nrow(mx)))
  
  return(hm)
}


############################### Main script ####################################
# read input file
input_file = as.data.frame(fread(file_path_input))

# create dataframes for pvalues and snp info to fill 
df_snp_p_heatmap <- data.frame()
vct_snp_info_heatmap = data.frame()

# Check if new columns exist. If not create them
if (("SNP_within_LD_block" %in% colnames(input_file)) == F) {
  input_file[,"SNP_within_LD_block"] <- NA
  input_file[,"Gene_name_of_SNP"] <- NA
  input_file[,"minP_within_LD_block"] <- NA
  input_file[,"Distance"] <- NA
  input_file[,"SNP_within_LD_block_region"] <- NA
  input_file[,"Gene_name_of_SNP_region"] <- NA
  input_file[,"minP_within_LD_block_region"] <- NA
  input_file[,"Distance_region"] <- NA
  input_file[,"SNP_within_gene"] <- NA
  input_file[,"minP_within_gene"] <- NA
  input_file[,"Comment"] <- NA
}


# loop over all genes in input file
for (i in 1:nrow(input_file)) {
  
  # Save the input gene to a variable
  input_gene = input_file$`Human gene names`[i]
  
  # Get index of input gene in input file
  input_gene_index = match(input_gene, input_file$`Human gene names`)
  
  # extract data from database file and catch error if gene is not in dataset
  # set error flag to false
  flag_error = TRUE
  
  tryCatch({
    # Create query with search window = 0 
    query_0 <- process_user_input(entry = input_gene,
                                  mapping.proteincoding = mapping.proteincoding,
                                  window = 0,
                                  subset = search_dataset,
                                  phenotype = phenotype,
                                  max_query = max_query,
                                  standalone = T)
    
    
    # Get index of input gene in input file
    input_gene_index = match(input_gene, input_file$`Human gene names`)# turn query into tabix query
    tabix_query_0 <- get_tabix_query(query_0,df.static.pos,df.static.rsid)
    
    if (search_dataset == "tophits"){
      start_end_pos <- extract_multiple_variants(tabix_query_0,
                                                 f.data_p=paste0(datadir,"/tophits_data/",phenotype,
                                                                 ".logP.outfile.tsv.gz.tophits.gz"),
                                                 f.data_beta=paste0(datadir,"/tophits_data/",phenotype,
                                                                    ".BETA.outfile.tsv.gz.tophits.gz"),
                                                 f.data_se=paste0(datadir,"/tophits_data/",phenotype,
                                                                  ".SE.outfile.tsv.gz.tophits.gz"),
                                                 f.data.index=paste0(datadir,"/tophits_data/",phenotype,
                                                                     ".logP.outfile.index.tsv.gz.tophits.gz"))
      # set error flag to false
      flag_error = FALSE
      
    } else if (search_dataset == "GW"){
      start_end_pos <- extract_multiple_variants(tabix_query_0,
                                                 f.data_p=paste0(datadir,"/full_data_combined/",phenotype,
                                                                 ".logP.outfile.tsv.gz"),
                                                 f.data_beta=paste0(datadir,"/full_data_combined/",phenotype,
                                                                    ".BETA.outfile.tsv.gz"),
                                                 f.data_se=paste0(datadir,"/full_data_combined/",phenotype,
                                                                  ".SE.outfile.tsv.gz"),
                                                 f.data.index=paste0(datadir,"/full_data_combined/",phenotype,
                                                                     ".logP.outfile.index.tsv.gz"))
      # set error flag to false
      flag_error = FALSE
      
    }
  }, error=function(e) {
    message(paste0(input_gene, " not found in ", phenotype, " dataset. Continue with next gene of input list."))
  })
  
  # check if error flag is TRUE
  if (flag_error == TRUE){
    
    # Get index of input gene in input file
    input_gene_index = match(input_gene, input_file$`Human gene names`)

    # print error in comment columns 
    input_file$Comment[input_gene_index] <- "Gene not in dataset"
    
    # add p values and snp_info to df_snp_p_heatmap and vct_snp_info_heatmap
    tmp <- matrix(ncol=500,nrow=1)
    colnames(tmp) <- colnames(start_end_pos[["df_snp_p"]]) # causes an error if the first gene in the input data is not present in the dataset
    df_snp_p_heatmap <- rbind(df_snp_p_heatmap, tmp)

    tmp <- matrix(ncol=15,nrow=1)
    colnames(tmp) <- colnames(cbind(start_end_pos[["df_snp_info"]], maxLogP_ecg_region=NA, minP_ecg_region=NA)) # causes an error if the first gene in the input data is not present in the dataset
    vct_snp_info_heatmap <- rbind(vct_snp_info_heatmap, tmp)
    
    # if no SNP for current gene was found skip to next one
    next
  }
  
  # Get smallest p-value
  min_p = min(abs(start_end_pos$df_snp_info$minP))
  
  # get index of min p-value
  idx_min_p = match(as.character(min_p), 
                    start_end_pos[["df_snp_info"]][["minP"]])
  
  # Get rsid of min p-value
  min_p_rsid = start_end_pos[["df_snp_info"]][["SNP"]][idx_min_p]
  
  # add snp and p value of most significant SNP_within_gene to input data frame
  input_file$SNP_within_gene[input_gene_index] <- min_p_rsid
  input_file$minP_within_gene[input_gene_index] <- min_p
  
  if (include_LD == TRUE){
    
    # catch error in case the max BP SNP is not in LDlink database  
    # and in case it isn't take the second to max and so on
    valid_snp = FALSE
    flag_no_ld = FALSE
    n = 0
    
    while (valid_snp == FALSE){
      
      tryCatch({
        if (length(sort(as.numeric(start_end_pos$df_snp_info$BP)))-n <= 0){
          input_file$Comment[input_gene_index] <- "No LD info"
          flag_no_ld = TRUE
          break
        }
        
        # Get snps with highest value in BP position
        highest_snp_dist = sort(as.numeric(start_end_pos$df_snp_info$BP), 
                                decreasing = FALSE)[length(sort(as.numeric(start_end_pos$df_snp_info$BP)))-n]
        
        # Get index of highest bp snp
        i_highest = match(as.character(highest_snp_dist), start_end_pos[["df_snp_info"]][["BP"]])
        
        # Get rsid's
        rs_highest = start_end_pos$df_snp_info$SNP[i_highest]
        
        # Create the two commands for the LD data retrieval
        ld_command_highest = paste0("https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=", 
                                    rs_highest,"&pop=",pops,"&r2_d=r2&token=70722d8a3958")
        
        # Create df's from the LD data
        df_ld_highest = as.data.frame(data.table::fread(ld_command_highest))
        
        # set valid_snp to true to break the loop and continue
        valid_snp = TRUE
        
      }, error=function(e) {
        
        # print error message
        message(paste(rs_highest, "not in LDlink database. Trying next SNP."))
        
      }, finally={
        
        # increase number of SNPs to skip by one
        n = n + 1
        
      })
    }
    
    # catch error in case the min BP SNP is not in LDlink database
    # and in case it isn't take the second to min and so on
    valid_snp = FALSE
    n = 0
    
    while (valid_snp == FALSE){
      
      tryCatch({
        if (length(sort(as.numeric(start_end_pos$df_snp_info$BP)))-n <= 0){
          input_file$Comment[input_gene_index] <- "No LD info"
          flag_no_ld = TRUE
          break
        }
        # Get snps with lowest and lowest value in BP'position
        lowest_snp_dist = sort(as.numeric(start_end_pos$df_snp_info$BP),
                               decreasing = TRUE)[length(sort(as.numeric(start_end_pos$df_snp_info$BP)))-n]
        
        # Get index of these snps
        i_lowest = match(as.character(lowest_snp_dist), start_end_pos[["df_snp_info"]][["BP"]])
        
        # Get rsid's
        rs_lowest = start_end_pos$df_snp_info$SNP[i_lowest]
        
        # Create the two commands for the LD data retrieval
        ld_command_lowest = paste0("https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=",
                                   rs_lowest,"&pop=",pops,"&r2_d=r2&token=70722d8a3958")
        
        # Create df's from the LD data
        df_ld_lowest = as.data.frame(data.table::fread(ld_command_lowest))
        
        # set valid_snp to true to break the loop and continue
        valid_snp = TRUE
        
      }, error=function(e) {
        
        # print error message
        message(paste(rs_lowest, "not in LDlink database. Trying next SNP."))
        
      }, finally={
        
        # increase number of SNPs to skip by one
        n = n + 1
      })
    }
    
    if (flag_no_ld == TRUE){
      # add p values and snp_info to df_snp_p_heatmap and vct_snp_info_heatmap
      df_snp_p_heatmap <- rbind(df_snp_p_heatmap,
                                start_end_pos[["df_snp_p"]][idx_min_p,])
      vct_snp_info_heatmap <- rbind(vct_snp_info_heatmap,
                                    cbind(start_end_pos[["df_snp_info"]][idx_min_p,],
                                          maxLogP_ecg_region=NA, minP_ecg_region=NA))
      
      # continue with next gene
      next
    }
    
    # Exclude all snps with MAF < 0.05
    df_ld_highest_MAF = subset(df_ld_highest, MAF > 0.05)
    df_ld_lowest_MAF = subset(df_ld_lowest, MAF > 0.05)
    
    # Only keep snps that have D' > dprime threshold (default is 0.5)
    df_ld_highest_Dp = subset(df_ld_highest_MAF, Dprime > dprime_thresh)
    df_ld_lowest_Dp = subset(df_ld_lowest_MAF, Dprime > dprime_thresh)
    
    # Only keep snps that have D' > dprime threshold (default is 0.5)
    df_ld_highest_Dp = subset(df_ld_highest_Dp, Dprime != 1)
    df_ld_lowest_Dp = subset(df_ld_lowest_Dp, Dprime != 1)
    
    # Keep only positive distances for highest and negative distances for lowest SNP
    df_ld_highest_Dp = subset(df_ld_highest_Dp, Distance > 0)
    df_ld_lowest_Dp = subset(df_ld_lowest_Dp, Distance < 0)
    
    # Get distance for the lowest and the highest with regards to the snp's in the
    # start_end_pos data.
    ld_highest_dist = max(df_ld_highest_Dp$Distance)
    ld_lowest_dist = min(df_ld_lowest_Dp$Distance)
    
    # check if any of the distances is -inf, i.e. no SNP outside of the gene is 
    # left. In case is is, set it to 0
      if (rs_lowest == rs_highest){

        # print error in comment columns
        input_file$Comment[input_gene_index] <- "First and last SNP in gene are the same after LD query"

        # add p values and snp_info to df_snp_p_heatmap and vct_snp_info_heatmap
        df_snp_p_heatmap <- rbind(df_snp_p_heatmap,
                                  start_end_pos[["df_snp_p"]][idx_min_p,])
        vct_snp_info_heatmap <- rbind(vct_snp_info_heatmap,
                                      cbind(start_end_pos[["df_snp_info"]][idx_min_p,],
                                            maxLogP_ecg_region=NA, minP_ecg_region=NA))

        # if no SNP for current gene was found skip to next one
        next
      }

      if (ld_highest_dist == -Inf){
        ld_highest_dist <- 0
        input_file$Comment[input_gene_index] <- "LD issue downstream"
      }
    
    if (ld_lowest_dist == Inf){
      ld_lowest_dist <- 0
      input_file$Comment[input_gene_index] <- "LD issue upstream"
    }
    
    # Write the distances (or BP localization) of the furthest snps from the ld
    # search as a region of the chromosome
    region_search = paste0(start_end_pos$df_snp_info$CHR[1], ":", 
                           as.character(lowest_snp_dist + ld_lowest_dist), "-", 
                           as.character(highest_snp_dist + ld_highest_dist))
    
    # Create a query with this region as input:
    query_region <- process_user_input(entry = region_search,
                                       mapping.proteincoding = mapping.proteincoding,
                                       window = 0,
                                       subset = search_dataset,
                                       phenotype = phenotype,
                                       max_query = max_query,
                                       standalone = T)
    
    # turn query into tabix query
    tabix_query_region <- get_tabix_query(query_region,df.static.pos,df.static.rsid)
    
    # extract data from database file
    if (search_dataset == "tophits"){
      data_region <- extract_multiple_variants(tabix_query_region,
                                               f.data_p=paste0(datadir,"/tophits_data/",phenotype,
                                                               ".logP.outfile.tsv.gz.tophits.gz"),
                                               f.data_beta=paste0(datadir,"/tophits_data/",phenotype,
                                                                  ".BETA.outfile.tsv.gz.tophits.gz"),
                                               f.data_se=paste0(datadir,"/tophits_data/",phenotype,
                                                                ".SE.outfile.tsv.gz.tophits.gz"),
                                               f.data.index=paste0(datadir,"/tophits_data/",phenotype,
                                                                   ".logP.outfile.index.tsv.gz.tophits.gz")
      )
    } else if (search_dataset == "GW"){
      data_region <- extract_multiple_variants(tabix_query_region,
                                               f.data_p=paste0(datadir,"/full_data_combined/",phenotype,
                                                               ".logP.outfile.tsv.gz"),
                                               f.data_beta=paste0(datadir,"/full_data_combined/",phenotype,
                                                                  ".BETA.outfile.tsv.gz"),
                                               f.data_se=paste0(datadir,"/full_data_combined/",phenotype,
                                                                ".SE.outfile.tsv.gz"),
                                               f.data.index=paste0(datadir,"/full_data_combined/",phenotype,
                                                                   ".logP.outfile.index.tsv.gz")
      )
    }
  } else{
    data_region <- start_end_pos
  }

################### Get most signif SNP using full ECG region ##################
  # Get smallest p-value
  smallest_p = min(abs(data_region$df_snp_info$minP))

  # get index of smalles p-value
  idx_smallest_p = match(as.character(smallest_p),
                         data_region[["df_snp_info"]][["minP"]])

  # Get rsid
  smallest_p_rsid = data_region[["df_snp_info"]][["SNP"]][idx_smallest_p]

  # Get distance of this snp and write tables
  smallest_p_bp = as.integer(data_region[["df_snp_info"]][["BP"]][idx_smallest_p])

  # get closest gene to most signif snp
  smallest_p_gene = data_region[["df_snp_info"]][["Gene"]][idx_smallest_p]

  # Determine if snp is inside the input gene (then get "0") or distance from
  # each end of gene
  if (smallest_p_bp < lowest_snp_dist) {
    smallest_p_dist = smallest_p_bp - lowest_snp_dist
  } else if (smallest_p_bp > highest_snp_dist) {
    smallest_p_dist = smallest_p_bp - highest_snp_dist
  } else {
    smallest_p_dist = 0
  }

  # Print rsid, p-value and distance to input file
  input_file$SNP_within_LD_block[input_gene_index] <- smallest_p_rsid
  input_file$Gene_name_of_SNP[input_gene_index] <- smallest_p_gene
  input_file$minP_within_LD_block[input_gene_index] <- smallest_p
  input_file$Distance[input_gene_index] <- smallest_p_dist

################# filtering for ECG region of interest #########################
  # Get smallest p-value of the selected region of the ECG
  list_p = c()
  
  for (snp in 1:nrow(data_region$df_snp_info)) {
    if (abs(min(as.numeric(unlist(data_region$df_snp_p[snp,ecg_region])))) > max(as.numeric(unlist(data_region$df_snp_p[snp,ecg_region])))){
      list_p[snp] = -max(abs(as.numeric(unlist(data_region$df_snp_p[snp,ecg_region]))))
      
    } else if(abs(min(as.numeric(unlist(data_region$df_snp_p[snp,ecg_region])))) <= max(as.numeric(unlist(data_region$df_snp_p[snp,ecg_region])))){
      list_p[snp] = max(abs(as.numeric(unlist(data_region$df_snp_p[snp,ecg_region]))))
    }
  }
  
  # Write these p-vaues to the data
  data_region$df_snp_info$maxLogP_ecg_region <- unlist(list_p)
  data_region$df_snp_info$minP_ecg_region <- (10^-abs(data_region$df_snp_info$maxLogP_ecg_region)) * sign(data_region$df_snp_info$maxLogP_ecg_region)

  # Get smallest p-value
  #mallest_p = min(abs(data_region$df_snp_info$minP_ecg_region))

  if (min(abs(subset(data_region$df_snp_info, minP_ecg_region < 0)[["minP_ecg_region"]])) < min(subset(data_region$df_snp_info, minP_ecg_region > 0)[["minP_ecg_region"]])){
    smallest_p = -min(abs(subset(data_region$df_snp_info, minP_ecg_region < 0)[["minP_ecg_region"]]))
    
  } else if (min(abs(subset(data_region$df_snp_info, minP_ecg_region < 0)[["minP_ecg_region"]])) >= min(subset(data_region$df_snp_info, minP_ecg_region > 0)[["minP_ecg_region"]])){
    smallest_p = min(abs(subset(data_region$df_snp_info, minP_ecg_region > 0)[["minP_ecg_region"]]))
  }
  
  # get index of smalles p-value
  idx_smallest_p = match(as.character(smallest_p),
                         data_region[["df_snp_info"]][["minP_ecg_region"]])

  # Get rsid
  smallest_p_rsid = data_region[["df_snp_info"]][["SNP"]][idx_smallest_p]

  # Get distance of this snp and write tables
  smallest_p_bp = as.integer(data_region[["df_snp_info"]][["BP"]][idx_smallest_p])

  # get closest gene to most signif snp
  smallest_p_gene = data_region[["df_snp_info"]][["Gene"]][idx_smallest_p]

  # Determine if snp is inside the input gene (then get "0") or distance from
  # each end of gene
  if (smallest_p_bp < lowest_snp_dist) {
    smallest_p_dist = smallest_p_bp - lowest_snp_dist
  } else if (smallest_p_bp > highest_snp_dist) {
    smallest_p_dist = smallest_p_bp - highest_snp_dist
  } else {
    smallest_p_dist = 0
  }

  # Print rsid, p-value and distance to input file
  input_file$SNP_within_LD_block_region[input_gene_index] <- smallest_p_rsid
  input_file$Gene_name_of_SNP_region[input_gene_index] <- smallest_p_gene
  input_file$minP_within_LD_block_region[input_gene_index] <- data_region[["df_snp_info"]][["minP_ecg_region"]][idx_smallest_p]
  input_file$Distance_region[input_gene_index] <- smallest_p_dist
  
  # add p values and snp_info to df_snp_p_heatmap and vct_snp_info_heatmap
  df_snp_p_heatmap <- rbind(df_snp_p_heatmap, 
                            data_region[["df_snp_p"]][idx_smallest_p,])
  vct_snp_info_heatmap <- rbind(vct_snp_info_heatmap, 
                                data_region[["df_snp_info"]][idx_smallest_p,])

################################## ECG plots ###################################
  # Check if plots are supposed to be created or not
  if (print_all_ecg_plots == "all"){
    
    # set correct ecg_stats depending on phenotype setting
    if (phenotype == "stretch"){
      ecg_stats=df_ecg_stretch
      
    } else if (phenotype == "unadjusted"){
      ecg_stats=df_ecg_unadjusted
    }
    
    # Create ECG for each snp in this region
    # Create a pdf file, in which all plots are saved
    filename = paste0("/media/sf_E_DRIVE/Software/ECG_plotter/ecgenetics/",
                      "plots/ECG_plots_", input_gene,".pdf")

    # create PDF file
    pdf(filename)

    # For each SNP create an ecg plot
    for (u in data_region$df_snp_info$SNP) {
      
      # Get index of the snp
      idx = match(u, data_region[["df_snp_info"]][["SNP"]])
      
      # Check if SNP is found in data
      if (is.na(idx) == F) {
        
        # create ECG profile plot
        ecg_plot <- make_ecg_plot(vct_snp_p=data_region$df_snp_p[idx,],
                                  vct_snp_beta=data_region$df_snp_beta[idx,],
                                  vct_snp_se=data_region$df_snp_se[idx,],
                                  vct_snp_info=data_region$df_snp_info[idx,],
                                  df_ecg_stats=ecg_stats,
                                  plot_adjusted_means = F,
                                  plot_beta = F,
                                  invert=FALSE)
        
        # save plot as pdf
        filename = paste0("/media/sf_E_DRIVE/Software/ECG_plotter/ecgenetics/",
                          "plots/ECG_plot_", input_gene,".pdf")
        ggsave(file=filename, plot=ecg_plot, height = 5.5 , width = 7)
      }
    }
    
  } else if (print_all_ecg_plots == "one"){
    
    # set correct ecg_stats depending on phenotype setting
    if (phenotype == "stretch"){
      ecg_stats=df_ecg_stretch
      
    } else if (phenotype == "unadjusted"){
      ecg_stats=df_ecg_unadjusted
    }
    
    # Create ECG for each snp in this region
    # Get smallest p-value
    smallest_p = min(abs(data_region$df_snp_info$minP_ecg_region))
    
    # Match p-value to rsid
    idx_smallest_p = match(as.character(smallest_p), 
                        data_region[["df_snp_info"]][["minP_ecg_region"]])
    
    # Check if SNP is found in data
    if (is.na(idx_smallest_p) == F) {
      
      # create ECG profile plot
      ecg_plot <- make_ecg_plot(vct_snp_p=data_region$df_snp_p[idx_smallest_p,],
                                vct_snp_beta=data_region$df_snp_beta[idx_smallest_p,],
                                vct_snp_se=data_region$df_snp_se[idx_smallest_p,],
                                vct_snp_info=data_region$df_snp_info[idx_smallest_p,],
                                df_ecg_stats=ecg_stats,
                                plot_adjusted_means = T,
                                plot_beta = F,
                                invert=FALSE)
      
      # save plot as pdf
      filename = paste0("/media/sf_E_DRIVE/Software/ECG_plotter/ecgenetics/",
                        "plots/ECG_plot_", input_gene,".pdf")
      ggsave(file=filename, plot=ecg_plot, height = 5.5 , width = 7)
    }

    
  } else{
    # do nothing
  }
  
  # print progress to console
  print(paste("Gene",input_gene, "done!"))
}

# write results to file
file_path_output = paste0(substr(file_path_input, 0, nchar(file_path_input)-4), 
                          "_output.tsv")

write.table(input_file, file=file_path_output, sep='\t',row.names =FALSE)

# rename rows
input_type = "Gene" # needs to be changed when SNPs are possible as input
rownames(df_snp_p_heatmap) <- input_file[[input_type]]
rownames(vct_snp_info_heatmap) <- input_file[[input_type]]

# add orca to r PATH 
#Sys.setenv(PATH=paste0(Sys.getenv("PATH"), ":/usr/bin/orca"))

# set path to save for heatmap
file_path_output = paste0(substr(file_path_input, 0, nchar(file_path_input)-4), 
                          "_heatmap.pdf")

# file_path_output = "heatmap.png"
                          
# create heatmap with most significant SNP of each gene in input file
hm <- make_heatmap_plot(df_snp_p = df_snp_p_heatmap,
                        vct_snp_info = vct_snp_info_heatmap)

# show heatmap
hm

# save heatmap
hm$width <- 1000
hm$height <- 90*length(input_file)
export(hm, file = file_path_output)
