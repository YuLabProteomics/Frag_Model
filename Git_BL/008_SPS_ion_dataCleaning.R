library(rawrr)# read raw ms data
library(PSMatch) # fragment the peptide
library(limma)
library(dplyr)
library(tidyr)
library(ggbreak)
library(funcTools)
library(RColorBrewer)
source('000_functions.R')



############### 1) filtering SPS-ions with specific threshold

## load all scan information
load('../002_Rdata_CSV/004_ScanInfor_all.Rdata')

fl = list.files('../002_Rdata_CSV/', pattern = 'run2F')
# f_numbers <- sub(".*(F[0-9]+)\\.Rdata", "\\1", fl)
numbers <- as.numeric(gsub(".*F([0-9]+).*", "\\1", fl))
# Sort the files based on the extracted numbers
fl <- fl[order(numbers)]

for(mi in seq_along(fl)){
  load(paste0('../002_Rdata_CSV/', fl[mi]))
  
  SPSions_data = lapply(matched_ion_all, function(x){
    df = x
    
    ### 1)add the signal-to-noise ratio
    df$SN_ratio = df$exp_intensity/df$exp_noises
    # normalize intensity to base peak within scannumber
    N_intensity = df%>%
      group_by(scanNumber) %>%
      reframe(N_intensity = exp_intensity/max(exp_intensity)) %>%
      ungroup()
    
    df$N_intensity = N_intensity$N_intensity
      
    # stat for total ion
     scan_total_ion = df %>% 
      group_by(scanNumber) %>%
      reframe(int_total_ion = sum(exp_intensity),
              total_ion_SN_ratio = sum(SN_ratio))
    
    ### 2) keep the top 200 (intensity) MS2 ions ################################# sort the intensity and take the top 200 rows
   # After the ion matching, we do not have any scan contains more than 200 spectra
    print(fl[[mi]])
    cat('the maximum number of matched ion is ', max(df$Num_matched_ion))
    
    ### 3)exclude b1/y1 ions
    exclude = c('b1', 'y1', 'y1_')
    tmp_NOb1y1 = subset(df, !(ion %in% exclude))
    
    ### 4) extract K-ending y-ions
    yions = subset(tmp_NOb1y1, grepl('^y[\\*\\_]?$', type))
    k_ending_yions <- subset(yions, grepl("K$", seq))
    
    # extract b-ions
    bions = subset(tmp_NOb1y1, grepl('^b[\\*\\_]?$', type))
    # combine b and y ions table
    b_y_ions = rbind(k_ending_yions, bions)
    
    
    ### 5)keep ions with exp_mz >= 400
    b_y_ions = b_y_ions[b_y_ions$exp_mz >= 400, ]
    
    
    
    ### 6)excluding ions m/z in the range of -50- +5 of the MS1 precursor m/z
    ## get scan information and precursor m/z from raw
    scan_raw = ScanInfor_all[[mi]]
    scan_raw = scan_raw[, c('scan','StartTime','charge',
                            'monoisotopicMz', 'type', 'scanEnergy')]
    scan_raw$precusorMZ_low = scan_raw$monoisotopicMz - 50
    scan_raw$precusorMZ_high = scan_raw$monoisotopicMz + 5
    # re-name of the scan raw
    colnames(scan_raw) = c("scan", 'StartTime','PrecusorCharge', "monoisotopicMz", 
                           "scantype","scanEnergy",     
                           "precusorMZ_low","precusorMZ_high")
    
    # Merge scan_raw into b_y_ions based on scanNumber
    merged_data <- merge(b_y_ions, scan_raw, by.x = "scanNumber", by.y = "scan")
    
    # Check if exp_mz is within the range
    index = merged_data$exp_mz >= merged_data$precusorMZ_low & merged_data$exp_mz <= merged_data$precusorMZ_high
    #
    filtered_b_y_ions = merged_data[!index, ]
    filtered_b_y_ions$NCE = paste0(filtered_b_y_ions$scantype,'_', filtered_b_y_ions$scanEnergy)
    
    # stat for total sps-ion
    scan_SPS_ion = filtered_b_y_ions %>%
      group_by(scanNumber) %>%
      reframe(Num_total_sps = length(scanNumber),
              int_total_sps = sum(exp_intensity),
              total_sps_SN_ratio = sum(SN_ratio))
    
    ### 7) stat for sps-ions m/z > precursor m/z
    index = filtered_b_y_ions$exp_mz > filtered_b_y_ions$monoisotopicMz
    df_tmp = filtered_b_y_ions[index, ]
    #
    bigger_sps = df_tmp %>% 
      group_by(scanNumber) %>%
      reframe(Num_bigger_sps = n(),
              int_bigger_sps = sum(exp_intensity), 
              bigger_sps_SN_ratio = sum(SN_ratio))
    
    
    ### 8)filter top10 strong ions in each scan
    top10_ions <- filtered_b_y_ions %>%
      group_by(scanNumber) %>%                  # Group by scanNumber
      arrange(scanNumber, desc(exp_intensity)) %>% # Sort by exp_intensity within each scanNumber group
      slice_head(n = 10)                        # Select the top 10 rows for each scanNumber group
    top10_ions = as.data.frame(top10_ions)
    # top10_ions$NCE = paste0(top10_ions$scantype,'_', top10_ions$scanEnergy)

    
    ### 9) do the statistics for the number of sps-ions and the total intensity/SN_ratio
    stat_SPSIon = top10_ions %>%
      group_by(scanNumber) %>%
      summarize( top10_SN_ratio = sum(SN_ratio),
                 PrecusorCharge = unique(PrecusorCharge),
                 PrecusorMZ = unique(monoisotopicMz),
                 Num_exp_ion = unique(Num_exp_ion),
                 Num_matched_ion = unique(Num_matched_ion),
                 RTtime = unique(StartTime),
                 NCE = unique(NCE),
                 Num_top10SPS_ions = n(),
                 int_top10_sps = sum(exp_intensity))
    stat_SPSIon = as.data.frame(stat_SPSIon)
    stat_SPSIon = merge(stat_SPSIon, scan_total_ion)
    stat_SPSIon = merge(stat_SPSIon, scan_SPS_ion)
    stat_SPSIon = merge(stat_SPSIon, bigger_sps, by = 'scanNumber', all.x = TRUE)
    
   
    return(list(top10_SPSions = top10_ions, stat_SPSIon = stat_SPSIon, total_SPSions = filtered_b_y_ions))
    
  })
  save(SPSions_data, file = paste0('../002_Rdata_CSV/008_SPS_ions_f', mi, '.Rdata'))
}




# ### 8) assign the first amino-acid a weight(frequency * 100)
# tmt_bearingAcid = strsplit2(top10_ions$seq, split = '')[,1]
# tmt_weight = (table(tmt_bearingAcid)/nrow(top10_ions) * 100)
# 
# 
# ### 9) encoding the ion-seq using one-hot encoding
# # Apply the encode_sequence function to each sequence in the seq column
# encoded_df <- t(sapply(top10_ions$seq, function(seq) {
#   encode_sequence(seq = seq, tmt_weight = tmt_weight)
# }, simplify = "data.frame"))
# rownames(encoded_df) = NULL
# 
# 
# # Combine the original data with the encoded sequences
# encoded_top10_ions <- cbind(top10_ions, encoded_df)
# 
# # View the result
# head(encoded_top10_ions)



############### 2) plot  of top10_SN_ratio in different NCE 

fl = list.files('../002_Rdata_CSV/', pattern = '008_SPS_ions')


numbers <- as.numeric(gsub(".*f([0-9]+).Rdata", "\\1", fl))
# Sort the files based on the extracted numbers
fl <- fl[order(numbers)]
n = c(3,9,16,20)
# for(ni in seq_along(fl)){
for(ni in n){
  load(paste0('../002_Rdata_CSV/',fl[ni]))
  stat_SPSIon = SPSions_data$matched_ion_IdenComet$stat_SPSIon
  SPSions = SPSions_data$matched_ion_IdenComet$top10_SPSions
  p = ggplot(data = stat_SPSIon, mapping = aes(NCE, top10_SN_ratio )) +
    # geom_point() +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA)+
    scale_y_continuous(limits = c(0,4500))+
    labs(
      y = 'Sum top10 SN_ratio',
      title = paste0('Sum SN_ratio of top10 SPS-ions f',ni ))+
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5)  # Center the title
    )
  print(p)
  ggsave(filename = paste0('../004_Figures/008_spsIons_SNRatio_f', ni, '.png'), plot = p,
         width = 5, height = 3)
}






############### 3) plot the number of SPS ions in different NCE 
fl = list.files('../002_Rdata_CSV/', pattern = '008_SPS_ions')

numbers <- as.numeric(gsub(".*f([0-9]+).Rdata", "\\1", fl))
# Sort the files based on the extracted numbers
fl <- fl[order(numbers)]

for(ni in seq_along(fl)){
  load(paste0('../002_Rdata_CSV/',fl[ni]))
  stat_SPSIon = SPSions_data$matched_ion_IdenComet$stat_SPSIon
  SPSions = SPSions_data$matched_ion_IdenComet$top10_SPSions
  
  # Create the SPSions_NCE list with adjusted tables
  SPSions_NCE <- sapply(unique(stat_SPSIon$NCE), function(x) {
    # Create a table for the current NCE value
    table_data <- table(subset(stat_SPSIon, NCE == x)$NumSPS_ions)
    
    # Ensure all numbers from 1 to 10 are represented, setting missing numbers as 0
    all_numbers <- 1:10  # Define the range of numbers you expect (1 to 10)
    
    # Create a result table with all numbers from 1 to 10, filling missing values with 0
    result_table <- setNames(as.integer(sapply(all_numbers, function(n) table_data[as.character(n)])), as.character(all_numbers))
    
    # Replace NA (missing numbers) with 0
    result_table[is.na(result_table)] <- 0
    
    return(result_table)
  })
  
  SPSions_NCE = as.data.frame(SPSions_NCE)
  SPSions_NCE$Num_sps_ions = 1:10
  df_long = reshape2::melt(data = SPSions_NCE, measure.vars = colnames(SPSions_NCE)[-ncol(SPSions_NCE)])
  df_long$variable = factor(df_long$variable, 
                            levels = c('cid_25', 'cid_27', 'cid_30', 'cid_32', 'cid_35',
                                       'hcd_25', 'hcd_27', 'hcd_30', 'hcd_32', 'hcd_35'))
  df_long$Num_sps_ions = factor(df_long$Num_sps_ions, levels = unique(df_long$Num_sps_ions))
  
  
  df_long$NCE = strsplit2(df_long$variable, split = '_')[,1]
  df_long$NCEValue = strsplit2(df_long$variable, split = '_')[,2]
  df_long$NCE_variable = paste0(df_long$NCEValue, "_", df_long$variable)
  df_long$NCE_variable = factor(df_long$NCE_variable, 
                                levels = c("25_cid_25", "27_cid_27", "30_cid_30", "32_cid_32", "35_cid_35" ,
                                           "25_hcd_25","27_hcd_27", "30_hcd_30", "32_hcd_32", "35_hcd_35"))
  

  p = ggplot(data = df_long, aes(Num_sps_ions, value, fill = NCE, group = NCE_variable)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Bars next to each other
    scale_fill_brewer(palette = "Set2")+  # Apply friendly colors
    scale_y_break(c(150, 1500))+ # , scales = 'free'
    labs(
      x = '# SPS ions',
      y = 'Number of scans',
      title = paste0('Number of SPS-ions f',ni))+
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      legend.title = element_text(size = 10),  # Adjust legend title size
      legend.text = element_text(size = 8)     # Adjust legend text size
    ) +
    facet_wrap(~NCEValue, nrow = 2, ncol = 5)  # Facet by NCEValue (each variable gets its own plot)
  print(p)
  ggsave(filename = paste0('../004_Figures/008_Number_spsIons_f', ni, '.png'), plot = p,
         width = 8, height = 4)
}





############### 4) plot the precursor charge 
fl = list.files('../002_Rdata_CSV/', pattern = '008_SPS_ions')
numbers <- as.numeric(gsub(".*f([0-9]+).Rdata", "\\1", fl))
# Sort the files based on the extracted numbers
fl <- fl[order(numbers)]


for(ni in seq_along(fl)){
  load(paste0('../002_Rdata_CSV/',fl[ni]))
  stat_SPSIon = SPSions_data$matched_ion_IdenComet$stat_SPSIon
  SPSions = SPSions_data$matched_ion_IdenComet$top10_SPSions
  
  # Create the PrecursorCharge list with adjusted tables
  PrecursorCharge <- sapply(unique(stat_SPSIon$NCE), function(x) {
    # Create a table for the current NCE value
    table_data <- table(subset(stat_SPSIon, NCE == x)$PrecusorCharge)
    # Ensure all numbers from 1 to 10 are represented, setting missing numbers as 0
    all_numbers <- 2:6  # 
    # Create a result table with all numbers from 1 to 10, filling missing values with 0
    result_table <- setNames(as.integer(sapply(all_numbers, function(n) table_data[as.character(n)])), as.character(all_numbers))
    # Replace NA (missing numbers) with 0
    result_table[is.na(result_table)] <- 0
    return(result_table)
  })
  
  PrecursorCharge = as.data.frame(PrecursorCharge)
  PrecursorCharge$Charge = rownames(PrecursorCharge)
  df_long = reshape2::melt(data = PrecursorCharge, measure.vars = colnames(PrecursorCharge)[-ncol(PrecursorCharge)])
  df_long$variable = factor(df_long$variable, 
                            levels = c('cid_25', 'cid_27', 'cid_30', 'cid_32', 'cid_35',
                                       'hcd_25', 'hcd_27', 'hcd_30', 'hcd_32', 'hcd_35'))
  df_long$Charge = factor(df_long$Charge, levels = unique(df_long$Charge))
  
  df_long$NCE = strsplit2(df_long$variable, split = '_')[,1]
  df_long$NCEValue = strsplit2(df_long$variable, split = '_')[,2]
  df_long$NCE_variable = paste0(df_long$NCEValue, "_", df_long$variable)
  df_long$NCE_variable = factor(df_long$NCE_variable, 
                                levels = c("25_cid_25", "27_cid_27", "30_cid_30", "32_cid_32", "35_cid_35" ,
                                           "25_hcd_25","27_hcd_27", "30_hcd_30", "32_hcd_32", "35_hcd_35"))
  
 
  p = ggplot(data = df_long, aes(Charge, value, fill = NCE, group = NCE_variable)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Bars next to each other
    scale_fill_brewer(palette = "Set2")+  # Apply friendly colors
    labs(
      x = 'Precursor charge',
      y = 'Number of scans',
      title = paste0('SPS-ions Precursor Charge f', ni)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      legend.title = element_text(size = 10),  # Adjust legend title size
      legend.text = element_text(size = 8)     # Adjust legend text size
    ) +
    facet_wrap(~NCEValue, nrow = 2, ncol = 5)  # Facet by NCEValue (each variable gets its own plot)
  print(p)
  
  ggsave(filename = paste0('../004_Figures/008_PreCharge_spsIons_f',ni, '.png'), plot = p,
         width = 6, height = 3)
}



