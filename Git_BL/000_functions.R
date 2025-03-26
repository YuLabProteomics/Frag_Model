library(rawrr)# read raw ms data
library(PSMatch) # fragment the peptide
library(dplyr)
library(ggplot2)
library(scales)

########### 1) generating the fragmentation for identified peptide
Fragment_generator = function(peptide= peptide,
                              charge = charge, 
                              tolerance = 0.01,
                              modification = c(Nterm = 295.1896, C = 57.02146)){
  
  pep = regmatches(peptide, regexpr("(?<=\\.).+?(?=\\.)", peptide, perl = TRUE))
  # check if the peptide ending with K or not
  if (substr(pep, nchar(pep), nchar(pep)) == "K") {
    modification = c(modification, Cterm = 295.1896)
  }
  # check position of Methionine and add the modification
  star_position = regexpr("M\\*", pep)
  if(star_position != -1){
    modification = c(modification, M = 15.994915)
    pep <- gsub("\\*", "", pep)
  }
  
  z = 1:(charge-1)
  FragIon = calculateFragments(pep, z = z,
                               type = c("y", "b"), 
                               modifications = modification,
                               neutralLoss = NULL)
  FragIon$mzRange_low =  FragIon$mz - tolerance
  FragIon$mzRange_high = FragIon$mz + tolerance
  return(FragIon)
}

##### consider neutral loss version
Fragment_generator_NeutralLoss = function(peptide= peptide,
                              charge = charge, 
                              tolerance = 0.01,
                              modification = c(Nterm = 295.1896, C = 57.02146)){
  
  pep = regmatches(peptide, regexpr("(?<=\\.).+?(?=\\.)", peptide, perl = TRUE))
  # check if the peptide ending with K or not
  if (substr(pep, nchar(pep), nchar(pep)) == "K") {
    modification = c(modification, Cterm = 295.1896)
  }
  # check position of Methionine and add the modification
  star_position = regexpr("M\\*", pep)
  if(star_position != -1){
    modification = c(modification, M = 15.994915)
    pep <- gsub("\\*", "", pep)
  }
  
  z = 1:(charge-1)
  FragIon = calculateFragments(pep, z = z,
                               type = c("y", "b"), 
                               modifications = modification,
                               neutralLoss = defaultNeutralLoss())
  FragIon$mzRange_low =  FragIon$mz - tolerance
  FragIon$mzRange_high = FragIon$mz + tolerance
  return(FragIon)
}




########### 2) spectra-spectrum matching function

Frag_match <- function(theoretical_ions = fragment, 
                       experimental_spectrum = experimental_spectrum) {
  matched_ions <- list()
  
  # Loop through each experimental ion
  for (i in 1:nrow(experimental_spectrum)) {
    exp_mz <- experimental_spectrum$mz[i]
    # print(exp_mz)
    exp_intensity <- experimental_spectrum$intensity[i]
    exp_noises = experimental_spectrum$noises[i]
    # Find matching theoretical ions within the tolerance range
    matches <- theoretical_ions %>%
     filter(mzRange_low <= exp_mz & exp_mz <= mzRange_high)
    
    if (nrow(matches) > 0) {
      matched_ions[[length(matched_ions) + 1]] <- list(
        experimental_mz = exp_mz,
        experimental_intensity = exp_intensity,
        experimental_noises = exp_noises,
        Num_exp_ion = nrow(experimental_spectrum),
        theoretical_matches = matches
      )
    }
  }
  
  return(matched_ions)
}



########### 3)one-hot encode ion sequences and assign the first amino-acid a weight
# encode_sequence <- function(seq = 'GNAGGLHH', 
#                             aa_list = c('A','C','D', 'E', 'F', 'G',
#                                         'H', 'I', 'K', 'L' ,
#                                         'M', 'N' , 'P', 'Q' , 'R' , 
#                                         'S', 'T' ,'V' , 'W' , 'Y'),
#                             tmt_weight) {
#   
#   # Create a vector of zeros for the one-hot encoding of the sequence
#   encoding <- rep(0, length(aa_list))
#   
#   # Loop through each amino acid in the sequence
#   for (i in seq_along(aa_list)) {
#     # Check if the amino acid is present in the sequence
#     if (any(strsplit(seq, '')[[1]] == aa_list[i])) {
#       # If it's the first amino acid in the sequence, use its weight from tmt_weight
#       if (substr(seq, 1, 1) == aa_list[i]) {
#         encoding[i] <- tmt_weight[aa_list[i]]  # Set the frequency-based weight for the first amino acid
#       } else {
#         encoding[i] <- 1  # Regular encoding for others
#       }
#     }
#   }
#   
#   names(encoding) = aa_list
#   return(encoding)
# }
# 
# 
# 
# 
# ########### 4)One-hot encoding function
# one_hot_encode_precCharge <- function(charge, max_charge = 6) {
#   encoding <- rep(0, max_charge - 1)
#   encoding[charge - 1] <- 1
#   return(encoding)
# }

########### 3) encode peptide 
# Define the alphabet mapping
ALPHABET <- c("A" = 1, "C" = 2, "D" = 3, "E" = 4, "F" = 5, 
              "G" = 6, "H" = 7, "I" = 8, "K" = 9, "L" = 10, 
              "M" = 11, "N" = 12, "P" = 13, "Q" = 14, "R" = 15, 
              "S" = 16, "T" = 17, "V" = 18, "W" = 19, "Y" = 20, 
              "m" = 21)
# "M(ox)" = 21)

# Function to encode a peptide sequence with correct handling of "M(ox)"
encode_peptide <- function(sequence, max_length = 20) {
  # Replace "M*" with "m" in the sequence
  sequence <- gsub("M\\*", "m", sequence)
  
  # Split the sequence into individual components
  peptide_parts <- unlist(strsplit(sequence, "(?=[A-Za-z()])", perl = TRUE))
  
  # Remove empty elements resulting from the split
  peptide_parts <- peptide_parts[peptide_parts != ""]
  
  # Map each component to its numeric value using ALPHABET
  encoded_peptide <- sapply(peptide_parts, function(x) ALPHABET[x])
  
  # Replace any invalid characters with 0
  encoded_peptide[is.na(encoded_peptide)] <- 0
  
  # # If the sequence is shorter than max_length, pad it with 0s at the end
  # if (length(encoded_peptide) < max_length) {
  #   encoded_peptide <- c(encoded_peptide, rep(0, max_length - length(encoded_peptide)))
  # }
  res = rep(0, max_length)
  res[1:length(encoded_peptide)] = encoded_peptide
  
  # # If the sequence is longer than max_length, truncate it
  # encoded_peptide <- encoded_peptide[1:max_length]
  
  return(unname(res))
}


########### 4) One-hot encoding function
one_hot_encode_precCharge <- function(charge, max_charge = 6) {
  encoding <- rep(0, max_charge - 1)
  encoding[charge - 1] <- 1
  return(encoding)
}



########### 5)plot for spectra and highlight the matched ions
specPlot = function(fraction = 'f5', scanNum = 29538){
  ## load raw data
  rawfile = list.files('../001_rawData/raw/', pattern = paste0(fraction, '.raw'))
  spec = readSpectrum(rawfile = paste0('../001_rawData/raw/', rawfile), scan = scanNum)
  scanSpectra = data.frame(mz = spec[[1]]$mZ, intensity = spec[[1]]$intensity, noises = spec[[1]]$noises) 
  scanSpectra$SN_ratio = scanSpectra$intensity/scanSpectra$noises
  
  ## load the identified data from the Comet to get the peptide information
  peptideFile = list.files('../001_rawData/raw/', pattern = paste0(fraction, '.csv'))
  peptideInfo = read.csv(paste0('../001_rawData/raw/', peptideFile), header = T,
                         stringsAsFactors = F, check.names = F)
  peptide = subset(peptideInfo, ScanF == scanNum)$Peptide
  pep = regmatches(peptide, regexpr("(?<=\\.).+?(?=\\.)", peptide, perl = TRUE))
  
  
  ## load matched ion data
  fl = list.files('../002_Rdata_CSV/', pattern = paste0('^007_matched_ion_run2', toupper(fraction), '\\.Rdata$'))
  load(paste0('../002_Rdata_CSV/', fl))
  # matchedIon =  matched_ion_all$matched_ion_IdenComet
  matchedIon = do.call(rbind, matched_ion_all)
  matchedIon = subset(matchedIon, scanNumber == scanNum)
  matchedIon$SN_ratio = matchedIon$exp_intensity/matchedIon$exp_noises
  matchedIon$z_label = paste0(matchedIon$ion,'+', matchedIon$z)
  
  ## load SPS-ions data
  load(paste0('../002_Rdata_CSV/008_SPS_ions_', fraction, '.Rdata'))
  top10_SPSIon = subset(rbind(SPSions_data$matched_ion_IdenComet$top10_SPSions, SPSions_data$matched_ion_Failed$top10_SPSions), scanNumber == scanNum)
  top10_SPSIon$SN_ratio = -top10_SPSIon$SN_ratio
  top10_SPSIon$z_label = paste0(top10_SPSIon$ion,  '+', top10_SPSIon$z)
  
  
  ## prepare data for plot
  # ymin <- min(c(scanSpectra$SN_ratio, top10_SPSIon$SN_ratio))
  ymax <- max(c(scanSpectra$SN_ratio, top10_SPSIon$SN_ratio))
  p = ggplot(data = scanSpectra, aes(x = mz, y = SN_ratio))+
    geom_segment(aes(x = mz, xend = mz, y = 0, yend = SN_ratio), color = "black") +  # Vertical lines from y = 0 to intensity
    scale_y_continuous(labels = scientific)+
    
    
    # add matched ions
    geom_segment(data = matchedIon, aes(x = mz, xend = exp_mz, y = 0, yend = SN_ratio), 
                 color = ifelse(grepl("^b", matchedIon$ion), "blue", "red")) +  # Vertical lines from y = 0 to intensity
    geom_text(data = matchedIon, aes(x = exp_mz, y = SN_ratio, label = z_label), 
              color = ifelse(grepl("^b", matchedIon$ion), "blue", "red"), 
              size = 3, hjust = 0.5, vjust = -0.5) +
    
    # add top10 sps ions
    geom_segment(data = top10_SPSIon, aes(x = mz, xend = exp_mz, y = 0, yend = SN_ratio), 
                 color = ifelse(grepl("^b", top10_SPSIon$ion), "blue", "red")) +  # Vertical lines from y = 0 to intensity
    # ylim(0, min(-scanSpectra$SN_ratio) - 1) +
    expand_limits(y = c(-ymax, ymax + 1)) +
    
    geom_text(data = top10_SPSIon, aes(x = exp_mz, y = SN_ratio, label = z_label), 
              color = ifelse(grepl("^b", top10_SPSIon$ion), "blue", "red"), 
              size = 4, hjust = 0.5, vjust = 1) +
    scale_y_continuous(labels = scientific)+
    
    # Add upper-right corner label
    annotate("text", x = Inf, y = Inf,
             label = paste0('ScanNum: ', scanNum), hjust = 2, vjust = 2, size = 4) +
    annotate("text", x = Inf, y = Inf,
             label = unique(top10_SPSIon$NCE), hjust = 4.5, vjust = 4, size = 4) +
    annotate("text", x = Inf, y = Inf,
             label = paste0('# total peaks: ', unique(top10_SPSIon$Num_exp_ion)), 
             hjust = 1.9, vjust = 6, size = 4) +
    annotate("text", x = Inf, y = Inf,
             label = paste0('# matched peaks: ', unique(top10_SPSIon$Num_matched_ion)), 
             hjust = 1.6, vjust = 8, size = 4) +
    
    
    # Add lower-right corner label
    annotate("text", x = Inf, y = Inf,
             label = "Top10 SPS ions", hjust = 1.7, vjust = 40, size = 4) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Sum SN_ratio: ", round(abs(sum(top10_SPSIon$SN_ratio)),digits = 2)), 
             hjust = 1.1, vjust = 42, size = 4) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Num SPS ions: ", nrow(top10_SPSIon)), 
             hjust = 1.1, vjust = 44, size = 4) +
    
    labs(title = paste0('Spectrum for ', pep),
         y = 'SN ratio') +
    
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5)# Center the title
    )
  print(p)
}






########### 6)plot for spectra and highlight the matched ions
specPlot_compare = function(fraction = 'f4', scanNum = 29538){
  ## load raw spectra
  rawfile = list.files('../001_rawData/raw/', pattern = paste0(fraction, '.raw'))
  spec = readSpectrum(rawfile = paste0('../001_rawData/raw/', rawfile), scan = scanNum)
  scanSpectra = data.frame(mz = spec[[1]]$mZ, intensity = spec[[1]]$intensity, noises = spec[[1]]$noises) 
  scanSpectra$N_intensity = scanSpectra$intensity/max(scanSpectra$intensity)
  
 
  ## load SPS-ions data
  load(paste0('../002_Rdata_CSV/008_SPS_ions_', fraction, '.Rdata'))
  total_SPSIon = subset(SPSions_data$matched_ion_IdenComet$total_SPSions, scanNumber == scanNum)
  # keep only b and y ion to plot
  total_SPSIon = subset(total_SPSIon, type %in% c('b', 'y'))
  # total_SPSIon$N_intensity = total_SPSIon$exp_intensity/max(total_SPSIon$exp_intensity)
  total_SPSIon$z_label = paste0(total_SPSIon$ion, '+', total_SPSIon$z)
  
  
  ## load the Prosit spectra
  load(paste0('../002_Rdata_CSV/0011_Prosit_Prediction_', fraction, '.Rdata'))
  prosit_spectra = subset(predictions, scanNumber == scanNum)
  prosit_spectra$intensities = -prosit_spectra$intensities
  
  # Process mz values based on ion type and peptide sequence ending
  prosit_spectra <- prosit_spectra %>%
    mutate(
      # Extract the peptide sequence by removing modifications (e.g., [UNIMOD:737]-LVTDLTK[UNIMOD:737] -> LVTDLTK)
      peptide_sequence = gsub("\\[UNIMOD:[0-9]+\\](-?)([A-Za-z]+)\\[UNIMOD:[0-9]+\\]", "\\2", peptide_sequences),
      
      # Extract the charge from the annotation column (e.g., y1+1 -> 1, y2+1 -> 2, etc.)
      charge = as.numeric(gsub(".*\\+([0-9]+)", "\\1", annotation)),
      
      # Determine if the ion is 'b' or 'y' and if it's a K-ending peptide
      ion_type = substr(annotation, 1, 1),  # Extracts 'b' or 'y' from annotation
      is_k_ending = grepl("K$", peptide_sequence),  # Checks if peptide ends with K
      
      # Adjust mz values based on ion type and whether it's a K-ending peptide
      mz_adjusted = case_when(
        ion_type == "b" ~ mz + (66 / charge),  # Adjust all b-ions
        ion_type == "y" & is_k_ending ~ mz + (66 / charge),  # Adjust y-ions if peptide ends with K
        TRUE ~ mz  # Do not adjust mz for other cases
      )
    )
  
  
  ## load the Prosit SPSion
  load(paste0('../002_Rdata_CSV/0011_Prosit_SPSIon_', fraction, '.Rdata'))
  peptide = unique(SPS_ions[SPS_ions$scanNumber == scanNum, ]$Peptide)
  prosit_SPSion = SPS_ions[SPS_ions$scanNumber == scanNum, ]
  prosit_SPSion$intensities = -prosit_SPSion$intensities
  
  
  
  ## prepare data for plot
  ggplot(data = scanSpectra, aes(x = mz, y = N_intensity))+
    geom_segment(aes(x = mz, xend = mz, y = 0, yend = N_intensity), color = "black") +  # Vertical lines from y = 0 to intensity
    # scale_y_continuous(labels = scientific)+
    
    # add sps ions
    geom_segment(data = total_SPSIon, aes(x = exp_mz , xend = exp_mz, y = 0, yend = N_intensity), 
                 color =  "blue", size = 1)+   # Vertical lines from y = 0 to intensity
    
    geom_text(data = total_SPSIon, aes(x = exp_mz, y = N_intensity, label = z_label),
              # color = ifelse(grepl("^b", matchedIon$ion), "blue", "red"),
              size = 3, hjust = 0.5, vjust = -0.5) +
    
    # add prosit ions
    geom_segment(data = prosit_spectra, aes(x = mz_adjusted, xend = mz_adjusted, y = 0, yend = intensities), 
                 color = 'black')+  # Vertical lines from y = 0 to intensity
    
    # add prosit SPS ions
      geom_segment(data = prosit_SPSion, aes(x = mz.y, xend = mz.y, y = 0, yend = intensities), 
                   color = 'orange', size = 1) +
  geom_text(data = prosit_SPSion, aes(x = mz.y, y = intensities, label = annotation ),
            # color = ifelse(grepl("^b", top10_SPSIon$ion), "blue", "red"),
            size = 3, hjust = 0.5, vjust = 1) +
  
    # Add upper-right corner label
    annotate("text", x = max(scanSpectra$mz), y = Inf,
             label = paste0('ScanNum: ', scanNum), hjust = 2, vjust = 2, size = 4) +
    annotate("text", x = max(scanSpectra$mz), y = Inf,
             label = unique(total_SPSIon$NCE), hjust = 4, vjust = 4, size = 4) +
    annotate("text", x = max(scanSpectra$mz), y = Inf,
             label = paste0('# total peaks: ', unique(total_SPSIon$Num_exp_ion)), 
             hjust = 2, vjust = 6, size = 4) +
    annotate("text", x = max(scanSpectra$mz), y = Inf,
             label = paste0('# SPS ion peaks: ', nrow(total_SPSIon)), 
             hjust = 1.6, vjust = 8, size = 4) +
    annotate("text", x = max(scanSpectra$mz), y = Inf,
             label = paste0('# Sum SPS ion intensity: ', format(round(sum(total_SPSIon$exp_intensity),1), scientific = TRUE, digits = 2)), 
             hjust = 1, vjust = 10, size = 4, color = 'blue') +
    
    
    # Add lower-right corner label
    annotate("text", x = Inf, y = Inf,
             label = paste0('# total predicted peaks: ', nrow(prosit_spectra)), 
             hjust = 1.2, vjust = 40, size = 4) +
    annotate("text", x = Inf, y = Inf,
             label = paste0('# predicted SPS ion peaks: ', nrow(prosit_SPSion)), 
             hjust = 1.2, vjust = 42, size = 4) +
    
    
    labs(title = paste0('Spectrum for ', unique(prosit_spectra$peptide_sequence)),
         y = 'Relative Abundance') +
    
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5)# Center the title
    )

}



########### 7) calculate the difference for neighbors whthin peptide
neighbor_diff <- function(data) {
  dist_all = c()
  for(ni in unique(data$pep_id)){
    abun_ni = subset(data, pep_id == ni)$pep_SN_ratio_1
    dist = c(0)
    for (i in 1:(length(abun_ni) - 1)) {
      
      dist_i = abs(abun_ni[i] - abun_ni[i + 1]) 
      dist = c(dist, dist_i)
      
    }
    dist_all = c(dist_all, dist)
  }
  return(dist_all)
}




############ 8) randomly remove rows for neighbor_dist < 0.03 & data_ni$neighbor_dist != 0
Random_remove_rows <- function(data, threshold = 0.03, remove_num = 2) {
  final_data <- list()
  # Loop through each unique pep_id
  # remove_rows = list()
  for(ni in unique(data$pep_id)) {
    # Subset the data for the current pep_id
    data_ni <- subset(data, pep_id == ni)
    # Identify rows with neighbor_dist < 0.02 & neighbor_dist != 0
    index <- data_ni$neighbor_dist < threshold & data_ni$neighbor_dist != 0
    
    # Determine how many rows to remove (2 if there are more than 2 such rows, otherwise 1)
    # rows_to_remove <- ifelse(sum(index) > 2, 2, 1)
    
    # If there are rows to remove, randomly sample those rows to remove
    if (sum(index) >=  remove_num) {
      rows_to_remove_indices <- sample(which(index), remove_num)
      
      # Remove the selected rows
      data_ni <- data_ni[-rows_to_remove_indices, ]
      
    }
    # if (sum(index) == 1) {
    #   rows_to_remove_indices <- which(index)
    #   
    #   # Remove the selected rows
    #   data_ni <- data_ni[-rows_to_remove_indices, ]
    #   
    # }
    # Append the modified data for this pep_id to the final list
    # remove_rows[[as.character(ni)]] = rows_to_remove_indices
    final_data[[ni]] <- data_ni
  }
  # Combine the list of data frames back into one data frame
  final_data_df <- do.call(rbind, final_data)
  return(final_data_df)
}



############ 9) plot to show the target density before and after down sampling
Density_target <- function(traindata, traindata_d, title = 'traindata') {
  plot(hist(traindata$pep_SN_ratio_1), col = "blue", lwd = 2, 
       main = paste0( "Density Plot before and after downsampling ", title),
       xlab = "Normalized SN_ratio", ylab = "Density")
  
  # Add the second density plot to the same figure
  lines(density(traindata_d$pep_SN_ratio_1), col = "red", lwd = 2)
 
  percent_removed <- (nrow(traindata) - nrow(traindata_d)) / nrow(traindata) * 100
  
  # Add the calculated percentage label to the plot (left-aligned text)
  text(x = max(traindata$pep_SN_ratio_1) * 0.1, y = max(density(traindata$pep_SN_ratio_1)$y) * 0.8,
       labels = paste("Removed:", round(percent_removed, 2), "%"),
       col = "black", cex = 1.2, adj = c(0, 0.5))  # Adjusted for left alignment
  
  # Optionally, add a legend to differentiate the two plots
  legend("topleft", legend = c("original", "downsampling"), 
         col = c("blue", "red"), lwd = 2)
}




########### 10) encode peptide (encode 'M*' as 'M')
# Define the alphabet mapping


# Function to encode a peptide sequence with correct handling of "M(ox)"
encode_peptide_M11 <- function(sequence, max_length = 20) {
  ALPHABET <- c("A" = 1, "C" = 2, "D" = 3, "E" = 4, "F" = 5, 
                "G" = 6, "H" = 7, "I" = 8, "K" = 9, "L" = 10, 
                "M" = 11, "N" = 12, "P" = 13, "Q" = 14, "R" = 15, 
                "S" = 16, "T" = 17, "V" = 18, "W" = 19, "Y" = 20, 
                "m" = 11)
  # "M(ox)" = 21)
  # Replace "M*" with "m" in the sequence
  sequence <- gsub("M\\*", "m", sequence)
  
  # Split the sequence into individual components
  peptide_parts <- unlist(strsplit(sequence, "(?=[A-Za-z()])", perl = TRUE))
  
  # Remove empty elements resulting from the split
  peptide_parts <- peptide_parts[peptide_parts != ""]
  
  # Map each component to its numeric value using ALPHABET
  encoded_peptide <- sapply(peptide_parts, function(x) ALPHABET[x])
  
  # Replace any invalid characters with 0
  encoded_peptide[is.na(encoded_peptide)] <- 0
  
  # # If the sequence is shorter than max_length, pad it with 0s at the end
  # if (length(encoded_peptide) < max_length) {
  #   encoded_peptide <- c(encoded_peptide, rep(0, max_length - length(encoded_peptide)))
  # }
  res = rep(0, max_length)
  res[1:length(encoded_peptide)] = encoded_peptide
  
  # # If the sequence is longer than max_length, truncate it
  # encoded_peptide <- encoded_peptide[1:max_length]
  
  return(unname(res))
}





















