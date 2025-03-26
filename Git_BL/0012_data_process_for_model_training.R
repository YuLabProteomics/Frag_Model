library(ggplot2)
library(dplyr)
###


## load the SPS_ion with peptide sequence and masterscan
load('../002_Rdata_CSV/008_01_SPSion_stat.Rdata')

# add the fraction for each data frame in the list and combine it in a dataframe
SPSion_stat <- lapply(names(SPSion_stat), function(name) {
  x <- SPSion_stat[[name]]  # Access the data frame using the list element name
  x$fraction <- rep(name, nrow(x))  # Add the fraction column with the name repeated for each row
  return(x)
})
SPSion_stat = do.call(rbind, SPSion_stat)

# subset one fraction to zixuan
# data_f1 = subset(SPSion_stat, fraction == 'f1')
# write.csv(SPSion_stat, file = '../002_Rdata_CSV/0012_spsion_zixuan.csv')

plot = FALSE
if(plot == TRUE){
## barplot show the frequency peptide length
pep_length =  table(unname(sapply(unique(SPSion_stat$Peptide), function(x) nchar(x))))
pep_length = data.frame(pep_length)

p = ggplot(pep_length, aes(x = Var1, y = Freq))+
  geom_bar(stat = "identity") + 
  theme_bw() +
  labs(x = 'Peptide length',
       y = 'Frequency',
       title = 'Length of peptide') +
  theme(
    plot.title = element_text(hjust = 0.5),# Center the title
    # axis.text.x = element_text(angle = -45, hjust = 1),
    axis.text.x.bottom = element_text(angle = 30, hjust = 1, size = 7)
  )
print(p)
ggsave(filename = '../004_Figures/0012_Peptide_length_allFraction.png', plot = p, width = 5,height = 4)
}



############# 1) calculate the average for the duplicated peptides

# Step 1: unique the peptide within fraction and masterScan 
#(some peptide has been fragmented several times in one fraction)
SPSion_stat$FM = paste0(SPSion_stat$fraction, '_', SPSion_stat$masterScan)
pep = SPSion_stat %>%
  group_by(masterScan, fraction) %>%
  reframe(Peptide = unique(Peptide))


#Step 2: filter the non duplicated peptides
non_dup_pep = SPSion_stat %>%
  filter(Peptide %in% pep$Peptide[!(duplicated(pep$Peptide) | duplicated(pep$Peptide, fromLast = TRUE))])

#Step 3: filter duplicated peptides by Peptide
duplicated_peptides <- SPSion_stat %>%
  filter(Peptide %in% pep$Peptide[duplicated(pep$Peptide) | duplicated(pep$Peptide, fromLast = TRUE)])

# check quality of scans (calculate the sum of top10_SN within masterscan)
check_dup_q  <- duplicated_peptides %>%
  group_by(Peptide, PrecusorCharge, FM) %>%
  mutate(total_top10_SN = sum(top10_SN_ratio)) %>% 
  ungroup()

# joint the table and sort the NCE
tmp <- check_dup_q %>%
  group_by(FM) %>%
  arrange(FM, NCE) %>% 
  ungroup()

# keep the duplicated peptide with precusorCharge unique to other precusorCharge
Keep_1 <- tmp %>% 
  group_by(Peptide, PrecusorCharge) %>% 
  filter(n_distinct(total_top10_SN) == 1) 

# calculate the correlation between reference masterscan and the rest masterscans
masterscan_cor <- tmp %>% 
  group_by(Peptide, PrecusorCharge) %>% 
  filter(n_distinct(total_top10_SN) > 1) %>%  # Ensure there is more than one unique total_top10_SN
  mutate(
    # Now calculate the correlation for each row based on the reference top10_SN_ratio
    cor_results = purrr::map_dbl(total_top10_SN, function(ni) {
      # Get all the top10_SN_ratio values for the reference (max total_top10_SN)
      ref_values = top10_SN_ratio[total_top10_SN == max(total_top10_SN)]
      
      # Get the current values for the current total_top10_SN (ni)
      current_values = top10_SN_ratio[total_top10_SN == ni]
      
      # Ensure the lengths of ref_values and current_values are the same before calculating correlation
      if(length(ref_values) == length(current_values)) {
        # Calculate the correlation (using 'complete.obs' to handle missing data)
        return(cor(ref_values, current_values, use = "complete.obs"))
      } else {
        # If lengths don't match, return NA for correlation
        return(NA_real_)
      }
    })
  ) %>% 
  ungroup()

# Keep the masterscan correlation bigger than 0.6
Keep_2 = masterscan_cor %>% 
  filter(cor_results > 0.6) 

# part 1, take the peptides have only 2 replicates
Keep_2_part1 =  Keep_2 %>% 
  group_by(Peptide, PrecusorCharge) %>% 
  filter(n_distinct(total_top10_SN) ==2) %>% 
  filter(total_top10_SN == max(total_top10_SN)) %>% 
  mutate(avg_top10_SN_ratio = top10_SN_ratio)

# part 1, take the peptides have more than 2 replicates
Keep_2_part2 =  Keep_2 %>% 
  group_by(Peptide, PrecusorCharge) %>% 
  filter(n_distinct(total_top10_SN) >2) %>%  
  group_by(Peptide, NCE, PrecusorCharge) %>%
  mutate(avg_top10_SN_ratio = mean(top10_SN_ratio, na.rm = TRUE)) %>% 
  filter(total_top10_SN == max(total_top10_SN))

# combine keep_1 and Keep_2_part1 and Keep_2_part2         
good_scan = rbind(Keep_2_part1,Keep_2_part2)
good_scan$top10_SN_ratio = good_scan$avg_top10_SN_ratio
good_scan = good_scan[, colnames(non_dup_pep)]

good_scan = rbind(good_scan, Keep_1[, colnames(non_dup_pep)])
  

# Step 6: rbind the non_dup and duplicated peptides
Uni_pep_SPSIons = rbind(non_dup_pep, good_scan)


# Step 7: keep peptides length <20
 pep_length = lapply(setNames((Uni_pep_SPSIons$Peptide), (Uni_pep_SPSIons$Peptide)), function(x) {nchar(x)})
 Uni_pep_SPSIons= Uni_pep_SPSIons[pep_length <=20, ]

 #save the dataset
 save(Uni_pep_SPSIons, file = '../002_Rdata_CSV/0012_Uni_pep_SPSIons.Rdata')
 
# Step 8: create indices for peptide
pep_indices = Uni_pep_SPSIons %>%
  group_by(Peptide, PrecusorCharge) %>%
  mutate(pep_id = cur_group_id()) %>%
  ungroup()

# Step 9: create SN_ratio and SN_rank within peptide 
pep_SN_sort = pep_indices %>%
  group_by(pep_id) %>%
  mutate(pep_SN_rank = rank(-top10_SN_ratio,ties.method = 'first'),
         pep_SN_ratio = max(top10_SN_ratio)/top10_SN_ratio,
         pep_SN_ratio_1 = top10_SN_ratio/max(top10_SN_ratio)) %>%
  # reframe(pep_SN_ratio = top10_SN_ratio/max(top10_SN_ratio)) %>%
  ungroup()

# Step 10: remove miss cleaved peptidts
filtered_peptides_KR <- grepl("([KR])$", pep_SN_sort$Peptide)
pep_SN_sort = pep_SN_sort[filtered_peptides_KR, ]

#save the dataset
save(pep_SN_sort, file = '../002_Rdata_CSV/0012_pep_SN_sort_SPSIons.Rdata')

#Step 12: keep peptides fragmented 10 times
indices = as.data.frame(table(pep_SN_sort$pep_id))
indices = indices[indices$Freq == 10, ]
pep_SN_sort = pep_SN_sort[pep_SN_sort$pep_id %in% indices$Var1,]
save(pep_SN_sort, file = '../002_Rdata_CSV/0012_pep_SN_sort_SPSIons_10Scans.Rdata')



############# 2)Split the data into Training, Validation, and Testing
load('../002_Rdata_CSV/0012_pep_SN_sort_SPSIons.Rdata')
  
########### 1) encode peptide 
# Define the alphabet mapping
ALPHABET <- c("A" = 1, "C" = 2, "D" = 3, "E" = 4, "F" = 5, 
              "G" = 6, "H" = 7, "I" = 8, "K" = 9, "L" = 10, 
              "M" = 11, "N" = 12, "P" = 13, "Q" = 14, "R" = 15, 
              "S" = 16, "T" = 17, "V" = 18, "W" = 19, "Y" = 20, 
              "M(ox)" = 21)

# Function to encode a peptide sequence with correct handling of "M(ox)"
encode_peptide <- function(sequence, max_length = 30) {
  # Replace "M*" with "M(ox)" in the sequence
  sequence <- gsub("M\\*", "M(ox)", sequence)
  
  # Split the sequence into individual components
  peptide_parts <- unlist(strsplit(sequence, "(?=[A-Za-z()])", perl = TRUE))
  
  # Remove empty elements resulting from the split
  peptide_parts <- peptide_parts[peptide_parts != ""]
  
  # Map each component to its numeric value using ALPHABET
  encoded_peptide <- sapply(peptide_parts, function(x) ALPHABET[x])
  
  # Replace any invalid characters with 0
  encoded_peptide[is.na(encoded_peptide)] <- 0
  
  # If the sequence is shorter than max_length, pad it with 0s at the end
  if (length(encoded_peptide) < max_length) {
    encoded_peptide <- c(encoded_peptide, rep(0, max_length - length(encoded_peptide)))
  }
  
  # If the sequence is longer than max_length, truncate it
  encoded_peptide <- encoded_peptide[1:max_length]
  
  return(unname(encoded_peptide))
}


########### 2) One-hot encoding function
one_hot_encode_precCharge <- function(charge, max_charge = 6) {
  encoding <- rep(0, max_charge - 1)
  encoding[charge - 1] <- 1
  return(encoding)
}



load('../002_Rdata_CSV/0012_pep_SN_sort_SPSIons.Rdata')

sub_data = pep_SN_sort

# sub_data = subset(pep_SN_sort, fraction == 'f1')
sub_data = as.data.frame(sub_data)
rownames(sub_data) = paste0('row_', 1:nrow(sub_data)) 

############### 2. index rows for data spiting

idx = unique(sub_data$pep_id)
total_size <- length(idx)

train_size <- floor(0.7 * total_size)  # 70% for training
val_size <- floor(0.15 * total_size)   # 15% for validation
test_size <- total_size - train_size - val_size  # # 15% for validatio

# Randomly shuffle the indices
shuffled_idx <- sample(idx)

# Split into training, validation, and testing
train_idx <- shuffled_idx[1:train_size]
val_idx <- shuffled_idx[(train_size + 1):(train_size + val_size)]
test_idx <- shuffled_idx[(train_size + val_size + 1):total_size]

# row indices for data set

train_row = rownames(sub_data[sub_data$pep_id %in% train_idx, ])
train_row = as.numeric(strsplit2(train_row, split = '_')[,2])

val_row = rownames(sub_data[sub_data$pep_id %in% val_idx, ])
val_row = as.numeric(strsplit2(val_row, split = '_')[,2])

test_row = rownames(sub_data[sub_data$pep_id %in% test_idx, ])
test_row = as.numeric(strsplit2(test_row, split = '_')[,2])



############### 3. Process and split the data

####### feature peptide

## 1)encode the peptide
pep_encode = sapply(sub_data$Peptide, function(x) encode_peptide(sequence = x, max_length = 30))
colnames(pep_encode) = NULL
pep_encode = t(as.matrix(pep_encode))

pep_train = pep_encode[train_row, ]
pep_val = pep_encode[val_row, ]
pep_test = pep_encode[test_row,]



####### meta features 

## 2)encode the precursor charge
pre_charge_encode = sapply(sub_data$PrecusorCharge, function(x) 
  one_hot_encode_precCharge(charge = x, max_charge = 6))

## 3)encode fragmentation type
fragType = limma::strsplit2(sub_data$NCE, split = '_')[,1]
df = data.frame(fragType)

fragType_encode <- df %>%
  mutate(fragType = case_when(
    fragType == "cid" ~ 1,
    fragType == "hcd" ~ 2
  ))

## 4)normalize NCE
NCE = limma::strsplit2(sub_data$NCE, split = '_')[,2]
normalized_NCE = as.numeric(NCE)/100

## 5)normalized to base_scan SN_ratio of top10 sps
SPS_SN_normalized = sub_data$pep_SN_ratio_1

## 6)eligible Number ratio (total number of sps vs matched number of ion) 
sps_ratio = sub_data$Num_total_sps/sub_data$Num_matched_ion

## 7) eligible sps intensity ratio (top10 sps int vs total(matched) ion int)
top10_sps_ratio = sub_data$int_top10_sps/sub_data$int_total_ion


## 8) combine all meta features together and spilit the dataset 
meta_all = cbind(t(pre_charge_encode), fragType_encode, normalized_NCE, SPS_SN_normalized, sps_ratio, top10_sps_ratio)
meta_all = as.matrix(meta_all)
rownames(meta_all) = NULL


meta_train = meta_all[train_row, ]
meta_val = meta_all[val_row, ]
meta_test = meta_all[test_row, ]

## 9) split target
target_train = sub_data$pep_SN_rank[train_row]
target_val = sub_data$pep_SN_rank[val_row]
target_test = sub_data$pep_SN_rank[test_row]

## save the dataset out
save(pep_train, meta_train, target_train, file = '../002_Rdata_CSV/0012_train.Rdata')
save(pep_val, meta_val, target_val, file = '../002_Rdata_CSV/0012_val.Rdata')
save(pep_test, meta_test, target_test, file = '../002_Rdata_CSV/0012_test.Rdata')

# save(pep_train, meta_train, target_train, file = '../002_Rdata_CSV/0012_train_f1.Rdata')
# save(pep_val, meta_val, target_val, file = '../002_Rdata_CSV/0012_val_f1.Rdata')
# save(pep_test, meta_test, target_test, file = '../002_Rdata_CSV/0012_test_f1.Rdata')









