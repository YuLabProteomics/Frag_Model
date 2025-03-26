source('000_functions.R')


############# Split the data into Training, Validation, and Testing
# load('../002_Rdata_CSV/0012_pep_SN_sort_SPSIons.Rdata') 
load('../002_Rdata_CSV/0012_pep_SN_sort_SPSIons_10Scans.Rdata') # keep peptides fragmented 10 times

########### 1. encoding functions (in 000_functions.R)




############### 2. index rows for data spiting
# load('../002_Rdata_CSV/0012_pep_SN_sort_SPSIons.Rdata')
load('../002_Rdata_CSV/0012_pep_SN_sort_SPSIons_10Scans.Rdata') # keep peptides fragmented 10 times

# save one fraction data to Zixuan
# sub_data = subset(pep_SN_sort, fraction == 'f1')
# write.csv(sub_data, file = '../002_Rdata_CSV/0012_spsion_f1.csv')
# pep_SN_sort = pep_SN_sort %>% arrange(pep_id)


sub_data = pep_SN_sort
sub_data = as.data.frame(sub_data)
rownames(sub_data) = paste0('row_', 1:nrow(sub_data)) 


##Min-Max normalization for each peptide's 10 target values
normalize_targets <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

sub_data = sub_data %>% group_by(pep_id) %>% 
  mutate(min_max_norm_SN = normalize_targets(top10_SN_ratio)) %>% 
  ungroup()
# replace the 0s by generate random values using normal distribution around the minimum of the peptideafter min-max normalization
set.seed(123)  # For reproducibility
options(scipen = 999)
# generate value
sub_data = sub_data %>% 
  group_by(pep_id) %>% 
  mutate(min_max_norm_SN_min =  min(min_max_norm_SN[min_max_norm_SN != 0]) ,
         min_max_norm_SN_min_normal = rnorm(1,mean = min_max_norm_SN_min, sd = min_max_norm_SN_min*0.02)) %>%
  ungroup()
sub_data = sub_data %>% 
  group_by(pep_id) %>% 
  arrange(pep_id, pep_SN_rank) %>% 
  ungroup()

#replace
sub_data[sub_data$min_max_norm_SN ==0, ]$min_max_norm_SN = unique(sub_data$min_max_norm_SN_min_normal)
sub_data = as.data.frame(sub_data)
rownames(sub_data) = paste0('row_', 1:nrow(sub_data)) 



## divide data into k-ending and R-ending peptides
K_pep = grepl("([K])$", sub_data$Peptide)
R_pep = grepl("([R])$", sub_data$Peptide)

split_ind = data.frame()
set.seed(seed = 123)
for(ni in list(K_pep, R_pep)){
  pep_ni = subset(sub_data, ni)
  idx = unique(pep_ni$pep_id)
  total_size <- length(idx)
  
  train_size <- floor(0.7 * total_size)  # 70% for training
  val_size <- floor(0.15 * total_size)   # 15% for validation
  test_size <- total_size - train_size - val_size  # # 15% for test
  
  # Randomly shuffle the indices
  shuffled_idx <- sample(idx)
  
  # Split into training, validation, and testing
  train_idx <- shuffled_idx[1:train_size]
  val_idx <- shuffled_idx[(train_size + 1):(train_size + val_size)]
  test_idx <- shuffled_idx[(train_size + val_size + 1):total_size]
  
  # row indices for data set
  train_row = rownames(pep_ni[pep_ni$pep_id %in% train_idx, ])
  
  val_row = rownames(pep_ni[pep_ni$pep_id %in% val_idx, ])
 
  test_row = rownames(pep_ni[pep_ni$pep_id %in% test_idx, ])

  tmp_df = data.frame(row_ind = c(train_row, val_row, test_row), 
                    train_ind =c(rep('train', length(train_row)), rep('val', length(val_row)), rep('test', length(test_row)))) 
  split_ind = rbind(split_ind, tmp_df)
  
}

############### 3. Process and split the data

####### feature peptide
train_row = as.numeric(limma::strsplit2(split_ind[split_ind$train_ind == 'train', ]$row_ind, split = '_')[,2])
val_row = as.numeric(limma::strsplit2(split_ind[split_ind$train_ind == 'val', ]$row_ind, split = '_')[,2])
test_row = as.numeric(limma::strsplit2(split_ind[split_ind$train_ind == 'test', ]$row_ind, split = '_')[,2])


## 1)encode the peptide
pep_encode = sapply(sub_data$Peptide, function(x) encode_peptide(sequence = x, max_length = 20))

pep_encode = t(as.matrix(pep_encode))
rownames(pep_encode) = NULL

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


## 8) combine all meta features together and spilit the dataset 
meta_all = cbind(t(pre_charge_encode), fragType_encode, normalized_NCE)
meta_all = as.matrix(meta_all)
rownames(meta_all) = NULL


meta_train = meta_all[train_row, ]
meta_val = meta_all[val_row, ]
meta_test = meta_all[test_row, ]

## 9) split target

# normalized to peptide max targer
target_train = sub_data$pep_SN_ratio_1[train_row]
target_val = sub_data$pep_SN_ratio_1[val_row]
target_test = sub_data$pep_SN_ratio_1[test_row]

## save the dataset out
save(pep_train, meta_train, target_train, file = '../002_Rdata_CSV/0012_01_train.Rdata')
save(pep_val, meta_val, target_val, file = '../002_Rdata_CSV/0012_01_val.Rdata')
save(pep_test, meta_test, target_test, file = '../002_Rdata_CSV/0012_01_test.Rdata')


# target_train = sub_data$pep_SN_rank[train_row]
# target_val = sub_data$pep_SN_rank[val_row]
# target_test = sub_data$pep_SN_rank[test_row]


# min-max target
# load('../002_Rdata_CSV/0012_01_wholeDataTable.Rdata')
# target_train = sub_data$min_max_norm_SN[train_row]
# target_val = sub_data$min_max_norm_SN[val_row]
# target_test = sub_data$min_max_norm_SN[test_row]
target_train = traindata$min_max_norm_SN
target_val = valdata$min_max_norm_SN
target_test = testdata$min_max_norm_SN

## save the dataset out
save(pep_train, meta_train, target_train, file = '../002_Rdata_CSV/0012_01_train_minMax.Rdata')
save(pep_val, meta_val, target_val, file = '../002_Rdata_CSV/0012_01_val_minMax.Rdata')
save(pep_test, meta_test, target_test, file = '../002_Rdata_CSV/0012_01_test_minMax.Rdata')






# save(pep_train, meta_train, target_train, file = '../002_Rdata_CSV/0012_train_f1.Rdata')
# save(pep_val, meta_val, target_val, file = '../002_Rdata_CSV/0012_val_f1.Rdata')
# save(pep_test, meta_test, target_test, file = '../002_Rdata_CSV/0012_test_f1.Rdata')

## save the whole table for datasets
traindata = sub_data[train_row,]
valdata = sub_data[val_row,]
testdata = sub_data[test_row, ]

save(traindata, valdata, testdata, file = '../002_Rdata_CSV/0012_01_wholeDataTable.Rdata')






######### format the target (sn_ratio for sps) table

# load('../002_Rdata_CSV/008_01_SPSions.Rdata')
# 
# SPSions <- lapply(names(SPSions), function(name) {
#   x <- SPSions[[name]]  # Access the data frame using the list element name
#   x$fraction <- rep(name, nrow(x))  # Add the fraction column with the name repeated for each row
#   return(x)
# })
# SPSions = do.call(rbind, SPSions)
