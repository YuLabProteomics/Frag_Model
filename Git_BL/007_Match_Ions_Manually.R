library(rawrr)# read raw ms data
library(PSMatch) # fragment the peptide
library(limma)
source('000_functions.R')

library(foreach)
library(doParallel)

# How many cores does your CPU have
n_cores <- detectCores()
# Register cluster
cluster <- makeCluster(3)
registerDoParallel(cluster)

# Export necessary functions and variables to the workers
clusterExport(cluster, list('Fragment_generator_NeutralLoss','calculateFragments',
                            'Frag_match', 'readSpectrum', 
                              'defaultNeutralLoss'))


########## 1) Generate fragment ions and do the match with spectrum in raw data
fl = list.files('../001_rawData/raw/')
fl = unique(strsplit2(fl, split = '\\.')[,1])


# for(mi in 1:2){
foreach(mi = 13:24, .packages = c("foreach", "doParallel", "dplyr")) %dopar% {
  flname = paste0('F_', mi)
  print(flname)
  ## load the data generate from code 006), 
  #including the identified scans (IdenScan_Comet)
  #scan need to find back (NeedScan) 
  #and the corresponding peptide information for fratmentation (PeptideForFrag)
  load(paste0('../002_Rdata_CSV/006_IdenScan_FailedScan_F', mi, '.Rdata'))
  raw_file = paste0('../001_rawData/raw/', fl[mi], '.raw')
  
  ### match failed scans
  matched_ion_Failed = lapply(NeedScan, function(x){
    ## localize the identified pep scan and extract peptide sequence and charge
    pepscan = x[length(x)]
    idenPep_table = subset(PeptideForFrag, ScanF == pepscan)
    idenPep = idenPep_table$Peptide
    charge = idenPep_table$z
    ## Fragment the peptide
    fragment = Fragment_generator_NeutralLoss(peptide= idenPep,
                                              charge = charge, 
                                              tolerance = 0.01,
                                              modification = c(Nterm = 295.1896, C = 57.02146))
    
    ## Extract the mz from raw file according to the scan number
    scans = x[-length(x)]
    spec = readSpectrum(rawfile = raw_file, scan = scans)
    
    scanMatched_ion = list()
    for(ni in 1:length(spec)){
      mz = spec[[ni]]$mZ
      intensity = spec[[ni]]$intensity
      noises = spec[[ni]]$noises
      experimental_spectrum = data.frame(mz = mz, 
                                         intensity = intensity,
                                         noises = noises)
      # print(dim(experimental_spectrum))
      
      ## do the match
      matched_ions = Frag_match(theoretical_ions = fragment,
                                experimental_spectrum = experimental_spectrum)
      
      scanMatched_ion[[as.character(spec[[ni]]$scan)]] = do.call(rbind, lapply(matched_ions, function(x) {
        # Add experimental mz and intensity to the theoretical matches dataframe
        theoretical_with_exp <- x$theoretical_matches
        theoretical_with_exp$exp_mz <- x$experimental_mz
        theoretical_with_exp$exp_intensity <- x$experimental_intensity
        theoretical_with_exp$exp_noises = x$experimental_noises
        theoretical_with_exp$Num_exp_ion = x$Num_exp_ion
        theoretical_with_exp$Num_matched_ion = length(matched_ions)
        theoretical_with_exp$scanNumber = spec[[ni]]$scan
        return(theoretical_with_exp)
      }))
    }
    return(scanMatched_ion)
  })
  delist =  unlist(matched_ion_Failed, recursive = FALSE)
  matched_ion_Failed = do.call(rbind, delist)
  rownames(matched_ion_Failed) = NULL
  
  ##### match the Identified scans
  matched_ion_IdenComet = lapply(IdenScan_Comet, function(x){
    ## localize the identified pep scan and extract peptide sequence and charge
    pepscan = x[1]
    idenPep_table = subset(PeptideForFrag, ScanF == pepscan)
    idenPep = idenPep_table$Peptide
    charge = idenPep_table$z
    ## Fragment the peptide
    fragment = Fragment_generator_NeutralLoss(peptide= idenPep,
                                              charge = charge, 
                                              tolerance = 0.01,
                                              modification = c(Nterm = 295.1896, C = 57.02146))
    
    ## Extract the mz from raw file according to the scan number
    scans = x
    spec = readSpectrum(rawfile = raw_file, scan = scans)
    
    scanMatched_ion = list()
    for(ni in 1:length(spec)){
      mz = spec[[ni]]$mZ
      intensity = spec[[ni]]$intensity
      noises = spec[[ni]]$noises
      experimental_spectrum = data.frame(mz = mz, 
                                         intensity = intensity,
                                         noises = noises)
      # print(dim(experimental_spectrum))
      
      ## do the match
      matched_ions = Frag_match(theoretical_ions = fragment,
                                experimental_spectrum = experimental_spectrum)
      
      scanMatched_ion[[as.character(spec[[ni]]$scan)]] = do.call(rbind, lapply(matched_ions, function(x) {
        # Add experimental mz and intensity to the theoretical matches dataframe
        theoretical_with_exp <- x$theoretical_matches
        theoretical_with_exp$exp_mz <- x$experimental_mz
        theoretical_with_exp$exp_intensity <- x$experimental_intensity
        theoretical_with_exp$exp_noises = x$experimental_noises
        theoretical_with_exp$Num_exp_ion = x$Num_exp_ion
        theoretical_with_exp$Num_matched_ion = length(matched_ions)
        theoretical_with_exp$scanNumber = spec[[ni]]$scan
        return(theoretical_with_exp)
      }))
    }
    return(scanMatched_ion)
  })
  flattened_list = unlist(matched_ion_IdenComet, recursive = FALSE)
  matched_ion_IdenComet = do.call(rbind, flattened_list)
  rownames(matched_ion_IdenComet) = NULL
  
  matched_ion_all = list(matched_ion_Failed = matched_ion_Failed, 
                         matched_ion_IdenComet = matched_ion_IdenComet)
  save(matched_ion_all,
       file = paste0('../002_Rdata_CSV/007_matched_ion_run2F', mi, '.Rdata'))
}
# Stop the cluster after processing is complete
stopCluster(cluster)
