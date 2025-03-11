library(keras3)
library(dplyr)
library(ggplot2)
library(limma)
library(reshape2)
library(caret)



########### 1) model use the normalized to the max (top10 SPS SN)

Model_009 <- function(shape_peptide = c(20), embedding_dim, 
                      unit_encoder1, dropout_1, 
                      att_num_heads, att_key_dim, 
                      shape_meta = c(7), dropout_meta, 
                      unit_decoder, lr_adam, 
                      save_path, save_base, ni, 
                      pep_train, meta_train, target_train,
                      epochs, batch_size,
                      pep_val, meta_val, target_val) {
  ############### Input peptide and encoders
  ## 1) input layerfor peptide
  input_peptides <- layer_input(shape = shape_peptide,dtype = 'int32', name = "peptides_in")  # peptides input layer
  
  ## 2) embedding 
  embedding_layer <- input_peptides %>%
    layer_embedding(input_dim = 22, output_dim = embedding_dim, dtype = 'float32', name = "embedding", mask_zero = TRUE) # dtype: float32
  
  ## 3) First Encoder GRU (Bidirectional)
  encoder1 <- embedding_layer %>%
    bidirectional(layer_gru(units = unit_encoder1, return_sequences = TRUE, activation = 'relu', name = "encoder1")) # units = 256
  
  ## 4) Dropout layer after first encoder
  dropout1 <- encoder1 %>%
    layer_dropout(rate = dropout_1, name = "dropout_1")
  
  ## 5) use multi head attention 
  # encoder_att <- layer_attention(list(dropout2, dropout2, dropout2))
  encoder_att <- layer_multi_head_attention(list(dropout1, dropout1, dropout1), 
                                            num_heads = att_num_heads, key_dim = att_key_dim)
  
  
  ############### Input meta (precursor charge, fragmentation method, NCE)
  ## 7) meta input layer
  input_meta = layer_input(shape = shape_meta, dtype = 'float32', name = 'meta_in')
  
  
  ## 8) Dense layer for meta processing
  meta_dense <- input_meta %>%
    layer_dense(units = 2*unit_encoder1, activation = "relu", name = "meta_dense") %>%
    layer_dropout(rate = dropout_meta, name = "dropout_meta")
  
  ## 9) reshape meta to the same shape as encoder_att for concatenation
  meta_dense_reshaped <- meta_dense %>%
    layer_repeat_vector(20)  # Repeat across 20 time steps
  
  ## 10) concatenate peptide and meta together by element-wise multiply
  add_meta = keras$layers$multiply(inputs = list(encoder_att, meta_dense_reshaped), name = 'add_meta')
  
  
  ## 12) flatten the layer after attention
  flattened <- add_meta %>%
    layer_flatten(name = "flattened")
  
  
  ## 13) Repeat the vector across time steps
  repeated <- flattened %>%
    layer_repeat_vector(20, name = "repeat")
  
  ## 14) Decoder GRU
  decoder <- repeated %>%
    layer_gru(units = unit_decoder, return_sequences = TRUE, name = "decoder", activation = 'relu')
  
  ## 15) Final dense layer for output
  final_output <- decoder %>%
    layer_flatten() %>%
    layer_dense(units = 1, name = "final_output")
  
  ## 16) Construct the model
  model <- keras_model(inputs = list(input_peptides, input_meta),
                       outputs = final_output)
  
  ## 17) Compile the model
  optimizer = optimizer_adam(learning_rate = lr_adam)
  
  model %>%
    compile(optimizer = optimizer, #"adam", # Optimizer
            loss = "mean_squared_error", # Loss function 
            metrics = c("mean_absolute_error")) # Metrics to track during training
  
  ## Model summary
  summary(model)
  
  
  ## 18) Define the early stopping callback
  early_stopping <- callback_early_stopping(
    monitor = 'val_loss',         # Metric to monitor (e.g., validation loss)
    patience = 5,                 # Number of epochs with no improvement to wait before stopping
    restore_best_weights = TRUE   # Restore the best model weights
  )
  
  ## 19) Define the model checkpoint callback to save the best model
  checkpoint <- callback_model_checkpoint(
    filepath = paste0(save_path,save_base, 'Models/best_model_', ni, '.keras'),   # File path to save the best model
    monitor = 'val_loss',         # Metric to monitor (usually validation loss)
    save_best_only = TRUE,        # Save only the best model
    mode = 'min',                 # Save model with minimum validation loss
    verbose = 1                   # Print messages when saving the model
  )
  
  
  ## 20) Define the learning rate scheduler
  lr_scheduler <- callback_reduce_lr_on_plateau(monitor = 'val_loss', 
                                                factor = 0.6, 
                                                patience = 5, 
                                                min_lr = 1e-7)
  
  ## 21) fit the model
  
  set.seed(123)
  history =  model  %>% 
    fit(
      x = list(pep_train,meta_train),    # List of inputs (train_peptides, train_meta)
      y = target_train,                      # Labels for training data
      epochs = epochs,                              # Number of epochs
      batch_size = batch_size,                          # Batch size
      validation_data = list(list(pep_val, meta_val), target_val),  # Validation data
      verbose = 1,                               # Verbosity (1 for progress bar)
      callbacks = list(lr_scheduler, early_stopping, checkpoint)  # List of callbacks
    )
}

