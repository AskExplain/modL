## 
setwd('~/Documents/main_files/GCP/explain-dpm/uploads/script_DRAFT_cell2pixel/')
image_list <- list.files('./shiny/v4/data/Bigger_image_new/')
image_list <- lapply(image_list,function(image_ID){
  png::readPNG(paste('./shiny/v4/data/Bigger_image_new/',image_ID,sep = ""))
})
saveRDS(object = image_list,file = "./shiny/v4/data/image_list.RDS")
