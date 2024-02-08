### This file is to generate and count simulated phage images by LoST for modeling the empirical bias correction method
### This file takes ~ 6h to run.
### Output file "result_model.csv" that will be used for fitting the model.

#get the current filepath
mainDir = getwd()
#create a subfolder to save simulated images for model estimation
dir.create(file.path(mainDir, "model_images/"))

#simulate phage images and save for model estimation
#replicate number for each fixed spot size and count
rep_num = 10
#different spot count to be modeled
ct_model = c(10,20,30,50,70,100,120,150,200,250,300)
#different spot diameter to be modeled
sz_model = c(20,30,40,50,60)
for (i in sz_model) {
  #create sub-folders for different spot size
  dir.create(paste(mainDir,"/model_images/model_images_",i,sep = ""))
  print(paste("Generating images with spot diameter",i))
  for (j in ct_model) {
    print(paste("Generating images with spot diameter",i, "and with spot count",j))
    for (k in 1:rep_num) {
      img0 <- PhageImage_gen(s0 = j, d2o = i,seed0 = k+100)
      writeImage(img0,paste(mainDir,"/model_images/model_images_",i,"/model_ct",j,"_sz",i,"_rep",k,".jpg",sep = ""))
    }
  }
}



## Obtain observed count and estimated spot size for each image and save results of images for modeling----

#for each sub-folder in the folder model_images, use LoST to obtain spot count and size, and combine the results
res_model = as.data.frame(NULL)
for (i in sz_model) {
  print(paste("Count images with spot diameter",i))
  # get a fake meta file about dilution information for LoST function
  info = fread(system.file("extdata", "metafile_example.csv",package="LoST"))
  names = list.files(paste(mainDir,"/model_images/model_images_",i,sep = ""))
  info <- do.call("rbind", replicate(length(names), info[1,], simplify = FALSE))
  info$image_name = names
  write.csv(info,paste(mainDir,"/model_images/info_model_",i,".csv",sep = ""))

  # use LoST function to obtain count and estimated spot size
  res0_model = LoST(paste(mainDir,"/model_images/model_images_",i,sep = ""),meta.file = paste(file.path(mainDir),"/model_images/info_model_",i,".csv",sep = ""),use_dilute = F, dilute.model = F)
  obs_count = as.numeric(res0_model@count)
  obs_size = unlist(lapply(res0_model@clus_size, function(t){median(sqrt(unlist(t)/pi)*2)}))
  true_count = as.numeric(gsub('[^[:alnum:] ]', '', substr(names,start = 9,stop = 11)))
  true_size = as.numeric(gsub('[^[:digit:] ]', '', substr(names,start = 14,stop = 16)))

  res_model = rbind(res_model, cbind(obs_count,obs_size,true_count,true_size))
  colnames(res_model) = c("obs_count","obs_size","true_count","true_size")

  #save the results for modeling
  write.csv(res_model,paste(mainDir,"/model_images/result_model.csv",sep = ""))
}
