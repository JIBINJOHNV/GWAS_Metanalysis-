library(UpSetR)

df=read.csv("For_UpSetR_Analysis.csv")
myNamedList <- as.list(df)
test<-lapply(myNamedList,function(v){v[!is.na(v)]})

listInput <- list(EDU_Inpt=EDU_INPUT,COG_Input=COG_INPUT,SCZ_Input=SCZ_INPUT,Metal=METAL,PLEIO_ls=PLEIO_LS,PLEIO_Pleio=PLEIO_PLEIO,Gpcacov=GpcaCov,Gpcacorr=GpcaCorr)

OutputFilename="Upset_output_venn"
tiff(file=file.path(getwd(),paste0(OutputFilename,"_venn.tiff")), units="in", width=12, height=12, res=100)
upset(fromList(listInput), order.by = "freq",mb.ratio = c(0.55, 0.45),)
dev.off()


OutputFilename="Upset_output_venn"
tiff(file=file.path(getwd(),paste0(OutputFilename,"_venn.tiff")), units="in", width=12, height=12, res=100)

upset(df, sets = c("SCZ","PLEIO_LS","METAL","PLEIO_PLEIO","EDU","COG","GpcaCov","GpcaCorr") ,
      mb.ratio = c(0.55, 0.45), order.by = "freq")

dev.off()
