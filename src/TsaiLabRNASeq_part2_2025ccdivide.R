###############
## functions for Rmarkdown

set_variables = function(){
  # use <<- to set global variables
  
  ## files in data directory
  #file_in <<-  "JGI_MDfloral_85_counts.tsv"
  #file_design <<- "MD_design_nov2024.txt"
  #file_compare <<- "compare2_MD.txt"
  file_in <<-  "pdel_txi_counts_genelevel_2025.tsv"
  file_design <<- "pdel_design.txt"
  file_compare <<- "pdel_compare.txt"
  ## files in db directory
  ## Changed to Pt-HAP1 2024
  ## file_gtf <<-  "h1v5.gtf"
  file_gtf <<-  "PdelWV94v2.1.gtf"
  ##file_annotation <<- "h1.gene_desc.txt"
  file_annotation <<- "pdel.gene_desc.2026.txt"
  go_map <<- "pdel.uni_jgi.GO.list"  # map file linking gene GO terms.
  ##go_map <<- "h1.uni_jgi.GO.list"  # map file linking gene GO terms.
  go_desc <<- "gene_ontology_ext_rename.tab"
}


## Keep the commands below not changed

######################################################
get_fpkms = function(){
  # Inputs
  file_in_full = paste("data//", file_in , sep="")
  file_design_full = paste("data//", file_design, sep="")
  file_gtf_full = paste("db//",  file_gtf, sep="")
  # Ouputs
  file_fpkm = paste("fpkm_",file_in, sep="")
  file_fpkm_full = paste("data//", file_fpkm, sep="")
  
  ## Functions
  if (!require("GenomicFeatures")) {
    ##install.packages("BiocManager",repos = "http://cran.us.r-project.org")
    #BiocManager::install("GenomicFeatures",force=TRUE)
  }
  if (!require("DESeq2")) {
    #BiocManager::install("DESeq2",force=TRUE)
    #biocLite("DESeq2")
    #biocLite("DESeq2")
  }
  library(GenomicFeatures)
  library(DESeq2)
  
  
  if(!file.exists(file_fpkm_full)){
    # gtf file
    txdb <- makeTxDbFromGFF(file_gtf_full, format = "gtf", circ_seqs = character())
    ebg <- exonsBy(txdb, by="gene")
    
    # design file
    sampleTable0 = read.table(file_design_full, head=TRUE)
    # data file
    data_in = read.table(file_in_full, head=TRUE,row.names =1)
    countData = as.data.frame(data_in)
    samples = names(data_in)
    # grouping
    group =  as.character(sampleTable0[["Group"]])
    file_id =  as.character(sampleTable0[["File"]])
    site = match(samples, file_id)
    group_real = group[site]
	group_real[is.na(group_real)] = "not"
    colData = as.data.frame(cbind(samples,group_real))
    names(colData)=c("samples","condition")
    ## load data
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                  colData = colData,
                                  design = ~ condition)
    rowRanges(dds) = ebg
    fpkm_out = fpkm(dds)
    write.table(as.data.frame(fpkm_out),file=file_fpkm_full,sep="\t",quote = FALSE)
  }# file exists
}


############################################
get_gene_num_above_cutoff = function(fpkm_cut_off_in){
  # Input
  file_fpkm = paste("data//", "fpkm_", file_in ,sep="")
  file_design_full = paste("data//", file_design, sep="")
  # Output
  file_gene_num = paste("data//", "expressed_num_",file_in,sep="")
  if(!file.exists(file_gene_num)){
    sampleTable0 = read.table(file_design_full, head=TRUE)
    data_in = read.table(file_fpkm, head=TRUE,row.names =1)
    group_uni = unique(as.character(sampleTable0[["Group"]]))
    out = c()
    for (group_in in group_uni ){
      sampleTable = sampleTable0[sampleTable0$Group == group_in,]
      samples_in = as.character(sampleTable$File)
      fpkm_test = data_in[, samples_in]
      if (is.null(dim(fpkm_test))){
        fpkm_min = fpkm_test
      } else{
        fpkm_min = apply(fpkm_test, 1,min)
      }
      fpkm_sel = fpkm_min[fpkm_min >= fpkm_cut_off]
      gene_num = length(fpkm_sel)
      out = c(out,gene_num)
    }
    out_final2 = cbind(group_uni , out)
    colnames(out_final2) = c("Samples", "Gene Number")
    write.table(out_final2, file=file_gene_num, row.names = FALSE, sep="\t")
  }
}


###################################################
get_cut_tag = function (fpkm_out_sel, real_group ){
  real_sample = names(fpkm_out_sel)
  unique_group = unique(real_group)
  #average 
  group_num = length(unique_group)
  data_sel = fpkm_out_sel
  #### get the sub set of the data 
  group_num
  for (group_index in c(1:group_num)){
    #group_index = 5
    group_id = unique_group[group_index]
    sample_in_group = real_sample[real_group == group_id]
    #first_sample = sample_in_group[1]
    data_in_group = data_sel[,sample_in_group]
    ## single replicate
    if(is.null(dim(data_in_group))){
      data_ave  = data_in_group
    }else{
      data_ave  = apply(data_in_group ,1,min)
    }
    
    if(group_index == 1){
      ave_name = c(group_id)
      ave_data_all = data_ave
    }else{
      ave_name = c(ave_name,group_id)
      ave_data_all = cbind(ave_data_all,data_ave)
    }
  }
  colnames(ave_data_all)= ave_name
  out_tag = apply(ave_data_all ,1,max)
  out_tag_F = rep(0,length(out_tag))
  out_tag_F[out_tag>=1] = 1
  out_tag_F[out_tag>=3] = 3
  out_tag_F[out_tag>=5] = 5
  out_tag_F[out_tag>=10] = 10
  return (out_tag_F)
}



########################################################
run_PCA_cluster = function(){
	## Inputs
	file_design_full = paste("data//", file_design,sep="")
	file_fpkm = paste("data//", "fpkm_", file_in,sep="")
	## output
	file_out_pcc = paste("data//", "pcc_", file_in,".csv",sep="")
	
	## read data
	sampleTable = read.table(file_design_full, head=TRUE)
	data_in = read.table(file_fpkm, head=TRUE,row.names =1)
	samples_in = as.character(sampleTable$File)
	fpkm_test = data_in[, samples_in]
	
	## get colors
	test_group = as.character(sampleTable$Group)
	test_unique = unique(test_group)
	unique_g_len = length(test_unique)
	palette(rainbow(unique_g_len))
	colors=palette()  
	
	## assign colors
	site = match(test_group,test_unique)
	color_in<-colors[site]
	
	#### 
	#filtered data
	fpkm_cut_tag = get_cut_tag(fpkm_test,test_group)
	affydata_sel = fpkm_test[fpkm_cut_tag>=5,]
	data.ma = as.matrix(log(affydata_sel+1) )
	
	#sample_name<-names(data0)
	#data.ma<-as.matrix(data0)
	pcdat = princomp(data.ma)
	load_data<-pcdat$loadings
	
	
	part<-load_data[,c(1,2)]
	point_name<-  as.character(sampleTable$sample_id)
	title = paste("PCA Analysis ",sep="")
	plot(part,main=title,type='n')
	points(part,col = color_in,pch = 16)
	part_ave = apply(part,2,mean)
	h_cut = round(part_ave[2],3)
	v_cut = round(part_ave[1],3)
	abline(col="red",h=h_cut)
	abline(col="red",v=v_cut)
	x<-part[,1]
	y<-part[,2]
	text(x, y+0.01, label=point_name,cex = 0.6)
	title = paste("Sample Clustering ",sep="")
	colnames(data.ma)= point_name
	plot(hclust(dist(t(data.ma))),main = title)
	
	##############
	### PCC calculation
	data = affydata_sel 
    names<- colnames(data)
    ### run PCC
    x <- data
    y <- data
    bl <- lapply(x, function(u){
			lapply(y, function(v){
						cor(u,v,method="pearson") # Function with column from x and column from y as inputs
					})
	})
    out = matrix(unlist(bl), ncol=ncol(y), byrow=T)

    ### change names
    rownames(out)= names
    colnames(out)= names
    write.csv(out,file_out_pcc)	
}



################################################################
run_DEseq_General = function(gene_list_for_go,qvalue_cut,fold_change_cut,fpkm_cut,pvalue_cut){
  #Input
  file_annotation_full = paste("db//",file_annotation,sep="")
  file_compare_full = paste("data//", file_compare ,sep="")
  file_design_full = paste("data//",file_design,sep="")
  file_counts = paste("data//",file_in ,sep="")
  file_fpkm  = paste("data//","fpkm_", file_in,sep="")
  # output
  file_deg_num = paste("data//","DE_",file_in,sep="")
  dir.create(file.path("data//", "GO"), showWarnings = FALSE) 
  file_final_csv  = paste("data//","DE_",file_in, "_Final_Out.csv",sep="")
  
  
  if (!require("DESeq2")) {
    #source("https://bioconductor.org/biocLite.R")
    #biocLite("DESeq2")
  }
  library(DESeq2)
  
  
  data_anno = read.table(file_annotation_full ,header=TRUE,sep="\t",comment.char="",quote = '')
  gene_id = as.character(data_anno[[1]])
  
  ### Run the DE analysis
  
  
  options(digits = 4)
  if(!file.exists(file_deg_num)){
    #data
    counts_in = read.table(file_counts, head=TRUE,row.names =1)
    fpkms_in =  read.table(file_fpkm, head=TRUE,row.names =1)
    # design file
    sampleTable0 = read.table(file_design_full, head=TRUE)
    # compare pairs
    compare = read.table(file_compare_full, head=TRUE)
    total_num = dim(compare)[1]
    tracking = 0
    condition_name = c()
    gene_num_out = c()
    for (index_num in c(1:total_num)){
      tracking = tracking + 1
      test_name = as.character(compare[index_num,1])
      group = as.character(compare[index_num,2])
      g2 = as.character(compare[index_num,3])
      g1 = as.character(compare[index_num,4])
      TestGroup = as.character(sampleTable0[[group]])
      sampleTable_add = cbind(sampleTable0,TestGroup)
      sampleTable1 = sampleTable_add[sampleTable_add$TestGroup == g1,]
      sampleTable2 = sampleTable_add[sampleTable_add$TestGroup == g2,]
      sampleTable = rbind(sampleTable1,sampleTable2)
      # finished design, get the data
      sampleTable_sel = sampleTable
      samples_in = as.character(sampleTable_sel$File)
      group_real = as.character(sampleTable_sel$TestGroup)
      fpkm_test = fpkms_in[, samples_in]
      countData = counts_in[, samples_in]
      colData = as.data.frame(cbind(samples_in,group_real))
      names(colData)=c("samples","condition")
      dds <- DESeqDataSetFromMatrix(countData = countData,
                                    colData = colData,
                                    design = ~ condition)
      # DESeq2
      dds <- DESeq(dds)
      resSFtreatment <- results(dds, cooksCutoff=FALSE, contrast=c("condition",g2,g1))
      out = as.data.frame(resSFtreatment)
      # Get FPKM average
      fpkm_cut_tag = get_cut_tag(fpkm_test,group_real)
      out_test = cbind(out,fpkm_cut_tag)
      
      final_each = cbind(out_test$log2FoldChange,out_test$pvalue,out_test$padj,out_test$fpkm_cut_tag)
      names = c('(logFC)','(Pvalue)','(Qvalue)','(FPKMcut)')
      final_name = paste(test_name,'_',names,sep="")
      colnames( final_each) = final_name
      
      if (tracking == 1){
        final_table = final_each
      }else{
        final_table = cbind(final_table, final_each)
      }
      # DESeq2 
      
      # filtering for lists and  numbers 
      if (pvalue_cut == 1 ){
        gene_sel = out_test[((!is.na(out_test$padj))&(!is.na(out_test$fpkm_cut_tag))&(!is.na(out_test$log2FoldChange)))& out_test$padj < qvalue_cut & out_test$fpkm_cut_tag >=fpkm_cut & abs(out_test$log2FoldChange) >= log2(fold_change_cut) ,]
        
      }else{
        gene_sel = out_test[((!is.na(out_test$pvalue))&(!is.na(out_test$fpkm_cut_tag))&(!is.na(out_test$log2FoldChange)))&out_test$pvalue < pvalue_cut & out_test$fpkm_cut_tag >=fpkm_cut & abs(out_test$log2FoldChange) >= log2(fold_change_cut) ,]
      }
	  cus_pval_cutoff=0.05
	  cus_qval_cutoff=0.05
	  gene_sel_up = out_test[which(out_test$log2FoldChange > 1.25 & out_test$pvalue < cus_pval_cutoff),]
      gene_sel_do = out_test[which(out_test$log2FoldChange < -1.25 & out_test$pvalue < cus_pval_cutoff),]
      file_out_up = paste("data//GO//","UP_",test_name,".Pvaluecut.txt",sep="")
      file_out_do = paste("data//GO//","DO_",test_name,".Pvaluecut.txt",sep="")
      file_out_do = paste("data//GO//","DO_",test_name,".Pvaluecut.txt",sep="")
      gene_list_up = rownames(gene_sel_up)
      gene_list_do = rownames(gene_sel_do)
      gene_num_up = length(gene_list_up)
      gene_num_do = length(gene_list_do)
      write.table(gene_list_up, file= file_out_up , row.names = FALSE,col.names = FALSE)
      write.table(gene_list_do, file= file_out_do , row.names = FALSE,col.names = FALSE)
      #Ouput DE number
      condition_name = c(condition_name, paste("UP_",test_name,spe=""),paste("DO_",test_name,spe=""))
      gene_num_out = c(gene_num_out,gene_num_up,gene_num_do)
    }# for loop of tests
    
    
    out_final2 = cbind(condition_name, gene_num_out)
    colnames(out_final2) = c("Tests", "DEG number")
    write.table(out_final2, file=file_deg_num, row.names = FALSE, sep="\t")
    
    ### final CSV
    gene_names = rownames(fpkms_in)
    final_table = cbind(final_table,gene_names)
    ### final_table = cbind(final_table,gene_names,fpkms_in)
    match_out = match(gene_names,gene_id)
    anno_reorder = data_anno[match_out,]
    final_table = cbind(gene_names,final_table,anno_reorder)
    ### final_table = cbind(gene_names,final_table)
    write.csv(final_table, file=file_final_csv, row.names = FALSE,quote=TRUE)
    
  } # requried
  
}

#########################################################

run_go_enrichment = function(data_folder){
  ## Inputs
  go_map_full = paste("db//",go_map,sep="")
  go_desc_full=paste("db//",go_desc,sep="")
  # Output
  

  if (!require("topGO")) {
    BiocManager::install("topGO")
  }
  

  ###########################
  ## run GO
  library(topGO)
  GO_space = c("BP","MF","CC")
  for (go_group in GO_space){
  #BP MF  CC
  #read map
  geneID2GO <- readMappings(file = go_map_full)
  geneNames <- names(geneID2GO)
  dir.create(file.path(data_folder, go_group), showWarnings = FALSE)
  outfolder=paste(data_folder,"/",go_group,sep="")
  my_list<-list.files(path = data_folder, pattern =".txt$", all.files = TRUE,
                      full.names = FALSE, recursive = FALSE)
  
  mygo<-function(gene_list, folder){
    file<-paste(folder,"/",gene_list,sep="")
    #read select genes
    read.table(file,header=FALSE,sep="\t")->data
    myInterestingGenes <-as.character(data[[1]])
    geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
    names(geneList) <- geneNames
    BPterms <- ls(GOBPTerm)
    # MF BP CC
    GOdata <- new("topGOdata", ontology = go_group, allGenes = geneList,nodeSize=3, annot = annFUN.gene2GO,
                  gene2GO = geneID2GO)
    resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    num<-length(usedGO(GOdata))
    allRes <- GenTable(GOdata, classic = resultFis,
                       ranksOf = "classic", topNodes = num)
	order_data = allRes[order(allRes$GO.ID),]
	order_data$Fold=order_data$Significant/order_data$Expected
    file_out= paste(outfolder,"/","out_",go_group,"_",gene_list,".tab",sep="")
    write.table(order_data,file=file_out,sep="\t",quote = FALSE,row.names = FALSE,col.names =TRUE)
	return(order_data)
  }


  file_trk=0;
  for(fff in my_list){
	if(file_trk ==0){
		des=mygo(fff, data_folder)
		prename=strsplit(fff,"\\.")[[1]][1]
		colnames(des) <-paste(prename,colnames(des),sep="_")
		 # get level
		 read.table(go_desc_full,header=FALSE,sep="\t",quote = "")->data_anno
		go_id = as.character(data_anno[[1]])
		go_name_anno = as.character(data_anno[[4]])
		go_level_anno = as.character(data_anno[[5]])
		go_list_final =  as.character(des[,1])
		match_out = match(go_list_final,go_id)
		go_name = go_name_anno[match_out]
		go_level = go_level_anno[match_out]
		des=cbind(des,go_name,go_level)
	}else{
	    bes=mygo(fff, data_folder)
		prename=strsplit(fff,"\\.")[[1]][1]
		colnames(bes) <-paste(prename,colnames(bes),sep="_")
		des=cbind(des,bes[,4:7])
	}
	file_trk=file_trk+1
  }
  file_out= paste(data_folder,"/","out_",go_group,"_new_merged.tab",sep="")
  write.table(des,file=file_out,sep="\t",quote = FALSE,row.names = FALSE,col.names =TRUE)
  ###########################################################################
  ## merge GO
  read.table(go_desc_full,header=FALSE,sep="\t",quote = "")->data_anno
  go_id = as.character(data_anno[[1]])
  go_name_anno = as.character(data_anno[[4]])
  go_level_anno = as.character(data_anno[[5]])
  
  my_list<-list.files(path = paste(data_folder,"/",go_group,sep=""), pattern ="out.+\\.tab$", all.files = TRUE,
                      full.names = FALSE, recursive = FALSE)
  get_sample = function (x,in_split="\\.txt",perl = T){
    strsplit(x, in_split)[[1]][1]
  }
  sample_id = as.character(sapply(my_list ,get_sample))
  
  num_tracking = 0
  for(fff in my_list){    
    num_tracking = num_tracking+ 1
    #fff = "out_BP_Geno3_FD55_WW.txt"
    fff_full = paste( data_folder,go_group,"/",fff,sep="")
    read.table(fff_full,header=TRUE,sep="\t",quote = "")->data_in
    order_data = data_in[order(data_in$GO.ID),]
    pvalue = as.character(order_data$classic)
    pvalue[pvalue == '< 1e-30'] = 1e-30
    #print(length(pvalue))
	#order_data$TYPE=go_group
	#order_data$filesource=fff
	#order_data$Fold=order_data$Significant/order_data$Expected
	
    convert = (-1)* log10(as.numeric(pvalue))	
    if (num_tracking == 1){
      annotation = order_data[,c(1:3)]
      out = convert
    }else{
      out= cbind(out,convert)
    }
  }
  
  colnames(out)= sample_id
  # get max
  out_max = apply(out,1,max)
  # get level 
  go_list_final =  as.character(annotation$GO.ID)
  match_out = match(go_list_final,go_id)
  go_name = go_name_anno[match_out]
  go_level = go_level_anno[match_out]
  gO_merge_file_out = paste("data","/GO_enrichment","_",go_group,"_","merge.tab",sep="")
  final = cbind(annotation,out_max , out, go_name, go_level )
  write.table(final,file=gO_merge_file_out,sep="\t",quote = FALSE,row.names = FALSE,col.names =TRUE)
  }
}
