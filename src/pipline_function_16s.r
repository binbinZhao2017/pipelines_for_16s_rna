# color set
mycols <<- c(
  "#e6194b", "#3cb44b", "#4363d8", "#ffe119", "#f58231",
  "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
  "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000",
  "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080",
  "#000000"
)
# some function
convertSeq2ASV<-function(phy){
  # rename sequence to ASV

  if(taxa_are_rows(phy) == F){
    phy <- t(phy)
    print("Converted to taxa name are rows")
  }

  if(all(rownames(otu_table(phy)) ==rownames(tax_table(phy)))){ #check names
    out <- cbind(paste0("ASV",seq(1:nrow(otu_table(phy)))),data.frame(refseq(phy))[,1],rownames(tax_table(phy)),tax_table(phy))
    write.csv(out,file=paste0("seq_asv_table.csv"),row.names = F,quote = F)
    taxa_names(phy) <- paste0("ASV",seq(1:nrow(otu_table(phy))))
  }else{
    message("rownames in otu table is not equal to rownames in tax table")
  }
  return(phy)
}


## group by other params

lefse_func1 <- function(group1, group2, phy_new, ref = "0"){
  # browser()
  drug_bash <- group1 %+% "vs" %+% group2 
  dosage_value <- strsplit(group1, "_")[[1]][1]
  batch_value <- strsplit(group1, "_")[[1]][2]
  input_level<-c(group1,group2)
  compari1 <- merge_phyloseq(subset_samples(phy_new, Dosage == "30" & Drug_batch == "201902"),subset_samples(phy_new, Groups == "Vehicle"))
  compari1_meta<- meta(compari1)
  metat_tmp_lefse <- compari1_meta %>% dplyr::select(Dosage,SampleID,Dosage)#u can change me!
  lefse_format_input<-paste0("output/structure/tables/lefse.input_",paste0(input_level,collapse = '_')) %+% drug_bash %+% ".txt" 
  prep_lefse(input_phy = compari1,outname=lefse_format_input ,meta_lefse_change=metat_tmp_lefse)
  system("sed -i 's/\\[//g' " %+% lefse_format_input)
  system("sed -i 's/\\]//g' " %+% lefse_format_input)
  system("sed -i 's/|/__/g' " %+% lefse_format_input)
  run_lefse_output<-'output/structure/tables/lefse.input_' %+% paste0(input_level,collapse = '_') %+% drug_bash %+% ".res"
  lefse_format_output<- ' output/structure/tables/lefse.input.in' %+% paste0(input_level,collapse = '_') %+% drug_bash
  str_bash<-'activate_local=/home/zhiyu/data/software/anaconda3/bin \n  source $activate_local/activate lefse2 \n  conda info --env \n'  %+%
    ' lefse-format_input.py ' %+% lefse_format_input  %+% lefse_format_output  %+% ' -c 1  -u 2 -o 1000000 \n' %+%
    ' run_lefse.py ' %+% lefse_format_output   %+% '  ' %+% run_lefse_output
  
  lefse_sh<-"lefse_sh" %+% paste0(input_level,collapse = '_') %+% drug_bash %+% '.sh'
  write_lines(str_bash,lefse_sh)
  system("bash -x " %+% lefse_sh)
  # system("rm" %+% lefse_format_output)
  run_lefse_output<-'output/structure/tables/lefse.input_' %+% paste0(input_level,collapse = '_') %+% drug_bash %+% ".res"
  dir.create("output/structure/figures/lefse/phy_subset")
  dir.create("output/structure/figures")
  system("sed -i 's/__/\\./g' " %+% run_lefse_output)
  plot_res_out<-"output/structure/figures/lefse/phy_subset/lefse_" %+% paste0(input_level,collapse = '_') %+% '' %+% drug_bash %+% ''  %+% '.svg' 
  # system("rm" %+% plot_res_out)
  plot_res(input_file=run_lefse_output,
           ref=ref,input_col=mycols[c(1,2)],dpi=300,width=25,height=30,ysize=20,lsize=25,
           outname= plot_res_out)
  knitr::include_graphics(plot_res_out)
}

spda_rere_plot<-function(dframe,Tax.filter,n,i){  
  # yy<<-dframe
  # uu<<-Tax.filter
  # dframe<-yy
  # Tax.filter<-uu
  # i<-"Genus"
  # n='n1'
  #dframe=res.splsda.plot_1
  dir.create("output/splsda",recursive = T)
  
  all_asv<-Tax.filter[, i]%>%data.frame()
  all_asv$SampleID<-rownames(all_asv)
  all_asv$SampleID2<-paste(all_asv[,i],"|",rownames(all_asv))
  dframe=merge(dframe,all_asv,by="SampleID",all.x = T)
  dframe=dframe[order(abs(dframe$importance),decreasing=T),]
  
  
  # browser()
  # dframe %>% select_if(function(importance){abs(importance)>=0.01})
  dframe$importance2<-dframe$importance %>% abs()
  dframe<-dframe %>% dplyr::filter(importance2>=0.01)
  tt<-as.numeric(as.character(dframe$importance) )
  
  
  col <- as.character(dframe$color)
  names(col) <- as.character(dframe$GroupContrib)
  
  
  p <- ggplot(data=dframe,aes(x= reorder(dframe$SampleID2, abs(tt) ) ,y= tt,fill = GroupContrib) )+ 
    #scale_y_discrete( limits = c( -1,1 ))+
    scale_fill_manual(values=col) +
    geom_bar(stat="identity",width = 0.5,position=position_dodge())  + 
    coord_flip() +theme_minimal() + 
    #scale_y_discrete( limits = c( -1,1 ))+
    labs(y="importance",x="") + 
    #theme(axis.text.y = element_text(angle = 0, size = 10, face = "bold"))
    theme(strip.text.x = element_text(size = 20, face = "bold"),
          text = element_text(size = 20),
          axis.text.x = element_text(angle = 0,size=15),
          axis.text.y =element_text(size = 10), 
          axis.title.y = element_blank() ,
          legend.text = element_text(size = 20) ,
          legend.title = element_text(size = 23) ,
          plot.title = element_text(size = 30, face="bold"))+
    ggtitle( paste0(" comp ",n," ",i) )
  
  ggsave(plot= p, filename = paste0("output/splsda/splda_test",".png"), dpi = 300, height = 10, width = 10, units = "in")
  
  p
  
  return(p)
}


## just group by groups
lefse_func <- function(group1, group2, phy_new, ref = group1, group = "Groups", batch){
  # browser()
  drug_bash <- batch %+% group1 %+% "vs" %+% group2 
  input_level<-c(group1,group2)
  rr <- phy_new@sam_data %>%data.frame()
  phy_new@sam_data$GROUPS <- rr[,group]
  s1 <- rownames(subset(sample_data(phy_new), GROUPS == group1))
  s2 <- rownames(subset(sample_data(phy_new), GROUPS == group2))
  compari1 <- merge_phyloseq(prune_samples(s1, phy_new), prune_samples(s2, phy_new))
  # compari1 <- merge_phyloseq(subset_samples(phy_new, Groups == group1),subset_samples(phy_tmp, Groups == group2))
  compari1_meta<- meta(compari1)
  metat_tmp_lefse <- compari1_meta %>% dplyr::select(GROUPS,SampleID,GROUPS)#u can change me!
  lefse_format_input<-paste0("output/structure/tables/lefse.input_",paste0(input_level,collapse = '_')) %+% drug_bash %+% ".txt" 
  prep_lefse(input_phy = compari1,outname=lefse_format_input ,meta_lefse_change=metat_tmp_lefse)
  system("sed -i 's/\\[//g' " %+% lefse_format_input)
  system("sed -i 's/\\]//g' " %+% lefse_format_input)
  system("sed -i 's/|/__/g' " %+% lefse_format_input)
  run_lefse_output<-'output/structure/tables/lefse.input_' %+% paste0(input_level,collapse = '_') %+% drug_bash %+% ".res"
  lefse_format_output<- ' output/structure/tables/lefse.input.in' %+% paste0(input_level,collapse = '_') %+% drug_bash
  str_bash<-'activate_local=/home/zhiyu/data/software/anaconda3/bin \n  source $activate_local/activate lefse2 \n  conda info --env \n'  %+%
    ' lefse-format_input.py ' %+% lefse_format_input  %+% lefse_format_output  %+% ' -c 1  -u 2 -o 1000000 \n' %+%
    ' run_lefse.py ' %+% lefse_format_output   %+% '  ' %+% run_lefse_output
  
  lefse_sh<-"lefse_sh" %+% paste0(input_level,collapse = '_') %+% "_" %+% drug_bash %+% '.sh'
  write_lines(str_bash,lefse_sh)
  system("bash -x " %+% lefse_sh)
  # system("rm" %+% lefse_format_output)
  run_lefse_output<-'output/structure/tables/lefse.input_' %+% paste0(input_level,collapse = '_') %+% drug_bash %+% ".res"
  dir.create("output/structure/figures/lefse/phy_subset")
  dir.create("output/structure/figures")
  system("sed -i 's/__/\\./g' " %+% run_lefse_output)
  plot_res_out<-"output/structure/figures/lefse/phy_subset/lefse_" %+% paste0(input_level,collapse = '_') %+% '' %+% drug_bash %+% ''  %+% '.svg' 
  # system("rm" %+% plot_res_out)
  plot_res(input_file=run_lefse_output,
           ref=ref,input_col=mycols[c(1,2)],dpi=300,width=25,height=30,ysize=20,lsize=25,
           outname= plot_res_out)
  knitr::include_graphics(plot_res_out)
}

ordered_phy<-function(phy,group){
  phy.meta<-meta(phy)
  phy.meta$gp <-as.factor(phy.meta[,group]) # add a new column gp
  sample_data(phy)<-phy.meta
  grp_caterogy<-levels(phy.meta$gp)
  obj_list<-list()

  for (cat_name in grp_caterogy){
    input_cat_name <<- cat_name
    phy_sub <- subset_samples(phy,gp == input_cat_name)
    obj_list<-c(obj_list,phy_sub)
  }
  order_phy <-  do.call(merge_phyloseq,obj_list) # merge phyloseq object
  return(order_phy)
}


plot_alpha <- function(input_phy,selected_col,input_levels,div_type,outname,col,legend=T){
  # add a new col
  # browser()
  sample_data(input_phy)$selected_col <- factor(unlist(sample_data(input_phy)[,colnames(meta(input_phy)) == selected_col]), levels = input_levels)
  #browser()

  # comparison
  allcomb <- combn(input_levels,m = 2)
  if(length(input_levels) >2){
    my_compare <- lapply(1:dim(allcomb)[2],function(i){ allcomb[,i]})
  }  else{
    my_compare <- list(input_levels)
  }

  # plot
  alpha <- plot_richness(input_phy, x = selected_col, measures = div_type, color = selected_col) +
    geom_boxplot(color=col) +
    scale_color_manual(values = col) +
    geom_jitter() +
    theme_bw() +
    xlab("") +
    stat_compare_means(comparisons = my_compare, method = "wilcox.test", label = "p.format",na.rm = T)+
    theme(strip.text.x = element_text(size = 15, face = "bold")) +
    theme(text = element_text(size = 15)) +
    theme(axis.text.x = element_text(angle = 90))

  if (legend != T ){
    alpha <- alpha + theme(legend.position="none")
  }
  ggsave(alpha, filename = outname, dpi = 300, height = 10, width = 4, units = "in")
  print(alpha)
}

prep_lefse <- function(input_phy,outname,meta_lefse_change){
  phyloseq <- input_phy
  # browser()
  phyloseq_otu<-as.data.frame(t(otu_table(phyloseq)))

  #transform to TSS
  phy.tss  = transform_sample_counts(phyloseq, function(x) x/sum(x) )

  # to data.frame
  ps.tax_level_to_df <- phyloseq_to_df(phy.tss)


  # form mega tax
  ps.tax_level_to_df <- ps.tax_level_to_df %>% unite("mega_tax","Phylum","Class","Order","Family","Genus",
                                                     "Species","OTU", sep = "|", remove = FALSE)

  ps.tax_level_to_df <- ps.tax_level_to_df[!is.na(ps.tax_level_to_df$Phylum),]
  ps.tax_level_to_df <- ps.tax_level_to_df %>%
    # dplyr::select(-Phylum,-Class,-Order,-Family,-Genus,-Species,-OTU, -Kingdom)
    dplyr::select(-Phylum,-Class,-Order,-Family,-Genus,-Species,-OTU)


  # prepare meta table
  t_meta <- as.data.frame(t(meta_lefse_change))
  # Meta_tab<-sample_data(phyloseq)
  # t_meta <- as.data.frame(t(Meta_tab))

  # merge table
  meta_tax_otu_df<- plyr::rbind.fill(t_meta,ps.tax_level_to_df)

  meta_tax_otu_df <-meta_tax_otu_df%>% dplyr::select("mega_tax", everything())
  # meta_tax_otu_df[3,] <- colnames(meta_tax_otu_df)
  # meta_tax_otu_df[4,] <- meta_tax_otu_df[13,]#add covariates cage number
  # meta_tax_otu_df <- meta_tax_otu_df[c(8,3,4,17:nrow(meta_tax_otu_df)),]
  # #head(meta_tax_otu_df)
  #meta_tax_otu_df[1:3,1] <- c("Class","Samples","Covariates")
    meta_tax_otu_df[1:2,1] <- c("Class","Samples")

  # write table
  write.table(meta_tax_otu_df, sep = "\t",file = outname,
              row.names = FALSE,col.names=FALSE,quote = FALSE)
  #write.table(meta_tax_otu_df, sep = "\t",file = paste0("output/30_07_2020/tables/lefse.input_",outname,".txt"),row.names = FALSE,col.names=FALSE,quote = FALSE)
}

plot_res <- function(input_file,ref,input_col,dpi=300,width=5,height=5,ysize=5,lsize=5,outname){
  tab <- read.delim(input_file,header = F)
  # browser()
  tab <- tab[!is.na(tab$V4),]
  tab$V3 <- gsub("01_","",tab$V3)
  tab$V3 <- gsub("02_","",tab$V3)
  tab$V3 <- gsub("03_","",tab$V3)
  try1 <- try(tab$V3 <- gsub("04_","",tab$V3))
  if('try-error' %in% class(try1)){
    next
  }else{
    tab$V3 <- gsub("04_","",tab$V3)
  }

  tab$V3 <- as.factor(tab$V3)

  # make rowname
  rname <- sapply(tab$V1, function(x){
    tax <- strsplit(as.character(x),"\\.")[[1]]
    if (length(tax) > 1){
          #out <- paste(paste0(c("p__","c__","o__","f__","g__","s__","")[1:length(tax)],strsplit(as.character(x),"\\.")[[1]])[(length(tax)-2):length(tax)],collapse = "|")
          out <- paste(paste0(c("p__","c__","o__","f__","g__","s__","")[1:length(tax)],strsplit(as.character(x),"\\.")[[1]])[c(4,5,7)],collapse = "|")
          return(out)
    }else{
      return(tax)
    }

    })
  rname <- gsub("NA","",rname)

  rownames(tab) <- rname
  for (i in 1:nrow(tab)){
    if( tab$V3[i] == ref){
      tab$V4[i] <- -1*tab$V4[i]
    }
  }
  tab <- tab[order(abs(tab$V4)),c(4,3)]
  tab <- cbind(rownames(tab),tab)
  colnames(tab) <- c("name","value","Group")

  p1 <- ggbarplot(tab, x = "name", y = "value",
          fill = "Group",
          color = "white",
          palette = input_col,
          sort.val = "desc",
          sort.by.groups = FALSE,
          x.text.angle = 0,
          ylab = "LDA score (log10)",
          rotate = TRUE,
          ggtheme = theme_pubr(),
          width=0.9
          ) +
    theme(axis.text.y=element_text(size=ysize),
          axis.text.x=element_text(size=ysize),
          legend.text=element_text(size=lsize),
          axis.title.x=element_text(size=(ysize+2)),
          axis.title.y=element_text(size=(ysize+2)),
          legend.title = element_text(size = lsize) )

  ggsave(p1,filename = outname,device = "svg",dpi = dpi,width = width,height = height)
}

plot_res_gene <- function(input_file,ref,input_col,dpi=300,width=5,height=5,ysize=5,lsize=5,outname, top =100){
  # 默认取前100行进行绘图，后面可以加if进行判断。
  tab <- read.delim(input_file,header = F)
  # browser()
  tab <- tab[!is.na(tab$V4),]
  tab$V3 <- gsub("01_","",tab$V3)
  tab$V3 <- gsub("02_","",tab$V3)
  tab$V3 <- gsub("03_","",tab$V3)
  try1 <- try(tab$V3 <- gsub("04_","",tab$V3))
  if('try-error' %in% class(try1)){
    next
  }else{
    tab$V3 <- gsub("04_","",tab$V3)
  }

  tab$V3 <- as.factor(tab$V3)

  # make rowname
  # rname <- sapply(tab$V1, function(x){
  #   tax <- strsplit(as.character(x),"\\.")[[1]]
  #   if (length(tax) > 1){
  #     #out <- paste(paste0(c("p__","c__","o__","f__","g__","s__","")[1:length(tax)],strsplit(as.character(x),"\\.")[[1]])[(length(tax)-2):length(tax)],collapse = "|")
  #     out <- paste(paste0(c("p__","c__","o__","f__","g__","s__","")[1:length(tax)],strsplit(as.character(x),"\\.")[[1]])[c(4,5,7)],collapse = "|")
  #     return(out)
  #   }else{
  #     return(tax)
  #   }
  #
  # })
  # rname <- gsub("NA","",rname)

  # rownames(tab) <- rname
  rownames(tab) <-tab$V1
  for (i in 1:nrow(tab)){
    if( tab$V3[i] == ref){
      tab$V4[i] <- -1*tab$V4[i]
    }
  }
  tab <- tab[order(abs(tab$V4)),c(4,3)][1:top,] # modified to V4  by binbinzhao@20200916
  tab <- cbind(rownames(tab),tab)
  colnames(tab) <- c("name","value","Group")

  p1 <- ggbarplot(tab, x = "name", y = "value",
                  fill = "Group",
                  color = "white",
                  palette = input_col,
                  sort.val = "desc",
                  sort.by.groups = FALSE,
                  x.text.angle = 0,
                  ylab = "LDA score (log10)",
                  rotate = TRUE,
                  ggtheme = theme_pubr(),
                  width=0.9
  ) +
    theme(axis.text.y=element_text(size=ysize),
          axis.text.x=element_text(size=ysize),
          legend.text=element_text(size=lsize),
          axis.title.x=element_text(size=(ysize+2)),
          axis.title.y=element_text(size=(ysize+2)),
          legend.title = element_text(size = lsize) )

  ggsave(p1,filename = outname,device = "svg",dpi = dpi,width = width,height = height)
}

NetCoMi_fc<-function(phyloseq_,target_column2,two_group_classfication2){
  library(NetCoMi)
  # browser()
  #group_column<-"disease"
  #group_columns_content<- c("Healthy Control", "Ulcerative Colitis")
  group_column<-target_column2
  group_columns_content<- two_group_classfication2
  # phyloseq_<<-phyloseq
  # phyloseq_<-phyloseq
  #file_meta<<-filter.meta

  sample_names(phyloseq_) <-  gsub("-",".",sample_names(phyloseq_))# maybe changed
  phyloseq_2 <- core(phyloseq_, detection = 0.0005, prevalence = 50/100)

  otuTable(phyloseq_2)<-otuTable(phyloseq_2) +1
  #sample_data(phyloseq_)<-file_meta
  otu.mt<-otu_table(phyloseq_2) %>% data.frame() %>% as.matrix() %>% t()
  otu_table(phyloseq_2) <-otu_table(otu.mt,taxa_are_rows=F)
  amgut2.filt.phy<-phyloseq_2
  amgut_split <- metagMisc::phyloseq_sep_variable(amgut2.filt.phy, group_column,drop_zeroes=FALSE)#"SEASONAL_ALLERGIES"
  meta_<-sample_data(amgut2.filt.phy) %>% data.frame()
  amgut_split
  ls
  phy1<-amgut_split[[1]]
  print(phy1)
  phy2<-amgut_split[[2]]
  print(phy2)

  net_season <- netConstruct(phy1, phy2, verbose = 2,
                             filtTax = "highestVar",filtTaxPar = list(highestVar = 50),measure = "spieceasi",
                             measurePar = list(method = "mb", nlambda=100, pulsar.params=list(rep.num=20)),
                             normMethod = "none", zeroMethod = "none",sparsMethod = "none", seed = 123456)

  props_season <- netAnalyze(net_season, clustMethod = "cluster_fast_greedy")
  labels_t<-paste0(props_season$input$adjaMat1 %>% rownames())

  #u can change the lables of nodes
  otu_tax<-phyloseq_to_df(phyloseq_)
  rownames(otu_tax)<-otu_tax$OTU
  otu_tax<-otu_tax[c(labels_t),]
  otu_tax$tax_otu<-paste0(otu_tax$Genus,'|',otu_tax$OTU)


  # netcomi_plot <-plot(props_season, sameLayout = TRUE, layoutGroup = 1,
  #                     nodeSize = "eigenvector", cexNodes = 1.5, cexLabels = 1.8,
  #                     groupNames = group_columns_content)
  netcomi_plot <-plot(props_season, sameLayout = TRUE, layoutGroup = 1,
                      nodeSize = "eigenvector", cexNodes = 1.5, cexLabels = 3,
                      groupNames = group_columns_content,labels=otu_tax$tax_otu)

  return(netcomi_plot)
}

