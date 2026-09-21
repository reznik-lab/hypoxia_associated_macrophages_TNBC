library(data.table)
library(ggplot2)

#####

filePath="/Users/lium2/work/Ed_lab/Perry_lab_project/co_register_intensity/without_spatial_smoothing/07_07_2024"
fileName="ad_mIF_metadata_df.csv"
fileName=file.path(filePath,fileName)

ad_mIF_metadata = fread(fileName, data.table=FALSE)

#####

filePath="/Users/lium2/work/Ed_lab/Perry_lab_project/co_register_intensity/without_spatial_smoothing/07_07_2024"
fileName="ad_mIF_data_df.csv"
fileName=file.path(filePath,fileName)

ad_mIF_data = fread(fileName, data.table=FALSE)

# cluster 0, fibroblast, vasculature
# cluster 1, macrophages, and DAPI
# cluster 2, macrophages
# cluster 3, unknown
# cluster 4, Ki67

#####

ad_mIF_df = cbind(ad_mIF_metadata,ad_mIF_data)

#####

filePath="/Users/lium2/work/Ed_lab/Perry_lab_project/co_register_intensity/without_spatial_smoothing/07_07_2024"
fileName="ad_MALDI_metadata_wo_scale_df.csv"
fileName=file.path(filePath,fileName)

ad_MALDI_metadata = fread(fileName, data.table=FALSE)

ad_MALDI_metadata$clusters=ad_mIF_metadata$clusters

#####

filePath="/Users/lium2/work/Ed_lab/Perry_lab_project/co_register_intensity/without_spatial_smoothing/07_07_2024"
fileName="ad_MALDI_data_wo_scale_df.csv"
fileName=file.path(filePath,fileName)

ad_MALDI_data = fread(fileName, data.table=FALSE)

ad_MALDI_data_norm = log10(ad_MALDI_data*10^15 + 1 )

#####

ad_MALDI_df = cbind(ad_MALDI_metadata,ad_MALDI_data_norm)

filePath="/Users/lium2/work/Ed_lab/Perry_lab_project/co_register_intensity/without_spatial_smoothing/07_07_2024"
fileName="ad_MALDI_combined_wo_scale_df.csv"
fileName=file.path(filePath,fileName)

write.table(ad_MALDI_df,file=fileName, sep="\t",row.names = FALSE, col.names = TRUE)

#####

picked_cluster=0
picked_cluster=1
picked_cluster=2
picked_cluster=3
picked_cluster=4

group_1 = "group_E0"
group_2 = "group_4T1"
day_1 = "day_14"
day_2 = "day_7"

criteria_1 = ad_MALDI_df$clusters %in% picked_cluster
criteria_2 = ad_MALDI_df$group_id %in% group_1
criteria_3 = ad_MALDI_df$group_id %in% group_2

#criteria_4 = ad_MALDI_df$day_id %in% day_1
#criteria_4 = ad_MALDI_df$day_id %in% c(day_1,day_2)


tmp_df_group_1 = ad_MALDI_df[criteria_1 & criteria_2, ]
tmp_df_group_2 = ad_MALDI_df[criteria_1 & criteria_3, ]

#tmp_df_group_1 = ad_MALDI_df[criteria_1 & criteria_2 & criteria_4, ]
#tmp_df_group_2 = ad_MALDI_df[criteria_1 & criteria_3 & criteria_4, ]


show_title="cluster_4_Ki67_day_14_group_E0_vs_group_4T1"
#show_title="cluster_3_unknown_day_14_group_E0_vs_group_4T1"
#show_title="cluster_0_Vasculature_and_Fibroblast_day_14_group_E0_vs_group_4T1"
#show_title="cluster_1_macrophages_and_DAPI_day_14_group_E0_vs_group_4T1"
#show_title="cluster_2_macrophages_day_14_group_E0_vs_group_4T1"



metabolite_names = colnames(ad_MALDI_data)

test_result = data.frame()

#picked_metabolite= metabolite_names[1]
    
for(idx in 1:length(metabolite_names)){    
  
    picked_metabolite= metabolite_names[idx]
    
    message(sprintf("Work on %s", picked_metabolite))
   
    temp1 = tmp_df_group_1[, colnames(tmp_df_group_1) %in% picked_metabolite]
    temp2 = tmp_df_group_2[, colnames(tmp_df_group_2) %in% picked_metabolite]
  
    meanValueGroup1<-mean(10^temp1)
    meanValueGroup2<-mean(10^temp2)
    
    test_result[picked_metabolite,'log2FC'] = log2( meanValueGroup2/meanValueGroup1 )
    test_result[picked_metabolite,'FC'] = meanValueGroup2/meanValueGroup1
    
    temptest = wilcox.test(temp1,temp2,alternative = "two.sided", exact=FALSE)
    
    test_result[picked_metabolite,'pval'] = temptest$p.value
    test_result[picked_metabolite,'UStat'] = temptest$statistic
    
}    
    
if( sum(test_result$pval == 0) != 0 ){
  
  smallest_pvalue<-5e-324
  test_result[test_result$pval == 0,]$pval<-smallest_pvalue
  test_result$padj<-p.adjust(test_result$pval,method="BH")
  
}  


#####

#library(dplyr)
library(plyr)

test_data<-tmp_df_group_1
sample_id_list<-unique(test_data$sample_id)

result<-lapply(1:length(sample_id_list), function(x){
  
  test_data_subset<-test_data[test_data$sample_id %in% sample_id_list[x],]
  test_data_subset_meta<-test_data_subset[1,4:7]
  
  test_data_subset_value<-test_data_subset[,10:ncol(test_data_subset)]
  test_data_subset_value_average<-data.frame(t(colSums(test_data_subset_value)/nrow(test_data_subset_value)) , check.names = FALSE)
  
  test_data_subset_row_average<-cbind(test_data_subset_meta,test_data_subset_value_average)
  
  return(test_data_subset_row_average)
})    


tmp_df_group_1_metabolite_average<-rbind.fill(result)

#####

test_data<-tmp_df_group_2
sample_id_list<-unique(test_data$sample_id)

result<-lapply(1:length(sample_id_list), function(x){
  
  test_data_subset<-test_data[test_data$sample_id %in% sample_id_list[x],]
  test_data_subset_meta<-test_data_subset[1,4:7]
  
  test_data_subset_value<-test_data_subset[,10:ncol(test_data_subset)]
  test_data_subset_value_average<-data.frame(t(colSums(test_data_subset_value)/nrow(test_data_subset_value)), check.names=FALSE )
  
  test_data_subset_row_average<-cbind(test_data_subset_meta,test_data_subset_value_average)
  
  return(test_data_subset_row_average)
})    


tmp_df_group_2_metabolite_average<-rbind.fill(result)

combined_df<-rbind(tmp_df_group_1_metabolite_average,tmp_df_group_2_metabolite_average)

#####

metabolite_names = colnames(ad_MALDI_data)

test_data<-combined_df

result_out<-list()

for(idx in 1:length(metabolite_names)){    
  
  #idx<-37
  picked_metabolite= metabolite_names[idx]
  
  message(sprintf("Work on %s", picked_metabolite))
  

  temp1 = data.frame(test_data[, colnames(test_data) %in% picked_metabolite],stringsAsFactors = FALSE)
  colnames(temp1)<-picked_metabolite
  temp2 = cbind(test_data[,c(1:4)],temp1)
  
  tmpData<-temp2
  
  ####
  criteria_1_group_1<-"group_E0"
  criteria_1_group_2<-"group_4T1"
  criteria_2<-"day_14"
  
  group_1_data<-tmpData[tmpData$group_id %in% criteria_1_group_1, 5]
  group_2_data<-tmpData[tmpData$group_id %in% criteria_1_group_2, 5]
  
  wilcox_test_result<-wilcox.test(group_1_data,
              group_2_data,
              alternative = "two.sided",
              exact=FALSE)
  
  
  group_1_data<-tmpData[tmpData$group_id %in% criteria_1_group_1 & tmpData$day_id %in% criteria_2, 5]
  group_2_data<-tmpData[tmpData$group_id %in% criteria_1_group_2 & tmpData$day_id %in% criteria_2, 5]  
  
  t_test_result<-t.test(group_1_data,
         group_2_data, 
         exact=FALSE)
  
  mean_value_group_1<-mean(10^group_1_data)
  mean_value_group_2<-mean(10^group_2_data)
  log2_fold_change<-log2(mean_value_group_2/mean_value_group_1)
  
  test_result_df<-data.frame("metabolite"=picked_metabolite,
                             "ref_group"=criteria_1_group_1,
                             "test_group"=criteria_1_group_2,
                             "log2FC"=log2_fold_change,
                             "wilcox_pval"=wilcox_test_result$p.value,
                             "t_test_pval"=t_test_result$p.value,
                             stringsAsFactors = FALSE)
  
  result_out[[idx]]<-test_result_df
  
  ####

  
}    

result_out<-rbind.fill(result_out)
result_out$wilcox_padj<-p.adjust(result_out$wilcox_pval,method="BH")
result_out$t_test_padj<-p.adjust(result_out$t_test_pval,method="BH")


#####

threshold_pvalue=0.05
threshold_log2fc=1
show_additional_labels=NULL
numOfCandidateLabels=6

makeVolcanoPlotPublication2 = function(res,
                                           threshold_pvalue=0.05,
                                           threshold_log2fc=1,
                                           numOfCandidateLabels=10,
                                           savepath = NULL,
                                           show_title = NULL,
                                           show_additional_labels = NULL,
                                           show_legend_position = "bottom"
                                           )
{

  #res=test_result
  
  res$log10pvalue = -log10(res$padj)
  
  res[which(res[,"padj"] < threshold_pvalue & res[,'log2FC'] > abs(threshold_log2fc)),'color'] = 'Increase'
  res[which(res[,"padj"] < threshold_pvalue & res[,'log2FC'] < -abs(threshold_log2fc)),'color'] = 'Decrease'
  res$label = rownames(res)
  res$name = res$label
  
  sigmets = rownames(res)[which(res$color != 'None')]
  
  #numOfCandidateLables<-30
  
  if (length(sigmets) > numOfCandidateLabels){
    # rank them, keep only top 50
    if(numOfCandidateLabels>0){
      
      topNumOfCandidates_up = sigmets[order(res[sigmets,'log2FC'],decreasing = TRUE)][1:round(numOfCandidateLabels/2)]
      topNumOfCandidates_down = sigmets[order(res[sigmets,'log2FC'],decreasing = FALSE)][1:round(numOfCandidateLabels/2)]
      
      topNumOfCandidates = c(topNumOfCandidates_up,topNumOfCandidates_down)
      
      #topNumOfCandidates<-unique(c(topNumOfCandidates,show_additional_labels))
      
    }else{
      topNumOfCandidates<-NULL
    }
    not_topNumOfCandidates = setdiff(sigmets,topNumOfCandidates)
    res[not_topNumOfCandidates,'label'] = ""
  }
  
  if(!is.null(show_additional_labels)){
    res[show_additional_labels,"label"]<-show_additional_labels
  }
  
  
  
  criteria1_label<-which(res[,"padj"] < threshold_pvalue & abs(res[,'log2FC']) > threshold_log2fc )
  criteria2_label<-which(res$name %in% show_additional_labels)
  
  #threshold_pvalue=0.05
  #threshold_log2fc=1
  
  foldChangeRange<-ceiling(max(abs(res$log2FC)))
  
  #foldChangeBreaks<-seq(-foldChangeRange,foldChangeRange,1)
  foldChangeBreaks<-c(-rev(seq(2,foldChangeRange+1,2)),0,seq(2,foldChangeRange+1,2))
  
  yAxisRange<-ceiling(max(abs(res$log10pvalue)))+1
  #yAxisBreaks<-seq(0,yAxisRange,1)
  yAxisBreaks<-seq(0,yAxisRange,2)
  
  #show_legend_position="bottom"
  
  p <- ggplot(res,aes(x=log2FC,y=log10pvalue, label=label)) 
  p <- p + geom_hline(yintercept = -log10(threshold_pvalue), linetype="dashed", color="grey65", size = 0.25)
  p <- p + geom_vline(xintercept = abs(threshold_log2fc), linetype="dashed", color="grey65", size = 0.25)
  p <- p + geom_vline(xintercept = -abs(threshold_log2fc), linetype="dashed", color="grey65", size = 0.25)
  p <- p + geom_point(aes(color = color), size = 0.5)
  p <- p + geom_text_repel(
    data=res[unique(c(criteria1_label,criteria2_label)),],
    size= 7 / .pt,
    min.segment.length = 0,
    box.padding = 0.2,
    max.overlaps = Inf,
    #nudge_x = 0.15,
    segment.size  = 0.2,
    segment.color = "grey50",
    nudge_y = 0.15,
    #segment.curvature = -0.1,
    #segment.ncp = 3,
    #segment.angle = 20
  )   
  
  
  
  
  #p <- p + geom_point()
  p <- p + labs(title=show_title)
  p <- p + xlab('Log2 Fold Change') 
  p <- p + ylab('-log10 (padj)') 
  p <- p + scale_color_manual(name = 'Differential Abundance',
                              values = c('Decrease' = 'blue','Increase' = 'red','None' = 'grey85','Black' = 'Black'))
  #p <- p + theme(legend.position = 'bottom')
  
  #p <- p + xlim(-foldChangeRange,foldChangeRange)
  p <- p + scale_x_continuous(breaks=foldChangeBreaks)
  #p <- p + ylim(0,max(res$log10pvalue)+5)
  p <- p + expand_limits(x=c(-foldChangeRange-1,foldChangeRange+1),y=c(0,yAxisRange))
  #p <- p + scale_y_continuous(expand = expansion(mult=c(0, 0),add=1),
  #                            breaks=yAxisBreaks)
  
  p <- p + scale_y_continuous(expand = expansion(mult=c(0,0.05)))
  
                            #                            breaks=yAxisBreaks)
  #p <- p + coord_trans(y="log10")
  #p <- p + scale_y_continuous(trans="log10")
  #p <- p + ylim(0,max(res$log10pvalue)+5)
  p <- p + theme_classic()
  p <- p + theme(axis.text = element_text(size = 7),
                 #axis.title = element_text(size = 5, face="bold"),
                 axis.title.x = element_text(size=7),
                 axis.title.y = element_text(size=7),
                 text = element_text(size=7),
                 axis.text.x=element_text(colour="black"),
                 axis.text.y=element_text(colour="black"),
                 axis.ticks=element_line(colour="black"),
                 #axis.text.x = element_text(angle=45,hjust=1), 
                 #axis.title.x = element_blank(),
                 #axis.title.y = element_blank(),
                 #panel.border = element_rect(linetype = "solid", colour = "black"),
                 panel.border = element_blank(),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line=element_line(colour="black"),
                 plot.margin= margin(10, 10, 5, 5, "pt"),
                 plot.title = element_text(lineheight=1.5, face="bold",size=5,hjust=0.5),
                 plot.subtitle = element_text(lineheight=1.5, face="bold",size=7,hjust=0.5),
                 legend.position=show_legend_position)
  print(p)

  return(p)

}


#graph<-makeVolcanoPlotPublication2(test_result, numOfCandidateLabels=6, show_title=show_title, show_legend_position="None")

show_additional_labels <-c("Arginine_M.H",
                           "5-Oxo-Proline/Pyroglutamate_M.H",
                           "Ribose 5-phosphate_M.H",
                           "Sedoheptulose 7-phosphate_M.H",
                           "Lactic acid_M.H",
                           "Lactic acid_M.Cl",
                           "Glutamine_M.H",
                           "Glutamine_M.Cl",
                           "Glutamic acid_M.H")

result_out_plot<-result_out
rownames(result_out_plot)<-result_out$metabolite
result_out_plot$padj<-result_out_plot$t_test_padj
graph<-makeVolcanoPlotPublication2(result_out_plot, 
                                   threshold_pvalue=0.05,
                                   numOfCandidateLabels=12, 
                                   show_title=show_title,
                                   show_additional_labels = show_additional_labels,
                                   show_legend_position="None")

filePath<-"/Users/lium2/work/Ed_lab/Perry_lab_project/results/diff_met/07_08_2024"
fileName<-paste0(show_title,".pdf")
fileName=file.path(filePath,fileName)

pdf(fileName,width=5,height=5)
print(graph)
dev.off()





    





