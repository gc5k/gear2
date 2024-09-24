library(ggplot2)
library(tidyr)
library(cowplot)
library(gridExtra)
library(grid)
library(dplyr)
library(viridis)
library(ggpointdensity)
library(stringr)
library(scales)
library(reshape2)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplotify)
library(parallel)

arg<-commandArgs(T)
#########*********PLOT FUNCTION*********##########
lm_eqn1 <- function(plotDT){
  m <- lm(rsq ~ snp_1, plotDT)
  eq <- substitute(italic(y) == a + b * italic(x), 
                   list(a = format(unname(coef(m)[1]), digits = 2), 
                        b = format(unname(coef(m)[2]), digits = 3)))
  as.character(as.expression(eq))
}
lm_eqn2 <- function(plotDT){
  m <- lm(rsq ~ snp_1, plotDT)
  eq <- substitute(italic(R)[I]^2~"="~r2~","~italic(p)~"="~pvalue,
                   list(r2 = formatC(summary(m)$r.squared,format = "f", digits = 2),
                        pvalue = format(summary(m)$coefficients[2,4],digits = 3,scientific =T)))
  
  as.character(as.expression(eq))
}

lm_eqn1_inter <- function(plotDT){
  m <- lm(rsq ~ snp_1-1, plotDT)
  eq <- substitute(italic(y) == a * italic(x)
                   #~","~italic(R)[II]^2~"="~r2
                   ,list(r2 = formatC(summary(m)$r.squared,format = "f", digits = 2),
                     a = format(unname(coef(m)[1]), digits = 2)))
  as.character(as.expression(eq))
}
lm_eqn2_inter <- function(plotDT){
  m <- lm(rsq ~ snp_1 + 0, plotDT)
  pval <- cor.test(plotDT$rsq,plotDT$snp_1)$p.value
  eq <- substitute(italic(rho)[2]~"="~r~","~italic(p)~"="~pvalue,
                   list(r = format(cor(plotDT$rsq,plotDT$snp_1), digits = 2),
                        pvalue = ifelse(pval == 0, "< 1e-324", format(pval, digits = 3, scientific = TRUE))))
  as.character(as.expression(eq))
}
#########*********PLOT FUNCTION*********##########

####High-resolution LD
input=arg[1]
num_cores <- 10
chrLD <- read.table(input,header=T)

xld<-0
if(length(strsplit(input,split = ".xld")[[1]])==2){
  xld<-2
}
#for blank
#chrLD[which(chrLD$rsq<0),5] <- min(chrLD[which(chrLD$rsq>0),5],na.rm = T)
chrLD$CHRi <- str_split_fixed(chrLD$Tagi, "_",2)[,1]
chrLD$CHRj <- str_split_fixed(chrLD$Tagj, "_",2)[,1]
chrLD$Tagi <- as.numeric(str_split_fixed(chrLD$Tagi, "_",2)[,1])*10000+as.numeric(str_split_fixed(chrLD$Tagi, "_",2)[,2])
chrLD$Tagj <- as.numeric(str_split_fixed(chrLD$Tagj, "_",2)[,1])*10000+as.numeric(str_split_fixed(chrLD$Tagj, "_",2)[,2])

chrLD <- mclapply(1:nrow(chrLD), function(x) {
  row <- chrLD[x, ]
  if (row[1] < row[3]) {
    tmp <- row[c(1:2,7+xld)]
    row[c(1:2,7+xld)] <- row[c(3:4,8+xld)]
    row[c(3:4,8+xld)] <- tmp
  }
  return(row) 
}, mc.cores = num_cores)

# 将结果合并为矩阵
chrLD <- do.call(rbind, chrLD)

chrLD <- chrLD[order(chrLD$Tagi,chrLD$Tagj),]
chromosome<-chrLD[which(chrLD$Tagj==chrLD[1,3]),]$CHRi
chrLD$rsq_log10 <- -log10(chrLD$rsq)
chrLD_Lower_Melt <- dcast(data = chrLD[,c(1,3,9+xld)],Tagi~Tagj)
chrLD_Lower_Melt <- chrLD_Lower_Melt[,-1]

###Only lower triangle
chrLD_Lower_Melt[upper.tri(chrLD_Lower_Melt,diag = F)] <- ""
# visulation
chrLD_Lower_Melt <- data.frame(chrLD_Lower_Melt)
colnames(chrLD_Lower_Melt) <- seq(1:ncol(chrLD_Lower_Melt))
chrLD_Lower_Melt=apply(chrLD_Lower_Melt,2,as.numeric)
max(chrLD$rsq_log10,na.rm = T)
min(chrLD$rsq_log10,na.rm = T)

col = c("#CB2A04FF","#F66B19FF","#FABA39FF","#C7EF34FF","#1AE4B6FF","#36AAF9FF","#4662D7FF","#30123BFF")
col_fun = colorRamp2(c(0, 1, 2, 4, 6, 8, 10), col[1:7])

P1 <- Heatmap(chrLD_Lower_Melt,
              column_title_side = "bottom",column_names_gp = gpar(fontsize = 12,fontfamily = "serif"),
              cluster_rows = FALSE, cluster_columns = FALSE,
              col = col_fun,show_row_names = FALSE,show_column_names = FALSE,na_col = "white",
              use_raster = F,
              heatmap_legend_param = list(
                title = expression(paste(-log[10],"LD")), at = c(0, 1, 2, 4, 6, 8 ,10),
                labels = c(0, 1, 2, 4, 6, 8 ,">=10")
              ),
              row_split = factor(chromosome,levels = c(1:length(unique(chrLD$Tagi)))), column_split = factor(chromosome,levels = c(1:length(unique(chrLD$Tagi))))
)
P2 <- as.ggplot(P1) +
  theme(text=element_text(size = 12,family = "serif"),
        plot.title = element_text(hjust=0.5),
        legend.title = element_blank(),
        legend.position = c(0.95,0.95),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

#}



####Norm I and II
List_Files=arg[2]
SNP_vs_Rsq <- matrix(NA,nrow = length(List_Files),ncol = 6)
norm_12<-NULL
for(i in 1:length(List_Files)){
  ##Read sample size
  n <- as.numeric(length(read.table(list.files(pattern=paste0(strsplit(List_Files[i], split = ".xld")[[1]][1],".fam$")) )[,1]))
  
  ##Collect eigenvalue files and record the origin chr-order in eigen_tag
  Eigenvalue_Files <- list.files(pattern=paste0(strsplit(List_Files[i], split = ".xld")[[1]][1],".*\\.evals.txt$"))
  eigen_tag <- gsub(paste0(strsplit(List_Files[i], split = ".xld")[[1]][1],"_"),"",gsub(".evals.txt","",Eigenvalue_Files))
  
  ##read xld files
  Data <- read.table(List_Files[i],header = T, fill = TRUE)[,c(1:5)]
  colnames(Data)[1]<-"CHRi"
  colnames(Data)[3]<-"CHRj"
  #subset inter and generate plot df
  inter <- Data[Data[,1] != Data[,3],]
  lij<-mean(Data$rsq)
  #subset intra
  Data <- Data[Data[,1] == Data[,3],]
  Data$snp_1 <- 1/Data$SNPi
  li<-mean(Data$rsq)
  ##correct the order of the chromosomes in eigenvalues adn .xld
  # 为 Data$Tagi 添加 'chr' 前缀
  chr_in_Data_Tagi <- Data$CHRi
  
  # 在 eigen_tag 中找到对应的位置
  position_in_eigen_tag <- match(chr_in_Data_Tagi, eigen_tag)
  print(position_in_eigen_tag)
  
  #the eigenvalue num
  k=1
  b0<-summary(lm(rsq ~ snp_1,Data))$coefficients[1,1]
  b1<-summary(lm(rsq ~ snp_1,Data))$coefficients[2,1]
  summary(lm(rsq ~ snp_1,Data))$r.squared
  SNP_vs_Rsq[i,2]<-summary(lm(rsq ~ snp_1,Data))$r.squared
  Rsq<-SNP_vs_Rsq[i,2]
  rho<-cor(Data$rsq,Data$snp_1)
  SNP_R<-cor(Data$snp_1,Data$rsq)
  SNP_vs_Rsq[i,1]<-cor(Data$rsq,Data$snp_1)
  post_cor=-1
  
  #stopping rule
  latent_thres<-1
  eigenvalue<-read.table(Eigenvalue_Files[1])
  while ( abs((SNP_vs_Rsq[i,2]-post_cor))>0.00001 & eigenvalue[k,1]>1 & k <=length(eigenvalue[,1]) &(latent_thres<eigenvalue[k,1])) {
    
    for (j in 1:length(Eigenvalue_Files)) {

      eigenvalue<-read.table(Eigenvalue_Files[position_in_eigen_tag[j]])
      beta=(eigenvalue[c(1:k),1]-latent_thres)^2
      
      if (j==1){
        ratio<-sum(beta)/((n+1)*n)
      }else{
        ratio<-c(ratio,sum(beta)/((n+1)*n))
      }
      
    }
    
    Data_peeling<- Data
    
    Data_peeling$rsq<-Data_peeling$rsq-ratio
    
    if(k==1){
      post_cor<--1
      rsq_record<-SNP_vs_Rsq[i,2]
      r_record<-SNP_R
    }else{
      post_cor<-SNP_vs_Rsq[i,2]
      rsq_record<-c(rsq_record,SNP_vs_Rsq[i,2])
      r_record<-c(r_record,SNP_R)
    }
    
    #SNP_vs_Rsq[i,2] <- cor(Data_peeling$rsq,Data_peeling$snp_1)
    SNP_vs_Rsq[i,2] <-summary(lm(rsq ~ snp_1,Data_peeling))$r.squared
    SNP_R<-cor(Data_peeling$rsq,Data_peeling$snp_1)
    #eigenvalue num +`1
    k<-k+1
  }
  rsq_record<-c(rsq_record,SNP_vs_Rsq[i,2])
  r_record<-c(r_record,SNP_R)
  
  rsq_record<-as.data.frame(cbind(rsq_record[-1],c(1:(length(rsq_record)-1)),r_record[-1]))
  colnames(rsq_record)<-c("R","Eigenvalue_num","R_cor")
  max(rsq_record$R)
  
  ##**********choose the max rsq******##
  
  max_rsq<-which(rsq_record$R==max(rsq_record$R))
  
  #set the best latent_thres
  legend_data <- data.frame(
    Eigenvalue_num = c(NA, NA),
    R = c(NA, NA),
    color = c("black", "#357EBD99")
  )
  
  for (j in 1:length(Eigenvalue_Files)) {
    
    eigenvalue<-read.table(Eigenvalue_Files[position_in_eigen_tag[j]])
    beta=(eigenvalue[c(1:max_rsq),1]-latent_thres)^2
    
    if (j==1){
      ratio<-sum(beta)/((n+1)*n)
    }else{
      ratio<-c(ratio,sum(beta)/((n+1)*n))
    }
    
  }
  Data_peeling<- Data
  
  Data_peeling$rsq<-Data_peeling$rsq-ratio
  SNP_vs_Rsq[i,2] <- cor(Data_peeling$rsq,Data_peeling$snp_1)
  SNP_vs_Rsq[i,3]<--log(sum(Data$SNPi^2*Data$rsq)/sum(Data$SNPi^2),10)
  SNP_vs_Rsq[i,4]<--log(sum(Data_peeling$SNPi^2*Data_peeling$rsq)/sum(Data_peeling$SNPi^2),10)
  ##**********choose the max rsq******##
  
  ##for inter##
  for (j in 1:length(Eigenvalue_Files)) {
    
    eigenvalue<-read.table(Eigenvalue_Files[j])
    
    if (j==1){
      top_eigenvalue<-cbind(eigenvalue[1:max_rsq,1],rep(eigen_tag[j],max_rsq))
    }else{
      top_eigenvalue<-rbind(top_eigenvalue,cbind(eigenvalue[1:max_rsq,1],rep(eigen_tag[j],max_rsq)))
    }
  }
  top_eigenvalue<-as.data.frame(top_eigenvalue)
  colnames(top_eigenvalue)<-c("value","chr")
  inter$lamda<-mapply(function(chr_i, chr_j) {
    sum(as.numeric(top_eigenvalue$value[which(top_eigenvalue$chr == chr_i)]) *
          as.numeric(top_eigenvalue$value[which(top_eigenvalue$chr == chr_j)]))
  }, inter$CHRi, inter$CHRj)
  inter$snp_1<-(inter$lamda-n)/((n+1)*n)
  inter$lamda<-(inter$lamda-n)/((n+1)*n)
  
  for (k in 1:(max_rsq+1)) {
    if(k == 1){
      ratio = 0  
    } else {
      for (j in 1:length(Eigenvalue_Files)) {
        eigenvalue <- read.table(Eigenvalue_Files[position_in_eigen_tag[j]])
        beta <- (eigenvalue[c(1:k-1),1]-1)^2
        
        if (j == 1){
          ratio <- sum(beta)/(n*n)
        } else {
          ratio <- c(ratio, sum(beta)/(n*n))
        }
      }
    }
    
    Data_peeling <- Data
    Data_peeling$rsq <- Data_peeling$rsq - ratio
    norm_2<-inter
    for (j in 1:length(Eigenvalue_Files)) {
      eigenvalue <- read.table(Eigenvalue_Files[j])
      
      if (j == 1){
        top_eigenvalue <- cbind(eigenvalue[k:max_rsq,1], rep(eigen_tag[j], max_rsq-k+1))
      } else {
        top_eigenvalue <- rbind(top_eigenvalue, cbind(eigenvalue[k:max_rsq,1], rep(eigen_tag[j], max_rsq-k+1)))
      }
    }
    
    if(k!=max_rsq+1){
      top_eigenvalue <- as.data.frame(top_eigenvalue)
      colnames(top_eigenvalue) <- c("value", "chr")
      norm_2$lamda <- mapply(function(chr_i, chr_j) {
        sum(as.numeric(top_eigenvalue$value[which(top_eigenvalue$chr == chr_i)]) *
              as.numeric(top_eigenvalue$value[which(top_eigenvalue$chr == chr_j)]))
      }, norm_2$CHRi, norm_2$CHRj)
      norm_2$snp_1 <- (norm_2$lamda - (max_rsq-k+1)) / ((n) * n)
      norm_2$lamda <- (norm_2$lamda - (max_rsq-k+1)) / ((n) * n)
      #norm_12<-rbind(norm_12,c(k-1,summary(lm(rsq ~ snp_1,Data_1KG_peeling))$r.squared ,summary(lm(rsq ~ snp_1-1, norm_2))$r.squared) )
      norm_12<-rbind(norm_12,c(k-1,cor(Data_peeling$rsq,Data_peeling$snp_1) ,cor(norm_2$rsq,norm_2$snp_1) ))
    }else{
      #norm_12<-rbind(norm_12,c(k-1,summary(lm(rsq ~ snp_1,Data_1KG_peeling))$r.squared ,0 ))
      norm_12<-rbind(norm_12,c(k-1,cor(Data_peeling$rsq,Data_peeling$snp_1) ,NA ))
    }
  }
  norm_12<-as.data.frame(norm_12)
  colnames(norm_12)<-c("Eigenvalue_num","Norm I","Norm II")
  
  norm_plot<-ggplot(norm_12, aes(x = Eigenvalue_num)) + 
    geom_point(aes(y = `Norm I`, color = "Norm I"), size = 3) +
    geom_point(aes(y = (`Norm II` + 1 )/2, color = "Norm II"), size = 3) +
    geom_line(aes(y = `Norm I`, color = "Norm I"), size = 1.5) +
    geom_line(aes(y = (`Norm II` + 1)/2, color = "Norm II"), size = 1.5) +
    theme_minimal() +
    theme(panel.grid = element_blank(),legend.position = c(0.9,0.5),
          text = element_text(family = 'serif', colour = 'black'),
          axis.ticks.x = element_line("black"), axis.line = element_line(color = "black"),
          axis.ticks.y = element_line("black"),legend.text = element_text(size=10)) +
    scale_y_continuous(limits = c(0,1),sec.axis = sec_axis(~ . * 2 - 1, name = as.expression(substitute(italic(rho)[2])))) +
    scale_x_continuous(breaks = pretty_breaks())+
    labs(x = "The number of Eigenvalue",
         y = as.expression(substitute(italic(rho)[1])),title = paste0(strsplit(List_Files[i], split = ".xld")[[1]][1]))+
    scale_color_manual(
      name = NULL,
      values = c("Norm I" = "black", "Norm II" = "#357EBD99"),
      labels = c(
        `Norm I` = as.expression(substitute(italic(rho)[1])),
        `Norm II` = as.expression(substitute(italic(rho)[2]))
      )
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  
  
  #plot for raw peeling inter
  
  y_ratio<-as.numeric(10^(SNP_vs_Rsq[i,4]-SNP_vs_Rsq[i,3]))
  Data$raw_rsq<-Data$rsq
  Data$rsq<-Data_peeling$rsq

  b0_peeling<-summary(lm(rsq ~ snp_1,Data_peeling))$coefficients[1,1]
  b1_peeling<-summary(lm(rsq ~ snp_1,Data_peeling))$coefficients[2,1]
  Rsq_peeling<-summary(lm(rsq ~ snp_1,Data_peeling))$r.squared
  rho_peeling<-cor(Data_peeling$rsq,Data_peeling$snp_1)

  scatter_plot_1 <- ggplot(Data, aes(x = snp_1)) +
    geom_point(aes(y = raw_rsq), color = "#ADB6B699", size = 8) +
    geom_smooth(method = "lm", se = TRUE, aes(y = raw_rsq), color = "#ADB6B699", fill = "#ADB6B699") +
    geom_text(aes(y=raw_rsq,label = CHRi,family='serif'),color="white") +
    geom_point(aes(y = rsq * y_ratio), color = ifelse(SNP_vs_Rsq[i, 2] > 0, "#357EBD99", "#FF9999"), size = 8) +
    geom_smooth(method = "lm", se = TRUE, aes(y = rsq * y_ratio), color = ifelse(SNP_vs_Rsq[i, 2] > 0, "#357EBD99", "#FF9999"), fill = ifelse(SNP_vs_Rsq[i, 2] > 0, "#357EBD99", "#FF9999")) +
    geom_text(aes(y = rsq * y_ratio, label = CHRi, family = 'serif')) +
    scale_y_continuous(
      labels = number_format(accuracy = 1e-4),
      sec.axis = sec_axis(~ . * (1 / y_ratio), name = "After peeling", labels = number_format(accuracy = 1e-5)),
      expand = c(0, 0)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    geom_text(x = median(Data$snp_1), y = (3 * max(Data$rsq) + median(Data$rsq)) / 4 * y_ratio, label = lm_eqn1(Data), parse = TRUE, size = 5, family = 'serif') +
    geom_text(x = median(Data$snp_1), y = (2 * max(Data$rsq) + median(Data$rsq)) / 3 * y_ratio, label = lm_eqn2(Data), parse = TRUE, size = 5, family = 'serif') +
    labs(x = "1/Number of SNPs", y = "Before peeling") +
    theme_minimal() +
    coord_cartesian(clip = 'off') +
    theme(
      plot.margin = margin(2, 10, 2, 2),
      panel.background = element_rect(fill = 'transparent', colour = NA),
      panel.grid = element_blank(),
      text = element_text(family = 'serif', colour = 'black'),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black")
    ) +
    labs(title = paste0(strsplit(List_Files[i], split = ".xld")[[1]][1],"-Norm I"))
  
  R_inter<-cor(inter$rsq,inter$lamda)
  SNP_vs_Rsq[i,6]<- -log(sum(as.numeric(inter$SNPi)*as.numeric(inter$SNPj)*as.numeric(inter$rsq))/sum(as.numeric(inter$SNPi)*as.numeric(inter$SNPj)),10)
  SNP_vs_Rsq[i,5]<-R_inter
  scatter_plot_3 <- ggplot(inter, aes(x = snp_1, y = rsq)) +
    geom_point(color=ifelse(R_inter>0.5,"#357EBD99","#FF9999"),size=8) +
    geom_smooth(method = "lm", se = TRUE, color = ifelse(R_inter>0.5,"#357EBD99","#FF9999"),fill = "#ADB6B699") +
    geom_text(aes(label = paste0(CHRi,"-",CHRj)),family='serif') +
    scale_y_continuous(expand = c(0,0))+
    scale_x_continuous(expand = c(0,0))+
    geom_text(x = (median(inter$lamda)+2*min(inter$lamda))/3, y = (3*max(inter$rsq)+ median(inter$rsq))/4, label = lm_eqn1_inter(inter), parse = TRUE, size = 5, family = 'serif',check_overlap = TRUE)+
    geom_text(x =  (median(inter$lamda)+2*min(inter$lamda))/3, y =  (2*max(inter$rsq)+ median(inter$rsq))/3, label = lm_eqn2_inter(inter), parse = TRUE, size = 5, family = 'serif',check_overlap = TRUE)+
    labs(x = "lamda", y = "Inter-c LD") +
    theme_minimal()+
    coord_cartesian(clip = 'off')+
    theme(plot.margin=margin(2,10,2,2),
          panel.background = element_rect(fill='transparent',colour = NA),
          panel.grid = element_blank(),
          text = element_text(family = 'serif', colour = 'black'),
          axis.line = element_line(colour = "black"),
          axis.ticks = element_line(colour = "black"))+
    labs(title =  paste0(strsplit(List_Files[i], split = ".xld")[[1]][1],"-Norm II") )

  
}

pdf(paste0(strsplit(List_Files[i], split = ".xld")[[1]][1],"_norm_high-resolutionLD.pdf"), width = 10, height = 10)
grid.newpage() 
pushViewport(viewport(layout = grid.layout(2, 2)))
print(P2, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(scatter_plot_1, vp = viewport(layout.pos.row = 1,layout.pos.col = 2))
print(scatter_plot_3, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(norm_plot, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
dev.off()

pdf(paste0("All_",strsplit(input,split = ".xld")[[1]][1],"_LD_block.pdf"), 20,19)
print(P2) 
dev.off()




