
X_LD_Visualization <- function(input=NULL,output=NULL,significant=FALSE){
    library(ggplot2)
    #  Check inputs
    if(!(file.exists(input)))
    stop("input file doesn't exist.")
    
    if(is.null(output))
      stop("output is NULL")
    
    chrLD <- na.omit(read.table(input,header=T))

    for(i in 1:nrow(chrLD)){
      if(chrLD[i,1] < chrLD[i,3]){
        tmp <- chrLD[i,1:2] 
        chrLD[i,1:2] <- chrLD[i,3:4]
        chrLD[i,3:4] <- tmp
      }
    }
    chrLD <- chrLD[order(chrLD$CHRi,chrLD$CHRj),]
    # Visulation
    # unscaled
    logchrLD <- data.frame(cbind(paste0("chr",chrLD[,1]),paste0("chr",chrLD[,3]),-log10(chrLD[,5]),chrLD[,7]))
    logchrLD[,3] <- as.numeric(logchrLD[,3])
    logchrLD[,4] <- as.numeric(logchrLD[,4])
    colnames(logchrLD) <- c("chri","chrj","-log10chrLD","p")
    logchrLD$chri <- factor(logchrLD$chri,levels = unique(logchrLD$chri))
    logchrLD$chrj <- factor(logchrLD$chrj,levels = unique(logchrLD$chrj))
    
    p0.01 <- 0.01/nrow(logchrLD)
    p0.05 <- 0.05/nrow(logchrLD)
    
    logchrLD[which(logchrLD$p <= p0.01),5] <- "**"
    logchrLD[which(logchrLD$p > p0.01 && logchrLD$p <= p0.05),5] <- "*"
    logchrLD[which(logchrLD$p > p0.05),5] <- ""
    
    colnames(logchrLD)[5] <- "significant"
    
    MAX <- max(logchrLD[,3])
    MIN <- min(logchrLD[,3])
    MID <- (MAX+MIN)/2
    
    if(significant != FALSE){
      label = logchrLD[,5]
    }else{
      label = ""
    }
    p <- ggplot(data = logchrLD, aes(logchrLD[,1], logchrLD[,2], fill = logchrLD[,3]))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "#FF0000", high = "#FFFFFF", mid = "#FF9E81",space = "Lab",
                           midpoint = MID, limit = c(MIN,MAX),name=expression(paste(-log[10],"LD"))) +
      geom_text(aes(label=label),col ="black",family="serif") +
      theme_minimal()+
      scale_y_discrete(position = "right") +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90,family="serif"),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(family="serif"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                   title.position = "top", title.hjust = 0.5)) +
      coord_fixed()
    ggsave(paste0(output,".pdf"),p,width=8,height=8)
    
    # scaled
    scaledchrLD <- data.frame(logchrLD[,c(1,2)],chrLD[,8],chrLD[,10])
    scaledchrLD[scaledchrLD[,1] == scaledchrLD[,2],4] <- NA
    scaledchrLD[,3] <- as.numeric(scaledchrLD[,3])
    scaledchrLD[,4] <- as.numeric(scaledchrLD[,4])
    colnames(scaledchrLD) <- c("chri","chrj","chrLD(scaled)","p")
    scaledchrLD$chri <- factor(scaledchrLD$chri,levels = unique(scaledchrLD$chri))
    scaledchrLD$chrj <- factor(scaledchrLD$chrj,levels = unique(scaledchrLD$chrj))
    
    
    p0.01 <- 0.01/nrow(na.omit(scaledchrLD))
    p0.05 <- 0.05/nrow(na.omit(scaledchrLD))
    
    scaledchrLD[which(scaledchrLD$p <= p0.01),5] <- "**"
    scaledchrLD[which(scaledchrLD$p > p0.01 && scaledchrLD$p <= p0.05),5] <- "*"
    scaledchrLD[which(scaledchrLD$p > p0.05 | is.na(scaledchrLD$p)),5] <- ""
    
    colnames(scaledchrLD)[5] <- "significant"
    
    if(significant != FALSE){
      label = scaledchrLD[,5]
    }else{
      label = ""
    }
    
    p <- ggplot(data = scaledchrLD, aes(scaledchrLD[,1], scaledchrLD[,2], fill = scaledchrLD[,3]))+
      geom_tile(color = "white")+
      scale_fill_gradient2(low = "#FFFFFF", high = "#FF0000",mid = "#FF9E81",
                           midpoint = 0.5, limit = c(0,1), space = "Lab",
                           name="LD (scaled)") +
      geom_text(aes(label=label),col ="black",family="serif") +
      theme_minimal()+
      scale_y_discrete(position = "right") +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90,family="serif"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(family="serif"),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(family="serif"),
        legend.title = element_text(family="serif"),
        legend.justification = c(1, 0),
        legend.position = c(0.6, 0.7),
        legend.direction = "horizontal")+
      guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                   title.position = "top", title.hjust = 0.5)) +
      coord_fixed()
    ggsave(paste0(output,".Scaled.pdf"),p,width=8,height=8)
  }
