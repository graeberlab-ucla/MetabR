#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param The data matrix to be used, the PC to be used for the x axis (default is 1), the PC to be used for the y axis (default is 2), and the cutoff value for minimum Correlation (default is 0.5).
#'
#' @return A pdf file with a scree plot, a pair plot with the first five PCs, and a PCA plot with PCs a and b.
#'
#' @examples make_PCA2_update(matrix = RelA)
#'
#' @export
#'

make_PCA2_update <- function(matrix, a=1, b=2, cutoff = 0.5)
{
  if (exists('Title')==F) stop('Title not specified')

  if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='Exp') {
    ext = 'Relative Amounts'
  } else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='MID') {
    ext = 'MIDs'
  } else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='FC') {
    ext = 'Fractional Contribution'
  } else if (gsub('(.)*_|[0-9]','',colnames(matrix))[1]=='Labeled') {
    ext = 'Percent Labeled'
  }

  pca <- prcomp(t(matrix), center=T, scale=T)
  var_PCs=round(summary(pca)$imp[2,]*100,1)
  CCP <- cor(scale(t(matrix), center=T), pca$x, use='pairwise') %>%
    data.frame() %>%
    mutate(Metabolite = rownames(.))

  CCP$Corr <- sqrt((CCP[,a])^2 + (CCP[,b])^2)

  PC <- pca$x %>%
    as.data.frame() %>%
    mutate(ROWNAMES = rownames(.),
           Sample.Name = gsub('Exp|MID|FC|Labeled','', ROWNAMES)) %>%
    left_join(., samples, by='Sample.Name')
  loadings <- data.frame(pca$rotation)
  loadings$Name <- rownames(loadings)
  loadings <- select(loadings, Name, everything())
  write.csv(loadings, file=paste0(Title,'-PC Loadings-',ext,'.csv'), row.names=T)
  scores=pca$x
  write.csv(scores, file=paste0(Title,'-PC Scores-',ext,'.csv'), row.names=T)

  PC.title=paste(Title,'-PC',a, ' vs. PC', b, '-PCA Plots2-', ext, '.pdf', sep='')
  pdf(file = PC.title, width=16, height=10)

  plot(var_PCs[1:min(10, length(pca$sdev))], type='b', pch=20, col='blue', ylab='Variance explained (%)', xlab='Principal Component', main='Screeplot')
  print(lattice::splom(data=PC, ~PC[,1:5],
                       groups = PC$Condition,
                       par.settings=list(superpose.symbol=list(col=colors, pch=19)),
                       auto.key=list(columns=4), pch=19))

  plot <- ggplot(PC, aes(PC[,a], PC[,b], label = Sample.Name),col=Condition)+
    geom_label_repel(color='black', fontface='bold', size=2.5, point.padding = .3,label.size = NA,segment.alpha = 0.6) +
    geom_point(size=3, shape=21, aes(fill=Condition))+
    labs(x=paste('PC',a, ' (',var_PCs[a],'%)',sep=''), y=paste('PC',b, ' (',var_PCs[b],'%)',sep=''), title=paste('PC',a,' vs. PC',b,': All samples', sep=''))+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          text = element_text(face='bold'))+
    scale_fill_manual('Condition', values=colors)

  CCP_plot <- filter(CCP, Corr >= cutoff) %>%
    ggplot(., aes(.[,a], .[,b], label=Metabolite))+
    geom_text_repel(aes(label=Metabolite, fontface=2),size=2.5,point.padding=.3,color="navy", min.segment.length= 0.1, segment.color = "grey", segment.alpha = 0.6)+
    geom_point(size=2, color='black')+
    #  geom_text(vjust=-1, color='navy', size=3)+
    labs(x=paste('PC',a, ' (',var_PCs[a],'%)',sep=''),
         y=paste('PC',b, ' (',var_PCs[b],'%)',sep=''),
         title = 'Correlation circle plot')+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          text = element_text(face='bold'))+
    xlim(min(CCP$PC1),max(CCP$PC1))+
    ylim(min(CCP$PC2),max(CCP$PC2))


  gridExtra::grid.arrange(plot, CCP_plot, nrow=1)
  dev.off()

  CCP1 <- suppressWarnings(CCP %>%
                             right_join(., Abbrev, by=c('Metabolite'='Abb')) %>%
                             select(KEGG.ID, paste0('PC',a)))
  names(CCP1)[2] <- 'Norm_Av'
  CCP1$Norm_Av[is.na(CCP1$Norm_Av)] <- 0
  write.csv(CCP1, paste0('CCP-PC', a, '-', ext,'.csv'), row.names=F)

  CCP2 <- suppressWarnings(CCP %>%
                             right_join(., Abbrev, by=c('Metabolite'='Abb')) %>%
                             select(KEGG.ID, paste0('PC',b)))
  names(CCP2)[2] <- 'Norm_Av'
  CCP2$Norm_Av[is.na(CCP2$Norm_Av)] <- 0
  write.csv(CCP2, paste0('CCP-PC', b, '-', ext,'.csv'), row.names=F)
}
