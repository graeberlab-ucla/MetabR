#' @title Functions to assist in the current UCLA Metabolomics Pipeline.
#'
#' @description This package contains functions to assist in the current UCLA Metabolomics Quality Control and Analysis Pipeline.
#'
#' @param the pathway of interest, a data frame containing the relative amounts of the different metabolites in the experiment, the number of conditions in the experiment, and the type of analysis or quality control being run.
#'
#' @return a page of plots of the relative amounts of the different metabolites that are part of the specified defined pathway ('metabolites' parameter)
#'
#' @examples bar_update_manual(metabolites = 'glycolysis', df = data4, n = num_conditions, type = "blank")
#'
#' @export
#'

bar_update_manual <- function(metabolites, df, repeats, n, type,index)
{
  if (length(metabolites) > 1)
  {
    ending <- 'select metabolites'
  }
  else if (metabolites == 'glycolysis')
  {
    metabolites = glycolysis
    ending='glycolytic metabolites'
  }
  else if (metabolites=='TCA')
  {
    metabolites = TCA
    ending <- 'TCA metabolites'
  }
  else if (metabolites=='PPP')
  {
    metabolites = PPP
    ending <- 'PPP metabolites'
  }
  else if (metabolites=='Curr')
  {
    metabolites <- Curr
    ending <- 'Currency metabolites'
  }
  else if (metabolites=='Cys')
  {
    metabolites <- Cys
    ending <- 'Cysteine metabolites'
  }
  else if (metabolites=='Adenine')
  {
    metabolites <- Adenine
    ending <- 'Adenosine derivatives'
  }
  else if (metabolites=='Cytosine')
  {
    metabolites <- Cytosine
    ending <- 'Cytidine derivatives'
  }
  else if (metabolites=='Guanine')
  {
    metabolites <- Guanine
    ending <- 'Guanine derivatives'
  }
  else if (metabolites=='Thymine')
  {
    metabolites <- Thymine
    ending <- 'Thymine derivatives'
  }
  else if (metabolites=='Uracil')
  {
    metabolites <- Uracil
    ending <- 'Uracil derivatives'
  }
  else if (metabolites=='AA')
  {
    metabolites <- AA
    ending <- 'Amino Acids'
  }
  else if (metabolites=='Hex')
  {
    metabolites <- Hex
    ending <- 'Hexosamine metabolites'
  }
  else if (metabolites=='FA')
  {
    metabolites <- FA
    ending <- 'Fatty Acids intermediates'
  }
  else if (metabolites=='Fru')
  {
    metabolites <- Fru
    ending <- 'Fructose Metabolism'
  }
  else if (metabolites=='CoAs')
  {
    metabolites <- CoAs
    ending <- 'CoA metabolism'
  }
  else if (metabolites=='Neurotrans')
  {
    metabolites <- Neurotrans
    ending <- 'Neurotransmitter levels'
  }
  else ending = ''

  if (sum(grepl('MID', names(df))) >= 1)
  {
    met = subset(df, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Iso=paste(Iso, Sig, sep='\n'),
                   Sig='') %>%
      mutate(Iso = factor(Iso, levels = paste(rep(paste('M', 0:50, sep=''), each=4),
                                              c('','*','**','***'), sep='\n')))

    Title = paste0("Isotopologue distribution of ",ending)
    x <- 'Isotopologue'
    y <- '% Labeled'
    a <-ggplot(met, aes(Iso, RelAmounts_Ave, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_text(size=11, face="bold")
    bar_plot_update_manual(a, met, Title, x, y, axis.text.x, scales='free', num_cond = n, type = type,index)
  }
  else if (sum(grepl('Exp', names(df))) >= 1)
  {
    if (type == "tf")
    {
      met = subset(df, Name %in% metabolites)
      stopifnot(length(unique(met$Name)) >= 1)
      met <-  mutate(met, Name=paste(Name, Sig, sep=' '),
                     Sig='')
      Title = paste0("Relative amounts of ",ending, " (not corrected for blank values)")
      x=''
      y='Relative Amounts'
      a <-ggplot(met, aes(Condition, RelAmounts_Ave, group=Condition, fill=Condition, label=Sig))
      axis.text.x=element_blank()
      bar_plot_update_manual(a, met, Title, x, y, axis.text.x, scales='free', num_cond = n, type = type,index)
    } else {
      met = subset(df, Name %in% metabolites)
      stopifnot(length(unique(met$Name)) >= 1)
      Title = paste0("Relative amounts of ",ending)
      x=''
      y='Relative Amounts'
      a <-ggplot(met, aes(Condition, RelAmounts_Ave, group=Condition, fill=Condition, label = Name))
      axis.text.x=element_blank()
      bar_plot_update_manual(a, met, Title, x, y, axis.text.x, scales='free', num_cond = n, type = type,index)
    }
  }

  else if (sum(grepl('Labeled', names(df))) >= 1)
  {
    met <- subset(df, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Name=paste(Name, Sig, sep=' '),
                   Sig='')
    Title <- paste0('Percent labeled in ', ending)
    x <- ''
    y <- '% Labeled'
    a <- ggplot(met, aes(Condition, RelAmounts_Ave, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_blank()
    bar_plot_update_manual(a, met, Title, x, y, axis.text.x, scales='fixed', num_cond = n, type = type,index)
  }

  else if (sum(grepl('FC', names(df))) >= 1)
  {
    met <- subset(df, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <-  mutate(met, Name=paste(Name, Sig, sep=' '),
                   Sig='')
    Title <- paste0('Fractional Contribution to ', ending)
    x <- ''
    y <- '% Fractional Contribution'
    a <- ggplot(met, aes(Condition, RelAmounts_Ave, group=Condition, fill=Condition, label=Sig))
    axis.text.x=element_blank()
    bar_plot_update_manual(a, met, Title, x, y, axis.text.x, scales='fixed', num_cond = n, type = type,index)
  }
}
