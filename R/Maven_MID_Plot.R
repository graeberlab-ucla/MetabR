#' Maven_MID_plot
#'
#' @param metabolites
#' @param df
#' @param repeats
#' @param n number of conditions
#' @param type
#' @param index
#' @param title_type
#' @param only_M0
#'
#' @return
#' @export
#'
#' @examples
Maven_MID_plot <- function(metabolites, df, repeats, n, type, index = NULL, title_type, only_M0 = NA)
{
  if (length(metabolites) > 1) {
    if (title_type == "nonpathway1") {
      ending <- "other metabolites - I"
    }
    else if (title_type == "nonpathway2") {
      ending <- "other metabolites - II"
    }
    else if (title_type == "nonpathway3") {
      ending <- "other metabolites - III"
    }
    else if (title_type == "medium") {
      ending <- "metabolites found in medium"
    }
    else {
      ending <- "select metabolites"
    }
  }
  else if (metabolites == "glycolysis") {
    metabolites = glycolysis
    ending = "glycolytic metabolites"
  }
  else if (metabolites == "TCA") {
    metabolites = TCA
    ending <- "TCA metabolites"
  }
  else if (metabolites == "PPP") {
    metabolites = PPP
    ending <- "PPP metabolites"
  }
  else if (metabolites == "Curr") {
    metabolites <- Curr
    ending <- "Currency metabolites"
  }
  else if (metabolites == "Cys") {
    metabolites <- Cys
    ending <- "Cysteine metabolites"
  }
  else if (metabolites == "Adenine") {
    metabolites <- Adenine
    ending <- "Adenosine derivatives"
  }
  else if (metabolites == "Cytosine") {
    metabolites <- Cytosine
    ending <- "Cytidine derivatives"
  }
  else if (metabolites == "Guanine") {
    metabolites <- Guanine
    ending <- "Guanine derivatives"
  }
  else if (metabolites == "Thymine") {
    metabolites <- Thymine
    ending <- "Thymine derivatives"
  }
  else if (metabolites == "Uracil") {
    metabolites <- Uracil
    ending <- "Uracil derivatives"
  }
  else if (metabolites == "AA") {
    metabolites <- AA
    ending <- "Amino Acids"
  }
  else if (metabolites == "Hex") {
    metabolites <- Hex
    ending <- "Hexosamine metabolites"
  }
  else if (metabolites == "FA") {
    metabolites <- FA
    ending <- "Fatty Acids intermediates"
  }
  else if (metabolites == "Fru") {
    metabolites <- Fru
    ending <- "Fructose Metabolism"
  }
  else if (metabolites == "CoAs") {
    metabolites <- CoAs
    ending <- "CoA metabolism"
  }
  else if (metabolites == "Neurotrans") {
    metabolites <- Neurotrans
    ending <- "Neurotransmitter levels"
  }
  else if (metabolites=='Custom')
  {
    metabolites = Custom
    ending <- 'Custom metabolites'
  }
  else ending = ""
  if (sum(grepl("MID", names(df))) >= 1) {
    met = subset(df, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met <- mutate(met, Iso = paste(Iso, Sig, sep = "\n"),
                  Sig = "")
    met <- met %>%
      mutate(Iso = factor(Iso, levels = sort(unique(Iso))))
    Title = paste0("Isotopologue distribution of ", ending, " (Corrected for natural abundance)")
    x <- "Isotopologue"
    y <- "% Labeled"
    a <- ggplot(met, aes(Iso, RelAmounts_Ave, group = Condition,
                         fill = Condition, label = Sig))
    temp <- a$data
    # if (sum(unique(temp$Name) %in% only_M0) > 0){
    #   temp <- temp[-which(temp$Name %in% only_M0), ]
    # }
    a$data <- temp
    axis.text.x = element_text(size = 7, face = "bold", angle = 90)
    if(dim(temp)[1]==0) return()
    bar_plot_update_manual(a, met, Title, x, y, axis.text.x,
                           scales = "free", num_cond = n, type = "tf")
  }
  else if (sum(grepl("Exp", names(df))) >= 1) {
    if (type == "tf") {
      met = subset(df, Name %in% metabolites)
      stopifnot(length(unique(met$Name)) >= 1)
      met <- mutate(met, Name = paste(Name, Sig, sep = " "),
                    Sig = "")
      Title = paste0("Relative amounts of ", ending, " (not corrected for blank values)")
      x = ""
      y = "Relative Amounts"
      a <- ggplot(met, aes(Condition, RelAmounts_Ave,
                           group = Condition, fill = Condition, label = Sig))
      axis.text.x = element_blank()
      bar_plot_update_manual(a, met, Title, x, y, axis.text.x,
                             scales = "free", num_cond = n, type = type,
                             index)
    }
    else {
      met = subset(df, Name %in% metabolites)
      stopifnot(length(unique(met$Name)) >= 1)
      Title = paste0("Relative amounts of ", ending)
      x = ""
      y = "Relative Amounts"
      a <- ggplot(met, aes(Condition, RelAmounts_Ave,
                           group = Condition, fill = Condition, label = Name))
      axis.text.x = element_blank()
      bar_plot_update_manual(a, met, Title, x, y, axis.text.x,
                             scales = "free", num_cond = n, type = type,
                             index)
    }
  }
  else if (sum(grepl("Labeled", names(df))) >= 1) {
    met <- subset(df, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met$Names <- met$Name
    met <- mutate(met, Name = paste(Name, Sig, sep = " "),
                  Sig = "")
    Title <- paste0("Percent labeled in ", ending)
    x <- ""
    y <- "% Labeled"
    a <- ggplot(met, aes(Condition, RelAmounts_Ave, group = Condition,
                         fill = Condition, label = Sig))
    temp <- a$data
    # if (sum(unique(temp$Names) %in% only_M0) > 0)
    #   temp <- temp[-which((temp$Names %in% only_M0)),]
    # }
    a$data <- temp
    axis.text.x = element_blank()
    bar_plot_update_manual(a, met, Title, x, y, axis.text.x,
                           scales = "fixed", num_cond = n, type = type, index)
  }
  else if (sum(grepl("FC", names(df))) >= 1) {
    met <- subset(df, Name %in% metabolites)
    stopifnot(length(unique(met$Name)) >= 1)
    met$Names <- met$Name
    met <- mutate(met, Name = paste(Name, Sig, sep = " "),
                  Sig = "")
    Title <- paste0("Fractional Contribution to ", ending)
    x <- ""
    y <- "% Fractional Contribution"
    a <- ggplot(met, aes(Condition, RelAmounts_Ave, group = Condition,
                         fill = Condition, label = Sig))
    temp <- a$data
    # if (sum(unique(temp$Names) %in% only_M0) > 0)
    #   temp <- temp[-which((temp$Names %in% only_M0)), ]
    # }
    a$data <- temp
    axis.text.x = element_blank()
    if(dim(temp)[1]==0) return()
    bar_plot_update_manual(a, met, Title, x, y, axis.text.x,
                           scales = "fixed", num_cond = n, type = type, index)
  }
}