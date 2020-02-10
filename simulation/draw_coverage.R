# parse filename
library(tidyverse)

parse_filename <- function(filename)
{
  Er <- str_match(filename, "er(.*?)p")[2] %>% as.numeric()
  beta <- str_match(filename, "b(.*?)e")[2] %>% as.numeric()
  pf1 <- str_match(filename, "p(.*?)d")[2] %>% as.numeric()
  delta <- str_match(filename, "d(.*?).r")[2] %>% as.numeric()
  
  return(c(Er, beta, pf1, delta))
}

merge_conditional_cp <- function(list_of_results)
{
  mat <- list_of_results[[1]]$Coverage_Probability
  for (i in 2:length(list_of_results))
  {
    mat <- mat + list_of_results[[i]]$Coverage_Probability
  }
  mat <- mat/length(list_of_results)
  
  cp_df <- reshape2::melt(mat)
  names(cp_df) <- c("Method", "Quantile", "Coverage")
  return(cp_df)
}

make_data_frame <- function(filename)
{
  list_of_results <- readRDS(filename)
  cp_df <- merge_conditional_cp(list_of_results)
  info_vec <- parse_filename(filename)
  
  Er <- info_vec[1]
  Beta <- info_vec[2]
  Pf1 <- info_vec[3]
  Delta <- info_vec[4]
  
  cp_df$Er <- Er
  cp_df$Beta <- Beta
  cp_df$Pf1 <- Pf1
  cp_df$Delta <- Delta
  
  LowerTail <- numeric(length(cp_df$Quantile))
  Nominal_Coverage <- numeric(length(LowerTail))
  for(i in 1:length(cp_df$Quantile))
  {
    LowerTail[i] <- substr(cp_df$Quantile[i], 1, 5)
    Nominal_Coverage[i] <- substr(cp_df$Quantile[i], 6, 7)
  }
  cp_df$LowerTail <- LowerTail
  cp_df$Nominal_Coverage <- Nominal_Coverage
  
  return(cp_df)
}

Pf1.lab <- c(0.01, 0.05, 0.1, 0.2)
names(Pf1.lab) <- c("Pf1 = 0.01", "Pf1 = 0.05", "Pf1 = 0.1", "Pf1 = 0.2")

plot.function <- function(cp_df, b, d)
{
  cp_df %>% filter( Beta == b, Delta == d) %>% 
    ggplot(aes(x = Er, y = Coverage, col = Method))+
    geom_point(aes(shape = Method), alpha = 0.75)+
    geom_line(alpha = 0.75)+
    geom_hline(aes(yintercept = as.numeric(Nominal_Coverage)/100),linetype="dashed")+
    facet_grid(Pf1~Quantile)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("Expected Number of Failures")+
    ggtitle(paste("Beta=",b,"; Delta=", d))
}

# 
setwd("./data/")
file_list <- list.files()
cp_df <- NULL
for (filename in file_list)
{
  cp_df <- rbind(cp_df, make_data_frame(filename))
}
cp_df$Pf1 <- factor(cp_df$Pf1, levels = c(0.01, 0.05, 0.1, 0.2), labels = c("Pf1 = 0.01","Pf1 = 0.05","Pf1 = 0.1","Pf1 = 0.2" ))

pdf("../four_methods_coverage.pdf",width = 11.69, height = 8.27)
plot.function(cp_df, 2, 0.1)
plot.function(cp_df, 2, 0.2)
plot.function(cp_df, 4, 0.1)
plot.function(cp_df, 4, 0.2)
dev.off()
