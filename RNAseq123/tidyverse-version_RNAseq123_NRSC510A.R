library(tidyverse)

# Set up a data frame for demonstration purposes (replace this with your actual data)
set.seed(123)

##LOADING IN THE DATA##
#this reads in the counts table
#row.names = 1 tells the code that the first row are just names

setwd("/Users/lhyannesoto/NRSC510_1/RNAseq_NRSC510/RNAseq123")
counts <- read.delim("R539_count.txt", row.names = 1)


#The original code is as follows
#genenames = cbind(counts$gene_name,rownames(counts))
#genenames = as.data.frame(genenames)


#Here is the tidyverse version
genenames <- counts %>%
  rownames_to_column(var = "ID") %>%
  select(gene_name, ID) %>%
  as_tibble()

d0 <- as.data.frame(matrix(rnorm(1000), ncol = 5))

# Define the color palette for samples
fill <- c("red", "green", "blue", "purple", "orange")

# Plotting using ggplot2 and tidyverse functions
pdf('FilteringCPM_plots_cutoff3.pdf')
par(mfrow = c(1, 2))

# Prefiltered data
lcpm <- cpm(d0, log = TRUE, prior.count = 2)
df_prefiltered <- data.frame(x = lcpm[, 1], sample = "Raw data")

for (i in 2:ncol(lcpm)) {
  df_prefiltered <- bind_rows(df_prefiltered, data.frame(x = lcpm[, i], sample = paste("Sample", i)))
}

ggplot(df_prefiltered, aes(x = x, fill = sample)) +
  geom_density(lwd = 2) +
  ylim(0, 0.5) +
  ggtitle("A. Raw data") +
  xlab("Log-cpm") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~sample) +
  theme_minimal()

# Filtered data
lcpm <- cpm(d0, log = TRUE, prior.count = 2)
df_filtered <- data.frame(x = lcpm[, 1], sample = "Filtered data")

for (i in 2:ncol(lcpm)) {
  df_filtered <- bind_rows(df_filtered, data.frame(x = lcpm[, i], sample = paste("Sample", i)))
}

ggplot(df_filtered, aes(x = x, fill = sample)) +
  geom_density(lwd = 2) +
  ylim(0, 0.5) +
  ggtitle("B. Filtered data") +
  xlab("Log-cpm") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~sample) +
  theme_minimal()

dev.off()
