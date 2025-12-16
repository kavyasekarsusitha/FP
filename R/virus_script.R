# Load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
#Read BLAST input file

blast <- read_tsv("blast_final.txt", 
                  col_names = TRUE)  # use the first row as column names


#Filter alignemnets <300 bp
blast <- blast %>% 
  filter(length >= 300)

# Add orientation
blast <- blast %>% 
  mutate(orientation = ifelse(send >= sstart, "forward", "reverse"))

library(dplyr)

# Unique viruses
viruses <- unique(blast$stitle)

for (v in viruses) {
  virus_data <- blast %>% filter(stitle == v)   # keeps stitle column
  filename <- paste0(output_dir, "/", v, ".csv")  # using stitle directly
  write.csv(virus_data, file = filename, row.names = FALSE)
}

#For CMV_RNA1

cmv_rna1 <- read.csv("NC_002034.1 Cucumber mosaic virus RNA 1, complete sequence.csv")

head(cmv_rna1)

cmv_rna1_plot <- ggplot(cmv_rna1, aes(x = pmin(sstart, send), xend = pmax(sstart, send), 
                  y = qseqid, yend = qseqid, color = orientation)) +
  geom_segment(size = 2) +
  labs(title = "Alignment coverage to CMV RNA 1 reference genome",
       x = "NC_002034.1 Cucumber mosaic virus RNA 1",
       y = "Query sequence") +
  theme_minimal()

ggsave("cmv_rna1.png", plot = cmv_rna1_plot, 
       width = 10, height = 6, dpi = 300)


#For CMV_RNA2
cmv_rna2 <- read.csv("NC_002035.1 Cucumber mosaic virus RNA 2, complete sequence.csv")
head(cmv_rna2)

cmv_rna2_plot <- ggplot(cmv_rna2, aes(x = pmin(sstart, send), xend = pmax(sstart, send), 
                                      y = qseqid, yend = qseqid, color = orientation)) +
  geom_segment(size = 2) +
  labs(title = "Alignment coverage to CMV RNA 2 reference genome",
       x = "NC_002035.1 Cucumber mosaic virus RNA 2",
       y = "Query sequence") +
  theme_minimal()

ggsave("cmv_rna2.png", plot = cmv_rna2_plot, 
       width = 10, height = 6, dpi = 300)

#For CMV_RNA3
cmv_rna3 <- read.csv("NC_001440.1 Cucumber mosaic virus RNA 3, complete sequence.csv")
head(cmv_rna3)

cmv_rna3_plot <- ggplot(cmv_rna3, aes(x = pmin(sstart, send), xend = pmax(sstart, send), 
                                      y = qseqid, yend = qseqid, color = orientation)) +
  geom_segment(size = 2) +
  labs(title = "Alignment coverage to CMV RNA 3 reference genome",
       x = "NC_001440.1 Cucumber mosaic virus RNA 3",
       y = "Query sequence") +
  theme_minimal()

ggsave("cmv_rna3.png", plot = cmv_rna3_plot, 
       width = 10, height = 6, dpi = 300)

#For Chilli vienal mottle virus
chvmv <- read.csv("NC_005778.1 Chilli veinal mottle virus, complete genome.csv")
head(chvmv)

chvmv_plot <- ggplot(chvmv, aes(x = pmin(sstart, send), xend = pmax(sstart, send), 
                                      y = qseqid, yend = qseqid, color = orientation)) +
  geom_segment(size = 2) +
  labs(title = "Alignment coverage to Chilli vienal mottle virus reference genome",
       x = "NC_005778.1 Chilli veinal mottle virus",
       y = "Query sequence") +
  theme_minimal()

ggsave("chvmv.png", plot = chvmv_plot, 
       width = 10, height = 6, dpi = 300)

#For Chilli vienal mottle virus
chvmv <- read.csv("NC_005778.1 Chilli veinal mottle virus, complete genome.csv")
head(chvmv)

chvmv_plot <- ggplot(chvmv, aes(x = pmin(sstart, send), xend = pmax(sstart, send), 
                                y = qseqid, yend = qseqid, color = orientation)) +
  geom_segment(size = 2) +
  labs(title = "Alignment coverage to Chilli vienal mottle virus reference genome",
       x = "NC_005778.1 Chilli veinal mottle virus",
       y = "Query sequence") +
  theme_minimal()

ggsave("chvmv.png", plot = chvmv_plot, 
       width = 10, height = 6, dpi = 300)


#For Watermelon bud necrosis segment S virus

wbnv_s <- read.csv("NC_038288.1 Watermelon bud necrosis virus strain JT segment S, complete sequence.csv")
head(wbnv_s)

wbnv_s_plot <- ggplot(wbnv_s, aes(x = pmin(sstart, send), xend = pmax(sstart, send), 
                                y = qseqid, yend = qseqid, color = orientation)) +
  geom_segment(size = 2) +
  labs(title = "Alignment coverage to Watermelon bud necrosis virus segment S reference genome",
       x = "NC_038288.1 Watermelon bud necrosis virus strain JT segment S",
       y = "Query sequence") +
  theme_minimal()

ggsave("wbnv_s.png", plot = wbnv_s_plot, 
       width = 10, height = 6, dpi = 300)

#For Watermelon bud necrosis segment M virus

wbnv_m <- read.csv("NC_038290.1 Watermelon bud necrosis virus strain JT segment M, complete sequence.csv")
head(wbnv_m)

wbnv_m_plot <- ggplot(wbnv_m, aes(x = pmin(sstart, send), xend = pmax(sstart, send), 
                                  y = qseqid, yend = qseqid, color = orientation)) +
  geom_segment(size = 2) +
  labs(title = "Alignment coverage to Watermelon bud necrosis virus segment M reference genome",
       x = "NC_038290.1 Watermelon bud necrosis virus strain JT segment M",
       y = "Query sequence") +
  theme_minimal()

ggsave("wbnv_m.png", plot = wbnv_m_plot, 
       width = 10, height = 6, dpi = 300)

#For Tomato leaf curl Gujarat virus
tlcgv <- read.csv("NC_004559.1 Tomato leaf curl Gujarat virus - [Varanasi] segment B, complete sequence.csv")
head(tlgv)

tlcgv_plot <- ggplot(tlcgv, aes(x = pmin(sstart, send), xend = pmax(sstart, send), 
                                  y = qseqid, yend = qseqid, color = orientation)) +
  geom_segment(size = 2) +
  labs(title = "Alignment coverage to tomato leaf curl gujarat virus reference genome",
       x = "NC_004559.1 Tomato leaf curl Gujarat virus",
       y = "Query sequence") +
  theme_minimal()

ggsave("tlcgv.png", plot = tlcgv_plot, 
       width = 10, height = 6, dpi = 300)

library(readxl)

#Proportion of viruses based on genetic material
gen_mat <- read_excel("Genetic_material.xls")

head(gen_mat)

# Summarize counts per genetic material
virus_counts <- gen_mat %>%
  group_by(genmat) %>%
  summarise(count = n()) %>%
  ungroup()

library(ggplot2)

# Add percentage labels
virus_counts <- virus_counts %>%
  mutate(percent = round(count / sum(count) * 100, 1),
         label = paste0(genmat, " (", percent, "%)"))

# Plot
ggplot(virus_counts, aes(x = "", y = count, fill = genmat)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
  labs(title = "Proportion of Virus Types") +
  theme_void() +
  theme(legend.position = "none")

library(dplyr)
library(ggplot2)

# Prepare angles and label positions
virus_counts <- virus_counts %>%
  arrange(desc(genmat)) %>%
  mutate(
    percent = round(count / sum(count) * 100, 1),
    label = paste0(genmat, " (", percent, "%)"),
    ypos = cumsum(count) - count/2,
    angle = 90 - 360 * ypos / sum(count)
  )

# Custom colors  
mycols <- c(
  "DNA" = "#377eb8",
  "RNA" = "#4daf4a",
  "satellite" = "#e41a1c"
)

# Plot
pie <- ggplot(virus_counts, aes(x = "", y = count, fill = genmat)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = mycols) +
  coord_polar(theta = "y") +
  geom_text(
    aes(y = ypos, label = label, angle = angle),
    hjust = 0.5,
    size = 4.5
  ) +
  labs(title = "Proportion of Virus Types") +
  theme_void() +
  theme(legend.position = "none")

ggsave("pie.png", plot = pie, 
       width = 10, height = 6, dpi = 300)




