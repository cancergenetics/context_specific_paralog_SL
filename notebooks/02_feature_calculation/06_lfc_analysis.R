## lfc_df in the following format:
## gene1, gene2, the cell lines 
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(dplyr)
library(ggpubr)
library(tidyr)
library(readxl)


LFC_analysis_script <- function(lfc_df, cell_lines, pc_genes, nc, data_name)
{
  all_genes = c(lfc_df %>% pull(gene1), lfc_df %>% pull(gene2)) %>% unique()
  pc_genes = pc_genes[pc_genes %in% all_genes]
  df_long <- pivot_longer(lfc_df, 
                          cols = cell_lines, 
                          names_to = "Cell_Line", 
                          values_to = "LFC")
  #df_long = df_long %>% filter(!gene1 %in% nc, !gene2 %in% nc)
  df_long = df_long %>% 
    #filter(!(gene1 %in% nc & gene2 %in% nc)) %>% 
    mutate(
      essential = case_when(
        (gene1 %in% pc_genes & !(gene2 %in% pc_genes)) | (!(gene1 %in% pc_genes) & gene2 %in% pc_genes) ~ "1 essential",
        gene1 %in% pc_genes & gene2 %in% pc_genes ~ "2 essential",
        gene1 %in% nc & gene2 %in% nc ~ "2 non-targetting",
        (gene1 %in% nc & !(gene2 %in% nc)) | (!(gene1 %in% nc) & gene2 %in% nc) ~ "SKO",
        
        TRUE ~ "Others"
      )
    )

  meds = df_long %>% group_by(essential) %>% summarise(m = median(LFC), n = n())
  print(meds)
  mid_point_log = (meds %>% filter(essential == '1 essential') %>% pull(m) + 
                     meds %>% filter(essential == 'Others') %>% pull(m)) /2
  print(mid_point_log)
#  mid_point_linear = log2((2^(meds %>% filter(essential == '1 essential') %>% pull(m)) + 
 #                                2^(meds %>% filter(essential == 'Others') %>% pull(m))) /2)
 # print(mid_point_linear)
  line_data <- data.frame(
    yintercept = c(mid_point_log),#, mid_point_linear),
    label = c("Mid-point")#, "Mid-point Linear")
  )
  # Create the box plot
  b = ggplot(df_long, aes(x = essential, y = LFC, fill = essential)) +
    geom_boxplot() +
    stat_compare_means(comparisons = list(c("1 essential", "Others")),method = "t.test",
                       label = "p.signif") + # Add significance labels
    stat_compare_means(label.y = 1.5) + # Overall p-value
    
    labs(
      title = data_name,
      x = "Group",
      y = "LFC"
    ) +
 #   stat_summary(
#      fun = median,  # Calculate the median
#      geom = "text", # Add it as text
#      aes(label = round(..y.., 2)), # Display the median, rounded to 2 decimal places
#      color = "red", # Color for the text
#      vjust = -0.5   # Adjust the vertical position of the label
#    ) +
    scale_fill_brewer(palette = "Set3") +  # Set3 palette for fill
    
    theme_minimal()+
    annotate("text", y = mid_point_log+ 0.2 , label = paste0("mid-point = ", round(mid_point_log,2)),x = 1.5, size = 3 )+
    geom_hline(data = line_data, aes(yintercept = yintercept, linetype = label), color = "black") +
    scale_linetype_manual(
      name = "Reference Lines",
      values = c("Mid-point" = 2, "Mid-point Linear" = 3)
    )
  
  print(b)
  ggsave(plot = b, paste0(getwd(),"\\LFC_ANALYSIS\\",data_name, ".jpeg"),width = 10, height = 10, units = "in")
  
 # scaled_df_long = df_long %>% 
 #   mutate(LFC = LFC - mid_point_log) 
  
#  df_long = df_long %>%
#    pivot_wider(
#      names_from = Cell_Line,  # Column to create new columns
#      values_from = LFC      # Column to populate the new columns
#    ) %>%
#    select( -essential)

 # return(scaled_lfc)
  return(b)
  
}


## ITO ##
Ito_LFC = read.csv("../data/input/LFC/ITO_LFC.csv") %>% separate(gene_pair, into = c("gene1", "gene2"), sep = "_")
Ito_NC = c("AAVS1")
Ito_PC = read.delim("CEGv2.txt", stringsAsFactors = F) %>% 
  pull(GENE)

Ito_cell_lines = Ito_LFC %>% select(-X,-gene1,-gene2) %>% names()
ito_genes = c(Ito_LFC %>% pull(gene1), Ito_LFC %>% pull(gene2)) %>% unique()
ito_pc_ = Ito_PC[which(Ito_PC %in% ito_genes)]
LFC_analysis_script(Ito_LFC, cell_lines = Ito_cell_lines,pc_genes = ito_pc_,nc = Ito_NC, "Ito Screen")

# Get gemini scores for ito
ito_gemini = read.csv("../data/input/GEMINI/Gemini_ITO_Sensitive_Lethality.csv") %>%
  separate(X, into = c("gene1", "gene2"), sep = ";") %>%
  pivot_longer(cols =3:13, names_to = "Cell_Line") %>%
  rename(Score = value) 

Ito_LFC = read.csv("../data/input/LFC/ITO_LFC.csv") %>% 
  separate(gene_pair, into = c("gene1", "gene2"), sep = "_") %>%
  pivot_longer(cols = 4:14,names_to = "Cell_Line") %>%
  filter(gene1 != "AAVS1" , gene2 != "AAVS1") %>%
  rename(LFC = value)%>%
  select(-X)

merged_ito = merge(ito_gemini, Ito_LFC,by = c("gene1", "gene2", "Cell_Line")) %>%
  drop_na() %>%
  mutate(top_10 = Score >= quantile(Score, 0.9)) %>%
  mutate(threshold = LFC < -0.51 & top_10) %>%
  mutate(status = ifelse(top_10, "Top 10", NA)) %>%
  mutate(status = ifelse(threshold, "Top 10 and LFC < xx", status)) %>%
  mutate(status = ifelse(top_10 & LFC >= -0.51, "Ambigious", status)) %>%
  mutate(status = ifelse(!top_10, "Not SL", status))

cor.test(merged_ito$Score, merged_ito$LFC)

ggplot(merged_ito %>% filter(Cell_Line == "A549"))+
  geom_point(aes(x = Score, y = LFC, color = status )) +
  geom_hline(aes(yintercept = -0.51 ))


# Now add gemini FDR

gemini_fdr = read_excel("../../Ito et al/ito_published_gemini_gis.xlsx", sheet = "Supplementary table 10",skip = 2)
gemini_fdr = gemini_fdr %>%
  rename(gene_pair = ...1) %>%
  separate(gene_pair, into = c("gene1", "gene2"), sep = ";") %>%
  rename_with(~ str_remove(., "_.*")) %>%
  pivot_longer(cols = 3:13,names_to = "Cell_Line") %>%
  mutate(Cell_Line = ifelse(Cell_Line == "GI1", "GI1_004",Cell_Line)) %>%
  mutate(Cell_Line = ifelse(Cell_Line == "MEL202", "MEL202_003",Cell_Line)) %>%
  rename(fdr = value)
  
  
ito_gemini_study = read_excel("../../Ito et al/ito_published_gemini_gis.xlsx", sheet = "Supplementary table 8",skip = 2) %>%
  rename(gene_pair = ...1) %>%
  separate(gene_pair, into = c("gene1", "gene2"), sep = ";") %>%
  rename_with(~ str_remove(., "_.*")) %>%
  pivot_longer(cols = 3:13,names_to = "Cell_Line") %>%
  mutate(Cell_Line = ifelse(Cell_Line == "GI1", "GI1_004",Cell_Line)) %>%
  mutate(Cell_Line = ifelse(Cell_Line == "MEL202", "MEL202_003",Cell_Line)) 


## Check if scores comply 
theirs_vs_mine = ito_gemini_study %>% 
  left_join(ito_gemini, by = c("gene1", "gene2", "Cell_Line"))
## Something going on here

merged_ito = merged_ito %>%
  left_join(gemini_fdr,by = c("gene1", "gene2", "Cell_Line")) 


ggplot(merged_ito %>% drop_na())+
  geom_point(aes(x = Score, y = LFC, color = fdr < 0.05 )) +
  geom_hline(aes(yintercept = -0.51 )) +
  geom_vline(aes(xintercept = quantile(Score, 0.9)))

## Klingbeil Screen 

kb = read_excel("cd-23-1529_supplementary_table_s2_suppst2.xlsx")

kb = kb %>%
  rename(gene1 = `GENE 1`,
         gene2 = `GENE 2`) %>%
  select(-domain_1, -domain_2) %>% 
  mutate(gene1 = str_replace(gene1, "POSITIVE_CONTROL_", "")) %>%
  mutate(gene2 = str_replace(gene2, "POSITIVE_CONTROL_", ""))


known_essentials = read.delim("CEGv2.txt", stringsAsFactors = F) %>% 
  pull(GENE)
kb_genes = c(kb %>% pull(gene1), kb %>% pull(gene2)) %>% unique()

#kb_PC = known_essentials[known_essentials %in% kb_genes]
kb_NC = c("CUTTING_CONTROL", "NONCUTTING")
kb_cell_lines = kb %>% select(-gene1,-gene2) %>% names()

# Positive control for Klienblin
kb_data = readxl::read_excel("../../Kleinbein et al/guides.xlsx", sheet = 1)
kb_pc = kb_data %>% filter(gene1 == "POS_CONTROL" | gene2 == "POS_CONTROL") %>%
  select(gene1, gene2, guideID1) %>%
  
  separate(guideID1, into = c("A", "B"), sep = "_") %>%
  pull(A) %>% unique()

kb_pc = c(known_essentials[which(known_essentials %in% kb_genes)], kb_pc)
LFC_analysis_script(kb, cell_lines = kb_cell_lines,pc_genes = kb_pc,nc = kb_NC, "Klienblien Screen")
rm(kb_data)

##Parrish Screen
Parrish_LFC = read.csv("../data/input/LFC/Parrish_LFC.csv") %>% separate(gene_pair, into = c("gene1", "gene2"), sep = "_")
Parrish_essentials = read.csv("../Parrish/AchillesCommonEssentialControls.csv")
Parrish_essentials <- Parrish_essentials %>%
  separate(Gene, into = c("gene_symbol", "entrez_id"), sep = "\\s*\\(", convert = TRUE) %>%
  mutate(entrez_id = str_replace(entrez_id, "\\)", ""))  %>% select(gene_symbol) %>% unlist() %>% unname()

parrish_genes = unique(c(Parrish_LFC %>% pull(gene1), Parrish_LFC %>% pull(gene2)))

Parrish_NC = NA#c("AAVS1")

Parrish_cell_lines = Parrish_LFC %>% select(-X,-gene1,-gene2) %>% names()
Parrish_essentials = c(known_essentials[which(known_essentials %in% parrish_genes)], Parrish_essentials)

LFC_analysis_script(Parrish_LFC, cell_lines = Parrish_cell_lines,pc_genes = c(known_essentials,Parrish_essentials),nc = NA, "Parrish Screen")


## For each screen, Get LFC , then see how it varies for pc_genes and nc_genes

Ito_LFC = read.csv("../data/input/LFC/ITO_LFC.csv") %>% separate(gene_pair, into = c("gene1", "gene2"), sep = "_")
Ito_NC = "AAVS1"
Ito_PC = read.delim("CEGv2.txt", stringsAsFactors = F) %>% 
  pull(GENE)
ito_genes = c(Ito_LFC %>% pull(gene1), Ito_LFC %>% pull(gene2)) %>% unique()
Ito_PC = Ito_PC[Ito_PC %in% ito_genes]

df_long <- pivot_longer(Ito_LFC, 
                        cols = c("Meljuso","GI1_004","MEL202_003","PK1","MEWO",
                                          "HS944T","IPC298","A549","HSC5","HS936T","PATU8988S"), 
                        names_to = "Cell_Line", 
                        values_to = "Value")

Hart_nonessential <- read.delim("NEGv1.txt", stringsAsFactors = F) %>% 
  pull(Gene)
Ito_NC = Hart_nonessential[Hart_nonessential %in% ito_genes]

df_long = df_long %>% 
  mutate(
    essential = case_when(
      (gene1 %in% Ito_PC & !(gene2 %in% Ito_PC)) | (!(gene1 %in% Ito_PC) & gene2 %in% Ito_PC) ~ "1 essential",
      gene1 %in% Ito_PC & gene2 %in% Ito_PC ~ "2 essential",
      #gene1 %in% Ito_NC &  gene2 %in% Ito_NC ~ "2 non-essential",
       
      TRUE ~ "Others"
    )
  ) 

df_long = df_long %>% filter(gene1 != "AAVS1", gene2 != "AAVS1")
# Create the strip plot
meds = df_long %>% group_by(essential) %>% summarise(m = median(Value))
a = ggplot(df_long , aes(x = Cell_Line, y = Value, color = essential)) +
  geom_jitter(
    
    width = 0.2, alpha = 0.5) + # Adds jitter to spread points and avoid overlap

  geom_hline(yintercept = -0.8) +
  annotate("text", x = 2, y = -0.8 + 0.1, label = "Threshold = -0.8", color = "black", size = 3, hjust = 0) + # Annotate the line
  geom_hline(yintercept = -0.459) +
  annotate("text", x = 2, y = -0.459 + 0.1, label = "Threshold = -0.459", color = "black", size = 3, hjust = 0) + # Annotate the line
  
  labs(title = "Ito cell lines",
       x = "Cell Line",
       y = "LFC")+
  scale_color_manual(
    values = c(
      "1 essential" = "#FFCDD2", # Light red
      "2 essential" = "#E53935",    # Dark red
      #"2 non-essential" = "#4CAF50", # Green
      "Others" = "lightgrey"         # Light grey
    ))+
  theme_minimal()

ggsave("Ito_strip.jpeg", plot = a, units = "in", width = 10, height = 5)

meds = df_long %>% group_by(essential) %>% summarise(m = median(Value))
mid_point_ito = (meds %>% filter(essential == '1 essential') %>% pull(m) + 
                   meds %>% filter(essential == 'Others') %>% pull(m)) /2

mid_point_ito_linear = log2((2^(meds %>% filter(essential == '1 essential') %>% pull(m)) + 
                      2^(meds %>% filter(essential == 'Others') %>% pull(m))) /2)
line_data <- data.frame(
  yintercept = c(mid_point_ito, min_point_ito_linear),
  label = c("Mid-point", "Mid-point Linear")
)
# Create the box plot
b = ggplot(df_long, aes(x = essential, y = Value, fill = essential)) +
  geom_boxplot() +
  labs(
    title = "Ito (All cell lines)",
    x = "Group",
    y = "LFC"
  ) +
  scale_fill_manual(
    values = c(
      "1 essential" = "#FFCDD2", # Light red
      "2 essential" = "#E53935",    # Dark red
      "Others" = "lightgrey"         # Light grey
    ))+
  theme_minimal()+
  geom_hline(data = line_data, aes(yintercept = yintercept, linetype = label), color = "black") +
  scale_linetype_manual(
    name = "Reference Lines",
    values = c("Mid-point" = 2, "Mid-point Linear" = 3)
  )

  #facet_wrap(~Cell_Line)
ggsave("Ito_box.jpeg", plot = b, units = "in", width = 10, height = 5)

c = ggplot(df_long, aes(x = essential, y = Value, fill = essential)) +
  geom_boxplot() +
  labs(
    title = "Ito ",
    x = "Group",
    y = "LFC"
  ) +
  scale_fill_manual(
    values = c(
      "1 essential" = "#FFCDD2", # Light red
      "2 essential" = "#E53935",    # Dark red
      "2 non-essential" = "#4CAF50", # Green
      "Others" = "lightgrey"         # Light grey
    ))+
  theme_minimal() +
  facet_wrap(~Cell_Line)

## Repeat for Klingbeil Screen

kb = xlsx::read.xlsx("cd-23-1529_supplementary_table_s2_suppst2.xlsx", sheetIndex = 1)

kb = kb %>%
  rename(gene1 = GENE.1,
         gene2 = GENE.2) %>%
  select(-domain_1, -domain_2) %>% 
  filter(!str_detect(gene1, "CONTROL")) 
known_essentials = read.delim("../../Ito et al/Data/Data/CEGv2.txt", stringsAsFactors = F) %>% 
  pull(GENE)
kb_genes = c(kb %>% pull(gene1), kb %>% pull(gene2)) %>% unique
kb_PC = known_essentials[known_essentials %in% kb_genes]


Hart_nonessential <- read.delim("../../Ito et al/Data/Data/NEGv1.txt", stringsAsFactors = F) %>% 
  pull(Gene)
kb_NC = Hart_nonessential[Hart_nonessential %in% kb_genes]

kb_long <- pivot_longer(kb, 
                        cols = c("A549","ASPC1","CORL311","H1048","H1299",
                                 "H1436","H1836","H209","HEL","HPAFII","K562",
                                 "MDAMB231","MOLM13","NOMO1","PATU8902",
                                 "RD", "RH30", "SET2", "T3M4", "THP1", "YAPC", "H211"   ), 
                        names_to = "Cell_Line", 
                        values_to = "Value")

kb_long = kb_long %>% 
  mutate(
    essential = case_when(
      (gene1 %in% kb_PC & !(gene2 %in% kb_PC)) | (!(gene1 %in% kb_PC) & gene2 %in% kb_PC) ~ "1 essential",
      gene1 %in% kb_PC & gene2 %in% kb_PC ~ "2 essential",
      gene1 %in% kb_NC &  gene2 %in% kb_NC ~ "2 non-essential",
      
      TRUE ~ "Others"
    )
  ) 

kb_long %>% group_by(essential) %>% summarise(med = median(Value))
k_a = ggplot(kb_long , aes(x = Cell_Line, y = Value, color = essential)) +
  geom_jitter(
    
    width = 0.2, alpha = 0.5) + # Adds jitter to spread points and avoid overlap
  geom_hline(yintercept =  -0.8) +
  
  labs(title = "Kleinblien cell lines",
       x = "Cell Line",
       y = "LFC")+
  scale_color_manual(
    values = c(
      "1 essential" = "#FFCDD2", # Light red
      "2 essential" = "#E53935",    # Dark red
      "2 non-essential" = "#4CAF50", # Green
      "Others" = "lightgrey"         # Light grey
    ))+
  annotate("text", x = 2, y = -0.8 + 0.2, label = "Threshold = -0.8", color = "black", size = 4, hjust = 0) + # Annotate the line
  
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) 
)

ggsave("KB_strip.jpeg", plot = k_a, units = "in", width = 10, height = 5)


k_c = ggplot(kb_long, aes(x = essential, y = Value, fill = essential)) +
  geom_boxplot() +
  labs(
    title = "Kleinblien ",
    x = "Group",
    y = "LFC"
  ) +
  scale_fill_manual(
    values = c(
      "1 essential" = "#FFCDD2", # Light red
      "2 essential" = "#E53935",    # Dark red
      "2 non-essential" = "#4CAF50", # Green
      "Others" = "lightgrey"         # Light grey
    ))+
    theme_minimal()# +
 # facet_wrap(~Cell_Line)
ggsave("KB_box.jpeg", plot = k_a, units = "in", width = 10, height = 5)



x <- 2^seq(-5, 5, length.out = 500)  # Values from 2^-5 to 2^5
midx = (max(x)+min(x))/2
y <- log2(x)  # Log base 2 of x
midy = (max(y)+min(y))/2

data <- data.frame(x = x, y = y)

# Create the plot
ggplot(data, aes(x = x, y = y)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = midy)+
  geom_vline(xintercept = midx)+
   labs(title = "Log2 Plot", x = "x (log scale)", y = "log2(x)") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray")