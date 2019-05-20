library(ggplot2)
library(stringr)
library(RColorBrewer)
library(cowplot)
library(ape)
library(ggtree)
library(gggenes)
library(Biostrings)
library(grid)
library(gridExtra)
library(dendextend)

#### create groups of colours

n <- 61
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]

cols <- c("Yersinia" = "#4F2F4F",
          "Klebsiella" = "#FF0000",
          "Saccharophagus" = "#00FF00",
          "Salmonella" = "#0000FF",
          "Teredinibacter" = "#FF00FF",
          "Pectobacterium" = "#00FFFF",
          "Erwinia" = "#FFFF00",
          "Dickeya" = "#000000",
          "Pantoea" = "#70DB93",
          "Plautia" = "#B5A642",
          "Rahnella" = "#5F9F9F",
          "Enterobacter" = "#B87333",
          "Serratia" = "#2F4F2F",
          "Shimwellia" = "#9932CD",
          "Caldisericum" = "#871F78",
          "Kosakonia" = "#855E42",
          "Bdellovibrio" = "#545454",
          "Citrobacter" = "#8E2323",
          "Providencia" = "#F5CCB0",
          "Pluralibacter" = "#238E23",
          "Morganella" = "#CD7F32",
          "Francisella" = "#DBDB70",
          "Phytobacter" = "#C0C0C0",
          "Cronobacter" = "#527F76",
          "Persicobacter" = "#9F9F5F",
          "Leclercia" = "#8E236B",
          "Lysinibacillus" = "#2F2F4F",
          "Proteus" = "#EBC79E",
          "Chryseobacterium" = "#CFB53B",
          "Gibbsiella" = "#FF7F00",
          "Escherichia" = "#DB70DB",
          "Enterobacteriaceae" = "#D9D9F3",
          "Rouxiella" = "#5959AB",
          "Lelliottia" = "#8C1717",
          "Oligella" = "#238E68",
          "Mesoplasma" = "#6B4226",
          "Buttiauxella" = "#8E6B23",
          "Metakosakonia" = "#007FFF",
          "Undibacterium" = "#00FF7F",
          "TssA" = "#3cb44b",
          "TssB" = "#ffe119",
          "TssC" = "#e6194b",
          "TssD" = "#4363d8",
          "TssE" = "#ff1493",
          "TssF" = "#911eb4",
          "TssG" = "#46f0f0",
          "TssH" = "#f032e6",
          "TssI" = "#bcf60c",
          "TssJ" = "#fabebe",
          "TssK" = "#008080",
          "TssL" = "#e6beff",
          "TssM" = "#9a6324",
          "DUF4150" = "#fffac8",
          "FHA" = "#800000",
          "ImpE" = "#aaffc3",
          "PAAR_motif" = "#808000",
          "PAAR motif" = "#808000",
          "Pkinase" = "#ffd8b1",
          "PP2C" = "#000075",
          "TagF_N" = "#808080",
          "TagF N" = "#808080",
          "X-unknown" = "#ffffff",
          "Unknown" = "#ffffff",
          "ssp6" = "#f58231",
          "1" = "#ffbfbf",
          "2" = "#ff9595",
          "3" = "#ff6a6a",
          "4" = "#ff5555",
          "5" = "#ff2b2b",
          "6" = "#ff0000")

####load in data
results <- read.csv(file = "../data/results.csv", comment.char = "", header = F, stringsAsFactors = F, quote="")
metadata <- read.csv(file="../data/NCBI_complete_genomes_only_published_strains_metadata.csv", comment.char = "", header = T, stringsAsFactors = F, quote="")
lengths <- read.csv(file="../data/lengths_formatted.csv",comment.char = "", header = F, stringsAsFactors = F, quote="")
#metadata$Short.species <- str_replace_all(metadata$Short.species, " ", "_")
metadata4later <- metadata
metadata4later$File <- str_replace_all(metadata4later$File ,'.gbff', '')


##merge data
complete <- merge(results, metadata, by.x = "V1", by.y = "Strain.name")
complete_2 <- merge(complete,lengths, by.x ="V2",by.y="V1")

#### plot hits in species:

#table for number of species:
species <- as.data.frame(table(complete_2$Short.species))

### colours for species
names(col_vector) <- levels(species$Var1)

##plots

ggplot(data = species, aes(x = reorder(Var1, -Freq), y= Freq,fill = Var1)) + geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, size = 11)) +
  xlab("Species") +
  ylab("Frequency") +
  scale_fill_manual(values=c(col_vector),name = "Species") +
  guides(fill=guide_legend(ncol =3,keywidth=0.1,
                            keyheight=0.1,
                            default.unit="inch")) +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=8, face="bold"),
        legend.justification=c(0,0),
        legend.position=c(0.45, 0.24),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  ggsave("species_with_ssp6_hits.png",dpi=300, width=9,height=7, units=c("in"))

##species in metadata:
species_of_genomes <- as.data.frame(table(metadata$Short.species))
filtered_species_of_genomes <- subset.data.frame(species_of_genomes, Freq >= 20 & Var1 != "Short species")
filtered_species_of_genomes_2 <- subset.data.frame(species_of_genomes, Freq < 20 & Freq > 7 & Var1 != "Short species")
filtered_species_of_genomes_3 <- subset.data.frame(species_of_genomes, Freq <= 7  & Var1 != "Short species")


ggplot(data = filtered_species_of_genomes, aes(x = reorder(Var1, -Freq), y= Freq)) + geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, size = 11))+
  xlab("Species") +
  ylab("Frequency") +
  ggsave("species_in_dataset_1.png",dpi=300, width=9,height=7, units=c("in"))

ggplot(data = filtered_species_of_genomes_2, aes(x = reorder(Var1, -Freq), y= Freq)) + geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, size = 11))+
  xlab("Species") +
  ylab("Frequency")+
  ggsave("species_in_dataset_2.png",dpi=300, width=9,height=7, units=c("in"))

ggplot(data = filtered_species_of_genomes_3, aes(x = reorder(Var1, -Freq), y= Freq)) + geom_col() +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90, size = 11))+
  xlab("Species") +
  ylab("Frequency")+
  ggsave("species_in_dataset_3.png",dpi=300, width=9,height=7, units=c("in"))

#### lengths of the hits against hmmer score:

# ggplot(data = complete_2, aes(V2.y,V7, colour = Short.species)) +
#   geom_count() +
#   theme_classic()

ggplot(data = complete_2, aes(V2.y,V7, colour = Short.species, alpha = 0.1)) +
  geom_count() +
  theme_classic() +
  scale_color_manual(values=c(col_vector), name = "Species") +
  theme(text = element_text(size = 16)) +
  guides(color=guide_legend(ncol =3,keywidth=0.1,
                            keyheight=0.1,
                            default.unit="inch")) +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=8, face="bold"),
        legend.justification=c(0,0),
        legend.position=c(0.45, 0.12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        #legend.position = "none",
        legend.box.background = element_rect(colour = "black")) +
  xlim(0,650) +
  xlab("Length of hmmer hit (AAs)") +
  ylab("Hmmer score") +
  geom_text(aes(label=ifelse(V1 =="DB11",as.character("Db11"),'')),hjust=1.2,vjust=0,colour='black',alpha = NA) +
  geom_hline(yintercept=20, linetype="dashed", color = "red") +
  ggsave("hit_length_vs_hmmer_score_coloured.png", dpi=300, width=9,height=7, units=c("in"))

###################  extracting names of the genomes for step 2:

files <- subset.data.frame(complete_2, select=c("File"))
files$File <- str_replace_all(files$File, ".gbff", "")
write.table(files, col.names = F, row.names = F, quote = F, file = "filenames_with_ssp6_hmmer_hits.txt")

##########

# ggplot(complete_2, aes(V7,V10)) +
#   geom_count() +
#   theme_classic()

ggplot(data = complete_2, aes(V7,V10, colour = Short.species, alpha = 0.1)) +
  geom_count() +
  theme_classic() +
  scale_color_manual(values=c(col_vector), name = "Species") +
  theme(text = element_text(size = 16)) +
  theme(legend.position = "none") +
  xlab("Complete sequence hmmer score") +
  ylab("Best domain hmmer score") +
  geom_text(aes(label=ifelse(V1 =="DB11",as.character("Db11"),'')),hjust=1.2,vjust=0,colour='black',alpha = NA) +
  geom_hline(yintercept=20, linetype="dashed", color = "red") +
  ggsave("complete_sequence_vs_best_domain_hmmer_scores.png", dpi=300, width=9,height=7, units=c("in"))


#### add in the T6SS data:

complete_2$file_prefix <- str_replace_all(files$File, ".gbff", "")

T6SS_operons <- read.csv("../data/strain_statistics.csv",comment.char = "", header = T, stringsAsFactors = F, quote="")
complete_with_T6SS <- merge(complete_2,T6SS_operons,by.x="file_prefix",by.y="strain")
least_1_T6SS <- subset.data.frame(complete_with_T6SS, number_of_operons >= 1)
no_T6SS <- subset.data.frame(complete_with_T6SS, number_of_operons == 0)

box_data <- data.frame(y=c(20,350), x=c(150,250))

### hmmer score without T6SS data:
ggplot(data = least_1_T6SS, aes(V2.y,V7, colour = Short.species, alpha = 0.1)) +
  geom_count() +
  #geom_count(data = no_T6SS, aes(V2.y,V7), colour=c("grey")) +
  theme_classic() +
  scale_color_manual(values=c(col_vector), name = "Species") +
  theme(text = element_text(size = 16)) +
  guides(color=guide_legend(ncol =3,keywidth=0.1,
                            keyheight=0.1,
                            default.unit="inch")) +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=8, face="bold"),
        legend.justification=c(0,0),
        legend.position=c(0.45, 0.12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  xlim(0,650) +
  xlab("Length of hmmer hit (AAs)") +
  ylab("Hmmer score") +
  geom_text(aes(label=ifelse(V1 =="DB11",as.character("Db11"),'')),hjust=1.2,vjust=0,colour='black',alpha = NA) +
  annotate("rect",xmin = 150, xmax = 250, ymin = 20, ymax = 350, alpha = 0.05, colour = "red", linetype = "dashed" ) +
  ggsave("hit_length_vs_hmmer_score_coloured.png", dpi=300, width=9,height=7, units=c("in"))

##hmmer score with T6SS data - shared as grey
ggplot(data = least_1_T6SS, aes(V2.y,V7, colour = Short.species, alpha = 0.1)) +
  geom_count() +
  geom_count(data = no_T6SS, aes(V2.y,V7), colour=c("grey")) +
  theme_classic() +
  scale_color_manual(values=c(col_vector), name = "Species") +
  theme(text = element_text(size = 16)) +
  guides(color=guide_legend(ncol =3,keywidth=0.1,
                            keyheight=0.1,
                            default.unit="inch")) +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=8, face="bold"),
        legend.justification=c(0,0),
        legend.position=c(0.45, 0.12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  xlim(0,650) +
  xlab("Length of hmmer hit (AAs)") +
  ylab("Hmmer score") +
  geom_text(aes(label=ifelse(V1 =="DB11",as.character("Db11"),'')),hjust=1.2,vjust=0,colour='black',alpha = NA) +
  annotate("rect",xmin = 150, xmax = 250, ymin = 20, ymax = 350, alpha = 0.05, colour = "red", linetype = "dashed" ) +
  ggsave("hit_length_vs_hmmer_score_coloured_with_T6SS_in_grey.png", dpi=300,width=9,height=7, units=c("in"))

### output list of strains with hmmer score >= 20 and at least 1 T6SS:

filtered_list <- subset.data.frame(complete_with_T6SS, number_of_operons > 0 & V7 >= 20 )
filtered_list_of_genomes <- subset.data.frame(filtered_list, select=c("file_prefix"))
#filtered_list_of_genomes <- filtered_list_of_genomes[duplicated(filtered_list_of_genomes),]
write.table(filtered_list_of_genomes, quote =F, col.names = F, row.names =F, file="hmmer_more_than_20_least_1_t6SS_list_of_genomes.txt")

write.table(subset.data.frame(subset.data.frame(least_1_T6SS, V7 >= 20 ), V2.y <250 & V2.y >150, select=c("file_prefix")), quote=F, row.names = F, file = "selecting_hmmer_more_than_20_length_correct_least_1_T6SS.txt")

#### species to use to import into itol :

headers_for_itol <- data.frame(as.list(col_vector))
species_colours <- data.frame(keyName=names(col_vector), value=col_vector, row.names=NULL)
with_colours <- merge(complete_with_T6SS,species_colours,by.x="Short.species",by.y = "keyName", all.x = T)

for_itol <- subset.data.frame(with_colours,select=c("file_prefix","value","Short.species"))
for_itol$Short.species <- str_replace_all(for_itol$Short.species, " ", "_")



for_itol$Short.species <- str_replace_all(for_itol$Short.species, " ", "_")
write.csv(for_itol,file="colours_species_for_itol.csv",quote = F,row.names = F)
write.csv(headers_for_itol,file="headers_for_itol.csv",quote = F,row.names = F)


### hamburger_stats:
hamburger_data <- read.csv("../data/1_ssp6_hamburger_10k_up_and_down_operon_stats.csv",
         quote = "",
         comment.char = "",
         stringsAsFactors = F,
         header = T)

merged_for_hamburger <- subset.data.frame(merge(hamburger_data,metadata4later, by.x = "strain", by.y = "File"), select = c("strain","operon_name","Short.species","Strain"))
merged_for_hamburger <- merged_for_hamburger[!duplicated(merged_for_hamburger$operon_name),]


#### change names for tip labels and operons:

merged_for_hamburger$number_for_name <- as.data.frame(matrix(unlist(noquote(str_split(merged_for_hamburger$operon_name, "_"))), nrow=103,byrow = T))$V6
merged_for_hamburger$new_name<- with(merged_for_hamburger, paste0(Short.species, sep=" ", Strain, sep=" ", number_for_name))

### ssp6 hmmer score for each operon - data from hamburger step 5 on each 10k fragment:
hmmer_step_5_ssp6_data <- read.csv("../data/5_hamburger_hmmer_scores.csv",header=T, quote = "", comment.char = "",stringsAsFactors = F)
hmmer_ssp6_data_pre <-read.csv("../data/5_hamburger_hmmer_scores.csv", header=T, comment.char = "", stringsAsFactors = F, quote="")
row.names(hmmer_ssp6_data_pre) <- hmmer_ssp6_data_pre$hmmer
hmmer_ssp6_data <- subset.data.frame(hmmer_ssp6_data_pre, select=c("score"))
colnames(hmmer_ssp6_data) <- "Hmmer score"
#colnames(hmmer_ssp6_data_pre) <- c("operon_name","Hmmer_score")

##create csv for outputting table for supplementary:

pre_pre_sup <- merge(hamburger_data, metadata4later, by.x = "strain", by.y="File")
pre_sup <- merge(merge(pre_pre_sup,hmmer_ssp6_data_pre, by.x = "operon_name", by.y = "hmmer"),merged_for_hamburger,by.x="operon_name",by.y="operon_name")
supplementary_pre <-  (subset.data.frame(pre_sup, select=c("Assembly","new_name","Short.species.y","Strain.y","contig","start","stop","score","evalue")))
supplementary <- supplementary_pre[!duplicated(supplementary_pre),]
write.csv(supplementary,quote = F, row.names = F, file = "hmmer_score_supplementary.csv")

## cat colour groups

new_cols <- c(col_vector, cols)

###Hmmer Evalues vs hmmer score:
ggplot(supplementary,aes(y = score,x = evalue, color = Short.species.y)) + geom_count() +
  scale_color_manual(values=c(new_cols)) +
  guides(color=guide_legend(ncol =3,keywidth=0.1,
                            keyheight=0.1,
                            default.unit="inch")) +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=8, face="bold"),
        legend.justification=c(0,0),
        legend.position=c(0.35, 0.12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  ylab("Hmmer score") +
  xlab("Hmmer E value") +
  geom_text(aes(label=ifelse(Strain.y =="DB11",as.character("Db11"),'')),hjust=1.2,vjust=0,colour='black',alpha = NA) +
  geom_hline(yintercept=20, linetype="dashed", color = "red") +
#  geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
  #annotate("rect",xmin = 150, xmax = 250, ymin = 20, ymax = 350, alpha = 0.05, colour = "red", linetype = "dashed" ) +
  ggsave("Hmmer_score_vs_Evalue.png",dpi=300, width=9,height=7, units=c("in"))

ggplot(supplementary,aes(y = score,x = log(evalue), color = Short.species.y)) + geom_count() +
  scale_color_manual(values=c(new_cols)) +
  guides(color=guide_legend(ncol =3,keywidth=0.1,
                            keyheight=0.1,
                            default.unit="inch")) +
  theme(legend.text = element_text(size=7),
        legend.title = element_text(size=8, face="bold"),
        legend.justification=c(0,0),
        legend.position="none",#c(0.35, 0.12),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  ylab("Hmmer score") +
  xlab("Hmmer E value") +
  geom_text(aes(label=ifelse(Strain.y =="DB11",as.character("Db11"),'')),hjust=1.2,vjust=0,colour='black',alpha = NA) +
  geom_hline(yintercept=20, linetype="dashed", color = "red") +
  #  geom_vline(xintercept=0.05, linetype="dashed", color = "red") +
  #annotate("rect",xmin = 150, xmax = 250, ymin = 20, ymax = 350, alpha = 0.05, colour = "red", linetype = "dashed" ) +
  ggsave("Hmmer_score_vs_Evalue_log.png",dpi=300, width=9,height=7, units=c("in"))

  ## species for tree with genes:

ssp6_tree_data <- subset.data.frame(merged_for_hamburger, select = c(Short.species))
row.names(ssp6_tree_data) <- merged_for_hamburger$operon_name

### sp6_tree from hamburger:

hamburger_ssp6_tree <- read.tree(file="../data/5_hamburger_ssp6_hits.fasta.treefile")
tip_labels = data.frame(label=merged_for_hamburger$operon_name, label2 = merged_for_hamburger$new_name)
t1 <- ggtree(hamburger_ssp6_tree)  %<+% tip_labels + geom_rootedge(0.2) + geom_text2(aes(subset = !isTip, label=label),size = 2, nudge_x = -0.04, nudge_y = 0.45) + geom_tiplab(aes(label=label2),size=2.8) #+ geom_tiplab(size=2) #+ xlim(0,10) #+ xlim_tree(16)
t8 <- t1 %>% gheatmap(hmmer_ssp6_data, color = NULL, offset = 1, width = 0.2, colnames_offset_y = -1) +
  scale_fill_gradient(low="#ffffff", high="#2980B9") +
  theme(legend.title=element_text()) +
  labs(fill="Hmmer score") +
  geom_treescale(width = 2, offset = -1.2) +
  ggsave("all_hits_with_hmmer_score.png",dpi=300, width=16.53,height=11.69, units=c("in")) #+ theme(legend.position = "none")

### load in gggenes

genes <- read.csv(file="../data/5_hamburger_same_direction_gggenes_input.csv", header = T, comment.char = "", quote = "")
genes2 <- read.csv(file="../data/5_hamburger_gggenes_input.csv", header = T, comment.char = "", quote = "")

new_data <- genes[,c(2,1,3:7)]
test_data <- genes2[,c(2,1,3:7)]





### plot tree with gggenes
####alter previous tree to get the bootstrap values in the right position
t1_species <- ggtree(hamburger_ssp6_tree)  %<+% tip_labels + geom_rootedge(0.2) + geom_text2(aes(subset = !isTip, label=label),size = 2, nudge_x = -0.08, nudge_y = 0.45) + geom_tiplab(aes(label=label2),size=2.8) #+ geom_tiplab(size=2) #+ xlim(0,10) #+ xlim_tree(16)
t2 <- t1_species %>% gheatmap(ssp6_tree_data, color = NULL, offset = 2, width = 0.2, colnames_offset_y = -1) + scale_fill_manual(values=c(col_vector)) + theme(legend.position = "none")

facet_plot(t2, panel='Operons',
           mapping = aes(xmin = start, xmax = end, fill = gene, x=start, forward = direction, y=y),
           data=new_data, geom=geom_gene_arrow, arrow_body_height = grid::unit(1.5, "mm"), arrowhead_height = grid::unit(2, "mm"), arrowhead_width = grid::unit(2,"mm")) +
  scale_fill_manual(values=c(new_cols)) +
  geom_treescale(width = 2, offset = -1.2) +
  ggsave("filtered_published_ssp6_with_operons.png", dpi=300, width=19,height=11.69, units=c("in"))


ggplot2::ggplot(genes, ggplot2::aes(xmin = start, xmax = end,y = operon, fill = gene)) +
  geom_gene_arrow() +
  scale_fill_manual(values=c(new_cols))


###select small number of isolates to look at
###load in tree of ssp6 from subset of isolates
small_hamburger_ssp6_tree <- read.tree(file="../data/5_hamburger_selected_ssp6_hits.fasta.treefile")
t3 <- ggtree(small_hamburger_ssp6_tree)  %<+% tip_labels + geom_text2(aes(subset = !isTip, label=label),size = 3, nudge_x = -0.18, nudge_y = 0.15) + geom_tiplab(aes(label=label2),size=4,hjust=-0.02) + geom_rootedge(0.2)#+ geom_tiplab(size=2) #+ xlim(0,10) #+ xlim_tree(16)
t4 <- t3 %>% gheatmap(ssp6_tree_data, color = NULL, offset = 5, width = 0.5, colnames_offset_y = -1) + scale_fill_manual(values=c(col_vector)) + theme(legend.position = "none")
t4
t9 <- t3 %>% gheatmap(hmmer_ssp6_data, color = NULL, offset = 2, width = 0.2, colnames_offset_y = -1) +
  scale_fill_gradient(low="#ffffff", high="#2980B9") +
  theme(legend.title=element_text()) +
  labs(fill="Hmmer score") +
  geom_treescale(width = 2, offset = -0.4) +
  ggsave("selected_hits_with_hmmer_score.png",dpi=300, width=16.53,height=11.69, units=c("in")) #+ theme(legend.position = "none")

### draw smaller set with genes alongside

f1 <- facet_plot(t4, panel='Operons',
           mapping = aes(xmin = start, xmax = end, fill = gene, x=start, forward = direction, y=y),
           data=new_data, geom=geom_gene_arrow, arrow_body_height = grid::unit(3, "mm"), arrowhead_height = grid::unit(4, "mm"), arrowhead_width = grid::unit(2,"mm")) +
  scale_fill_manual(values=c(new_cols)) +
  geom_treescale(width = 2, offset = -0.4) +
  ggsave("smaller_set_ssp6with_operons.png",dpi=300, width=16.53,height=9, units=c("in"))


### legends for the species and genes separately:

#data frame for the operons in the tree - for the genes
operon_names <- as.data.frame(small_hamburger_ssp6_tree$tip.label)
colnames(operon_names) <- "name"
for_legend_plotting <- merge(operon_names,genes, by.x = "name", by.y = "operon")
for_legend_plotting$gene <- str_replace_all(for_legend_plotting$gene, "_", " ")

#plot species figure legend
l2 <- t3 %>% gheatmap(ssp6_tree_data, color = NULL, offset = 5, width = 0.5, colnames_offset_y = -1) +
  scale_fill_manual(values=c(col_vector),name = "Species") +
  guides(fill=guide_legend(override.aes=list(colour=NA),ncol=1,keywidth=0.3,
                           keyheight=0.3,
                           title="Species",
                           default.unit="inch")) +
  theme(legend.position = "right",
        legend.text = element_text(face = "italic"),
        legend.title = element_text(),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'))

l2



####species legend for bigger tree

l3 <- t2 +
  scale_fill_manual(values=c(col_vector),name = "Species") +
  guides(fill=guide_legend(override.aes=list(colour=NA),ncol=2,keywidth=0.3,
                           keyheight=0.3,
                           title="Species",
                           default.unit="inch")) +
  theme(legend.position = "right",
        legend.title = element_text(),
        legend.text = element_text(face = "italic"),
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'))
l3


#plot genes figure legend
l1 <- ggplot(for_legend_plotting, aes(xmin = start, xmax = end, y = name, fill = gene, forward = direction)) +
  geom_gene_arrow() +
  scale_fill_manual(values=c(cols)) +
  guides(fill=guide_legend(override.aes=list(colour=NA),ncol=2,keywidth=0.3,
                           keyheight=0.3,
                           title = "Gene",
                           default.unit="inch")) +
  theme(legend.position = "right",
        legend.title = element_text("Gene")) +
  theme_genes()



legend_1 <- cowplot::get_legend(l1)
legend_2 <- cowplot::get_legend(l2)
legend_3 <- cowplot::get_legend(l3)
legend_4 <- cowplot::get_legend(t8)

ggdraw() +
   draw_plot(legend_1, x = -0.1, y = 0.45, width = .5, height = .5) +
   draw_plot(legend_2, x = .15, y = 0.45, width = .5, height = .5) +
   draw_plot(legend_3, x = .5, y = 0.25, width = .5, height = .5) +
   draw_plot(legend_4, x = -0.1, y = 0, width = .5, height = .5) +
  ggsave("legends.png",dpi=300, width=10,height=8, units=c("in"))
