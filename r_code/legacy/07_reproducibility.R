library(tidyverse)
library(ggplot2)
library(ggthemes)
library(pheatmap)
library(stringr)
theme_set(
    theme_light(base_size = 11L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")

counts = .  %>% 
    filter(iso  %in% c("shift3p", "shift5p",  "snp", "reference"),
           pct>5) %>% 
    group_by(protocol) %>% 
    mutate(max_n = length(unique(replicate))) %>% 
    group_by(mi_rna, iso_detail, iso, protocol, max_n) %>%
    summarise(n_rep=length(unique(replicate)),
              counts=sum(normalized)) %>% 
    group_by(protocol, iso, max_n) %>%
    summarise(detect90=quantile(n_rep, .10),
              detect75=quantile(n_rep, .25),
              detect50=quantile(n_rep, .50),
              detect25=quantile(n_rep, .75),
              detect99=quantile(n_rep, .01)) %>% 
    ungroup()

cols = c("#dc3545",#red
         "#28a745",#green
         "#343a40",#black
         "#ffc107",#y
         "orange", 
         "#007bff") #blue
bind_rows(
    equimolar_razer3 %>% 
        filter(protocol!="clean") %>%
        mutate(protocol = paste0("tewari_synthetic_", protocol),
               replicate = lab) %>% 
        counts,
    dsrg %>% 
        mutate(protocol = paste0("dsrg_synthetic_", protocol),
               replicate = paste0(
                   str_match(short, "_([0-9]+)$")[,2],
               "_",
               lab)) %>%
        counts,
    narrykim %>% 
        mutate(protocol = paste0("narrykim_synthetic", protocol),
               replicate = str_match(short, "_([0-9]+)$")[,2]) %>% 
        counts,
    vandijk  %>% 
        mutate(protocol = paste0("vandijk_synthetic_", protocol),
               replicate = str_match(short, "_rep([0-9]+)_")[,2]) %>% 
        counts,
    plasma %>%
        filter(protocol!="clean") %>%
        mutate(protocol = paste0("tewari_human_", protocol, "_", lab),
               replicate = str_match(short, "_([0-9]+)$")[,2]) %>%
        counts
) %>% mutate(protocol=factor(protocol, levels = unique(protocol))) -> full
axiscol <- ifelse(grepl("syn", levels(full$protocol)), "black", "blue")
ggplot(full) +
    geom_bar(aes(protocol, max_n, fill = "0%"), stat = "identity") +
    geom_bar(aes(protocol, detect25, fill = "<25%"), stat = "identity") +
    geom_bar(aes(protocol, detect50, fill = "25%-50%"), stat = "identity") +
    geom_bar(aes(protocol, detect75, fill = "50%-75%"),stat = "identity") +
    geom_bar(aes(protocol, detect90, fill = "75%-90%"), stat = "identity") +
    geom_bar(aes(protocol, detect99, fill = ">99%"),stat = "identity") +
    # geom_point(aes(protocol, max_n), shape = 8, color = "black") +
    facet_wrap(~iso, nrow = 4) +
    ylab("# replicates") +
    theme_few() +
    scale_y_continuous(breaks = c(1,3,5,7, 9)) +
    scale_fill_manual(values=cols, 
                      name="sequences",
                      breaks=c("0%",
                               "<25%", "25%-50%",
                               "50%-75%", "75%-90%",
                               ">99%"))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color=axiscol),
          panel.grid.minor = element_line(color="lightgrey")) +
    ggsave("results/07_reproducibility/07_across_study_g5.pdf", height = 9, width = 12)
    


# equimolar_razer3 %>% filter(!grepl("x4n_[abcd]", sample),
#                             pct > 1) %>% 
#     distinct(read, sample, pct) %>% 
#     #mutate(pct=ifelse(pct>0,1,pct)) %>% 
#     spread(sample, pct, fill=-50) %>%
#     column_to_rownames("read") %>% 
#     as.data.frame() %>% 
#     pheatmap::pheatmap(., show_rownames = F, 
#                        color = c("white", RColorBrewer::brewer.pal(7,"YlOrRd")),
#                        breaks = c(-100,0,1,10,20,30,40,50,101))
# 

