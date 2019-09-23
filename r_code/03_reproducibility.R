library(tidyverse)
library(ggplot2)
library(ggthemes)
theme_update(
    legend.justification = "center",
    legend.position = "bottom")
synthetic = readRDS("data/synthetic_2019_srr_mirgff1.2.rds")


synthetic %>% 
    group_by(sample_n) %>% 
    summarize(minn=min(normalized)) %>% 
    as.data.frame() %>% 
    pull(minn) %>% 
    max ->cuto


synthetic %>% 
    filter(normalized>=cuto) %>% 
    distinct(read) -> seen
synthetic %>% 
    filter(normalized>=cuto, pct>=5) %>% 
    distinct(read) -> seen5

cols = c("#343a40",#black
         "red4",#red
         "#dc3545",#red
         "#ffc107",#y
         "orange", 
         "#007bff", #blue
         "#28a745")#green

# synthetic %>% 
#     filter(read %in% seen$read) %>% 
#     group_by(protocol) %>% 
#     mutate(maxsamples=length(unique(sample_n))) %>% 
#     filter(maxsamples>3) %>% 
#     group_by(read, iso, protocol, pct_cat, maxsamples) %>% 
#     summarise(n_rep=length(unique(sample_n))) %>% 
#     group_by(protocol, iso, maxsamples, pct_cat) %>%
#     summarise(detect99=quantile(n_rep, .01),
#               detect90=quantile(n_rep, .10),
#               detect75=quantile(n_rep, .25),
#               detect50=quantile(n_rep, .50),
#               detect25=quantile(n_rep, .75),
#               detect10=quantile(n_rep, .90),
#               detect01=quantile(n_rep, .99)
#               ) %>% 
#     ungroup() -> parsed
# 
# ggplot(parsed) +
# geom_bar(aes(pct_cat, maxsamples, fill = "1%"), stat = "identity") +
#     geom_bar(aes(pct_cat, detect10, fill = "10%"), stat = "identity") +
#     geom_bar(aes(pct_cat, detect25, fill = "25%"), stat = "identity") +
#     geom_bar(aes(pct_cat, detect50, fill = "50%"), stat = "identity") +
#     geom_bar(aes(pct_cat, detect75, fill = "75%"),stat = "identity") +
#     geom_bar(aes(pct_cat, detect90, fill = "90%"), stat = "identity") +
#     geom_bar(aes(pct_cat, detect99, fill = "99%"),stat = "identity") +
#     # geom_point(aes(protocol, max_n), shape = 8, color = "black") +
#     facet_grid(protocol~iso, scales="free_y") +
#     ylab("# replicates") +
#     theme_clean() +
#     scale_y_continuous(breaks = c(1,3,5,7, 9,11,13)) +
#     scale_fill_manual(values=cols, 
#                       name="sequences")+
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           panel.grid.minor = element_line(color="lightgrey")) 

fct_nrep_pct=c("<25", "25-50", "50-75", ">75")

analyze = . %>% 
    group_by(protocol) %>% 
    mutate(maxsamples=length(unique(sample_n))) %>% 
    mutate(nmirs=length(unique(mi_rna[iso=="reference"]))) %>% 
    filter(maxsamples>3) %>% 
    group_by(read, iso, protocol, pct_cat, maxsamples, nmirs) %>% 
    summarise(n_rep=length(unique(sample_n))) %>% 
    mutate(n_rep_cat=cut(n_rep/maxsamples,
                         c(-1,0.25,0.5,0.75,1.1),
                         fct_nrep_pct),
           n_rep_cat=factor(n_rep_cat, levels=fct_nrep_pct)) %>%
    group_by(iso, protocol, pct_cat, n_rep_cat, maxsamples, nmirs) %>%
    summarise(nreads=n(), nreads_mi = n()/unique(nmirs))
    
synthetic %>% 
    filter(read %in% seen$read) %>% 
    analyze() %>% 
    ggplot(aes(pct_cat, y=nreads_mi, fill=n_rep_cat)) +
    geom_bar(stat = "identity") +
    facet_grid(iso~protocol, scale="free_y") +
    theme_clean() +
    ylab("# of miRNAs") +
    xlab("Importance") +
    scale_fill_manual("% of replicates",values=cols[c(2,4,6,7)]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggsave("results/03_reproducibility/all_fixed_coveraged.pdf", width=12, height=10)


synthetic %>% 
    filter(read %in% seen5$read) %>% 
    analyze() %>% 
    ggplot(aes(pct_cat, y=nreads_mi, fill=n_rep_cat)) +
    geom_bar(stat = "identity") +
    facet_grid(iso~protocol, scale="free_y") +
    theme_clean() +
    scale_fill_manual(values=cols[c(2,4,6,7)]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
    ggsave("results/03_reproducibility/all_fixed_coveraged_5pct.pdf", width=12, height=10)


# illumina only across diff study to show similar noise results
synthetic %>% 
    filter(study %in% c("nkim", "tewari", "dsrg"),
           protocol  %in% c("aq", "ill")) %>% 
    mutate(protocol = paste(protocol, study)) %>% 
    group_by(protocol) %>% 
    mutate(maxsamples=length(unique(sample_n))) %>% 
    mutate(nmirs=length(unique(mi_rna[iso=="reference"]))) %>% 
    filter(maxsamples>2) %>% 
    group_by(read, iso, protocol, pct_cat, maxsamples, nmirs) %>% 
    summarise(n_rep=length(unique(sample_n))) %>% 
    mutate(n_rep_cat=cut(n_rep/maxsamples, c(-1,0.25,0.5,0.75,1.1))) %>%
    group_by(iso, protocol, pct_cat, n_rep_cat, maxsamples, nmirs) %>%
    summarise(nreads=n(), nreads_mi = n()/unique(nmirs)) %>% 
    ggplot(aes(pct_cat, y=nreads_mi, fill=n_rep_cat)) +
    geom_bar(stat = "identity") +
    facet_grid(iso~protocol, scale="free_y") +
    theme_clean() +
    scale_fill_manual(values=cols[c(2,4,6,7)]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
