library(tidyverse)
library(ggplot2)
library(ggthemes)
theme_update(
    legend.justification = "center",
    legend.position = "bottom")
synthetic = readRDS("data/synthetic_2019_mirgff1.2.rds")

prepare =  . %>% filter(ref_is_1 == 1) %>%
    dplyr::count(short, pct_cat, iso, protocol, study, libsize) %>% 
    group_by(short, protocol, study, libsize) %>% 
    mutate(pct_total = n/sum(n)*100)

synthetic %>% 
    group_by(short, study, protocol) %>% 
    mutate(libsize=log10(sum(value)),
           iso=relevel(as.factor(iso), "reference")) %>% 
    filter(total>500) %>% 
    prepare() %>%
    ggplot(aes(x = paste(study, short),
               fill = pct_cat, y = pct_total,
               alpha=libsize>5)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso, ncol = 1, scales = "free_y") +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    scale_alpha_discrete(range = c(0.3,1)) +
    ggtitle("isomiRs importance by type") +
    ylab("% of sequences") +
    ggsave("results/02_importance/all_samples_barplot.pdf", width=12, height = 9)


synthetic %>% 
    group_by(short, study, protocol) %>% 
    mutate(libsize=log10(sum(value)),
           libsizen=n(),
           iso=relevel(as.factor(iso), "reference")) %>% 
    filter(total>500) %>% 
    filter(ref_is_1 == 1) %>%
    group_by(short, protocol, study, libsize, iso, pct_cat) %>% 
    summarise(pct_total=sum(value)/unique(10^libsize)*100,
              pctn_total=n()/unique(libsizen)*100) %>% 
    ungroup() %>% 
    ggplot(aes(x = paste(study, short),
               y = pct_total,
               alpha=pctn_total,
               fill=pct_cat)) +
    geom_bar(stat="identity") +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    facet_wrap(~iso, ncol = 1, scales = "free_y") +
    scale_alpha_continuous(range = c(0.5,1)) +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggsave("results/02_importance/all_samples_barplot_abundance.pdf", width=12, height = 9)


synthetic %>% 
    group_by(short, study, protocol) %>%
    mutate(libsize=log10(sum(value)),
           libsizen=n(),
           libsizem=length(unique(mi_rna)),
           iso=relevel(as.factor(iso), "reference")) %>% 
    filter(total>500) %>% 
    filter(ref_is_1 == 1) %>%
    group_by(short, protocol, study, libsize, mi_rna, iso) %>% 
    arrange(short, protocol, study, libsize, mi_rna, iso, desc(pct)) %>% 
    top_n(1, wt = pct) %>% 
    group_by(short, protocol, study, libsize, iso, pct_cat) %>% 
    summarise(pct_total=sum(value)/unique(10^libsize)*100,
              mirn_total=length(unique(mi_rna))) %>% 
    ungroup() %>% 
    mutate(protocol=ifelse(protocol=="tru","ill",protocol),
           protocol=ifelse(protocol=="ilmn","ill",protocol),
           protocol=ifelse(grepl("nextf", short),"nex",protocol)) %>% 
    ggplot(aes(x = paste(study, short),
               y = mirn_total,
               alpha=pct_total,
               fill=pct_cat)) +
    geom_bar(stat="identity") +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    facet_wrap(~iso, ncol = 1) +
    scale_alpha_continuous(range = c(0.5,1)) +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggsave("results/02_importance/all_samples_barplot_num_mirna.pdf", width=12, height = 9)    

synthetic %>% 
    group_by(short, study, protocol) %>% 
    mutate(libsize=log10(sum(value)),
           libsizen=n(),
           libsizem=length(unique(mi_rna)),
           iso=relevel(as.factor(iso), "reference")) %>% 
    filter(total>500, pct>5) %>% 
    filter(ref_is_1 == 1) %>%
    group_by(short, protocol, study, libsize, mi_rna, iso) %>% 
    arrange(short, protocol, study, libsize, mi_rna, iso, desc(pct)) %>% 
    top_n(1,wt = pct) %>% 
    group_by(short, protocol, study, libsize, iso, pct_cat) %>% 
    summarise(pct_total=sum(value)/unique(10^libsize)*100,
              mirn_total=length(unique(mi_rna))) %>% 
    ungroup() %>% 
    mutate(protocol=ifelse(protocol=="tru","ill",protocol),
           protocol=ifelse(protocol=="ilmn","ill",protocol),
           protocol=ifelse(grepl("nextf", short),"nex",protocol)) %>% 
    ggplot(aes(x = paste(study, short),
               y = mirn_total,
               alpha=pct_total,
               fill=pct_cat)) +
    geom_bar(stat="identity") +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")[4:7]) +
    facet_wrap(~iso, ncol = 1) +
    scale_alpha_continuous(range = c(0.5,1)) +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggsave("results/02_importance/all_samples_barplot_num_mirna_5pct.pdf", width=12, height = 9)

prepare =  . %>% filter(ref_is_1 == 1) %>%
    dplyr::count(short, pct_cat, iso, protocol, study, libsize, iso_n) %>% 
    group_by(short, protocol, study, libsize) %>% 
    mutate(pct_total = n/sum(n)*100)
    
synthetic %>% 
    group_by(short, study, protocol) %>% 
    mutate(libsize=log10(sum(value))) %>% 
    filter(total>500, iso=="iso_3p") %>% 
    prepare() %>%
    ggplot(aes(x = paste(study, short),
               fill = pct_cat, y = n,
               alpha=libsize>5)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso_n, ncol = 1, scales = "free_y") +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    scale_alpha_discrete(range = c(0.3,1)) +
    ggtitle("isomiRs importance by size") +
    ylab("% of sequences") 

synthetic %>% 
    group_by(short, study, protocol) %>% 
    mutate(libsize=log10(sum(value))) %>% 
    filter(total>500, iso=="iso_5p") %>% 
    prepare() %>%
    ggplot(aes(x = paste(study, short),
               fill = pct_cat, y = n,
               alpha=libsize>5)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso_n, ncol = 1, scales = "free_y") +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    scale_alpha_discrete(range = c(0.3,1)) +
    ggtitle("isomiRs importance by size") +
    ylab("% of sequences") 

synthetic %>% 
    group_by(short, study, protocol) %>% 
    mutate(libsize=log10(sum(value))) %>% 
    filter(total>500, iso=="iso_snp") %>% 
    prepare() %>%
    ggplot(aes(x = paste(study, short),
               fill = pct_cat, y = n,
               alpha=libsize>5)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso_n, ncol = 1, scales = "free_y") +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    scale_alpha_discrete(range = c(0.3,1)) +
    ggtitle("isomiRs importance by size") +
    ylab("% of sequences") 
#    ggsave("results/04_other_data/04_other_data_pct_g5.pdf", height = 7)

synthetic %>% 
    group_by(short, study, protocol) %>% 
    mutate(libsize=log10(sum(value))) %>% 
    filter(total>500, iso=="iso_snp") %>% 
    prepare() %>%
    filter(!(pct_cat %in% c("<0.1", "0.1-1", "1-5"))) %>% 
    ggplot(aes(x = paste(study, short),
               fill = pct_cat, y = n)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso_n, ncol = 1, scales = "free_y") +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    scale_alpha_discrete(range = c(0.3,1)) +
    ggtitle("isomiRs importance by size") +
    ylab("% of sequences") 
#    ggsave("results/04_other_data/04_other_data_pct_g5.pdf", height = 7)

synthetic %>% 
    group_by(short, study, protocol) %>% 
    mutate(libsize=log10(sum(value))) %>% 
    filter(total>500, iso=="iso_add3p") %>% 
    prepare() %>%
    ggplot(aes(x = paste(study, short),
               fill = pct_cat, y = n,
               alpha=libsize>5)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso_n, ncol = 1, scales = "free_y") +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    scale_alpha_discrete(range = c(0.3,1)) +
    ggtitle("isomiRs importance by size") +
    ylab("% of sequences")

synthetic %>% 
    group_by(short, study, protocol) %>% 
    mutate(libsize=log10(sum(value))) %>% 
    filter(total>500, iso=="iso_add5p") %>% 
    prepare() %>%
    ggplot(aes(x = paste(study, short),
               fill = pct_cat, y = n,
               alpha=libsize>5)) +
    geom_bar(stat = "identity") +
    facet_wrap(~iso_n, ncol = 1, scales = "free_y") +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")) +
    scale_alpha_discrete(range = c(0.3,1)) +
    ggtitle("isomiRs importance by size") +
    ylab("% of sequences")
