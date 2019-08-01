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
    ylab("% of sequences") 


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
