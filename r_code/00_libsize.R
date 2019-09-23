synthetic %>% 
    group_by(short, study, protocol) %>% 
    summarise(n=n(),size=sum(value)) %>% 
    mutate(ytext = size + size*0.1,
           text = paste0(round(n/1000, digits = 0), "K")) %>% 
    ggplot(aes(paste(study, short), size, fill = protocol)) +
    geom_bar(stat = "identity") +
    ylab("reads") +
    geom_text(aes(paste(study, short), ytext, label = text), size = 3) +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggsave("results/01_library/all_samples_barplot.pdf", width=12, height = 9)


#    ggsave("results/00_library_size/00_libsize_diversity_bar.pdf", width = 9)
