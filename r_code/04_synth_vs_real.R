library(tidyverse)
library(ggplot2)
library(ggthemes)
theme_update(
  legend.justification = "center",
  legend.position = "bottom")

synthetic = readRDS("data/synthetic_2019_srr_mirgff1.2.rds")
real = readRDS("data/real_mirgff1.2.rds")

data = bind_rows(real %>% mutate(type="real"),
      synthetic %>% filter(study=="cwrigth" | study=="nkim" | study=="fratta") %>% mutate(type="synth")
      )

data_sm = data %>% 
    filter(total>100, pct>5) %>%
    group_by(sample_n, study, protocol) %>% 
    mutate(libsize=log10(sum(value)),
           libsizen=n(),
           libsizem=length(unique(mi_rna)),
           iso=relevel(as.factor(iso), "reference")) %>% 
    group_by(sample_n, protocol, type, study, libsize, mi_rna, iso) %>% 
    arrange(sample_n, protocol, type, study, libsize, mi_rna, iso, desc(pct)) %>% 
    top_n(1,wt = pct) %>% 
    group_by(sample_n, protocol, type, study, libsize, iso, pct_cat) %>% 
    summarise(pct_total=sum(value)/unique(10^libsize)*100,
              mirn_total=length(unique(mi_rna))/unique(libsizem)*100) %>% 
    ungroup() %>% 
    mutate(x=paste(protocol, type, study, sample_n))

axis_face = ifelse(grepl("real", sort(unique(data_sm$x))), "bold.italic", "plain")

cols = RColorBrewer::brewer.pal(6, "Dark2")
names(cols) = sort(unique(data_sm$protocol))
axis_col = sapply(sort(unique(data_sm$x)), function(x){
    cols[strsplit(x, " ")[[1]][1]]
})

ggplot(data_sm, aes(x = x,
               y = mirn_total,
               alpha=pct_total,
               fill=pct_cat)) +
    geom_bar(stat="identity") +
    scale_fill_manual("IMPORTANCE", values = RColorBrewer::brewer.pal(7, "Dark2")[4:7]) +
    # facet_wrap(~iso+protocol ,ncol=7, scales = "free_x") +
    facet_wrap(~iso, ncol=1) + 
    scale_alpha_continuous(range = c(0.5,1)) +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,face = axis_face,  color = axis_col)) +
    ggsave("results/04_synth_vs_real/barplot_total100_pct5.pdf", width=12, height = 10)
