library(tidyverse)

label_iso = function(read, mix, snp, i3p, i5p, ia3p, ia5p){
    if (mix == -1)
        return("reference")
    if (mix > 1)
        return("mix")
    if (i3p!=0)
        return("iso_3p")
    if (i5p!=0)
        return("iso_5p")
    if (ia3p!=0)
        return("iso_add3p")
    if (ia5p!=0)
        return("iso_add5p")
    if (snp!=0)
        return("iso_snp")
    stop(c(read, mix, snp, i3p, i5p, ia3p, ia5p))
}

annotate = . %>% 
#    gff %>%
    group_by(read) %>% 
    mutate(hits=n()) %>%
    # only uniquely mapped
    filter(hits==1) %>% 
    mutate( # consider additiongs always non-template additions
           iso_add3p = ifelse(iso_3p > 0, iso_3p, iso_add3p),
           iso_3p = ifelse(iso_3p > 0, 0, iso_3p),
           iso_add5p = ifelse(iso_5p < 0, abs(iso_5p), 0),
           iso_5p = ifelse(iso_5p < 0, 0, iso_5p),
           iso_mix = str_count(variant,  "iso"),
           iso_mix = ifelse(is.na(iso_mix), -1, iso_mix)) %>%
    ungroup() %>% 
    # this will counts how many changes in total
    mutate(iso_snp= sapply(str_match_all(iso_snp_nt,"[0-9][ACTG-]"), function(n){nrow(n)}),
           iso_n =  iso_snp + abs(iso_5p) + abs(iso_3p) + abs(iso_add3p)  + abs(iso_add5p)) %>% 
    # and keep only sequences with up to 4 changes
    filter(iso_n < 5) %>% 
    select(-uid, -variant, -hits ) %>%
    unite("iso_text", iso_add3p_nt, iso_3p_nt, iso_5p_nt, iso_snp_nt, sep =  " ") %>%
    ungroup %>% 
    rowwise %>% 
    mutate(iso = label_iso(read, iso_mix, iso_snp, iso_3p, iso_5p, iso_add3p, iso_add5p)) %>% 
    gather(sample, value, -mi_rna,
           -iso, -iso_mix,
           -iso_3p,-iso_5p,
           -iso_add3p, -iso_add5p, -iso_snp, -iso_n, -read, -iso_text) %>% 
    filter(value>0)

analysis = . %>% 
    group_by(sample, mi_rna) %>% 
    filter(value > 0) %>% 
    arrange(sample, mi_rna, desc(value)) %>%
    mutate(rank = 1:n(),
           total = sum(value),
           pct = value / total * 100) %>% 
    ungroup() %>% 
    group_by(read) %>% 
    mutate(reproducible_protocol = paste(unique(protocol), collapse = " "),
           reproducible_lab = length(unique(lab))) %>%
    group_by(mi_rna, sample) %>%
    mutate(ref_is_1 = length(mi_rna[rank == 1 & iso == "reference"])) %>% 
    ungroup %>% 
    mutate(pct_cat = cut(pct,
                         breaks = c(-1, 0.1, 1, 5, 10, 20, 50, 101),
                         labels = c("<0.1", "0.1-1", "1-5", "5-10", "10-20", "20-50", ">50")),
           mir_coverage = cut(total,
                              breaks = c(-1, 1, 10, 100, 1000, 1e30),
                              labels = c("<1", "1-10", "10-100", "100-1000", ">1000"))) %>% 
    ungroup()

library(edgeR)
norm = function(fn, cache=NULL){
    counts = fn[,13:ncol(fn)] %>% 
        as.matrix()
    counts = counts[,colSums(counts) > 10]
    # if (!is.null(cache)){
    #     if (file.exists(cache)){
    #         load(cache)
    #     }else{
    #         t = apply(counts, 2, function(c){
    #             thresholdSeq(c[c>0])[["threshold"]]
    #         })
    #         threshold = data.frame(sample = names(t), 
    #                                threshold = t, stringsAsFactors = F)
    #     }
    #     
    # }
    # save(threshold, file = cache)
    dds = DGEList(counts)
    dds = calcNormFactors(dds)
    counts = cpm(dds, normalized.lib.sizes = FALSE)
    normalized = cbind(fn[,"read"],
                       counts) %>% 
        gather(sample, normalized,  -read) 
    # %>% 
    #     left_join(threshold, by = "sample")
    
}

# tewari
gff = read_tsv("input/synthetic_tewari_mirtop.tsv") %>% 
    janitor::clean_names()
normalized = norm(gff)
# get all sequence that are hsa annotated or not in mirx(remove potential crossmapping)
tewari = gff %>% 
    annotate %>%
    left_join(normalized, by = c("sample", "read")) %>% 
    mutate(lab=stringr::str_extract(sample,"lab[0-9]"),
           protocol=stringr::str_remove_all(sample, "_.*$"),
           index = as.numeric(as.factor(sample))) %>% 
    unite("short", c("protocol", "lab", "index"), remove = FALSE) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis

tewari = tewari %>% 
    mutate(short = ifelse(sample=="x4n_nex_tflex_lab8_synth_eq",
                          "nextf_lab8_23",
                          short))

#carrie
gff = read_tsv("input/synthetic_carrie_mirtop.tsv") %>% 
    janitor::clean_names()
normalized = norm(gff)
carrie = gff %>% 
    annotate %>%
    left_join(normalized, by = c("sample", "read")) %>% 
    mutate(lab=stringr::str_extract(sample,"lab[0-9]"),
           protocol=stringr::str_remove_all(sample, "_.*$"),
           index = as.numeric(as.factor(sample))) %>% 
    unite("short", c("protocol", "lab", "index"), remove = FALSE) %>% 
    select(-index) %>% 
    distinct() %>% 
    analysis

# dsrg
gff = read_tsv("input/synthetic_dsrg_mirtop.tsv") %>% 
    janitor::clean_names() %>% 
    set_names(gsub("_mur_d", "_murd", names(.))) %>%
    .[,c(1:12, grep("_mur_", names(.)))] %>% 
    .[rowSums(as.matrix(.[,13:ncol(.)]))>0,]
normalized = norm(gff)
dsrg = gff %>% 
    annotate %>%
    left_join(normalized, by = c("sample", "read")) %>% 
    separate(sample, remove = F, into = c("protocol", "source", 
                                          "lab", "index", "snumber")) %>% 
    unite("short", c("protocol", "lab", "index"),
          remove = FALSE) %>% 
    select(-index,-source) %>% 
    distinct() %>% 
    analysis

synthetic = bind_rows(
    tewari %>% mutate(study="tewari"),
    carrie %>% mutate(study="carrie"),
    dsrg %>% mutate(study="dsrg")
)


saveRDS(synthetic, "data/synthetic_2019_mirgff1.2.rds")