source("source_me.R")

# line 1 = Sl
# line 2 = St
# line 3 = Nb
# line 4 = Sl St Nb
# line 5 = mix

# plot depth

depths <- lapply(setNames(nm = list.files(path = "depths/", pattern = ".*mapped_3.*depth.*")),
                 function(x) {
                   y <- read.table(paste0("depths/", x), col.names = c("ref", "position", "depth"))
                   y$line.rep.time <- x %>% gsub(".*be_|.*po_", "", .) %>% gsub(".*ori", "0", .) %>% gsub("_mapped.*", "", .) %>%
                     gsub("^1", "Sl.L", .) %>%
                     gsub("^2", "St.L", .) %>%
                     gsub("^3", "Nb.L", .) %>%
                     gsub("^4", "CTF.L", .) %>%
                     gsub("^5", "MIX.L", .) %>%
                     
                     gsub("L(\\d)", "L\\1 [", .) %>%
                     gsub("$", "th]", .) %>%
                     gsub("1th]", "1st]", .) %>%
                     gsub("2th]", "2nd]", .) %>%
                     gsub("3th]", "3rd]", .)

                   y$ref <- y$ref %>% gsub("_.*", "", .) %>%
                     gsub("PVYbe", "PVYNb", .) %>%
                     gsub("PVYpo", "PVYSt", .)
                     
                   y$vir.line.rep.time <- paste(y$ref, y$line.rep.time, sep = ".") %>%
                     gsub("\\.0th]", "", .)
                   y
                 })

depths.gg <- depths %>% bind_rows() %>%
  ggplot(aes(position, depth, fill = "grey30")) + geom_line() + facet_wrap(~vir.line.rep.time, scales = "free", ncol = 4) +
  theme(axis.text = element_text(size = 6),
        strip.text = element_text(size = 8))

ggsave("cov.pdf", depths.gg, width = 6.85, height = 9.21)

# cistrons on NC_001616

cistron <- data.frame(x1 = c(185, 1037, 2405, 3500, 3656, 5558, 5714, 6278, 7010, 8573, 2920),
                      x2 = c(1036, 2404, 3499, 3655, 5557, 5713, 6277, 7009, 8572, 9373, 3144),
                      y1= c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                      y2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                      cistron = c("P1", "HC-Pro", "P3", "6K1", "Cl", "6K2", "Nla-VPg", "Nla-Pro", "Nlb", "CP", "PIPO"),
                      size = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .75))

cistron <- dplyr::mutate(cistron, x1 := x1 - 12,
                         x2 := x2 - 12,
                         cistron := factor(cistron,
                                           levels = c("P1", "HC-Pro", "P3", "6K1", "Cl", "6K2", "Nla-VPg", "Nla-Pro", "Nlb", "CP", "PIPO")))

# levels for plots

sample_levels <- c("PVYNb", "PVYSt",
                   "Sl.L1 [1st]", "Sl.L2 [1st]",
                   "St.L1 [4th]", "St.L1 [6th]", "St.L2 [4th]", "St.L2 [10th]",
                   "Nb.L1 [4th]", "Nb.L1 [10th]",
                   "Nb.L2 [4th]", "Nb.L2 [10th]",
                   "CTF.L1 [4th]", "CTF.L2 [1st]", "CTF.L2 [4th]",
                   "MIX.L1 [2nd]", "MIX.L1 [3rd]", "MIX.L2 [4th]", "MIX.L2 [9th]")
  

# filter matrices based on coverage

mean.depth <- depths %>% lapply(function(x) mean(x$depth)) %>% unlist()

# match file names

names(mean.depth) <- names(mean.depth) %>% gsub("sorted.*", "recal_lofreq.vcf", .)

be.samples.discarded <- lofreq.be.vcfs[!names(lofreq.be.vcfs) %in% names(mean.depth[mean.depth > 50])] %>% names()
po.samples.discarded <- lofreq.po.vcfs[!names(lofreq.po.vcfs) %in% names(mean.depth[mean.depth > 50])] %>% names()

lofreq.be.vcfs <- lofreq.be.vcfs[names(lofreq.be.vcfs) %in% names(mean.depth[mean.depth > 50])]
lofreq.po.vcfs <- lofreq.po.vcfs[names(lofreq.po.vcfs) %in% names(mean.depth[mean.depth > 50])]

# benthamiana

be.freqs <- mapply(function(y, x) {
  line.rep.time <- y %>% gsub(".*be_", "", .) %>% gsub(".*ori", "PVYNb", .) %>% gsub("_mapped.*", "", .)
  mut <- paste0(x$REF, x$POS, x$ALT)
  freq <- x$INFO %>% gsub(".*AF=", "", .) %>% gsub(";.*", "", .) %>% as.numeric()
  df <- data.frame(mut = mut, freq = freq, line.rep.time = line.rep.time)
  df
},
names(lofreq.be.vcfs),
lofreq.be.vcfs, SIMPLIFY = FALSE) %>% bind_rows()

be.mtx <- lapply(list(`Sl.L1 [4th]` = c("PVYNb", "111"),
                      `Sl.L1 [4th]` = c("PVYNb", "121"),
                      `St.L1 [4th]` = c("PVYNb", "214"),
                      #line2pop1time5 = c("0", "215"),
                      `St.L2 [4th]` = c("PVYNb", "224"),
                      #line2pop2time5 = c("0", "225"),
                      `Nb.L2 [10th]` = c("PVYNb", "3110"),
                      `Nb.L2 [10th]` = c("PVYNb", "3210"),
                      #line4pop2time8 = c("0", "428"),
                      `MIX.L2 [9th]` = c("PVYNb", "513"),
                      `MIX.L2 [9th]` = c("PVYNb", "529"),
                      `Nb.L2 [4th]` = c("PVYNb", "314"),
                      `Nb.L2 [4th]` = c("PVYNb", "324"),
                      `CTF.L1 [4th]` = c("PVYNb", "414"),
                      `CTF.L2 [4th]` = c("PVYNb", "424"),
                      `MIX.L2 [4th]` = c("PVYNb", "524")),
                 function(x) {
                   
                   df <- be.freqs %>% filter(line.rep.time %in% x)
                   df <- df %>% pivot_wider(names_from = mut, values_from = freq) %>% as.data.frame()
                   rownames(df) <- df$line.rep.time
                   df <- df[-1] %>% t()
                   df[is.na(df)] <- 0
                   df <- df[,x] %>% as.data.frame()
                   
                   colnames(df) <- colnames(df) %>%
                     gsub("^1", "Sl.L", .) %>%
                     gsub("^2", "St.L", .) %>%
                     gsub("^3", "Nb.L", .) %>%
                     gsub("^4", "CTF.L", .) %>%
                     gsub("^5", "MIX.L", .) %>%
                     
                     gsub("L(\\d)", "L\\1 [", .) %>%
                     gsub("\\[(\\d+)", "[\\1th]", .) %>%
                     gsub("1th]", "1st]", .) %>%
                     gsub("2th]", "2nd]", .) %>%
                     gsub("3th]", "3rd]", .)
                   
                   df
                 }) 

# potato

po.freqs <- mapply(function(y, x) {
  line.rep.time <- y %>% gsub(".*po_", "", .) %>% gsub(".*ori", "PVYSt", .) %>% gsub("_mapped.*", "", .)
  mut <- paste0(x$REF, x$POS, x$ALT)
  freq <- x$INFO %>% gsub(".*AF=", "", .) %>% gsub(";.*", "", .) %>% as.numeric()
  df <- data.frame(mut = mut, freq = freq, line.rep.time = line.rep.time)
  df
},
names(lofreq.po.vcfs),
lofreq.po.vcfs, SIMPLIFY = FALSE) %>% bind_rows()

po.mtx <- lapply(list(#line1pop1time4 = c("0", "114"),
                      #line1pop1time5 = c("0", "115"),
                      `Sl.L2 [1st]` = c("PVYSt", "121"),
                      `St.L1 [6th]` = c("PVYSt", "216"),
                      `St.L2 [10th]` = c("PVYSt", "2210"),
                      `Nb.L1 [10th]` = c("PVYSt", "3110"),
                      `Nb.L2 [10th]` = c("PVYSt", "3210"),
                      `CTF.L2 [1st]` = c("PVYSt", "421"),
                      `MIX.L1 [2nd]` = c("PVYSt", "512"),
                      `St.L1 [4th]` = c("PVYSt", "214"),
                      `St.L2 [4th]` = c("PVYSt", "224"),
                      `Nb.L1 [4th]` = c("PVYSt", "314"),
                      `Nb.L2 [4th]` = c("PVYSt", "324")),
                 function(x) {
                   
                   df <- po.freqs %>% filter(line.rep.time %in% x)
                   df <- df %>% pivot_wider(names_from = mut, values_from = freq) %>% as.data.frame()
                   rownames(df) <- df$line.rep.time
                   df <- df[-1] %>% t()
                   df[is.na(df)] <- 0
                   df <- df[,x] %>% as.data.frame()
                   
                   colnames(df) <- colnames(df) %>%
                     gsub("^1", "Sl.L", .) %>%
                     gsub("^2", "St.L", .) %>%
                     gsub("^3", "Nb.L", .) %>%
                     gsub("^4", "CTF.L", .) %>%
                     gsub("^5", "MIX.L", .) %>%
                     
                     gsub("L(\\d)", "L\\1 [", .) %>%
                     gsub("\\[(\\d+)", "[\\1th]", .) %>%
                     gsub("1th]", "1st]", .) %>%
                     gsub("2th]", "2nd]", .) %>%
                     gsub("3th]", "3rd]", .)
                   
                   df
                 })

# proportion of mutations in inocule that are present at time 4

be.mut_in_4 <- lapply(be.mtx, function(x) {
  mut.0 <- be.freqs %>% filter(line.rep.time == 0) %>% nrow()
  mut.0_in_4 <- x[x[1] > 0 & x[2] > 0,] %>% nrow()
  mut.0_in_4 / mut.0
})

po.mut_in_4 <- lapply(po.mtx, function(x) {
  mut.0 <- po.freqs %>% filter(line.rep.time == 0) %>% nrow()
  mut.0_in_4 <- x[x[1] > 0 & x[2] > 0,] %>% nrow()
  mut.0_in_4 / mut.0
})

# syn and nonsyn

be.syn_nonsyn <- lapply(lofreq.be.vcfs, function(x) {
  if(nrow(x) != 0) check_snp(ref$PVYbe_final_cons, x, 173, 9358)
})

po.syn_nonsyn <- lapply(lofreq.po.vcfs, function(x) {
  if(nrow(x) != 0) check_snp(ref$PVYpo_final_cons, x, 173, 9358)
})

# mutations fixed at time 4

be.fix_in_4 <- lapply(be.mtx, function(x) {

  y <- be.syn_nonsyn %>% bind_rows() %>% apply(2, function(y) y[!is.na(y)] %>% unique())
  syn <- y[y == "syn"] %>% names()
  nonsyn <- y[y == "nonsyn"] %>% names
  
  genome.binary <- rep(0, length(ref$PVYbe_final_cons))
  
  fix <- x[x[2] == 1,] %>% row.names()
  
  fix.pos <- fix %>% gsub("[ACTG]", "", .) %>% as.numeric()
  fix.pos.bin <- genome.binary
  fix.pos.bin[fix.pos] <- 1

  fix.syn <- fix[fix %in% syn]
  fix.syn.pos <- fix.syn %>% gsub("[ACTG]", "", .) %>% as.numeric()
  fix.syn.pos.bin <- genome.binary
  fix.syn.pos.bin[fix.syn.pos] <- 1
  
  fix.nonsyn <- fix[fix %in% nonsyn]
  fix.nonsyn.pos <- fix.nonsyn %>% gsub("[ACTG]", "", .) %>% as.numeric()
  fix.nonsyn.pos.bin <- genome.binary
  fix.nonsyn.pos.bin[fix.nonsyn.pos] <- 1
  
  list(fix = fix,
       fix.pos.bin = fix.pos.bin,
       fix.syn = fix.syn,
       fix.syn.pos.bin = fix.syn.pos.bin,
       fix.nonsyn = fix.nonsyn,
       fix.nonsyn.pos.bin = fix.nonsyn.pos.bin)
})

po.fix_in_4 <- lapply(po.mtx, function(x) {
  
  y <- po.syn_nonsyn %>% bind_rows() %>% apply(2, function(y) y[!is.na(y)] %>% unique())
  syn <- y[y == "syn"] %>% names()
  nonsyn <- y[y == "nonsyn"] %>% names
  
  genome.binary <- rep(0, length(ref$PVYpo_final_cons))
  
  fix <- x[x[2] == 1,] %>% row.names()
  
  fix.pos <- fix %>% gsub("[ACTG]", "", .) %>% as.numeric()
  fix.pos.bin <- genome.binary
  fix.pos.bin[fix.pos] <- 1
  
  fix.syn <- fix[fix %in% syn]
  fix.syn.pos <- fix.syn %>% gsub("[ACTG]", "", .) %>% as.numeric()
  fix.syn.pos.bin <- genome.binary
  fix.syn.pos.bin[fix.syn.pos] <- 1
  
  fix.nonsyn <- fix[fix %in% nonsyn]
  fix.nonsyn.pos <- fix.nonsyn %>% gsub("[ACTG]", "", .) %>% as.numeric()
  fix.nonsyn.pos.bin <- genome.binary
  fix.nonsyn.pos.bin[fix.nonsyn.pos] <- 1
  
  list(fix = fix,
       fix.pos.bin = fix.pos.bin,
       fix.syn = fix.syn,
       fix.syn.pos.bin = fix.syn.pos.bin,
       fix.nonsyn = fix.nonsyn,
       fix.nonsyn.pos.bin = fix.nonsyn.pos.bin)
})

# syn and nonsyn mutations

be.mtx.cat <- be.mtx %>% lapply(function(x) x[2] %>% t() %>% as.data.frame()) %>%
  c(list(`0` = be.mtx[[1]][1] %>% t() %>% as.data.frame())) %>%
  bind_rows() %>% t()

be.mtx.cat[is.na(be.mtx.cat)] <- 0
be.mtx.cat <- be.mtx.cat %>% as.data.frame()

be.common.mut <- be.syn_nonsyn %>% bind_rows(.id = "data") %>% as.data.frame()
row.names(be.common.mut) <- be.common.mut$data %>% gsub(".*be_", "", .) %>% gsub(".*ori", "0", .) %>% gsub("_mapped.*", "", .)
be.common.mut <- be.common.mut[-1] %>% t() %>% as.data.frame()

be.mtx.cat.orf <- be.mtx.cat[row.names(be.mtx.cat) %in% row.names(be.common.mut),] # all mutations in coding region so far

be.common.nonsyn <- be.mtx.cat.orf[apply(be.common.mut, 1, function(x) any(grepl("nonsyn", x))),]
be.common.syn <- be.mtx.cat.orf[apply(be.common.mut, 1, function(x) any(grepl("^syn", x))),]

# potato

po.mtx.cat <- po.mtx %>% lapply(function(x) x[2] %>% t() %>% as.data.frame()) %>%
  c(list(`0` = po.mtx[[1]][1] %>% t() %>% as.data.frame())) %>%
  bind_rows() %>% t()

po.mtx.cat[is.na(po.mtx.cat)] <- 0
po.mtx.cat <- po.mtx.cat %>% as.data.frame()

po.common.mut <- po.syn_nonsyn %>% bind_rows(.id = "data") %>% as.data.frame()
row.names(po.common.mut) <- po.common.mut$data %>% gsub(".*po_", "", .) %>% gsub(".*ori", "0", .) %>% gsub("_mapped.*", "", .)
po.common.mut <- po.common.mut[-1] %>% t() %>% as.data.frame()

po.mtx.cat.orf <- po.mtx.cat[row.names(po.mtx.cat) %in% row.names(po.common.mut),] # all mutations in coding region so far

po.common.nonsyn <- po.mtx.cat.orf[apply(po.common.mut, 1, function(x) any(grepl("nonsyn", x))),]
po.common.syn <- po.mtx.cat.orf[apply(po.common.mut, 1, function(x) any(grepl("^syn", x))),]

# Shannon entropy

be.entropy.mtx <- be.mtx.cat %>%
  apply(2, calc.entropy)

be.entropy.mtx <- be.entropy.mtx/length(ref$PVYbe_final_cons)

po.entropy.mtx <- po.mtx.cat %>%
  apply(2, calc.entropy)

po.entropy.mtx <- po.entropy.mtx/length(ref$PVYpo_final_cons)

# entropy per position

be.entropy.pos.mtx <- be.mtx.cat %>%
  apply(2, calc.entropy, return.sum = FALSE)

be.entropy.pos.mtx <- be.entropy.pos.mtx %>% lapply(function(x) x %>% as.data.frame %>% t() %>% as.data.frame() %>%
                                                      dplyr::mutate(., position = rownames(.) %>%
                                                                      gsub("X", "", .) %>% as.numeric())) %>%
  bind_rows(.id = "sample") %>% rename(entropy = V1)

po.entropy.pos.mtx <- po.mtx.cat %>%
  apply(2, calc.entropy, return.sum = FALSE)

po.entropy.pos.mtx <- po.entropy.pos.mtx %>% lapply(function(x) x %>% as.data.frame %>% t() %>% as.data.frame() %>%
                                                      dplyr::mutate(., position = rownames(.) %>%
                                                                      gsub("X", "", .) %>% as.numeric())) %>%
  bind_rows(.id = "sample") %>% rename(entropy = V1)

# also plot fixed mutations

be.fixed.tb <- list(nonsyn = be.common.nonsyn,
                    syn = be.common.syn,
                    UTR = be.mtx.cat[!rownames(be.mtx.cat) %in% rownames(be.common.mut),]) %>% bind_rows(.id = "type") %>%
  mutate(., position = rownames(.) %>% gsub("[A-Z]", "", .) %>% as.numeric()) %>%
  pivot_longer(-c(position, type), names_to = "sample", values_to = "frequency") %>% dplyr::filter(frequency > 0)

po.fixed.tb <- list(nonsyn = po.common.nonsyn,
                    syn = po.common.syn,
                    UTR = po.mtx.cat[!rownames(po.mtx.cat) %in% rownames(po.common.mut),]) %>% bind_rows(.id = "type") %>%
  mutate(., position = rownames(.) %>% gsub("[A-Z]", "", .) %>% as.numeric()) %>%
  pivot_longer(-c(position, type), names_to = "sample", values_to = "frequency") %>% dplyr::filter(frequency > 0)

entropy.pos.gg <- ggarrange(ggplot() +
                              geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2, color = cistron),
                                           cistron %>% dplyr::mutate(y1 := max(be.entropy.pos.mtx$entropy)/2,
                                                                     y2 := max(be.entropy.pos.mtx$entropy)/2),
                                           linewidth = 3) +
                              
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = be.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "syn") %>%
                                         dplyr::mutate(frequency := max(be.entropy.pos.mtx$entropy)),
                                       color = "grey40", size = .05) +
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = be.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "nonsyn") %>%
                                         dplyr::mutate(frequency := max(be.entropy.pos.mtx$entropy)),
                                       color = "red", size = .05) +
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = be.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "UTR") %>%
                                         dplyr::mutate(frequency := max(be.entropy.pos.mtx$entropy)),
                                       color = "grey70", size = .05) +
                              
                              geom_point(mapping = aes(x = position, y = entropy),
                                       data = be.entropy.pos.mtx,
                                       fill = "black", size = .05) +
                              
                              facet_wrap(~sample %>% factor(levels = sample_levels),
                                         ncol = 1, strip.position = "left") +
                              theme(axis.title = element_blank(),
                                    strip.text = element_text(size = 5),
                                    strip.background = element_blank(),
                                    strip.placement = "outside",
                                    title = element_text(size = 6),
                                    axis.text = element_text(size = 5),
                                    legend.text = element_text(size = 5),
                                    legend.title = element_blank(),
                                    legend.key.size = unit(1, "mm")) +
                              scale_color_viridis(discrete = TRUE) +
                              ggtitle("PVYNb"),
                            
                            ggplot() +
                              geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2, color = cistron),
                                           cistron %>% dplyr::mutate(y1 := max(po.entropy.pos.mtx$entropy)/2,
                                                                     y2 := max(po.entropy.pos.mtx$entropy)/2),
                                           linewidth = 3) +
                              
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = po.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "syn") %>%
                                         dplyr::mutate(frequency := max(po.entropy.pos.mtx$entropy)),
                                       color = "grey40", size = .05) +
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = po.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "nonsyn") %>%
                                         dplyr::mutate(frequency := max(po.entropy.pos.mtx$entropy)),
                                       color = "red", size = .05) +
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = po.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "UTR") %>%
                                         dplyr::mutate(frequency := max(po.entropy.pos.mtx$entropy)),
                                       color = "grey70", size = .05) +
                              
                              geom_point(mapping = aes(x = position, y = entropy),
                                       data = po.entropy.pos.mtx,
                                       color = "black", size = .05) +
                              facet_wrap(~sample %>% factor(levels = sample_levels),
                                         ncol = 1, strip.position = "left") +
                              theme(axis.title = element_blank(),
                                    strip.text = element_text(size = 5),
                                    strip.background = element_blank(),
                                    strip.placement = "outside",
                                    title = element_text(size = 6),
                                    axis.text = element_text(size = 5))  +
                              scale_color_viridis(discrete = TRUE) +
                              ggtitle("PVYSt"),
                            common.legend = TRUE, legend = "right") %>%
  annotate_figure(left = text_grob("Entropy", size = 6, rot = 90),
                  bottom = text_grob("Position", size = 6))

ggsave("entropy.pdf", entropy.pos.gg, width = 6.85, height = 9.21)

# APD

be.apd <-  calc.apd(lofreq.be.vcfs, 9686) %>% setNames(., nm = names(.) %>%
                                                         gsub(".*be_", "", .) %>%
                                                         gsub(".*ori", "PVYNb", .) %>%
                                                         gsub("_mapped.*", "", .) %>%
                                                         
                                                         gsub("^1", "Sl.L", .) %>%
                                                         gsub("^2", "St.L", .) %>%
                                                         gsub("^3", "Nb.L", .) %>%
                                                         gsub("^4", "CTF.L", .) %>%
                                                         gsub("^5", "MIX.L", .) %>%
                                                         
                                                         gsub("L(\\d)", "L\\1 [", .) %>%
                                                         gsub("\\[(\\d+)", "[\\1th]", .) %>%
                                                         gsub("1th]", "1st]", .) %>%
                                                         gsub("2th]", "2nd]", .) %>%
                                                         gsub("3th]", "3rd]", .)) %>%
  as.data.frame() %>% t() %>% as.data.frame() %>%
  dplyr::mutate(., t = rownames(.) %>% gsub(".*\\.\\.", "", .) %>%
                  gsub("[a-z]*\\.$", "", .) %>%
                  as.numeric(),
                treatment = rownames(.) %>% gsub("\\..*", "", .),
                Line = rownames(.) %>% gsub(".*L(\\d)", "\\1", .) %>% gsub("\\.\\..*", "", .)) %>%
  rename(APD = V1)

po.apd <-  calc.apd(lofreq.po.vcfs, 9686) %>% setNames(., nm = names(.) %>%
                                                         gsub(".*po_", "", .) %>%
                                                         gsub(".*ori", "PVYSt", .) %>%
                                                         gsub("_mapped.*", "", .) %>%
                                                         
                                                         gsub("^1", "Sl.L", .) %>%
                                                         gsub("^2", "St.L", .) %>%
                                                         gsub("^3", "Nb.L", .) %>%
                                                         gsub("^4", "CTF.L", .) %>%
                                                         gsub("^5", "MIX.L", .) %>%
                                                         
                                                         gsub("L(\\d)", "L\\1 [", .) %>%
                                                         gsub("\\[(\\d+)", "[\\1th]", .) %>%
                                                         gsub("1th]", "1st]", .) %>%
                                                         gsub("2th]", "2nd]", .) %>%
                                                         gsub("3th]", "3rd]", .)) %>%
  as.data.frame() %>% t() %>% as.data.frame() %>%
  dplyr::mutate(., t = rownames(.) %>% gsub(".*\\.\\.", "", .) %>%
                  gsub("[a-z]*\\.$", "", .) %>%
                  as.numeric(),
                treatment = rownames(.) %>% gsub("\\..*", "", .),
                Line = rownames(.) %>% gsub(".*L(\\d)", "\\1", .) %>% gsub("\\.\\..*", "", .)) %>%
  rename(APD = V1)

be.apd <- rbind(be.apd,
                data.frame(APD = be.apd["PVYNb", "APD"],
                           t = 0,
                           treatment = be.apd$treatment,
                           Line = be.apd$Line)) %>%
  dplyr::filter(treatment != "PVYNb")

po.apd <- rbind(po.apd,
                data.frame(APD = po.apd["PVYSt", "APD"],
                           t = 0,
                           treatment = po.apd$treatment,
                           Line = po.apd$Line)) %>%
  dplyr::filter(treatment != "PVYSt")

# save figures

formatter <- function(...){
  function(x) format(x, ..., scientific = TRUE, digit = 2)
}

be.apd.gg <- be.apd %>% ggplot(aes(t, APD, group = Line, linetype = Line)) +
  geom_line() + facet_wrap(~treatment %>% factor(levels = c("Sl", "St", "Nb", "CTF", "MIX")),
                           scales = "free", nrow = 1) +
  ggtitle("PVYNb") +
  scale_y_continuous(labels = formatter()) +
  theme(axis.text = element_text(size = 5),
        title = element_text(size = 8))

po.apd.gg <- po.apd %>% ggplot(aes(t, APD, group = Line, linetype = Line)) +
  geom_line() + facet_wrap(~treatment %>% factor(levels = c("Sl", "St", "Nb", "CTF", "MIX")),
                           scales = "free", nrow = 1) +
  ggtitle("PVYSt") +
  scale_y_continuous(labels = formatter()) +
  theme(axis.text = element_text(size = 5),
        title = element_text(size = 8))

apd.gg <- ggarrange(be.apd.gg +
                      theme(axis.title = element_blank()) +
                      scale_x_continuous(breaks = seq(0, 10, 2), expand = c(0, 1)),
                    po.apd.gg +
                      theme(axis.title = element_blank()) +
                      scale_x_continuous(breaks = seq(0, 10, 2), expand = c(0, 1)),
                    ncol = 1,
                    common.legend = TRUE, legend = "right") %>%
  annotate_figure(left = "APD", bottom = "Passage", fig.lab = "b")

apd.diff <- list(PVYNb = be.apd %>% dplyr::filter(t != 0) %>%
                   dplyr::mutate(., diff = .$APD - dplyr::filter(be.apd, t == 0)$APD),
                 PVYSt = po.apd %>% dplyr::filter(t != 0) %>%
                   dplyr::mutate(., diff = .$APD - dplyr::filter(po.apd, t == 0)$APD)) %>%
  bind_rows(.id = "virus") %>% dplyr::mutate(treatment := as.factor(treatment), virus := as.factor(virus))

apd.lm <- lm(APD ~ treatment * virus + t, apd.diff)
apd.aov <- aov(apd.lm)

ggsave("apd.pdf", apd.gg, width = 6.85, height = 3.6)

# AFD

afd.pipe <- function(df,
                     rownm = NULL,
                     sum = TRUE,
                     limits = NULL) {
  
  if (!is.null(limits)) {
    df <- df %>% dplyr::filter(., as.numeric(gsub("[ACTG]", "", rownames(.))) %in% limits)
  }
  
  afd <- lapply(setNames(nm = colnames(df)), function(x) {
    lapply(setNames(colnames(df),
                    nm = colnames(df) %>% paste0("vs.", .)), function(y) {
      
      if (!is.null(rownm)) df <- df[rownames(df) %in% rownm,]
      
      calc.afd(df[c(x, y)], sum)
    })
  })
  
  if (!sum) return(afd)
  
  afd <- afd %>% unlist() %>% as.data.frame()
    
  colnames(afd) <- "total AFD"
  afd$sample1 <- rownames(afd) %>% gsub(".*\\.vs\\.", "", .)
  afd$sample2 <- rownames(afd) %>% gsub("\\.vs\\..*", "", .)
  
  return(afd)
}

be.afd <- afd.pipe(be.mtx.cat)
po.afd <- afd.pipe(po.mtx.cat)

afd.tb <- bind_rows(list(PVYNb = be.afd, PVYSt = po.afd), .id = "isolate")
afd.tb$sample1 <- factor(afd.tb$sample1, levels = sample_levels)

afd.tb$sample2 <- factor(afd.tb$sample2, levels = sample_levels)

afd.gg <- afd.tb %>%
  ggplot(aes(sample1, sample2, fill = `total AFD`)) +
  geom_tile() +
  scale_fill_viridis() +
  facet_wrap(~isolate, scales = "free") +
  scale_y_discrete(limits = rev) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))

ggsave("afd.pdf", afd.gg, width = 6.85, height = 3.6)

ggsave("afd_apd.pdf",
       ggarrange(afd.gg %>% annotate_figure(fig.lab = "a"),
                 apd.gg,
                 ncol = 1),
       width = 6.85, height = 7)

# save results

write_xlsx(list(`pct of snps from 0` = bind_rows(list(PVYNb = be.mut_in_4,
                                                      PVYSt = po.mut_in_4),
                                                 .id = "isolate"),
                
                `n snps` = data.frame(`n snps` = c(lofreq.be.vcfs %>% lapply(nrow) %>% unlist(),
                                                   lofreq.po.vcfs %>% lapply(nrow) %>% unlist()),
                                      sample = c(lofreq.be.vcfs %>% lapply(nrow) %>% unlist(),
                                                 lofreq.po.vcfs %>% lapply(nrow) %>% unlist()) %>% names()),
                
                `n snps fixed` = bind_rows(list(be = be.fix_in_4 %>% lapply(function(x) x$fix %>% length()),
                                                po = po.fix_in_4 %>% lapply(function(x) x$fix %>% length())),
                                           .id = "isolate"),
                
                `n snps fix syn` = bind_rows(list(be = be.fix_in_4 %>% lapply(function(x) x$fix.syn %>% length()),
                                                  po = po.fix_in_4 %>% lapply(function(x) x$fix.syn %>% length())),
                                             .id = "isolate"),
                
                `n snps fix nonsyn` = bind_rows(list(be = be.fix_in_4 %>% lapply(function(x) x$fix.nonsyn %>% length()),
                                                     po = po.fix_in_4 %>% lapply(function(x) x$fix.nonsyn %>% length())),
                                                .id = "isolate"),
                
                `diversity (Shannon entropy)` = bind_rows(list(be = be.entropy.mtx,
                                                               po = po.entropy.mtx),
                                                          .id = "isolate"),
                `AFD` = afd.tb),
           "results.xlsx")

# mean nonsyn AFD / mean syn AFD

afd.pipe.time <- function(df,
                          rownm = NULL,
                          sum = TRUE,
                          limits = NULL) {
  
  if (!is.null(limits)) {
    df <- df %>% dplyr::filter(., as.numeric(gsub("[ACTG]", "", rownames(.))) %in% limits)
  }
  
  afd <- lapply(setNames(2:ncol(df),
                         nm = paste0(colnames(df[1:(ncol(df) - 1)]),
                                     ".vs.",
                                     colnames(df[2:ncol(df)]))),
                function(x) {
                  
                  if (!is.null(rownm)) df <- df[rownames(df) %in% rownm,]
                  
                  df <- df[c(x - 1, x)]
                  df <- df[rowSums(df) > 0,]
                  
                  afd.res <- calc.afd(df, sum)
                  
                  if (sum) return(list(AFD = afd.res,
                                       n.variable.snps = nrow(df)))
                  
                  return(list(AFD = afd.res,
                              n.variable.snps = nrow(df),
                              Position = names(afd.res) %>% gsub("\\..*", "", .) %>% as.numeric()))
  }) %>% bind_rows(.id = "name")
  
  if (!sum) return(afd)

  afd$sample1 <- afd$name %>% gsub("\\.vs\\..*", "", .)
  afd$sample2 <- afd$name %>% gsub(".*\\.vs\\.", "", .)
  
  return(afd %>% dplyr::select(-name))
}

rownames(cistron) <- cistron$cistron

be.afd.syn_nonsyn.per_cistron <- apply(cistron %>% filter(cistron != "PIPO"), 1, function(x) {
  
  limits <- x["x1"]:x["x2"]
  
  lapply(be.mtx.cat %>% colnames() %>%
           gsub(" \\[.*", "", .) %>%
           grep("PVY", ., value = TRUE, invert = TRUE) %>%
           unique() %>% setNames(nm = .),
         function(y) {
           
           # sort samples in line by time and include original isolate
           
           samples <- c("PVYNb", grep(y, colnames(be.mtx.cat), value = TRUE))
           samples <- samples[c(1,
                                order(grep(y, colnames(be.mtx.cat), value = TRUE) %>%
                                        gsub(".* \\[", "", .) %>%
                                        gsub("[a-z]*\\]", "", .) %>%
                                        as.numeric()) + 1)]
           
           df.line.time <- be.mtx.cat[samples]
           
           list(nonsyn = afd.pipe.time(df.line.time, rownames(be.common.nonsyn), TRUE, limits),
                syn = afd.pipe.time(df.line.time, rownames(be.common.syn), TRUE, 173:9361)) %>% bind_rows(.id = "mutation") %>%
             dplyr::mutate(Passage = sample2 %>% gsub(".* \\[", "", .) %>%
                             gsub("[a-z]*\\]", "", .) %>%
                             as.numeric(),
                           `Mean AFD` = AFD / n.variable.snps)
         }) %>% bind_rows(.id = "Line")
  
}) %>% bind_rows(.id = "Cistron")

po.afd.syn_nonsyn.per_cistron <- apply(cistron %>% filter(cistron != "PIPO"), 1, function(x) {
  
  limits <- x["x1"]:x["x2"]
  
  lapply(po.mtx.cat %>% colnames() %>%
           gsub(" \\[.*", "", .) %>%
           grep("PVY", ., value = TRUE, invert = TRUE) %>%
           unique() %>% setNames(nm = .),
         function(y) {
           
           # sort samples in line by time and include original isolate
           
           samples <- c("PVYSt", grep(y, colnames(po.mtx.cat), value = TRUE))
           samples <- samples[c(1,
                                order(grep(y, colnames(po.mtx.cat), value = TRUE) %>%
                                        gsub(".* \\[", "", .) %>%
                                        gsub("[a-z]*\\]", "", .) %>%
                                        as.numeric()) + 1)]
           
           df.line.time <- po.mtx.cat[samples]
           
           list(nonsyn = afd.pipe.time(df.line.time, rownames(po.common.nonsyn), TRUE, limits),
                syn = afd.pipe.time(df.line.time, rownames(po.common.syn), TRUE, 173:9361)) %>% bind_rows(.id = "mutation") %>%
             dplyr::mutate(Passage = sample2 %>% gsub(".* \\[", "", .) %>%
                             gsub("[a-z]*\\]", "", .) %>%
                             as.numeric(),
                           `Mean AFD` = AFD / n.variable.snps)
         }) %>% bind_rows(.id = "Line")
}) %>% bind_rows(.id = "Cistron")

mean.afd.syn_nonsyn <- list(PVYNb = be.afd.syn_nonsyn.per_cistron,
                            PVYSt = po.afd.syn_nonsyn.per_cistron) %>%
  lapply(function(x) {
    
    x %>%
      dplyr::select(Cistron, Line, mutation, Passage, `Mean AFD`) %>%
      pivot_wider(names_from = mutation, values_from = `Mean AFD`)
    
  }) %>%
  bind_rows(.id = "Virus") %>%
  dplyr::select(Cistron, Line, Passage, Virus, nonsyn, syn) %>%
  dplyr::mutate(`Mean AFD nonsyn / mean AFD syn` = nonsyn / syn) %>%
  dplyr::rename(`Mean AFD syn` = syn,
                `Mean AFD nonsyn` = nonsyn)

mean.afd.syn_nonsyn.cistron.gg <- mean.afd.syn_nonsyn %>%
  dplyr::mutate(Treatment = Line %>% gsub("\\..*", "", .),
                Sample = paste0(Line, "[", Passage, "]")) %>%
  ggplot(aes(Cistron %>% factor(levels = cistron$cistron), `Mean AFD nonsyn / mean AFD syn`)) +
  geom_boxplot() +
  geom_point(aes(color = Treatment)) +
  geom_text_repel(aes(label = Sample), size = 1.5, seed = 100) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~Virus) +
  xlab("Cistron") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

ggsave("mean.afd.syn_nonsyn.cistron.pdf", mean.afd.syn_nonsyn.cistron.gg, width = 6.85, height = 5)

# AFD by position, then for nonsyn, calculate (AFD - (mean AFD syn)) / (SD syn)
# not in final manuscript

be.afd.syn_nonsyn.per_site <- apply(cistron %>% filter(cistron != "PIPO"), 1, function(x) {
  
  limits <- x["x1"]:x["x2"]
  
  lapply(be.mtx.cat %>% colnames() %>%
           gsub(" \\[.*", "", .) %>%
           grep("PVY", ., value = TRUE, invert = TRUE) %>%
           unique() %>% setNames(nm = .),
         function(y) {
           
           # sort samples in line by time and include original isolate
           
           samples <- c("PVYNb", grep(y, colnames(be.mtx.cat), value = TRUE))
           samples <- samples[c(1,
                                order(grep(y, colnames(be.mtx.cat), value = TRUE) %>%
                                        gsub(".* \\[", "", .) %>%
                                        gsub("[a-z]*\\]", "", .) %>%
                                        as.numeric()) + 1)]
           
           df.line.time <- be.mtx.cat[samples]
           
           syn <- afd.pipe.time(df.line.time, rownames(be.common.syn), FALSE, 173:9361)
           
           summarise.syn <- syn %>% group_by(name) %>% summarise(`Mean syn AFD` = mean(AFD), `SD syn AFD` = sd(AFD))
           
           afd.tb <- list(nonsyn = afd.pipe.time(df.line.time, rownames(be.common.nonsyn), FALSE, limits) %>%
                  filter(n.variable.snps > 0),
                syn = syn) %>% bind_rows(.id = "mutation") %>%
             mutate(.,
                    Passage = name %>% gsub(".* \\[", "", .) %>%
                      gsub("[a-z]*\\]", "", .) %>%
                      as.numeric(),
                    Sample = name %>% gsub(".*vs\\.", "", .))
           
           mean.syn.afd <- summarise.syn[match(afd.tb$name, summarise.syn$name), "Mean syn AFD"] %>% unlist()
           sd.syn.afd <- summarise.syn[match(afd.tb$name, summarise.syn$name), "SD syn AFD"] %>% unlist()
           
           afd.tb %>% mutate(`Standardized AFD` = (AFD - mean.syn.afd) / sd.syn.afd,
                             `Mean syn AFD` = mean.syn.afd,
                             `SD syn AFD` = sd.syn.afd)
           
         }) %>% bind_rows(.id = "Line")
  
}) %>% bind_rows(.id = "Cistron")

po.afd.syn_nonsyn.per_site <- apply(cistron %>% filter(cistron != "PIPO"), 1, function(x) {
  
  limits <- x["x1"]:x["x2"]
  
  lapply(po.mtx.cat %>% colnames() %>%
           gsub(" \\[.*", "", .) %>%
           grep("PVY", ., value = TRUE, invert = TRUE) %>%
           unique() %>% setNames(nm = .),
         function(y) {
           
           # sort samples in line by time and include original isolate
           
           samples <- c("PVYSt", grep(y, colnames(po.mtx.cat), value = TRUE))
           samples <- samples[c(1,
                                order(grep(y, colnames(po.mtx.cat), value = TRUE) %>%
                                        gsub(".* \\[", "", .) %>%
                                        gsub("[a-z]*\\]", "", .) %>%
                                        as.numeric()) + 1)]
           
           df.line.time <- po.mtx.cat[samples]
           
           syn <- afd.pipe.time(df.line.time, rownames(po.common.syn), FALSE, 173:9361)
           
           summarise.syn <- syn %>% group_by(name) %>% summarise(`Mean syn AFD` = mean(AFD), `SD syn AFD` = sd(AFD))
           
           afd.tb <- list(nonsyn = afd.pipe.time(df.line.time, rownames(po.common.nonsyn), FALSE, limits) %>%
                            filter(n.variable.snps > 0),
                          syn = syn) %>% bind_rows(.id = "mutation") %>%
             mutate(.,
                    Passage = name %>% gsub(".* \\[", "", .) %>%
                      gsub("[a-z]*\\]", "", .) %>%
                      as.numeric(),
                    Sample = name %>% gsub(".*vs\\.", "", .))
           
           mean.syn.afd <- summarise.syn[match(afd.tb$name, summarise.syn$name), "Mean syn AFD"] %>% unlist()
           sd.syn.afd <- summarise.syn[match(afd.tb$name, summarise.syn$name), "SD syn AFD"] %>% unlist()
           
           afd.tb %>% mutate(`Standardized AFD` = (AFD - mean.syn.afd) / sd.syn.afd,
                             `Mean syn AFD` = mean.syn.afd,
                             `SD syn AFD` = sd.syn.afd)
           
         }) %>% bind_rows(.id = "Line")
}) %>% bind_rows(.id = "Cistron")

afd.syn_nonsyn.per_site <- list(PVYNb = be.afd.syn_nonsyn.per_site,
                                PVYSt = po.afd.syn_nonsyn.per_site) %>%
  bind_rows(.id = "Virus") %>%
  dplyr::mutate(Treatment = Line %>% gsub("\\..*", "", .))

standardized.afd.gg <- ggarrange(be.afd.syn_nonsyn.per_site %>%
                                   ggplot(aes(x = Position, y = `Standardized AFD`)) +
                                   
                                   geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2, color = cistron),
                                                cistron %>% dplyr::mutate(y1 := 0,
                                                                          y2 := 0),
                                                linewidth = 3) +
                                   
                                   geom_point(data = be.afd.syn_nonsyn.per_site %>%
                                                dplyr::filter(mutation == "syn"),
                                              color = "black", size = .025) +
                                   geom_point(data = be.afd.syn_nonsyn.per_site %>%
                                                dplyr::filter(mutation == "nonsyn"),
                                              color = "red", size = .025) +
                                   
                                   facet_wrap(~Sample %>% factor(levels = sample_levels),
                                              ncol = 1, scales = "free_y", strip.position = "left") +
                                   theme(axis.title = element_blank(),
                                         strip.text = element_text(size = 5),
                                         strip.background = element_blank(),
                                         strip.placement = "outside",
                                         title = element_text(size = 6),
                                         axis.text = element_text(size = 5),
                                         legend.text = element_text(size = 5),
                                         legend.title = element_blank(),
                                         legend.key.size = unit(1, "mm")) +
                                   scale_color_viridis(discrete = TRUE) +
                                   ggtitle("PVYNb"),
                                 
                                 po.afd.syn_nonsyn.per_site %>%
                                   ggplot(aes(x = Position, y = `Standardized AFD`)) +
                                   geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2, color = cistron),
                                                cistron %>% dplyr::mutate(y1 := 0,
                                                                          y2 := 0),
                                                linewidth = 3) +
                                   geom_point(data = po.afd.syn_nonsyn.per_site %>%
                                                dplyr::filter(mutation == "syn"),
                                              color = "black", size = .025) +
                                   geom_point(data = po.afd.syn_nonsyn.per_site %>%
                                                dplyr::filter(mutation == "nonsyn"),
                                              color = "red", size = .025) +
                                   facet_wrap(~Sample %>% factor(levels = sample_levels),
                                              ncol = 1, scales = "free_y", strip.position = "left") +
                                   theme(axis.title = element_blank(),
                                         strip.text = element_text(size = 5),
                                         strip.background = element_blank(),
                                         strip.placement = "outside",
                                         title = element_text(size = 6),
                                         axis.text = element_text(size = 5),
                                         legend.text = element_text(size = 5),
                                         legend.title = element_blank(),
                                         legend.key.size = unit(1, "mm")) +
                                   scale_color_viridis(discrete = TRUE) +
                                   ggtitle("PVYSt"),
                                 common.legend = TRUE, legend = "right") %>%
  annotate_figure(bottom = text_grob("Position", size = 6),
                  left = text_grob("Standardized AFD", size = 6, rot = 90))

ggsave("standardized_afd.pdf", standardized.afd.gg, width = 6.85, height = 9.21)

afd.syn_nonsyn.per_site %>% filter(mutation == "nonsyn" & `Standardized AFD` > 1.5) %>% print(n = 70)
