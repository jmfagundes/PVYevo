source("../PVYevo/source_me.R")

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

mean.depth <- depths %>% lapply(function(x) mean(x$depth)) %>% unlist()

depths.gg <- depths %>% bind_rows() %>%
  ggplot(aes(position, depth)) + geom_col(color = "grey40") +
  theme(axis.text = element_text(size = 6),
        strip.text = element_text(size = 8)) +
  geom_text(aes(label = paste0("Mean: ", round(mean_depth, digits = 2)), y = max/2, x = 5000),
            data = depths %>% bind_rows() %>%
              group_by(vir.line.rep.time) %>%
              summarise(mean_depth = mean(depth),
                        sd_depth = sd(depth),
                        max = max(depth)),
            color = "black", size = 2) +
  geom_text(aes(label = paste0("SD: ", round(sd_depth, digits = 2)), y = max/4, x = 5000),
            data = depths %>% bind_rows() %>%
              group_by(vir.line.rep.time) %>%
              summarise(mean_depth = mean(depth),
                        sd_depth = sd(depth),
                        max = max(depth)),
            color = "black", size = 2) +
  facet_wrap(~vir.line.rep.time, scales = "free", ncol = 4) +
  ylab("Depth") + xlab("Position")

ggsave("cov.pdf", depths.gg, width = 6.85, height = 9.21)

# compare with deltaQc

R_qPCR_data <- read_xlsx("R_qPCR_data.xlsx")
R_qPCR_data$Sample <- paste0(R_qPCR_data$Virus %>%
                               gsub("PVYNb", "PVYbe", .) %>%
                               gsub("PVYSt", "PVYpo", .),
                             "_",
                             R_qPCR_data$Group %>%
                               gsub("Sl", "1", .) %>%
                               gsub("St", "2", .) %>%
                               gsub("Nb", "3", .) %>%
                               gsub("CTF", "4", .) %>%
                               gsub("MIX", "5", .),
                             R_qPCR_data$Line,
                             R_qPCR_data$Passage) %>% gsub("None.*", "ori", .)

mean.depth.corr.names <- setNames(mean.depth, nm = names(mean.depth) %>% gsub("_mapped.*", "", .))

R_qPCR_data$Mean_cov <- mean.depth.corr.names[match(R_qPCR_data$Sample,
                                                    names(mean.depth.corr.names))]

R_qPCR_data %>% filter(!is.na(Mean_cov)) %>%
  ggplot(aes(`-dCq`, Mean_cov)) + geom_point() + geom_smooth(method = "lm")

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

# export GTF

GTF <- data.frame(source = "lofreq",
                  feature = "CDS",
                  start = cistron$x1,
                  end = cistron$x2,
                  score = ".",
                  strand = "+",
                  frame = "0",
                  attributes = paste0("gene_id",
                                      " ",
                                      cistron$cistron %>%
                                        gsub("^", "\"", .) %>%
                                        gsub("$", "\";", .)))

write.table(bind_cols(seqname = "PVYbe_final_cons",
                      GTF), file = "PVYbe.gtf", col.names = FALSE, quote = FALSE, row.names = FALSE,
            sep = "\t")

write.table(bind_cols(seqname = "PVYpo_final_cons",
                      GTF), file = "PVYpo.gtf", col.names = FALSE, quote = FALSE, row.names = FALSE,
            sep = "\t")

# levels for plots

sample_levels <- c("PVYNb", "PVYSt",
                   "Sl.L1 [1st]", "Sl.L2 [1st]",
                   "St.L1 [4th]", "St.L1 [6th]", "St.L2 [4th]", "St.L2 [10th]",
                   "Nb.L1 [4th]", "Nb.L1 [10th]",
                   "Nb.L2 [4th]", "Nb.L2 [10th]",
                   "CTF.L1 [4th]", "CTF.L2 [1st]", "CTF.L2 [4th]",
                   "MIX.L1 [2nd]", "MIX.L1 [3rd]", "MIX.L2 [4th]", "MIX.L2 [9th]")

# filter matrices based on coverage

# match file names

names(mean.depth) <- names(mean.depth) %>% gsub("sorted.*", "recal_lofreq.vcf", .)

be.samples.discarded <- lofreq.be.vcfs[!names(lofreq.be.vcfs) %in% names(mean.depth[mean.depth > 100])] %>% names()
po.samples.discarded <- lofreq.po.vcfs[!names(lofreq.po.vcfs) %in% names(mean.depth[mean.depth > 100])] %>% names()

lofreq.be.vcfs <- lofreq.be.vcfs[names(lofreq.be.vcfs) %in% names(mean.depth[mean.depth > 100])]
lofreq.po.vcfs <- lofreq.po.vcfs[names(lofreq.po.vcfs) %in% names(mean.depth[mean.depth > 100])]

# eliminate samples with no SNPs after filtering

be.samples.discarded <- c(be.samples.discarded, names(lofreq.be.vcfs[lofreq.be.vcfs %>% lapply(function(x) nrow(x) == 0) %>% unlist()]))
po.samples.discarded <- c(po.samples.discarded, names(lofreq.po.vcfs[lofreq.po.vcfs %>% lapply(function(x) nrow(x) == 0) %>% unlist()]))

lofreq.be.vcfs <- lofreq.be.vcfs[lofreq.be.vcfs %>% lapply(function(x) nrow(x) > 0) %>% unlist()]
lofreq.po.vcfs <- lofreq.po.vcfs[lofreq.po.vcfs %>% lapply(function(x) nrow(x) > 0) %>% unlist()]

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

# filter alleles based on frequency

#be.freqs <- be.freqs[be.freqs$freq >= freq.thresh,]

be.mtx <- lapply(list(#`Sl.L1 [1st]` = c("PVYNb", "111"),
                      `Sl.L2 [1st]` = c("PVYNb", "121"),
                      `St.L1 [4th]` = c("PVYNb", "214"),
                      #line2pop1time5 = c("0", "215"),
                      `St.L2 [4th]` = c("PVYNb", "224"),
                      #line2pop2time5 = c("0", "225"),
                      `Nb.L1 [10th]` = c("PVYNb", "3110"),
                      `Nb.L2 [10th]` = c("PVYNb", "3210"),
                      #line4pop2time8 = c("0", "428"),
                      `MIX.L1 [3rd]` = c("PVYNb", "513"),
                      `MIX.L2 [9th]` = c("PVYNb", "529"),
                      `Nb.L1 [4th]` = c("PVYNb", "314"),
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

# filter alleles based on frequency

#po.freqs <- po.freqs[po.freqs$freq >= freq.thresh,]

po.mtx <- lapply(list(#line1pop1time4 = c("0", "114"),
                      #line1pop1time5 = c("0", "115"),
                      #`Sl.L2 [1st]` = c("PVYSt", "121"),
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
  mut.0 <- be.freqs %>% filter(line.rep.time == "PVYNb") %>% nrow()
  mut.0_in_4 <- x[x[1] > 0 & x[2] > 0,] %>% nrow()
  mut.0_in_4 / mut.0
})

po.mut_in_4 <- lapply(po.mtx, function(x) {
  mut.0 <- po.freqs %>% filter(line.rep.time == "PVYSt") %>% nrow()
  mut.0_in_4 <- x[x[1] > 0 & x[2] > 0,] %>% nrow()
  mut.0_in_4 / mut.0
})

# syn and nonsyn

be.syn_nonsyn <- lapply(lofreq.be.vcfs, function(x) {
  if(nrow(x) != 0) check_snp.fix(ref$PVYbe_final_cons, x, 173, 9355)
})

po.syn_nonsyn <- lapply(lofreq.po.vcfs, function(x) {
  if(nrow(x) != 0) check_snp.fix(ref$PVYpo_final_cons, x, 173, 9355)
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
row.names(be.common.mut) <- be.common.mut$data %>% gsub(".*be_", "", .) %>% gsub(".*ori", "0", .) %>%
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
  gsub("3th]", "3rd]", .) %>%
  gsub("^0", "PVYNb", .)
be.common.mut <- be.common.mut[-1] %>% t() %>% as.data.frame()

be.mtx.cat.orf <- be.mtx.cat[row.names(be.mtx.cat) %in% row.names(be.common.mut),] # all mutations in coding region so far

be.common.mut <- be.common.mut[be.mtx.cat.orf %>% colnames()]
be.common.mut <- be.common.mut[be.mtx.cat.orf %>% rownames(),]

substitute.freq.syn_nonsyn <- function(x, y, mut) {
  
  mtx <- x
  
  for (i in 1:nrow(x)) {
    
    row <- x[i,]
    
    for (j in 1:length(row)) {
      
      if (y[i, j] != mut | is.na(y[i, j])) mtx[i, j] <- NA
      
    }
  }
  return(mtx)
}

be.common.nonsyn <- substitute.freq.syn_nonsyn(be.mtx.cat.orf, be.common.mut, "nonsyn")
be.common.syn <- substitute.freq.syn_nonsyn(be.mtx.cat.orf, be.common.mut, "syn")
be.common.stop <- substitute.freq.syn_nonsyn(be.mtx.cat.orf, be.common.mut, "stop")

# potato

po.mtx.cat <- po.mtx %>% lapply(function(x) x[2] %>% t() %>% as.data.frame()) %>%
  c(list(`0` = po.mtx[[1]][1] %>% t() %>% as.data.frame())) %>%
  bind_rows() %>% t()

po.mtx.cat[is.na(po.mtx.cat)] <- 0
po.mtx.cat <- po.mtx.cat %>% as.data.frame()

po.common.mut <- po.syn_nonsyn %>% bind_rows(.id = "data") %>% as.data.frame()
row.names(po.common.mut) <- po.common.mut$data %>% gsub(".*po_", "", .) %>% gsub(".*ori", "0", .) %>%
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
  gsub("3th]", "3rd]", .) %>%
  gsub("^0", "PVYSt", .)
po.common.mut <- po.common.mut[-1] %>% t() %>% as.data.frame()

po.mtx.cat.orf <- po.mtx.cat[row.names(po.mtx.cat) %in% row.names(po.common.mut),] # all mutations in coding region so far

po.common.mut <- po.common.mut[po.mtx.cat.orf %>% colnames()]
po.common.mut <- po.common.mut[po.mtx.cat.orf %>% rownames(),]

po.common.nonsyn <- substitute.freq.syn_nonsyn(po.mtx.cat.orf, po.common.mut, "nonsyn")
po.common.syn <- substitute.freq.syn_nonsyn(po.mtx.cat.orf, po.common.mut, "syn")
po.common.stop <- substitute.freq.syn_nonsyn(po.mtx.cat.orf, po.common.mut, "stop")

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
                    stop = be.common.stop,
                    UTR = be.mtx.cat[!rownames(be.mtx.cat) %in% rownames(be.common.mut),]) %>% bind_rows(.id = "type") %>%
  mutate(., position = rownames(.) %>% gsub("[A-Z]", "", .) %>% gsub("\\..*[0-9]$", "", .) %>% as.numeric()) %>%
  pivot_longer(-c(position, type), names_to = "sample", values_to = "frequency") %>% dplyr::filter(frequency > 0)

po.fixed.tb <- list(nonsyn = po.common.nonsyn,
                    syn = po.common.syn,
                    stop = po.common.stop,
                    UTR = po.mtx.cat[!rownames(po.mtx.cat) %in% rownames(po.common.mut),]) %>% bind_rows(.id = "type") %>%
  mutate(., position = rownames(.) %>% gsub("[A-Z]", "", .) %>% gsub("\\..*[0-9]$", "", .) %>% as.numeric()) %>%
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
                                       color = "grey40", linewidth = .05) +
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = be.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "nonsyn") %>%
                                         dplyr::mutate(frequency := max(be.entropy.pos.mtx$entropy)),
                                       color = "red", linewidth = .05) +
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = be.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "stop") %>%
                                         dplyr::mutate(frequency := max(be.entropy.pos.mtx$entropy)),
                                       color = "grey40", linewidth = .05) +
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = be.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "UTR") %>%
                                         dplyr::mutate(frequency := max(be.entropy.pos.mtx$entropy)),
                                       color = "grey70", linewidth = .05) +
                              
                              geom_point(mapping = aes(x = position, y = entropy),
                                       data = be.entropy.pos.mtx,
                                       fill = "black", size = .01) +
                              
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
                                       color = "grey40", linewidth = .05) +
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = po.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "nonsyn") %>%
                                         dplyr::mutate(frequency := max(po.entropy.pos.mtx$entropy)),
                                       color = "red", linewidth = .05) +
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = po.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "stop") %>%
                                         dplyr::mutate(frequency := max(po.entropy.pos.mtx$entropy)),
                                       color = "grey40", linewidth = .05) +
                              geom_col(mapping = aes(x = position, y = frequency),
                                       data = po.fixed.tb %>%
                                         dplyr::filter(frequency == 1 & type == "UTR") %>%
                                         dplyr::mutate(frequency := max(po.entropy.pos.mtx$entropy)),
                                       color = "grey70", linewidth = .05) +
                              
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

# FST 

fst <- lapply(setNames(c("pvybe", "pvypo"), nm = c("PVYNb", "PVYSt")), function(x) {
  list(fst = read.table(paste0("sam_aln/", x, "_fullgenome.fst.tb"), col.names = c("sample1", "sample2", "fst")),
       sample.names = read.table("sam_aln/sample_order.txt", col.names = "sample") %>%
         filter(grepl(x, sample, ignore.case = TRUE)))
})

fst$PVYNb$sample.names <- fst$PVYNb$sample.names %>%
  mutate(sample := sample %>% gsub(".*be_", "", .) %>%
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
           gsub("3th]", "3rd]", .))

fst$PVYSt$sample.names <- fst$PVYSt$sample.names %>%
  mutate(sample := sample %>%
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
           gsub("3th]", "3rd]", .))

correct.fst.sample.names <- function() {
  lapply(fst, function(x) {
    
    out <- data.frame(sample1 = character(),
                      sample2 = character(),
                      fst = numeric())
    
    for (i in 1:nrow(x$fst)) {
      
      row <- x$fst[i,]
      sample1 <- x$sample.names[row$sample1, "sample"]
      sample2 <- x$sample.names[row$sample2, "sample"]
      
      out <- out %>% rbind(list(sample1 = sample1, sample2 = sample2, fst = row$fst))
    }
    samples1 <- out$sample1
    out <- out %>% rbind(out %>% mutate(sample1 := sample2,
                                        sample2 := samples1))
  })
}

fst.tb <- correct.fst.sample.names() %>%
  bind_rows(.id = "isolate") %>%
  mutate(sample1 := factor(sample1, levels = sample_levels),
         sample2 := factor(sample2, levels = sample_levels))

fst.gg <- fst.tb %>%
  ggplot(aes(sample1, sample2, fill = `fst`)) +
  geom_tile() +
  scale_fill_viridis(name = "F<sub>ST</sub>") +
  facet_wrap(~isolate, scales = "free") +
  scale_y_discrete(limits = rev) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.title = ggtext::element_markdown(size = 8),
        legend.text = element_text(size = 6))

ggsave("fst_apd.pdf",
       ggarrange(fst.gg %>% annotate_figure(fig.lab = "a"),
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
                `FST` = fst.tb),
           "results.xlsx")

# results

# number of fixed SNPs

be.fixed.tb %>% filter(frequency == 1) %>% dplyr::select(sample) %>% table()
po.fixed.tb %>% filter(frequency == 1) %>% dplyr::select(sample) %>% table()

# syn/nonsyn fixed SNPs

be.fixed.tb$sample %>% unique() %>% setNames(nm = .) %>%
  lapply(function(x) {
    
    tb <- be.fixed.tb %>% dplyr::filter(sample == x)
    syn <- tb %>% filter(frequency == 1 & type == "syn") %>% nrow()
    nonsyn <- tb %>% filter(frequency == 1 & type == "nonsyn") %>% nrow()
    stop <- tb %>% filter(frequency == 1 & type == "stop") %>% nrow()
    list(syn = syn, nonsyn = nonsyn, stop = stop)
    
  })

po.fixed.tb$sample %>% unique() %>% setNames(nm = .) %>%
  lapply(function(x) {
    
    tb <- po.fixed.tb %>% dplyr::filter(sample == x)
    syn <- tb %>% filter(frequency == 1 & type == "syn") %>% nrow()
    nonsyn <- tb %>% filter(frequency == 1 & type == "nonsyn") %>% nrow()
    stop <- tb %>% filter(frequency == 1 & type == "stop") %>% nrow()
    list(syn = syn, nonsyn = nonsyn, stop = stop)
    
  })

# SNPgenie

snpgenie.products <- lapply(paste0("lofreqout_vslow/products_results/", list.files("lofreqout_vslow/products_results/", "product_results.txt")),
                            function(x) {
                              read.table(x, header = TRUE) %>%
                                mutate(mean_gdiv_polymorphic := mean_gdiv_polymorphic %>%
                                         gsub("\\*", "0", .) %>% as.numeric(),
                                       mean_N_gdiv := mean_N_gdiv %>%
                                         gsub("\\*", "0", .) %>% as.numeric(),
                                       mean_S_gdiv := mean_S_gdiv %>%
                                         gsub("\\*", "0", .) %>% as.numeric())
                            }) %>% bind_rows() %>%
  mutate(Virus = file %>%
           gsub("_.*", "", .) %>%
           gsub("PVYbe", "PVYNb", .) %>%
           gsub("PVYpo", "PVYSt", .),
         Treatment = file %>%
           gsub("PVYbe_ori.*", "PVYNb", .) %>%
           gsub("PVYpo_ori.*", "PVYSt", .) %>%
           gsub("PVY[b,p]._", "", .) %>%
           gsub("^1", "Sl.L", .) %>%
           gsub("^2", "St.L", .) %>%
           gsub("^3", "Nb.L", .) %>%
           gsub("^4", "CTF.L", .) %>%
           gsub("^5", "MIX.L", .) %>%
           gsub("\\.L.*", "", .),
         Line = file %>%
           gsub("PVY.._[0-9]", "", .) %>%
           gsub("[0-9]_.*|10_.*", "", .),
         Passage = file %>%
           gsub("PVY.._[0-9][0-9]", "", .) %>%
           gsub("_.*", "", .))

snpgenie.products.tb <- snpgenie.products %>%
  dplyr::mutate(`piN/piS` = piN/piS,
                `piN/piS` := replace(`piN/piS`, `piN/piS` == Inf | is.nan(`piN/piS`) | `piN/piS` == 0, NA)) #%>%
  #rename(`Mean Genetic Diversity` = mean_gdiv_polymorphic)

snpgenie.products.lm <- lm(`piN/piS` ~ Virus + Treatment + product, snpgenie.products.tb %>%
                             dplyr::select(product, `piN/piS`, #`Mean Genetic Diversity`,
                                           Virus, Treatment, Line, Passage))

snpgenie.gg <- snpgenie.products.tb %>%
  dplyr::select(product, `piN/piS`, #`Mean Genetic Diversity`, 
                Virus, Treatment, Line, Passage) %>%
  pivot_longer(-c(product, Virus, Treatment, Line, Passage), names_to = "metric") %>%
  dplyr::mutate(Sample = paste0(Treatment, ".L", Line, "[", Passage, "]") %>%
                  gsub("\\.LPVY.*", "", .),
                Treatment := Treatment %>% gsub("PVY.*", "Original", .) %>% factor(levels = c("Original", "Sl", "St", "Nb", "CTF", "MIX")),
                metric := metric %>% gsub("piN/piS", "pi*N/pi*S", .) %>% gsub(" ", "~", .)) %>%
  rename(Cistron = product) %>%
  dplyr::filter(!is.na(value)) %>%
  ggplot(aes(Cistron %>% factor(levels = cistron$cistron), value)) +
  geom_boxplot() +
  geom_point(aes(color = Treatment)) +
  geom_text_repel(aes(label = Sample), size = 1.5, seed = 100) +
  facet_grid(rows = vars(metric), cols = vars(Virus), scales = "free", switch = "y", labeller = label_parsed) +
  xlab("Cistron") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.y = element_blank(),
        strip.background.y = element_blank(),
        strip.placement = "outside")

snpgenie.sup.gg <- snpgenie.products.tb %>%
  dplyr::select(product, mean_gdiv_polymorphic, mean_N_gdiv, mean_S_gdiv, piN, piS,
                Virus, Treatment, Line, Passage) %>%
  pivot_longer(-c(product, Virus, Treatment, Line, Passage), names_to = "metric") %>%
  dplyr::mutate(Sample = paste0(Treatment, ".L", Line, "[", Passage, "]") %>%
                  gsub("\\.LPVY.*", "", .),
                Treatment := Treatment %>% gsub("PVY.*", "Ancestral", .) %>% factor(levels = c("Ancestral", "Sl", "St", "Nb", "CTF", "MIX")),
                metric := metric %>% gsub("piN/piS", "pi*N/pi*S", .) %>% gsub(" ", "~", .)) %>%
  rename(Cistron = product) %>%
  dplyr::filter(!is.na(value)) %>%
  ggplot(aes(Cistron %>% factor(levels = cistron$cistron), value)) +
  geom_boxplot() +
  geom_point(aes(color = Treatment)) +
  geom_text_repel(aes(label = Sample), size = 1.5, seed = 100) +
  facet_grid(rows = vars(metric), cols = vars(Virus), scales = "free", switch = "y", labeller = label_parsed) +
  xlab("Cistron") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title.y = element_blank(),
        strip.background.y = element_blank(),
        strip.placement = "outside")

ggsave("snpgenie_results.pdf", snpgenie.gg, width = 6.85, height = 3)
ggsave("snpgenie_results_sup.pdf", snpgenie.sup.gg, width = 6.85, height = 9.21)
