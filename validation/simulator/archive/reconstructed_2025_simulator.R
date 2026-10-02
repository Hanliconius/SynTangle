# Reconstructed 2025 chromosome mutation simulator
#
# Historical reference only. This file is reconstructed from surviving chat
# transcript code blocks and is not claimed to be a byte-for-byte snapshot of
# one original file. The mutation logic below is preserved as recovered so that
# future validation work can distinguish historical behavior from later fixes.
#
# Known historical rough edges are intentionally left visible here; use
# ../chromosome_simulator.R for current validation work.
#
# This archived function expects dplyr/tibble verbs to be available.

mutate_genome_structure <- function(genome) {
  genes <- genome$genes
  chroms <- genome$chroms
  used_chroms <- character()
  mutation_log <- list()

  mutate_id <- 1

  available_chroms <- unique(chroms$chrom)
  n_mutations <- sample(2:4)

  for (i in seq_len(n_mutations)) {
    possible_chroms <- setdiff(available_chroms, used_chroms)
    if (length(possible_chroms) < 1) break

    action <- sample(c("fusion", "fission", "inversion"), 1)

    if (action == "fusion" && length(possible_chroms) >= 2) {
      chr_pair <- sample(possible_chroms, 2)
      used_chroms <- c(used_chroms, chr_pair)

      genes1 <- genes %>% filter(chrom == chr_pair[1])
      genes2 <- genes %>% filter(chrom == chr_pair[2])

      offset <- max(genes1$end) + 10000
      genes2 <- genes2 %>%
        mutate(start = start + offset, end = end + offset)

      fused_genes <- bind_rows(genes1, genes2) %>%
        mutate(chrom = paste0("fusion_", chr_pair[1], "_", chr_pair[2], "_", mutate_id))

      genes <- genes %>% filter(!chrom %in% chr_pair)
      genes <- bind_rows(genes, fused_genes)

      new_length <- max(fused_genes$end)
      chroms <- chroms %>% filter(!chrom %in% chr_pair)
      chroms <- bind_rows(chroms, tibble(
        chrom = fused_genes$chrom[1],
        length = new_length,
        species = unique(fused_genes$species)
      ))

      mutation_log[[mutate_id]] <- paste("Fusion:", chr_pair[1], "+", chr_pair[2], "→", fused_genes$chrom[1])
    }

    if (action == "fission" && length(possible_chroms) >= 1) {
      chr_split <- sample(possible_chroms, 1)
      used_chroms <- c(used_chroms, chr_split)

      genes_split <- genes %>% filter(chrom == chr_split)
      split_point <- median(genes_split$start)

      f1 <- genes_split %>% filter(start <= split_point) %>%
        mutate(chrom = paste0(chr_split, "_f1_", mutate_id))
      f2 <- genes_split %>% filter(start > split_point) %>%
        mutate(
          chrom = paste0(chr_split, "_f2_", mutate_id),
          start = start - min(start) + 1,
          end = end - min(start) + 1
        )

      genes <- genes %>% filter(chrom != chr_split)
      genes <- bind_rows(genes, f1, f2)

      chroms <- chroms %>% filter(chrom != chr_split)
      chroms <- bind_rows(chroms,
                          tibble(chrom = f1$chrom[1], length = max(f1$end), species = unique(f1$species)),
                          tibble(chrom = f2$chrom[1], length = max(f2$end), species = unique(f2$species)))

      mutation_log[[mutate_id]] <- paste("Fission:", chr_split, "→", f1$chrom[1], "+", f2$chrom[1])
    }

    if (action == "inversion" && length(possible_chroms) >= 1) {
      chr_inv <- sample(possible_chroms, 1)
      used_chroms <- c(used_chroms, chr_inv)

      genes_chr <- genes %>% filter(chrom == chr_inv)
      chr_len <- max(genes_chr$end)

      min_size <- as.integer(0.4 * chr_len)
      max_size <- as.integer(0.9 * chr_len)
      inv_size <- sample(min_size:max_size, 1)
      inv_start <- sample(1:(chr_len - inv_size), 1)
      inv_end <- inv_start + inv_size

      inv_genes <- genes_chr %>%
        filter(start >= inv_start & end <= inv_end) %>%
        mutate(
          new_start = inv_end - (end - inv_start),
          new_end = inv_end - (start - inv_start),
          start = new_start,
          end = new_end,
          strand = ifelse(strand == "+", "-", "+")
        ) %>%
        select(-new_start, -new_end)

      genes_chr <- genes_chr %>%
        filter(start < inv_start | end > inv_end) %>%
        bind_rows(inv_genes) %>%
        arrange(start)

      genes <- genes %>% filter(chrom != chr_inv)
      genes <- bind_rows(genes, genes_chr)

      mutation_log[[mutate_id]] <- sprintf("Inversion: %s (%0.0f-%0.0f)", chr_inv, inv_start, inv_end)
    }

    mutate_id <- mutate_id + 1
  }

  list(genes = genes, chroms = chroms, log = mutation_log)
}
