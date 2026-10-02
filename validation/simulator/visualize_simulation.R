# Lightweight base-R visual audit for synthetic validation cases.

map_display_positions <- function(genome, display_state, gap = 0.08) {
  sp <- unique(genome$chroms$species)[1]
  state <- display_state[display_state$species == sp, , drop = FALSE]
  state <- state[order(state$display_rank), , drop = FALSE]
  lengths <- genome$chroms$length[match(state$chrom, genome$chroms$chrom)]
  widths <- lengths / max(lengths)
  starts <- numeric(length(widths))
  if (length(widths) > 1L) {
    for (i in 2:length(widths)) starts[i] <- starts[i - 1L] + widths[i - 1L] + gap
  }
  data.frame(
    species = sp,
    chrom = state$chrom,
    orientation = state$orientation,
    x0 = starts,
    width = widths,
    length = lengths,
    stringsAsFactors = FALSE
  )
}

anchor_x <- function(chrom, pos, map) {
  i <- match(chrom, map$chrom)
  if (is.na(i)) return(NA_real_)
  frac <- pos / map$length[i]
  if (map$orientation[i] < 0) frac <- 1 - frac
  map$x0[i] + frac * map$width[i]
}

plot_layout_panel <- function(genomes, display_state, main = "layout", link_alpha = 0.18) {
  species_names <- names(genomes)
  maps <- lapply(genomes, map_display_positions, display_state = display_state)
  names(maps) <- species_names
  xmax <- max(vapply(maps, function(m) max(m$x0 + m$width), numeric(1)))
  y <- rev(seq_along(species_names))

  plot(NA, xlim = c(0, xmax), ylim = c(0.5, length(species_names) + 0.5),
       xlab = "", ylab = "", yaxt = "n", xaxt = "n", bty = "n", main = main)
  axis(2, at = y, labels = species_names, las = 1, tick = FALSE)

  # Homology links are drawn first so chromosome bars remain legible.
  if (length(species_names) >= 2L) {
    for (i in seq_len(length(species_names) - 1L)) {
      g1 <- genomes[[species_names[i]]]
      g2 <- genomes[[species_names[i + 1L]]]
      shared <- intersect(g1$genes$name, g2$genes$name)
      a <- g1$genes[match(shared, g1$genes$name), , drop = FALSE]
      b <- g2$genes[match(shared, g2$genes$name), , drop = FALSE]
      for (j in seq_along(shared)) {
        x1 <- anchor_x(a$chrom[j], (a$start[j] + a$end[j]) / 2, maps[[i]])
        x2 <- anchor_x(b$chrom[j], (b$start[j] + b$end[j]) / 2, maps[[i + 1L]])
        segments(x1, y[i], x2, y[i + 1L],
                 col = grDevices::adjustcolor("grey30", alpha.f = link_alpha))
      }
    }
  }

  for (i in seq_along(species_names)) {
    m <- maps[[i]]
    for (j in seq_len(nrow(m))) {
      segments(m$x0[j], y[i], m$x0[j] + m$width[j], y[i], lwd = 7, lend = 2)
      label <- if (m$orientation[j] > 0) m$chrom[j] else paste0(m$chrom[j], " (rev)")
      text(m$x0[j] + m$width[j] / 2, y[i] + 0.12, label, cex = 0.55, srt = 35, adj = 0)
    }
  }
}

write_simulation_audit_pdf <- function(genomes, native_state, tangled_state, path) {
  grDevices::pdf(path, width = 16, height = 9, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mfrow = c(1, 2), mar = c(1, 5, 3, 1))
  plot_layout_panel(genomes, native_state, main = "Extant chromosomes: native presentation")
  plot_layout_panel(genomes, tangled_state, main = "Same biology: deliberately tangled presentation")
  invisible(path)
}
