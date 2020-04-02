library(sgtr)
setwd("/home/romain/work/code/sgtr/test/")

###################################
########## DATA LOADING ###########
###################################

# Load chromosome names from file
chromosomes = load_chromosome_names("chromosomes.tsv")

# Load psass window output
psass_window_chr = load_genome_metrics("psass_window.tsv", chromosomes)
psass_window_chr_detect = load_genome_metrics("psass_window.tsv", chromosomes=NULL, detect_chromosomes = TRUE)
psass_window_no_chr = load_genome_metrics("psass_window.tsv", detect_chromosomes = FALSE)
psass_window_no_unplaced = load_genome_metrics("psass_window_no_unplaced.tsv")

# Load psass snp output
psass_snp_chr = load_genome_metrics("psass_snps.tsv", chromosomes)
psass_snp_chr_detect = load_genome_metrics("psass_snps.tsv", chromosomes=NULL, detect_chromosomes = TRUE)
psass_snp_no_chr = load_genome_metrics("psass_snps.tsv", detect_chromosomes = FALSE)
psass_snp_no_unplaced = load_genome_metrics("psass_snps_no_unplaced.tsv")

# Load psass fst output
psass_fst_chr = load_genome_metrics("psass_fst.tsv", chromosomes)
psass_fst_chr_detect = load_genome_metrics("psass_fst.tsv", chromosomes=NULL, detect_chromosomes = TRUE)
psass_fst_no_chr = load_genome_metrics("psass_fst.tsv", detect_chromosomes = FALSE)
psass_fst_no_unplaced = load_genome_metrics("psass_fst_no_unplaced.tsv")

# Draw circos plot
draw_circos(psass_window_chr$data, psass_window_chr$lengths,
            tracks = list(circos_track("Fst", label = expression("F"["ST"])),
                          circos_track(c("Snps_females", "Snps_males"), label = "SNPs", color = c("firebrick2", "dodgerblue3")),
                          circos_track(c("Abs_depth_females", "Abs_depth_males"), label = "Depth", color = c("firebrick2", "dodgerblue3"))),
            output.file = "circos.png")


# Plot circos
plot_circos("psass_window.tsv",
            tracks = list(track("Fst", label = expression("F"["ST"])),
                          track(c("Snps_females", "Snps_males"), label = "SNPs", colors = c("firebrick2", "dodgerblue3"))),
            chromosomes_file = "chromosomes.tsv",
            output_file = "circos2.png")

# Plot circos multiple tracks
plot_circos("psass_window.tsv",
            tracks = list(track("Fst", label = expression(bold("F"["ST"])),
                                major_lines_y = TRUE),
                          track("Snps_males", label = "M. SNPs", colors = c("dodgerblue3"),
                                major_lines_x = TRUE),
                          track("Snps_females", label = "F. SNPs", colors = c("firebrick2"),
                                major_lines_y = FALSE)),
            chromosomes_file = "chromosomes.tsv",
            output_file = "circos3.png",
            highlight = 'Chr24',
            highlight_bg_color = "thistle1",
            default_bg_colors = "white")


# Plot circos
plot_circos("psass_window.tsv",
            tracks = list(circos_track("Fst", label = expression("F"["ST"])),
                          circos_track(c("Snps_females", "Snps_males"), label = "SNPs", colors = c("firebrick2", "dodgerblue3")),
                          circos_track(c("Abs_depth_females", "Abs_depth_males"), label = "Depth", color = c("firebrick2", "dodgerblue3"))),
            chromosomes.file.path = "chromosomes.tsv",
            output.file = "circos2.png")

# Draw region plot
region = draw_region(psass_window_chr$data, psass_window_chr$lengths, "Chr24",
                     tracks = list(track("Fst", label = expression("F"["ST"])),
                                   track(c("Snps_females", "Snps_males"), label = "Pool-specific SNPs", color = c("firebrick2", "dodgerblue3"), alpha=0.6),
                                   track(c("Abs_depth_females", "Abs_depth_males"), label = "Absolute depth", color = c("firebrick2", "dodgerblue3"), alpha=c(0.4, 0.4))),
                     output.file = "region.png", width = 12, track.height = 4, res = 300)

# Plot region
plot_region("psass_window.tsv", "Chr24:0-6000000",
            tracks = list(track("Fst", label = expression("F"["ST"])),
                          track(c("Snps_females", "Snps_males"), label = "Pool-specific SNPs", color = c("firebrick2", "dodgerblue3"), alpha=0.6),
                          track(c("Abs_depth_females", "Abs_depth_males"), label = "Absolute depth", color = c("firebrick2", "dodgerblue3"), alpha=c(0.4, 0.4))),
            chromosomes.file.path = "chromosomes.tsv",
            output.file = "region2.png", width = 12, track.height = 4, res = 300)

# Draw manhattan
draw_manhattan_plot(psass_window_chr$data, psass_window_chr$lengths,
                    tracks = list(manhattan_track("Fst", label = expression("F"["ST"])),
                                  manhattan_track("Snps_females", label = "SNPs", point.color = c("firebrick1", "firebrick4"))),
                    output.file = "manhattan.png")

# Plot manhattan
plot_manhattan("psass_window.tsv",
               tracks = list(manhattan_track("Fst", label = expression("F"["ST"])),
                             manhattan_track("Snps_females", label = "SNPs", point.color = c("firebrick1", "firebrick4"))),
               output.file = "manhattan.png")
