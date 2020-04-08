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

# Plot circos
plot_circos("psass_window.tsv",
            tracks = list(single_metric_track("Fst", type = "ribbon"),
                          multi_metrics_track(c("Snps_females", "Snps_males"),
                                              label = "SNPs",
                                              colors = c("red", "blue"),
                                              type = "points")),
            chromosomes_file = "chromosomes.tsv",
            highlight = 'Chr24',
            highlight_bg_color = "thistle1",
            output_file = "circos.png")


# Plot circos
plot_manhattan("psass_window.tsv",
               tracks = list(single_metric_track("Fst", point_size = 0.25, colors = c("grey10", "grey50")),
                             single_metric_track("Snps_females", color = c("red", "firebrick3")),
                             single_metric_track("Snps_males", color = c("blue", "dodgerblue3"))),
               chromosomes_file = "chromosomes.tsv",
               output_file = "manhattan.png")

# Plot circos
manhattan <- plot_manhattan("psass_window.tsv",
               tracks = list(single_metric_track("Snps_males", label = "Male-specific SNPs"),
                             single_metric_track("Snps_females", label = "Female-specific SNPs",
                                                 alpha = 0.5, point_size = 1.5)),
               chromosomes_file = "chromosomes.tsv",
               chromosomes_as_numbers = TRUE,
               output_file = "manhattan.png")

# Plot region
region <- plot_region("psass_window.tsv",
                      "Chr24:0-6000000",
                      tracks = list(single_metric_track("Snps_males", label = "Male-specific SNPs"),
                                    single_metric_track("Snps_females", label = "Female-specific SNPs",
                                                        alpha = 0.5, point_size = 1.5, type = "points")),
                      chromosomes_file = "chromosomes.tsv",
                      output_file = "region.png",
                      width = 12,
                      track_height = 4,
                      res = 300)

region <- plot_region("psass_window.tsv",
                      "Chr24:0-3000000",
                      tracks = list(single_metric_track("Fst", label = "Fst", color = "darkgreen", alpha = 1,
                                                        type = "points"),
                                    multi_metrics_track(c("Snps_females", "Snps_males"),
                                                        metric_labels = c("F. SNPs", "M. SNPs"),
                                                        label = "SNPs",
                                                        colors = c("red", "blue"),
                                                        type = "ribbon",
                                                        legend_position = c(0.8, 0.5))),
                      chromosomes_file = "chromosomes.tsv",
                      output_file = "region2.png",
                      width = 12,
                      track_height = 4,
                      res = 300)

