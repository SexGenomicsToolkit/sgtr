library(sgtr)

setwd("examples")

###################################
########## CIRCOS PLOTS ###########
###################################

# Circos plot showing a Fst ribbon track with alternating colors and a combined
# SNPs points track with colors associated to groups, highlighting Chr24.
plot_circos("psass/psass_window.tsv",
            tracks = list(single_metric_track("Fst", type = "ribbon",
                                              colors = c("green", "purple")),
                          multi_metrics_track(c("Snps_females", "Snps_males"),
                                              label = "SNPs",
                                              colors = c("red", "blue"),
                                              type = "points")),
            chromosomes_file = "psass/chromosomes.tsv",
            highlight = 'Chr24',
            highlight_bg_color = "thistle1",
            output_file = "circos.png")


# Circos plot showing a Fst ribbon track with alternating colors and a combined
# SNPs points track with colors associated to groups, highlighting Chr24.
plot_circos("psass/psass_window.tsv",
            tracks = list(single_metric_track("Fst",
                                              type = "points",
                                              colors = c("green", "purple")),
                          multi_metrics_track(c("Snps_females", "Snps_males"),
                                              label = "SNPs",
                                              colors = c("red", "blue"),
                                              type = "points")),
            chromosomes_file = "psass/chromosomes.tsv",
            highlight = 'Chr24',
            highlight_bg_color = "thistle1",
            output_file = "circos.png")

###################################
######### MANHATTAN PLOTS #########
###################################

# Manhattan plot showing a Fst track with alternating grey colors, a females
# SNPs track with alternating red colors, and a male SNPs track with alternating
# blue colors
plot_manhattan("psass/psass_window.tsv",
               tracks = list(single_metric_track("Fst",
                                                 point_size = 0.25,
                                                 colors = c("grey10", "grey50")),
                             single_metric_track("Snps_females",
                                                 color = c("red", "firebrick3")),
                             single_metric_track("Snps_males",
                                                 color = c("blue", "dodgerblue3"),
                                                 h_lines = list(h_line(y = 75,
                                                                       label = "test",
                                                                       color = "red",
                                                                       type = 2)
                                                                )
                                                 )
                             ),
               chromosomes_file = "psass/chromosomes.tsv",
               output_file = "manhattan.png")

# Manhattan plot showing example of customization for a males and females SNPs
# tracks
manhattan <- plot_manhattan("psass/psass_window.tsv",
               tracks = list(single_metric_track("Snps_males",
                                                 label = "Male-specific SNPs"),
                             single_metric_track("Snps_females",
                                                 label = "Female-specific SNPs",
                                                 alpha = 0.5,
                                                 point_size = 1.5)),
               chromosomes_file = "psass/chromosomes.tsv",
               chromosomes_as_numbers = TRUE,
               output_file = "manhattan.png")

###################################
########## REGION PLOTS ###########
###################################

# Region plot showing a males SNPs ribbon track and a females SNPs points track
# for the first 6Mbp of Chr24
region <- plot_region("psass/psass_window.tsv",
                      "Chr24:0-6000000",
                      tracks = list(single_metric_track("Snps_males",
                                                        label = "Male-specific SNPs",
                                                        h_lines = list(h_line(y=50,
                                                                              label = "yo"))),
                                    single_metric_track("Snps_females",
                                                        label = "Female-specific SNPs",
                                                        colors = "purple",
                                                        alpha = 0.5,
                                                        point_size = 1.5,
                                                        type = "points")),
                      chromosomes_file = "psass/chromosomes.tsv",
                      output_file = "region.png",
                      width = 12,
                      track_height = 4,
                      res = 300)

# Region plot showing a Fst points track and a combined SNPs ribbon track with
# labels for the entire Chr12.
region <- plot_region("psass/psass_window.tsv",
                      "Chr12",
                      tracks = list(single_metric_track("Fst",
                                                        label = "Fst",
                                                        color = "darkgreen",
                                                        alpha = 1,
                                                        type = "points"),
                                    multi_metrics_track(c("Snps_females", "Snps_males"),
                                                        metric_labels = c("F. SNPs", "M. SNPs"),
                                                        label = "SNPs",
                                                        colors = c("red", "blue"),
                                                        type = "ribbon",
                                                        legend_position = c(0.8, 0.5))),
                      chromosomes_file = "psass/chromosomes.tsv",
                      output_file = "region2.png",
                      width = 12,
                      track_height = 4,
                      res = 300)


###################################
########## DATA LOADING ###########
###################################

# Load chromosome names from file
chromosomes = load_chromosome_names("psass/chromosomes.tsv")

# Load psass window output
psass_window_chr = load_genome_metrics("psass/psass_window.tsv",
                                       chromosomes)
psass_window_chr_detect = load_genome_metrics("psass/psass_window.tsv",
                                              chromosomes=NULL,
                                              detect_chromosomes = TRUE)
psass_window_no_chr = load_genome_metrics("psass/psass_window.tsv",
                                          detect_chromosomes = FALSE)
psass_window_no_unplaced = load_genome_metrics("psass/psass_window_no_unplaced.tsv")

# Load psass snp output
psass_snp_chr = load_genome_metrics("psass/psass_snps.tsv",
                                    chromosomes)
psass_snp_chr_detect = load_genome_metrics("psass/psass_snps.tsv",
                                           chromosomes=NULL,
                                           detect_chromosomes = TRUE)
psass_snp_no_chr = load_genome_metrics("psass/psass_snps.tsv",
                                       detect_chromosomes = FALSE)
psass_snp_no_unplaced = load_genome_metrics("psass/psass_snps_no_unplaced.tsv")

# Load psass fst output
psass_fst_chr = load_genome_metrics("psass/psass_fst.tsv", chromosomes)
psass_fst_chr_detect = load_genome_metrics("psass/psass_fst.tsv", chromosomes=NULL,
                                           detect_chromosomes = TRUE)
psass_fst_no_chr = load_genome_metrics("psass/psass_fst.tsv",
                                       detect_chromosomes = FALSE)
psass_fst_no_unplaced = load_genome_metrics("psass/psass_fst_no_unplaced.tsv")

