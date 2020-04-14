library(sgtr)

# Default distrib plot
radsex_distrib("radsex/distrib.tsv",
               output_file = "distrib.png")

# Distrib plot with specified groups and labels
radsex_distrib("radsex/distrib.tsv",
               groups = c("M", "F"),
               group_labels = c("Males", "Females"),
               output_file = "distrib2.png")

# Default marker depths plot
radsex_marker_depths("radsex/signif.tsv",
                     output_file = "signif.png")

# Marker depths plot with specified groups, labels, and group info to color
# individual labels
radsex_marker_depths("radsex/signif.tsv",
                     group_info_file = "popmap.tsv",
                     group_labels = c("Females", "Males"),
                     label_colors = c("firebrick2", "dodgerblue3"),
                     output_file = "signif2.png")

# Default circos plot
radsex_map_circos("radsex/map.tsv",
                  chromosomes_file = "chromosomes.tsv",
                  output_file = "circos.png")

# Default region plot for Chr01
radsex_map_region("radsex/map.tsv",
                  "Chr01",
                  chromosomes_file = "chromosomes.tsv",
                  output_file = "region.png")

# Default manhattan plot
radsex_map_manhattan("radsex/map.tsv",
                     chromosomes_file = "chromosomes.tsv",
                     output_file = "manhattan.png")

# Depth boxplot for number of reads with groups, labels, and colors
radsex_depth("radsex/depths.tsv",
             groups = c("M", "F"),
             group_labels = c("Males", "Females"),
             group_colors = c("green", "purple"),
             metric = "Reads",
             ylab = "Number of reads")

# Depth barplot for median depth with labels and colors
radsex_depth("radsex/depths.tsv",
             type = "barplot",
             group_colors = c("dodgerblue3", "red2"),
             metric = "Median_depth",
             ylab = "Median read depth")

# Frequencies plot
radsex_freq("radsex/freq.tsv",
            color = "purple",
            alpha = 0.5,
            output_file = "freq.png")
