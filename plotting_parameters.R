### Script containing plotting parameters  
#############################################################################

plot_width = 11
rank_plot_width = 13
plot_height = 8
### Response plotting colours

FF_colour = "black"
PBE_colour =  "#009E73"
PBE0_colour = "#0072B2"
PBE_subcolour = "black"
PBE0_subcolour = "black"

### Plots with predicted/true values, uncertainly bars, and highlighted structures
pred_point_type = 4
ff_rank_point_type = 3

pred_point_size = 4
pred_rank_point_size = 3
ff_rank_point_size = 3
abline_width = 2
interval_bar_width = 0.5
label_text_size = 8

pred_point_colour = "#999999"
ff_rank_colour = "#009E73"
abline_colour = "#56B4E9"
label_colour = "#D55E00"

pred_ylab = "Predicted response (kJ/mol)\n"
true_xlab = "\nTrue response (kJ/mol)"

title_size = 28
subtitle_size = 24
axis_title_size = 24
axis_tick_label_size = 22
legend_title_size = 22
legend_text_size = 22
label_size = 20   # for alpha and beta labels

### Barcharts
barplot_width = 13
barplot_height = 8

## Colourblind-friendly palettes
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette_13 = c("#000000", RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"), "#999999")

