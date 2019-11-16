# Preamble ----------------------------------------------------------------

# Load libraries
packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter",
           "RColorBrewer",  "mixtools","matrixStats", "sva", "ggplot2", "sva")
loaded_libs <- sapply(packs, library, character.only = T, warn.conflicts = F, quietly = F)



# Multiplot  --------------------------------------------------------------

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# Plotting ----------------------------------------------------------------

MixA <- c("Jurkat" = 2.5, "IM-9" = 1.25, "Raji" =  2.5,"THP-1" = 3.75)/10
MixB <-  c("Jurkat" = 0.5, "IM-9" = 3.17, "Raji" = 4.75, "THP-1" = 1.58)/10
MixC <-  c("Jurkat" = 0.1, "IM-9" = 4.95, "Raji" = 1.65, "THP-1" = 3.3)/10
MixD  <- c("Jurkat" = 0.02, "IM-9" = 3.33, "Raji" = 3.33, "THP-1" = 3.33)/10
Mix <- rbind(MixA, MixB, MixC, MixD)
Mix <- data.table("Mix" = paste0("Mix", LETTERS[1:4]), Mix)
print(Mix)



plot_dt_cibersort <- fread(file.path(".", "sim_results", "cibersort_abbas.txt"))
plot_dt <- fread(file.path(".", "sim_results", "plot_dt_music_abbas_v2.txt"))
p_sims <- fread(file.path(".", "sim_results", "p_sims_abbas.txt"))
plot_dt_LM1 <- fread(file.path(".", "sim_results", "plot_dt_LM1_abbas.txt"))
p_sims_LM1 <- fread(file.path(".", "sim_results", "p_sims_LM1_abbas.txt"))
plot_dt_ts <- fread(file.path(".", "sim_results", "plot_dt_TS_abbas_v2.txt"))
plot_dt_ts_ave <- fread(file.path(".", "sim_results", "plot_dt_TS_ave_abbas.txt"))

p_sims_ts_sva <- fread(file.path(".", "sim_results", "fullTS_v2.txt"))
p_sims_ts <- fread(file.path(".", "sim_results", "p_sims_fullTS_with_filter_abbas.txt"))
p_sims_ts_with_sva_ave <- fread(file.path(".", "sim_results", "fullTS_v2_ave.txt"))

cols <- RColorBrewer::brewer.pal(length(unique(plot_dt$type)), "Dark2")
names(cols) <- unique(plot_dt$type)

r_plot_dt <- with(plot_dt, cor(x = Truth, y = Estimate))
r_plot_dt_LM1 <- with(plot_dt_LM1, cor(x = Truth, y = Estimate))
r_plot_dt_cibersort <- with(plot_dt_cibersort, cor(x = Truth, y = Estimate))
r_plot_dt_ts <- with(plot_dt_ts, cor(x = Truth, y = Estimate))
r_plot_dt_ts_ave <- with(plot_dt_ts_ave, cor(x = Truth, y = Estimate))
r_plot_dt <- round(r_plot_dt, 2)
r_plot_dt_LM1 <- round(r_plot_dt_LM1, 2)
r_plot_dt_cibersort <- round(r_plot_dt_cibersort, 2)
r_plot_dt_ts <- round(r_plot_dt_ts, 2)
r_plot_dt_ts_ave <- round(r_plot_dt_ts_ave, 2)

p1 <- ggplot(data = plot_dt, aes(x = Truth, y = Estimate, color = type, shape = type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x = 0.05, y = .4, label = paste("r", "=", r_plot_dt)) + 
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ggtitle("Model estimates vs Truth \n(MuSiC)") + 
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 8),
        legend.position = c(0.9,0.3),
        legend.title = element_blank())

p2 <- ggplot(data = plot_dt_cibersort, aes(x = Truth, y = Estimate, color = type, shape = type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x = 0.05, y = .4, label = paste("r", "=", r_plot_dt_cibersort)) + 
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ggtitle("Model estimates vs Truth \n(Cibersort)") + 
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 8),
        legend.position = c(0.9,0.3),
        legend.title = element_blank())

p3 <- ggplot(data = plot_dt_LM1, aes(x = Truth, y = Estimate, color = type, shape = type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x = 0.05, y = .4, label = paste("r", "=", r_plot_dt_LM1)) + 
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ggtitle("Model estimates vs Truth \n(LM; p > .75)") + 
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 8),
        legend.position = c(0.9,0.3),
        legend.title = element_blank())

p4 <- ggplot(data = plot_dt_ts, aes(x = Truth, y = Estimate, color = type, shape = type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x = 0.05, y = .4, label = paste("r", "=", r_plot_dt_ts)) + 
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ggtitle("Model estimates vs Truth \n(2-Step; p > .75)") + 
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 8),
        legend.position = c(0.9,0.3),
        legend.title = element_blank())
p5 <- ggplot(data = plot_dt_ts_ave, aes(x = Truth, y = Estimate, color = type, shape = type)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) +
  annotate("text", x = 0.05, y = .4, label = paste("r", "=", r_plot_dt_ts_ave)) + 
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ggtitle("Model estimates vs Truth \n(2-Step; Model averaged)") + 
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 8),
        legend.position = c(0.9,0.3),
        legend.title = element_blank())

png(file.path(".", "figs", "c1.png"), res = 200, width = 1920, height = 1080)
multiplot(p1, p2, p3, p4, p5, cols=3)
dev.off()

cols <- RColorBrewer::brewer.pal(length(unique(p_sims$type)), "Dark2")
names(cols) <- unique(p_sims$type)

g <- ggplot(data = p_sims, aes(x = p, y = ngenes, color = cols[1])) + 
  geom_point() + 
  geom_smooth(se = F) +
  transparent_legend + remove_grid +
  scale_y_continuous(limits = c(50, 54700), breaks = seq(50, 54700, by = 1e4)) +
  ylab("# Genes in Reference") + xlab("CDF") +
  ggtitle("Relationship between \ncount of genes and cutoff") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = "none",
        legend.title = element_blank())

png(file.path(".", "figs", "c3.png"), res = 200, width = 1280, height = 720)
g
dev.off()

g1 <- ggplot(data = p_sims_LM1, aes(x = p, y = pearson, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ylim(.6, 1) +
  ylab("R") + xlab("CDF") +
  ggtitle("Relationship between \ncutoff and model fit (LM)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.75,0.25),
        legend.title = element_blank())

l <- min(min(p_sims_LM1$rmse), min(p_sims_ts_sva$rmse))
l <- round(l, 2)
u <- max(max(p_sims_LM1$rmse), max(p_sims_ts_sva$rmse))
u <- round(u, 2)
my_breaks <- seq(l, u+5e-3, by = 5e-3)

g2 <- ggplot(data = p_sims_LM1, aes(x = p, y = rmse, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  scale_y_continuous(breaks = my_breaks ) +
  # ylab(expression(sqrt(over(sum( (p[i] - hat(p)[i])^2, i==1, k),k)))) +
  ylab("RMSE") +
  xlab("CDF") +
  ggtitle("Relationship between \ncutoff and model fit (LM)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.5,0.23),
        legend.direction="horizontal",
        legend.title = element_blank())

l_r <- min(min(p_sims_LM1$mean_ratio), min(p_sims_ts_sva$mean_ratio))
l_r <- round(log(l_r), 2)
u_r <- max(max(p_sims_LM1$mean_ratio), max(p_sims_ts_sva$mean_ratio))
u_r <- round(log(u_r), 2)
my_breaks_r <- seq(l_r, u_r+5e-1, by = 5e-1)

g3 <- ggplot(data = p_sims_LM1, aes(x = p, y = log(mean_ratio), color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  # ylim(0, NA) +
  scale_y_continuous(breaks = my_breaks_r ) +
  geom_hline(yintercept = 0) +
  ylab("Mean Ratio") + xlab("CDF") +
  ggtitle("Relationship between \ncutoff and model fit (LM)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.80,0.75),
        legend.title = element_blank())

# ggplot(data = p_sims_LM1, aes(x = pearson^2, y = rmse, color = type, shape = type)) + 
#   geom_point() + 
#   geom_smooth(se = F) +
#   scale_color_manual(values = cols) +
#   transparent_legend + remove_grid +
#   scale_y_continuous(breaks = my_breaks ) +
#   ylab(expression(sqrt(sum( (p[i] - hat(p)[i])^2, i==1, k)))) +
#   xlab("Pearson correlation coefficient") +
#   ggtitle("Relationship between \ncutoff and model fit (LM)") +
#   theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
#         text = element_text(size = 12),
#         axis.title = element_text(face="bold", size = 9),
#         axis.text.y=element_text(size = 8, face="bold"),
#         axis.text.x=element_text(size = 8, face="bold"),
#         legend.position = c(0.75,0.25),
#         legend.direction="horizontal",
#         legend.title = element_blank())

g4 <- ggplot(data = p_sims_ts_sva, aes(x = p, y = pearson, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ylim(.6, 1) +
  ylab("R") + xlab("CDF") +
  ggtitle("Relationship between \ncutoff and model fit (2-Step)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.75,0.25),
        legend.title = element_blank())


g5 <- ggplot(data = p_sims_ts_sva, aes(x = p, y = rmse, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  scale_y_continuous(breaks = my_breaks ) +
  # ylab(expression(sqrt(over(sum( (p[i] - hat(p)[i])^2, i==1, k),k)))) +
  ylab("RMSE") +
  xlab("CDF") +
  ggtitle("Relationship between \ncutoff and model fit (2-Step)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.5,0.23),
        legend.direction="horizontal",
        legend.title = element_blank())

g6 <- ggplot(data = p_sims_ts_sva, aes(x = p, y = log(mean_ratio), color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  # ylim(0, NA) +
  scale_y_continuous(breaks = my_breaks_r ) +
  geom_hline(yintercept = 0) +
  ylab(expression(log~Mean~Ratio)) + xlab("CDF") +
  ggtitle("Relationship between \ncutoff and model fit (2-Step)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.80,0.75),
        legend.title = element_blank())

png(file.path(".", "figs", "c2.png"), res = 200, width = 1920, height = 1080)
multiplot(g1, g4, g2, g5, g3, g6, cols=3)
dev.off()



k1 <- ggplot(data = p_sims_ts, aes(x = p, y = pearson, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ylim(.6, 1) +
  ylab("R") + xlab("CDF") +
  ggtitle("Relationship between \ncutoff and model fit (2-Step no sva)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.75,0.25),
        legend.title = element_blank())


k2 <- ggplot(data = p_sims_ts, aes(x = p, y = rmse, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  scale_y_continuous(breaks = my_breaks ) +
  # ylab(expression(sqrt(over(sum( (p[i] - hat(p)[i])^2, i==1, k),k)))) +
  ylab("RMSE") +
  xlab("CDF") +
  ggtitle("Relationship between \ncutoff and model fit (2-Step no sva)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.5,0.23),
        legend.direction="horizontal",
        legend.title = element_blank())

k3 <- ggplot(data = p_sims_ts, aes(x = p, y = log(mean_ratio), color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  # ylim(0, NA) +
  scale_y_continuous(breaks = my_breaks_r ) +
  geom_hline(yintercept = 0) +
  ylab(expression(log~Mean~Ratio)) + xlab("CDF") +
  ggtitle("Relationship between \ncutoff and model fit (2-Step no sva)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.65,0.55),
        legend.title = element_blank())

png(file.path(".", "figs", "sva_vs_nsva.png"), res = 200, width = 1920, height = 1080)
multiplot(k1, g4, k2, g5, k3, g6, cols=3)
dev.off()

l1 <- ggplot(data = p_sims_ts_with_sva_ave, aes(x = B, y = pearson, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  ylim(.6, 1) +
  ylab("R") + xlab("B") +
  ggtitle("Relationship between \ncutoff and model fit (2-Step)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.75,0.25),
        legend.title = element_blank())


l2 <- ggplot(data = p_sims_ts_with_sva_ave, aes(x = B, y = rmse, color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  scale_y_continuous(breaks = my_breaks ) +
  # ylab(expression(sqrt(over(sum( (p[i] - hat(p)[i])^2, i==1, k),k)))) +
  ylab("RMSE") +
  xlab("B") +
  ggtitle("Relationship between \ncutoff and model fit (2-Step)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.5,0.23),
        legend.direction="horizontal",
        legend.title = element_blank())

l3 <- ggplot(data = p_sims_ts_with_sva_ave, aes(x = B, y = log(mean_ratio), color = type, shape = type)) + 
  geom_point() + 
  geom_smooth(se = F) +
  scale_color_manual(values = cols) +
  transparent_legend + remove_grid +
  # ylim(0, NA) +
  scale_y_continuous(breaks = my_breaks_r ) +
  geom_hline(yintercept = 0) +
  ylab(expression(log~Mean~Ratio)) + xlab("B") +
  ggtitle("Relationship between \ncutoff and model fit (2-Step)") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold", size = 9),
        axis.text.y=element_text(size = 8, face="bold"),
        axis.text.x=element_text(size = 8, face="bold"),
        legend.position = c(0.65,0.35),
        legend.title = element_blank())

png(file.path(".", "figs", "c4.png"), res = 200, width = 1920, height = 1080)
multiplot(l1,l2, l3, cols=3)
dev.off()

