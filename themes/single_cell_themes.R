# Custom themes single cell

#custom themes
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(2.5, "cm")
)


mini_umap_theme <- list(
  #ggtitle(element_blank()),
  guides(x = axis, y = axis),
  labs(x = "UMAP1",y = "UMAP2"),
  theme(plot.title = element_text(hjust = 0.5, size = 5),
        axis.title = element_text(hjust = 0), 
        legend.text = element_text(size = 5),
        axis.title.x  = element_text(size = 5),
        axis.title.y = element_text(size = 5),
        legend.position.inside = c(0.015,0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.background = element_blank(),
        panel.background = element_blank(), #transparent panel bg
        plot.background = element_blank(), #transparent plot bg
        panel.border = element_blank(),
        axis.line = element_line(arrow = arrow(angle = 18,
                                               length = unit(0.5,"cm"), 
                                               type = "closed"),
                                 linewidth = 0.5)
  )
)


mini_tsne_theme <- list(
  #ggtitle(element_blank()),
  guides(x = axis, y = axis),
  labs(x = "tSNE1",y = "tSNE2"),
  theme(plot.title = element_text(hjust = 0.5, size = 5),
        axis.title = element_text(hjust = 0), 
        legend.text = element_text(size = 5),
        axis.title.x  = element_text(size = 5),
        axis.title.y = element_text(size = 5),
        legend.position.inside = c(0.015,0.5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.background = element_blank(),
        panel.background = element_blank(), #transparent panel bg
        plot.background = element_blank(), #transparent plot bg
        panel.border = element_blank(),
        axis.line = element_line(arrow = arrow(angle = 18,
                                               length = unit(0.5,"cm"), 
                                               type = "closed"),
                                 linewidth = 0.5)
  )
)

fp_theme <- list(
  guides(x = axis, y = axis),
  labs(x = "UMAP1",y = "UMAP2"),
  theme(text = element_text(size = 5),
    plot.title = element_text(hjust = 0.5, size = 5),
        axis.title = element_text(hjust = 0), 
        legend.text = element_text(size = 5),
        axis.title.x  = element_text(size = 5),
        axis.title.y = element_text(size = 5),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.background = element_blank(),
        panel.background = element_blank(), #transparent panel bg
        plot.background = element_blank(), #transparent plot bg
        panel.border = element_blank(),
        axis.line = element_line(arrow = arrow(angle = 18,
                                               length = unit(0.5,"cm"), 
                                               type = "closed"),
                                 linewidth = 0.5)
  )
)

barplot_theme <- list(
  labs(x = NULL, y = NULL, title = NULL),
  theme_classic(),
  theme(text = element_text(size = 7),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(vjust = 5), 
        axis.line.x.bottom = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  )
)

