ggplot(df) +
  geom_histogram(aes({aes}),{opts.binwidth}) +
  labs(x="{opts.log_x_open}{var_x}{opts.log_x_close}", y="Count", title="{opts.title}") +
  theme(title=element_text(size=20), 
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) {opts.add}