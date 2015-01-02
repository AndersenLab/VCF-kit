ggplot(df) +
  geom_point(aes({aes}),{opts.geom}) +
  labs(x="{xlab}", y="{ylab}", title="{opts.title}") +
  theme(title=element_text(size=20), 
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) {opts.add}