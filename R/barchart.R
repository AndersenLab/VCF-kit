ggplot(df) +
  geom_bar(aes(x=df${var1}),{opts.binwidth}) +
  labs(x="{var1}", y="Count", title="{opts.title}") +
  theme(title=element_text(size=20), 
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) {opts.add}