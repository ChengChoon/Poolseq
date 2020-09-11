
## define a data set as gwscan contained: (scaffold, pos, snp.FST )

load librarys= ggplot2, cowplot

plot.gwscan <- function (gwscan, size = 0.2) {

  # Add a column with the marker index.
  n      <- nrow(gwscan)
  gwscan <- cbind(gwscan,marker = 1:n)

  # Convert the p-values to the -log10 scale.
  #gwscan <- transform(gwscan,snp.FST = -log10(snp.FST))

  # Add column "odd.scaffold" to the table, and find the positions of the
  # scaffoldomosomes along the x-axis.
  gwscan <- transform(gwscan,odd.scaffold = (scaffold %% 2) == 1)
  x.scaffold  <- tapply(gwscan$marker,gwscan$scaffold,mean)

  # Create the genome-wide scan ("Manhattan plot").
  return(ggplot(gwscan,aes(x = marker,y = snp.FST,color = odd.scaffold)) +
           geom_point(size = size,shape = 1) +
           scale_x_continuous(breaks = x.scaffold) +
           scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
           scale_color_manual(values = c("skyblue","darkblue"),guide = "none") +
           labs(x = "Scaffolds",y = "SNP specific FST") +
           theme_cowplot(font_size = 10) +
           theme(axis.line = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.text.x = element_text(size = 3.5)))
}

