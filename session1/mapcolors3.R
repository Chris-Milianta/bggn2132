map.colors3 <- function (x,
                         high.low = range 
                         ,palette) {
  # determine percent vlaues of the 'high.low' range
  percent <- ((x - high.low[1])/(high.low[2]-high.low[1]))
  
  #find corresponding index position in the color 'palette'
  # note catch for 0 percent calues to 1 
  index <- round ( (length(palette)-1) * percent ) +1
  return (palette[index])
} 