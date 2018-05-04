
rescale <-  function(x) {
rng <-  range(x)
x <- (x-rng[1])/(rng[2]-rng[1])
return(x)
}



