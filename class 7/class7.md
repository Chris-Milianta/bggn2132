---
title: "Bioinformatics Class 7"
output: 
  html_document: 
    keep_md: yes

---



#functions again
we can source any file or F code with the `source`() function. 

```r
source("http://tinyurl.com/rescale-R")
```
lets make sure things are here


```r
ls()
```

```
##  [1] "both_na"         "both_na2"        "both_na3"       
##  [4] "df1"             "df2"             "df3"            
##  [7] "gene_intersect"  "gene_intersect2" "gene_intersect3"
## [10] "gene_intersect4" "rescale"         "rescale2"
```

check `rescale()` function is working


```r
rescale(1:10)
```

```
##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
##  [8] 0.7777778 0.8888889 1.0000000
```


break 

```r
rescale( c(1:10,"string"))
```

Lets check if `rescale2()` does better 


```r
rescale2( c(1:10,"string"))
```
## Function for finding missing values in two datasets

write a `both_na()` function to do this 

```r
x <-  c(1, 2, NA, 3, NA)
y <-  c(NA, 3, NA, 3, 4)
is.na(x)
```

```
## [1] FALSE FALSE  TRUE FALSE  TRUE
```

```r
which( is.na(x) )
```

```
## [1] 3 5
```

```r
sum( is.na(x) )
```

```
## [1] 2
```

using `&` operator and `sum()` to find how many things are missing in both experiments

```r
is.na(x) & is.na(y)
```

```
## [1] FALSE FALSE  TRUE FALSE FALSE
```

```r
sum(is.na(x) & is.na(y))
```

```
## [1] 1
```

now make function from this snippet. 


```r
both_na <- function(x,y) {
  sum(is.na(x) & is.na(y) ) 
}
```

"eejit proof" testing 

```r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)

#both_na(x,y2)
#gives error!
```
using `both_na2()` with != operator


```r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
both_na2(x,y2)
```



