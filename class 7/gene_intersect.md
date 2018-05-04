---
title: "Finding intersection"
output: 
  html_document: 
    keep_md: yes
---



# a function that is actually useful. 

```r
source("http://tinyurl.com/rescale-R")

x <- df1$IDs
y <- df2$IDs
x
```

```
## [1] "gene1" "gene2" "gene3"
```

```r
y
```

```
## [1] "gene2" "gene4" "gene3" "gene5"
```

```r
## finding intersect

intersect(x,y)
```

```
## [1] "gene2" "gene3"
```

```r
## or `%in%` works
x %in% y
```

```
## [1] FALSE  TRUE  TRUE
```

we can use the logical output in `%in%` to get at our data 

```r
##will output gene names
x[x %in% y]
```

```
## [1] "gene2" "gene3"
```

```r
y[y %in% x]
```

```
## [1] "gene2" "gene3"
```

lets put together as columns of a matrix


```r
cbind( x[x %in% y], y[y %in% x] )
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```

turn into function

```r
gene_intersect <- function(x, y) { 
   cbind( x[ x %in% y ], y[ y %in% x ] )
}
## test it 
gene_intersect(x,y)
```

```
##      [,1]    [,2]   
## [1,] "gene2" "gene2"
## [2,] "gene3" "gene3"
```

trying with data.frame input rather than vectors


```r
gene_intersect2(df1,df2)
```

```
##     IDs exp df2[df2$IDs %in% df1$IDs, "exp"]
## 2 gene2   1                               -2
## 3 gene3   1                                1
```

looks good but its skateboard not a car. 
Adding some flexibility for col name to match by 

```r
gene_intersect3(df1,df2)
```

```
##     IDs exp exp2
## 2 gene2   1   -2
## 3 gene3   1    1
```

making it more human friendly 

```r
(gene_intersect4(df1,df2))
```

```
##     IDs exp exp2
## 2 gene2   1   -2
## 3 gene3   1    1
```

`merge()` also works and is built in. 

```r
merge(df1, df2, by="IDs")
```

```
##     IDs exp.x exp.y
## 1 gene2     1    -2
## 2 gene3     1     1
```

