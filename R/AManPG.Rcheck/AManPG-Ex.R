pkgname <- "AManPG"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('AManPG')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("hello")
### * hello

flush(stderr()); flush(stdout())

### Name: hello
### Title: Hello, World!
### Aliases: hello

### ** Examples

hello()



cleanEx()
nameEx("prox.l1")
### * prox.l1

flush(stderr()); flush(stdout())

### Name: prox_l1
### Title: Proximal L1 Mapping
### Aliases: prox_l1

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }



cleanEx()
nameEx("spca.amanpg")
### * spca.amanpg

flush(stderr()); flush(stdout())

### Name: spca_amanpg
### Title: Alternativing Manifold Proximal Gradient Method for Sparse PCA
### Aliases: spca_amanpg

### ** Examples

##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }



cleanEx()
nameEx("svd.econ")
### * svd.econ

flush(stderr()); flush(stdout())

### Name: svd.econ
### Title: Economy-Size Singular Value Decomposition
### Aliases: svd.econ

### ** Examples

set.seed(10)
x <- matrix(rnorm(3 * 4), 3, 4)
udv <- svd.econ(x)

udv.0 <- svd.econ(x, 0)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
