library(testthat)
## R 4.3.2 on some machines: waldo only loads if glue is already loaded
if (requireNamespace("glue", quietly = TRUE)) loadNamespace("glue")
if (requireNamespace("waldo", quietly = TRUE)) loadNamespace("waldo")
library(impliedIV)
test_check("impliedIV")
