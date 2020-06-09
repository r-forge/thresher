### this may work to hold the results from one file
setClass("LGF1",
         slots = c(
           raw = "data.frame",   # rows = clones, columns = LGF-cytobands
           frequency = "data.frame",   # rows = cytobands, columns = L,G,F
           source = "character", # file name
           name = "character",   # user-defined data set name
           CL = "data.frame"     # cytoband locations
         ))
## Should the LGF frequency slot be a matrix instead of a data frem?

### Suppose we read multiple JSON files. Do we create a list of LGF
### objects, or do we instead create a classwhere the entries themselves
### are lists instead of data frames? Part of the sticking point is the
### CL slot, which must be the same for all LGF objects in a set.

### Second try. Here we make things into lists.
setClass("LGF2",
         slots = c(
           raw = "list",  # of data frame, rows = clones, columns = LGF-cytobands
           frequency = "data.frame", # rows = cytobands, columns = L,G,F (repeat)
           source = "character", # vector of file names
           CL = "data.frame"     # cytoband locations
         ))
### Note that the "names" slot goes away, ### since we can assign names to
### the elements in a list.
### Note also that if
###    all( names(raw) == c("ALL", "AML", "CLL")
### then 
###    all(colnames(frequency) == c("ALL.Loss", "ALL.Gain", "ALL.Fusion"
###                                 "AML.Loss", "AML.Gain", "AML.Fusion",
###                                 "CLL.Loss", "CLL.Gain", "CLL.Fusion") )
