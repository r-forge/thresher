### get idioGPS.csv and save it as cytobandLocations
### make these into Rdata files
idiocolors <- c(acen='#aa0000', gneg='white', 
                gpos100='black', gpos25='gray75',
                gpos50='gray50', gpos75='gray25',
                gvar='darkblue', stalk='lightblue')

cytobandLocations <- read.csv("idioGPS.csv", row.names=1)
cytobandLocations$Chromosome <- factor(cytobandLocations$Chromosome,
                                       levels = paste("chr", c(1:22,"X", "Y"),
                                                      sep=""))
cytobandLocations <- cytobandLocations[order(cytobandLocations$Chromosome,
                                             cytobandLocations$loc.start),]

save(cytobandLocations, idiocolors,
     file = file.path("..", "..", "data", "cytobandLocations.rda"))
