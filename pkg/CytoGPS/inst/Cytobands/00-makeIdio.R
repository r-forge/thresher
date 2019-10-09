### get idioGPSRaw.csv and save it as idio
### make these into Rdata files
idiocolors <- c(acen='#aa0000', gneg='white', 
                gpos100='black', gpos25='gray75',
                gpos50='gray50', gpos75='gray25',
                gvar='darkblue', stalk='lightblue')

cytobandLocations <- read.csv("idioGPS.csv", row.names=1)

save(cytobandLocations, idiocolors,
     file = file.path("..", "..", "data", "idiogram.rda"))
