bands <- read.table("ChrBands.txt", as.is = TRUE)[,1]

partial <- strsplit(bands, "(p|q)")
chr <- factor(unlist(sapply(partial, function(x) x[[1]])),
              levels = c(1:22, "X", "Y"))
summary(chr)

nc <- nchar(as.character(chr))
arm <- factor(substring(bands, 1, nc + 1))

type <- rep(c("Loss", "Gain", "Fusion"), 916)
label <- paste(type, rep(bands, 3), sep = "_")

length(label)
length(type)

lgfFeatures <- data.frame(Label = label,
                          Type = type,
                          Band = rep(bands, each = 3),
                          Chr = rep(chr, each = 3),
                          Arm = rep(arm, each = 3),
                          Index = 1:(3*916))
rownames(lgfFeatures) <- as.character(lgfFeatures$Label)
rm(arm, bands, chr, label, nc, partial, type)
save(lgfFeatures, file = file.path("..", "..", "data", "lgfFeatures.rda"))
