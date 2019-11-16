
# Load Libs ---------------------------------------------------------------

packs <- c("data.table", "dplyr", "SummarizedExperiment", "recount", "genefilter", "RColorBrewer", 
           "mixtools","matrixStats", "MuSiC", "convert", "xbioc", "ggplot2", "sva", "plotly",
           "readxl", "sp", "rgdal")
libs_loaded <- sapply(packs, library, character.only = T)



# Load Data ---------------------------------------------------------------

dapi <- read_xlsx("stats_byChannel_RNAscope.xlsx", sheet = "DAPI") %>% data.table
ubc <- read_xlsx("stats_byChannel_RNAscope.xlsx", sheet = "UBC") %>% data.table
snap25 <- read_xlsx("stats_byChannel_RNAscope.xlsx", sheet = "SNAP25") %>% data.table



# quick plot --------------------------------------------------------------


# mybb <- cbind(x=c(363498.5, 480497.5, 480497.5, 363498.5), y=c(5894630, 5894630, 5806524, 5806524))
mybb <- cbind(x=c(363498.5, 480497.5, 480497.5, 363498.5), y=c(5894630, 5894630, 5806524, 5806524))
crs <-CRS("+proj=utm +zone=12 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
mybb <- SpatialPolygons(list(Polygons(list(Polygon(mybb)),"1")), proj4string=crs)
plot(mybb)


plot(dapi[name == "BR5180_Slide 15_Section A_63x z-stack_gain 750_405 to 3_488 to 40_561 to pt2_637 to 40_UBC-690_SNAP25-570_PPIB-520_Image 10_Linear unmixing", .(Centroid_1,Centroid_2 )])
plot(dapi[name == "BR5180_Slide 15_Section A_63x z-stack_gain 750_405 to 3_488 to 40_561 to pt2_637 to 40_UBC-690_SNAP25-570_PPIB-520_Image 10_Linear unmixing", .(BoundingBox_1,BoundingBox_2 )])

summary(dapi$Area)
summary(snap25$Area)

unique(dapi$name)
unique(snap25$name)


tapply(snap25$Area, snap25$name, var)
tapply(dapi$Area, dapi$name, var)