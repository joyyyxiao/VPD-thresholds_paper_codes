library(raster)
# ESA land cover map downloaded from https://maps.elie.ucl.ac.be/CCI/viewer/download.php
lc2000 <- raster("ESACCI/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2000-v2.0.7.tif")
lc2022 <- raster("ESACCI/C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.nc")

###########
lc2000[lc2000 == 0] <- NA
lc2022[lc2022 == 0] <- NA

## Code info
## 1 = Cropland / Cropland-mosaic
## 2 = Forest (all tree cover,includding flooded)
## 3 = Shrubland (including mosaic shrub/tree/herb)
## 4 = Grassland
## 5 = Sparse / Bare / Lichens & mosses
## 6 = Urban
## 7 = Water
## 8 = Snow/Ice

crop_codes   <- c(10,11,12,20,30,40)
forest_codes <- c(50,60,61,62,70,71,72,80,81,82,90,160,170)
shrub_codes  <- c(100,110,120,121,122,180)
grass_codes  <- c(130)
sparse_codes <- c(140,150,151,152,153,200,201,202)
urban_codes  <- c(190)
water_codes  <- c(210)
snow_codes   <- c(220)

lookup <- data.frame(
  code = c(crop_codes,
           forest_codes,
           shrub_codes,
           grass_codes,
           sparse_codes,
           urban_codes,
           water_codes,
           snow_codes),
  big  = c(rep(1, length(crop_codes)),
           rep(2, length(forest_codes)),
           rep(3, length(shrub_codes)),
           rep(4, length(grass_codes)),
           rep(5, length(sparse_codes)),
           rep(6, length(urban_codes)),
           rep(7, length(water_codes)),
           rep(8, length(snow_codes)))
)


lc2000_big <- subs(lc2000, lookup, by = "code", which = "big",
                   filename = "LC_bigClass_2000_300m.tif",
                   overwrite = TRUE)

lc2022_big <- subs(lc2022, lookup, by = "code", which = "big",
                   filename = "LC_bigClass_2022_300m.tif",
                   overwrite = TRUE)

change_0022 <- overlay(
  lc2000_big, lc2022_big,
  fun = function(a, b) {
    
    out <- ifelse(is.na(a) | is.na(b), NA,
                  ifelse(a != b, 1, 0))  
    return(out)
  },
  filename  = "LC_change_2000_2022_bigclass.tif",
  overwrite = TRUE
)


res_deg <- res(change_0022)        
fact    <- round(0.5 / res_deg[1])    

change_frac_0022_05 <- aggregate(
  change_0022,
  fact = fact,
  fun  = function(x, ...) mean(x, na.rm = TRUE),
  filename = "LC_change_frac_0022_05deg.tif",
  overwrite = TRUE
)


mask_50 <- change_frac_0022_05 >= 0.50
mask_40 <- change_frac_0022_05 >= 0.40
mask_30 <- change_frac_0022_05 >= 0.30
mask_20 <- change_frac_0022_05 >= 0.20

writeRaster(mask_20,"./ESACCI/LULCC/lulcc_0022_20perc.tif")
writeRaster(mask_30,"./ESACCI/LULCC/lulcc_0022_30perc.tif")
writeRaster(mask_40,"./ESACCI/LULCC/lulcc_0022_40perc.tif")
writeRaster(mask_50,"./ESACCI/LULCC/lulcc_0022_50perc.tif")
