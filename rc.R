#' Convert a ratified raster into polygons by finding connected clumps by value
#'
#' Clumps with values equal to each value are detected and returned as individual polygons. Current only works for discrete
#' values so classify your raster first if using continuous data. 
#' 
#' Input raster must have a raster attribute table, and column "type" is assumed as the identifying label. 
#' @param r RasterLayer with levels
#' @param vals set of values to restrict polygonization to (optional)
#' @param typeColumn character name of type column, defaults to "type"
#' @param patchColumn character name of patch column, defaults to "patch"
#' @param factorV character name of factorValue column
#' @param directions integer 8 (Queen) or 4 (rook) see \code{\link{clump}}
#' 
#' @return SpatialPolygonsDataFrame with columns for type and patch
#' @export
#'
ratToPolygons <- function(r, vals = NULL, typeColumn = "type", patchColumn = "patch", factorV = "type", directions = 8) {
  if (!is.null(vals)) {
    uval <- vals 
  }   else  {
    uval <- sort(na.omit(unique(values(r))))
  }
  l <- vector("list", length(uval))
  for (i in seq_along(l)) {
    cl <- clump(r == uval[i])
    names(cl) <- patchColumn
    l[[i]] <- rasterToPolygons(cl, dissolve = TRUE)
    l[[i]][[typeColumn]] <- factorValues(r, uval[i])[, factorV]
    print(sprintf("%i of %i", i, length(l)))
  }
  poly <- l[[1]]
  for (i in seq_along(l)[-1]) {
    x <- spChFIDs(l[[i]], as.character(seq(nrow(l[[i]])) + max(as.integer(row.names(poly)))))
    poly <- spRbind(poly, x)
  } 
  ## makes patches unique
  poly[[patchColumn]] <- seq(nrow(poly))
  poly
}

#' Ratify a raster by prepending a prefix to pixel values
#'
#' The value understood by factorValues can be specified. 
#' @param r a RasterLayer, should contain only discrete values
#' @param prefix character string to prepend to pixel values
#' @param factorV name of factorValue 
#'
#' @return ratified raster
#' @export
#'
ratifyIndex <- function(r, lab = NULL, prefix = "prefix_", factorV = "type") {
  levs <- levels(r)[[1]]
  if (is.null(lab)) lab <- levs$ID
  levs[[factorV]] <- paste0(prefix, lab)
  levels(r) <- levs
  r
}

## -- smoothing applied to topography
## set to 1 for no smoothing
reductionfactor <- 14
## smooth nxn kernel applied to get *nice* contour boundaries
toposmooth <- TRUE
smoothingwindow <- 9
smoothingmodel <- "Gauss"
## --------------------


## -- intervals applied to classify the raster
intervals <- c(-4000, -3500, -2500, -2000, -1500, -1000, -500, 0)
## --------------------
library(raster)
library(rgdal)
library(maptools)
topo <- raster(file.path("data", "ETOPO2v2c_f4.nc"))
if (reductionfactor > 1) {
  topo <- aggregate(topo, fact = reductionfactor, fun = mean, na.rm = TRUE)
}
if (toposmooth) {
  ## smoothwindow is the neighbourhood of pixels to smooth with  
  gf <- focalWeight(topo, res(topo) * smoothingwindow, smoothingmodel) 
  ## focal applys the gf  kernel to apply a weighting to each pixel from its neighbourhood bathy 
  ## padding was not done in the draft FSA report 2015-09-24
  #bathy <- focal(bathy1, gf, na.rm = TRUE) ## pad = TRUE, padValue = 0)
  topo <- focal(topo, gf, na.rm = TRUE, pad = TRUE, padValue = 0)
}

dd <- c( min(c(floor(cellStats(topo, min)), head(intervals, 1))), tail(intervals, -1))
topo <- clamp(topo, lower = min(dd), upper = max(dd)-1)
topo_def <- data.frame(from = dd[-length(dd)], to = dd[-1], seq(0, length(dd)-2))
topo_rc <- reclassify(topo, topo_def)

topo_rat0 <- ratify(topo_rc)

topo_rat <- ratifyIndex(topo_rat0, prefix = "", lab = sprintf("<%gm", abs(dd)))
bathyPal <- mkPal(setNames(bathycols, levels(bathy_rat)[[1]]$type))
bathy_poly <- ratToPolygons(bathy_rat)



# ---- reporting-a ----


#' Convert  colours to grey scale 
#'
#' @param x character colour
#'
#' @return grey scale colours
#' @export
#'
togrey <- function(x) {
  grey(apply(col2rgb(x)/256, 2, mean))
}


#' Add vertices to a linear path to ensure curvature 
#'
#' This function takes a matrix of coordinates in longlat that define a pathy, and returns a denser matrix with more vertices based on a minimum allowed distance. 
#' 
#' The function assumes each segment is a great circle, on the Vincenty ellipsoid \code{\link[raster]{distVincentyEllipsoid}}. 
#' @param x matrix of coordinates, long/lat
#' @param minm default distance in metres on a great cicle
#'
#' @return matrix
#' @export
#'
#' @importFrom raster distVincentyEllipsoid
densifyMat <- function(x, minm = 1852 * 10) {
  ## assume this is the closed polygon
  l <- vector("list", nrow(x) - 1)
  for (i in seq_along(l)) {
    dist <- distVincentyEllipsoid(x[i, ], x[i + 1, ])
    nsegs <- max(c(dist %/% minm, 3))
    l[[i]] <- cbind(seq(x[i,1], x[i+1,1], length = nsegs), seq(x[i,2], x[i+1,2], length = nsegs))
  }
  do.call(rbind, l)
}


#' Select the line with the most vertices from a multi-line Lines object. 
#'
#' This function assumes a single-row SpatialLinesDataFrame, and drops all but the most complex line part. 
#' @param x single row SpatialLinesDataFrame
#'
#' @return SpatialLinesDataFrame
#' @export
#'
keepOnlyMostComplexLine <- function(x) {
  for (iObj in seq_len(nrow(x))) {
    if (inherits(x, "SpatialLinesDataFrame")) {
      wmax <- which.max(sapply(x[iObj, ]@lines[[1]]@Lines, function(x)
        nrow(x@coords)))
      x@lines[[iObj]]@Lines <- x@lines[[iObj]]@Lines[wmax]
    }
  }
  x
}


#' Identify point-in-triangle by conversion to polygons
#'
#' @param tri list P n*2 coordinates and T matrix of n*3 indices defining triangles
#' @param pts 
#' @importFrom sp Polygon Polygons SpatialPolygons CRS proj4string over
tri_pip <- function(tri, pts) {
  ps <- lapply(split(tri$T, seq(nrow(tri$T))), function(x) Polygon(tri$P[c(x, x[1]), ]))
  sp <- lapply(seq_along(ps), function(x) Polygons(ps[x], x))
  spp <- SpatialPolygons(sp, proj4string = CRS(proj4string(pts)))
  over(pts, spp)
}


#' Convert a ratified raster into polygons by finding connected clumps by value
#'
#' Clumps with values equal to each value are detected and returned as individual polygons. Current only works for discrete
#' values so classify your raster first if using continuous data. 
#' 
#' Input raster must have a raster attribute table, and column "type" is assumed as the identifying label. 
#' @param r RasterLayer with levels
#' @param vals set of values to restrict polygonization to (optional)
#' @param typeColumn character name of type column, defaults to "type"
#' @param patchColumn character name of patch column, defaults to "patch"
#' @param factorV character name of factorValue column
#' @param directions integer 8 (Queen) or 4 (rook) see \code{\link{clump}}
#' 
#' @return SpatialPolygonsDataFrame with columns for type and patch
#' @export
#'
ratToPolygons <- function(r, vals = NULL, typeColumn = "type", patchColumn = "patch", factorV = "type", directions = 8) {
  if (!is.null(vals)) {
    uval <- vals 
  }   else  {
    uval <- sort(na.omit(unique(values(r))))
  }
  l <- vector("list", length(uval))
  for (i in seq_along(l)) {
    cl <- clump(r == uval[i])
    names(cl) <- patchColumn
    l[[i]] <- rasterToPolygons(cl, dissolve = TRUE)
    l[[i]][[typeColumn]] <- factorValues(r, uval[i])[, factorV]
   print(sprintf("%i of %i", i, length(l)))
  }
  poly <- l[[1]]
  for (i in seq_along(l)[-1]) {
    x <- spChFIDs(l[[i]], as.character(seq(nrow(l[[i]])) + max(as.integer(row.names(poly)))))
    poly <- spRbind(poly, x)
  } 
  ## makes patches unique
  poly[[patchColumn]] <- seq(nrow(poly))
  poly
}

#' Ratify a raster by prepending a prefix to pixel values
#'
#' The value understood by factorValues can be specified. 
#' @param r a RasterLayer, should contain only discrete values
#' @param prefix character string to prepend to pixel values
#' @param factorV name of factorValue 
#'
#' @return ratified raster
#' @export
#'
ratifyIndex <- function(r, lab = NULL, prefix = "prefix_", factorV = "type") {
  levs <- levels(r)[[1]]
  if (is.null(lab)) lab <- levs$ID
  levs[[factorV]] <- paste0(prefix, lab)
  levels(r) <- levs
  r
}
#' Colour palette for Southern Ocean Bioregionalisation
#'
#' There are 20 colours, this function returns them all or matches values 1:20 to the right colour. 
#' @param x (optional) values for which colours are required
#' @param alpha optional opacity value (not yet implemented)
#' @return hex character string
#' @export
#'
#' @examples
pelregPal <- function(x, alpha = 1) {
  cc <- c("#BF9830", "#A64A00", "#E64400", "#EC6059", "#A60400", "#FF010A",
          "#FFC790", "#9ECAE1", "#6BAED6", "#2171B5", "#08306B", "#EEC6EF",
          "#B40097", "#750062", "#C0F400", "#9AB72E", "#00C90D", "#008209",
          "#FFBA00", "#FFFD00")
  
  if (alpha < 1) warning("alpha not yet implemented")
  if (missing(x)) cc else cc[x]
}

mkPal <- function(namedcols) {
  cols <- namedcols
  function(x) {
    cols[x]
  }
} 

  
#' Convert a raster into polygons by finding connected clumps by value
#'
#' Clumps with values equal to each value are detected and returned as individual polygons. Current only works for discrete
#' values so classify your raster first if using continuous data. 
#' @param r RasterLayer
#' @param vals set of values
#'
#' @return SpatialPolygonsDataFrame with columns for value ID and clump within value
#' @export
#'
#' @examples
ras2poly <- function(r, vals, valColumn = "val", clumpColumn = "clumpID") {
  uval <- sort(na.omit(unique(values(r))))
  l <- vector("list", length(uval))
  for (i in seq_along(l)) {
    cl <- clump(r == uval[i])
    names(cl) <- clumpColumn
    l[[i]] <- rasterToPolygons(cl, dissolve = TRUE)
    l[[i]][[valColumn]] <- uval[i]
  }
  poly <- c_poly(l) 
  poly
}

#' Combine a list of PolygonDataFrame objects into a single object. 
#' 
#' Objects are combined by modifying sequential IDs to increase. 
c_poly <- function(x) {
  poly <- x[[1]]
  if (length(x) == 1L) {
    warning("input of length one, returning first element")
    return(poly)
  }
  for (i in seq_along(x)[-1]) {
    x0 <- spChFIDs(x[[i]], as.character(seq(nrow(x[[i]])) + max(as.integer(row.names(poly)))))
    poly <- maptools::spRbind(poly, x0)
  }
  poly
}
#' Extract the longest bounded polygon ring. 
largestbound <- function(x, id = "1") {
  ## expecting a single row polygons object
  x1 <- x@polygons[[1]]@Polygons[[which.max(sapply(x@polygons[[1]]@Polygons, function(x) nrow(x@coords)))]]
  SpatialPolygons(list(Polygons(list(x1), id)), proj4string = CRS(proj4string(x)))
}
