readLeafTrace <- function(file) {
  x <- read.table(file, stringsAsFactors=FALSE, header=FALSE)
  names(x)[1:3] <- c("chrom", "chromStart", "chromEnd")
  x
}

# Read the lookup table
readLookupTable <- function(file) {
  lookup <- read.table(file, stringsAsFactors=FALSE, header=TRUE)
  names(lookup) <- c("leafName", "value", "include")
  return(lookup)
}

# Map the leaf values to colors
getLeafColors <- function(lookup, color_map=NULL) {
  # Create a color palette (you can use any method to generate colors, here we'll use a simple scale)
  if (is.null(color_map)) {
    #color_map <- c("AZ" = "red", "AZ1" = "orange", "CAR" = "blue", "CS" = "limegreen", "Pho" = "darkgray", "Tm" = "brown", "Tu" ="deeppink", "unk" = "cyan") 
    color_map <- c("Fb"="red", "Flou" = "orange", "Fger"= "limegreen", "NA1" = "blue", "NA2"="deeppink", "NA3"="cyan")
    }
  
  # Match leaf names in lookup with the color map
  leaf_colors <- setNames(color_map[lookup$value], lookup$leafName)
  
  return(leaf_colors)
}

# Modify the plotLeafTrace function to use the color map
plotLeafTrace <- function(file=NULL, x=NULL, col="black", xlim=NULL, ylim=NULL, xlab="Coordinate", ylab="", add=FALSE, subset=NULL, lwd=1, lookup=NULL, color_map=NULL, ...) {
  if (is.null(file) && is.null(x)) {
    stop("Need to provide either a layout file (file) or a layout table (x)")
  }
  if (!is.null(file))
    x <- readLeafTrace(file)
  
  # If lookup table is provided, get the corresponding colors
  if (!is.null(lookup)) {
    leaf_colors <- getLeafColors(lookup, color_map)
  } else {
    leaf_colors <- col  # Default to a single color
  }
  
  if (is.null(xlim))
    xlim <- range(c(x$chromStart, x$chromEnd))
  
  if (!is.null(subset))
    subset <- c(subset, paste0(subset, "_1"), paste0(subset, "_2"))
  
  nameCols <- seq(from=4, to=ncol(x), by=2)
  yCols <- nameCols+1
  if (is.null(ylim))
    ylim <- range(x[,yCols])
  
  if (!add)
    plot(c(0),c(0), type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
  
  for (i in nameCols) {
    ind <- unique(x[,i])
    if (length(ind) != 1) stop("got different inds in same col")
    if (!is.null(subset) && !is.element(ind, subset)) next
    
    # Use the color from the lookup table based on leaf name
    color <- as.character(leaf_colors[ind])
    
    segments(x0=x$chromStart, x1=x$chromEnd, y0=x[,i+1], col=color, lwd=lwd)
  }
  
  invisible(NULL)
}

# Read in the leaf trace data
x <- readLeafTrace("~/FgramARGchr3_1-1000000.txt_145118-999584.199.layout")

# Read in the lookup table
lookup <- readLookupTable("~/FgramARGstrainsPlus.txt")

# Plot with colors from the lookup table
plotLeafTrace(x=x, col="black", xlim=c(250000, 1000000), lookup=lookup, xlab="Coordinate", ylab="Leaf Trace")


