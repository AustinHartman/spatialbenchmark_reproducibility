# This file contains load functions for loading in datasets from sources not supported by Seurat
# 
# LoadMERFISH
# LoadResolve
# LoadOsmFISH
# LoadSTARmapPLUS
# LoadEELFISH
#
# 10x Xenium, Vizgen MERSCOPE, and Nanostring CosMx SMI data are also analyzed
# in the associated manuscript and can be loaded in Seurat using the LoadXenium,
# LoadMERSCOPE, and LoadNanostring functions within Seurat.
#

ReadMERFISH <- function(
  counts.h5ad.file,
  segmentations.file,
  molecules.file
) {
  # load up and format the segmentations matrix for loading into a FOV object
  # also generate some centroids from the segs
  segs.raw <- read.csv(segmentations.file, numerals = "no.loss")
  segs.raw$X <- as.character(segs.raw$X)
  cell.segs.list <- list()
  cell.cents.list <- list()
  for (r in 1:nrow(segs.raw)) {
    row <- segs.raw[r,]
    x <- as.numeric(strsplit(row$boundaryX, ", ")[[1]])
    y <- as.numeric(strsplit(row$boundaryY, ", ")[[1]])
    cell.segs <- data.frame("cell" = rep(row$X, length(x)), "x" = x, "y" = y)
    cell.cents <- data.frame("cell" = row$X, "x" = mean(x), "y" = mean(y))
    cell.segs.list[[row$X]] <- cell.segs
    cell.cents.list[[row$X]] <- cell.cents
  }
  segs <- bind_rows(cell.segs.list)
  cents <- bind_rows(cell.cents.list)
  rm(cell.segs.list)
  rm(segs.raw)
  
  # load in the counts - for some reason I can't find the counts for lots of the cells in the dataset
  file.h5 <- H5File$new(counts.h5ad.file, mode = "r")
  cellnames <- file.h5[["obs"]][]
  cellnames$idxs <- rownames(cellnames)
  colnames(cellnames) <- c("cellname", "idx")
  cell.idxs <- cellnames[cellnames$cellname %in% unique(cents$cell),]$idx
  cell.names <- cellnames[cellnames$cellname %in% unique(cents$cell),]$cellname
  counts <- file.h5[["X"]][, as.numeric(cell.idxs)]
  colnames(counts) <- cell.names
  rownames(counts) <- file.h5[["var"]][]$index
  
  # load in the moleculese
  mols <- read.csv(molecules.file)
  mols$X <- NULL
  mols$global_z <- NULL
  colnames(mols) <- c("x", "y", "gene")
  
  return(list("cents" = cents, "segs" = segs, "counts" = counts, "mols" = mols))
}

#### Read in the Resolve Mouse Brain Dataset
#' Read Resolve Data
#'
#' Read in ...
#'
#' @inheritParams ReadVitessce
#' @param data.dir Path to directory with Resolve data; requires at least one
#' file present with the following extensions:
#' \itemize{
#'  \item \dQuote{\code{expression_matrix.tsv}}: used for reading count matrix
#'  \item \dQuote{\code{raw_locations.txt}}: used for reading molecule spatial
#'  coordinate matrix
#' }
#' Files in \code{data.dir} (eg. \dQuote{\code{Resolve_Data/}}) should either be
#' named as extension above
#' (eg. \dQuote{\code{Resolve_Data/expression_matrix.tsv}}) or as
#' \dQuote{\code{basename(data.dir))_extension}}
#' (eg. \dQuote{\code{Resolve_Data/ResolveData_expression_matrix.tsv}})
#' @param roi.dir Path to directory with ROI files; requires at least on ROI
#' file that matches the pattern \dQuote{\code{^Cell_\\\\d+\\\.$}} present; pass
#' \code{NULL} to skip
#'
#' @return A list with some combination of the following value
#' \itemize{
#'  \item \dQuote{\code{transcripts}}: ...
#'  \item \dQuote{\code{molecules}}: ...
#'  \item \dQuote{\code{...}}: ...
#' }
#'
#' @export
#' @keywords internal
#'
#' @concept input
#'
#' @note This function requires the
#' \href{https://CRAN.R-project.org/package=RImageJROI}{RImageJROI} and
#' \href{https://cran.r-project.org/package=data.table}{\pkg{data.table}}
#' packages to be installed
#'
ReadResolve <- function(
  data.dir,
  roi.dir = file.path(data.dir, 'roi'),
  filter = NA_character_,
  verbose = TRUE
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install 'data.table' for this function")
  }
  imagej <- requireNamespace("RImageJROI", quietly = TRUE)
  roi.pattern <- 'Cell_\\d+\\.roi$'
  my.lapply <- ifelse(test = isTRUE(x = verbose), yes = pblapply, no = lapply)
  if (!is.null(x = data.dir) && !dir.exists(paths = data.dir)) {
    stop("Cannot find input data directory")
  }
  files <- vapply(
    X = c(
      transcripts = 'expression_matrix.tsv',
      molecules = 'raw_locations.txt'
    ),
    FUN = function(ext, bname) {
      if (is.null(x = data.dir)) {
        return(NA_character_)
      }
      fname <- file.path(data.dir, paste(bname, ext, sep = '_'))
      if (!file.exists(fname)) {
        fname <- file.path(data.dir, ext)
      }
      return(fname)
    },
    FUN.VALUE = character(length = 1L),
    bname = basename(path = data.dir),
    USE.NAMES = TRUE
  )
  files[!file.exists(files)] <- NA_character_
  if (!is.null(x = roi.dir)) {
    if (isFALSE(x = imagej)) {
      warning(
        "Cannot find RImageJROI; unable to load polygon coordiantes",
        immediate. = TRUE
      )
    } else if (!dir.exists(paths = roi.dir)) {
      warning("Cannot find directory with ROI files", immediate. = TRUE)
    } else if (!length(x = list.files(path = roi.dir, pattern = roi.pattern))) {
      warning("Cannot find ROI files in ", roi.dir, immediate. = TRUE)
    } else {
      files[['segmentations']] <- roi.dir
    }
  }
  if (all(is.na(x = files))) {
    stop("Cannot find Resolve input files in ", data.dir)
  }
  files <- files[!is.na(x = files)]
  outs <- vector(mode = 'list', length = length(x = files))
  names(x = outs) <- names(x = files)
  for (otype in names(x = outs)) {
    outs[[otype]] <- switch(
      EXPR = otype,
      'transcripts' = {
        if (isTRUE(x = verbose)) {
          message("Reading transcripts")
        }
        cts <- data.table::fread(
          file = files[[otype]],
          sep = '\t',
          header = TRUE,
          data.table = FALSE
        )
        rownames(x = cts) <- cts[, 'GeneID']
        cts <- cts[, -which(x = colnames(x = cts) == 'GeneID')]
        if (!is.na(x = filter)) {
          cts <- cts[!grepl(pattern = filter, x = rownames(x = cts)), , drop = FALSE]
        }
        idx <- unlist(x = lapply(
          X = seq_len(length.out = ncol(x = cts)),
          FUN = function(i) {
            if (!any(is.na(x = cts[, i]))) {
              return(i)
            }
            return(NULL)
          }
        ))
        cts <- as.matrix(x = cts[, idx, drop = FALSE])
        if ((sum(cts == 0) / length(x = cts)) > 0.4) {
          cts <- as.sparse(x = cts)
        }
        cts
      },
      'segmentations' = {
        if (isTRUE(x = verbose)) {
          message("Reading segmentations")
        }
        segs <- my.lapply(
          X = list.files(
            path = files[[otype]],
            pattern = roi.pattern,
            full.names = TRUE
          ),
          FUN = function(roi) {
            dat <- RImageJROI::read.ijroi(file = roi)
            df <- as.data.frame(x = dat$coords)
            df$cell <- dat$name
            return(df)
          }
        )
        do.call(what = 'rbind', segs)
      },
      'molecules' = {
        if (isTRUE(x = verbose)) {
          message("Reading molecules")
        }
        mols <- data.table::fread(
          file = files[[otype]],
          sep = '\t',
          header = FALSE,
          data.table = FALSE
        )
        colnames(x = mols) <- c('x', 'y', 'UNKNOWN', 'gene')
        mols
      },
      stop("Unknown Resolve input type: ", otype)
    )
  }
  return(outs)
}

LoadMERFISH <- function(
  data.dir,
  assay="RNA",
  fov="fov"
) {
  dfs <- ReadMERFISH(data.dir$counts, data.dir$segs, data.dir$mols)
  centroids.obj <- CreateCentroids(dfs$cents)
  fov.obj <- CreateFOV(
    coords = list("centroids" = centroids.obj),
    type = c("centroids")
  )
  obj <- CreateSeuratObject(counts = dfs$counts, assay = assay)
  fov.obj <- subset(fov.obj, cells = Cells(obj))
  obj[[fov]] <- fov.obj
  return(obj)
}

LoadEELFISH <- function(h5.path, assay="RNA", fov="fov") {
  h5 <- H5File$new(filename = h5.path, mode = 'r')
  print(h5[["X"]])
  mtx <- as.sparse(h5[["X"]][,])
  x.coords <- h5[["obs"]][["X_um"]][]
  y.coords <- h5[["obs"]][["Y_um"]][]
  colnames(mtx) <- seq(1:length(x.coords))
  rownames(mtx) <- h5[["var"]][["_index"]][]
  centroids.df <- data.frame(x=x.coords, y = y.coords, cell = seq(1:length(x.coords)))
  cents <- CreateCentroids(centroids.df)
  segmentations.data <- list(
    "centroids" = cents
  )
  coords <- CreateFOV(
    coords = segmentations.data,
    type = c("centroids")
  )
  object <- CreateSeuratObject(counts = mtx, assay = assay)
  object[[fov]] <- coords
  return(object)
}

LoadResolve <- function(data.dir, assay = "RNA", fov = "fov") {
  resolve.coords <- ReadResolve(data.dir)
  molecules.obj <- CreateMolecules(resolve.coords$molecules)
  segmentations.obj <- CreateSegmentation(resolve.coords$segmentations)
  centroids.obj <- as.Centroids(segmentations.obj, radius = 1)
  fov.obj <- CreateFOV(
    coords = list("segmentation" = segmentations.obj, "centroids" = centroids.obj),
    type = c("segmentation", "centroids"),
    molecules = molecules.obj
  )
  resolve.obj <- CreateSeuratObject(resolve.coords$transcripts, assay = assay)
  resolve.obj[[fov]] <- fov.obj
  return(resolve.obj)
}

LoadOsmFISH <- function(
  data.dir,
  assay = "RNA",
  fov = "fov"
) {
  osmFISH.coords <- ReadVitessce(
    counts = data.dir$counts,
    coords = data.dir$coords,
    molecules = data.dir$molecules,
    type = "segmentations")
  molecules.obj <- CreateMolecules(osmFISH.coords$molecules)
  segmentations.obj <- CreateSegmentation(osmFISH.coords$segmentations)
  centroids.obj <- as.Centroids(segmentations.obj)
  fov.obj <- CreateFOV(
    coords = list("segmentation" = segmentations.obj, "centroids" = centroids.obj),
    type = c("segmentation", "centroids"),
    molecules = molecules.obj
  )
  osmFISH.obj <- CreateSeuratObject(osmFISH.coords$counts, assay = assay)
  osmFISH.obj[[fov]] <- fov.obj
  return(osmFISH.obj)
}

LoadSTARmapPlus <- function(data.dir, fov = "fov", assay = "RNA") {
  files <- list.files(data.dir)
  counts.path <- file.path(data.dir, files[grepl("*raw_expression_pd.csv", files)])
  cell.positions.path <- file.path(data.dir, files[grepl("*_spatial.csv", files)])
  counts <- data.table::fread(counts.path)
  genes <- tolower(counts$GENE)
  substr(genes, 1, 1) <- toupper(substr(genes, 1, 1))
  rownames(counts) <- genes
  counts$GENE <- NULL

  cell.positions <- data.table::fread(cell.positions.path)
  cell.positions <- cell.positions[2:nrow(cell.positions),]

  # hold on to everything as metadata which is added to the object
  metadata.table <- cell.positions
  rownames(metadata.table) <- metadata.table$NAME
  metadata.table$NAME <- NULL

  cell.positions <- cell.positions[, c("NAME", "X", "Y")]
  colnames(cell.positions) <- c("cell", "x", "y")
  rownames(cell.positions) <- 1:nrow(cell.positions)
  cell.positions$x <- as.numeric(cell.positions$x)
  cell.positions$y <- as.numeric(cell.positions$y)
  object <- CreateSeuratObject(counts = counts, assay = assay)
  cents <- CreateCentroids(cell.positions)
  segmentations.data <- list("centroids" = cents)

  coords <- CreateFOV(
    coords = segmentations.data,
    type = "centroids",
    assay = fov
  )

  # subset both object and coords based on the cells shared by both
  cells <- intersect(Cells(object), Cells(x = coords, boundary = "centroids"))
  coords <- subset(x = coords, cells = cells)
  object[[fov]] <- coords
  object <- AddMetaData(object, metadata = metadata.table)
  return(object)
}

LoadMERFISH23 <- function(h5, s = 21) {
  print(h5)
  file.h5 <- H5File$new(h5, mode = "r")
  
  # first save all of the matrices
  full.mtx <- file.h5[["X"]][,]
  sample.ids <- file.h5[["obs"]][["sample_id"]][["codes"]][]
  gene.names <- file.h5[["var"]][["_index"]][]
  cellnames <- file.h5[["obs"]][["cell_id"]][]

  message(paste("Starting with", s))
  cell.idxs <- 1:length(cellnames)
  cell.idxs <- cell.idxs[sample.ids == s]

  mtx <- full.mtx[, cell.idxs]
  rownames(mtx) <- gene.names
  colnames(mtx) <- cellnames[sample.ids == s]

  # first save all of the matrices
  idxs <- 1:length(sample.ids)
  x.positions <- file.h5[["obs"]][["center_x"]][]
  y.positions <- file.h5[["obs"]][["center_y"]][]

  sample.idxs <- idxs[sample.ids == s]
  sample.cells <- cellnames[sample.ids == s]
  pos <- data.frame(
    x = x.positions[sample.idxs],
    y = y.positions[sample.idxs]
  )
  rownames(pos) <- sample.cells
  pos$cell <- rownames(pos)
  fov <- CreateFOV(
    coords = list(
      "centroids" = CreateCentroids(coords = pos)
    ), assay = "RNA"
  )
  obj <- CreateSeuratObject(counts = as.sparse(mtx))
  print(fov)
  print(obj)
  print(head(Cells(fov)))
  print(head(Cells(obj)))
  fov <- subset(fov, cells = intersect(Cells(obj), Cells(fov)))
  obj[["fov"]] <- fov

  return(obj)
}