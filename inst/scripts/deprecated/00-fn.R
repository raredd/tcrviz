### utils for tracking plots
# %||%, %inside%, cbindx, get_data, get_track, int, oig, rescaler, rm_alpha,
# tcol, write_htmlTable
###


`%||%` <- function(x, y) {
  if (is.null(x))
    y else x
}

`%inside%` <- function(x, interval) {
  interval <- sort(interval)
  x >= interval[1L] & x <= interval[2L]
}

cbindx <- function(..., fill = NA) {
  ## cbind objects with varying numbers of rows
  ## cbindx(1, 1:2, 1, 1:5)
  l <- list(...)
  l <- Filter(Negate(is.null), l)
  n <- max(sapply(l, NROW))
  l <- lapply(l, function(x) {
    if (is.null(dim(x))) {
      y <- `length<-`(x, n)
      replace(y, seq_along(y) > length(x), fill)
    } else {
      if (nrow(x) == n)
        x else rbind(x, matrix(fill, n - nrow(x), ncol(x)))
    }
  })
  do.call('cbind', l)
}

diversity <- function(x, which = c('shannon', 'simpson', 'invsimpson')) {
  x <- x / sum(x)
  
  ## shannon, simpson
  sh <- -x * log(x)
  si <- x * x
  
  div <- c(shannon = sum(sh, na.rm = TRUE),
           simpson = 1 - sum(si),
           invsimpson = 1 / sum(si))
  div[match.arg(which, several.ok = TRUE)]
}

#' Read tracking data
#'
#' @param path,pattern,sheets path of data, pattern of file name, and sheets
#'   to be used
#'
#' @return
#' A named list of data with each element containing the data for a single
#'   time point.
#'
#' @export

get_data <- function(path, pattern, sheets = NULL) {
  lf <- list.files(path, pattern, full.names = TRUE)

  ## if sheets is null, treat each file as time point
  ## otherwise, treat each sheet as separate time point
  ll <- if (is.null(sheets))
    lapply(lf, function(x) {
      as.data.frame(readxl::read_excel(x, sheet = 'Output'))
    })
  else lapply(sheets, function(x) {
    as.data.frame(readxl::read_excel(lf, sheet = x))
  })

  names(ll) <- if (is.null(sheets))
    gsub('\\.xlsx?$', '', basename(lf))
  else as.character(sheets)

  invisible(ll)
}

#' Get tracking data
#'
#' Get all data for figure, returns list with \code{"track"} and
#' \code{"first"}.
#'
#' @param path,pattern,sheets passed to \code{\link{get_data}}
#' @param inter intermediate value used instead of 0s sandwiched between
#'   non zeros
#'
#' @export

get_track <- function(path, pattern, sheets = NULL, inter = 0) {
  suppressMessages({
    dat <- get_data(path, pattern, sheets)
  })

  ## get all unique sequences and convert to factor
  clo <- lapply(dat, function(x) x$Clonotype)
  tbl <- sort(table(unlist(clo)), decreasing = FALSE)
  lvl <- names(tbl)
  # lbl <- seq_along(lvl)
  lbl <- lvl

  dat <- lapply(dat, function(x) {
    within(x, cl_track <- factor(Clonotype, lvl, lbl))
  })
  
  ## insert "inter" between non "value"-sandwiched sequences
  track <- sapply(dat, function(x) table(x$cl_track))
  
  ## apply fudge for tracking over time of 0s/undetected
  if (inter != 0)
    track[] <- t(apply(track, 1L, int, inter = inter))

  first <- max.col(+!!track, ties.method = 'first')

  invisible(
    structure(list(track = track, first = first), class = 'tcr')
  )
}

int <- function(x, value = 0, inter = 0.1, dups = FALSE) {
  ## insert "inter" between non "value"-sandwiched sequences
  ## x <- c(0, 1, 0, 0, 1, 0, 1); int(x, inter = 0.1)
  
  ## if inserted values should be unique (ascending) or dups are ok
  FUN <- if (dups)
    identity else cumsum
  
  rl <- rle(x)
  sa <- seq_along(rl$values)
  rp <- rl$values == value & sa > 1L & sa < length(sa)
  ## divide by 100 to avoid getting close to 1
  rl$values[rp] <- inter + FUN(rep(inter / 100, sum(rp)))
  res <- inverse.rle(rl)
  
  rp <- res < 1 & res > 0
  replace(res, rp, FUN(res[rp] - inter))
}

oig <- function(i, g, which = c('first', 'last')) {
  ## get index of firsts or lasts by a grouping variable
  ## g <- c(1, 2, 2, 2, 3, 3, 3)
  ## sapply(seq_along(g), function(x) oig(x, g, 'first'))
  ## sapply(seq_along(g), function(x) oig(x, g, 'last'))
  stopifnot(
    length(i) == 1L,
    length(g) > 0L,
    i %in% seq_along(g)
  )
  which <- match.arg(which)

  any(i == seq_along(g) & !duplicated(g, fromLast = which == 'last'))
}

rescaler <- function (x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
  ## rawr::rescaler
  ## rescaler(1:5)
  (x - from[1L]) / diff(from) * diff(to) + to[1L]
}

rm_alpha <- function(x) {
  ## remove alpha from hex string
  gsub('(^#[A-F0-9]{6})[A-F0-9]{2}$', '\\1', x, ignore.case = TRUE)
}

tcol <- function(colors, alpha = 1) {
  ## add alpha to colors
  v <- Vectorize(adjustcolor)
  res <- v(colors, replace(alpha, is.na(alpha) | !alpha %inside% 0:1, 1))
  unname(replace(res, is.na(colors), NA))
}

write_htmlTable <- function(x, file = '', attributes = TRUE) {
  ## rawr::write_htmlTable
  ## write an htmlTable object to file (or console)
  if (attributes) {
    x <- paste('<!DOCTYPE html>\n<html>\n<body>',
               gsub('gr[ea]y\\s*(?=;)', '#bebebe', x, perl = TRUE),
               '</body>\n</html>', sep = '\n')
    x <- structure(x, class = 'htmlTable', html = TRUE)
  }
  
  if (!is.character(file))
    x else cat(x, file = file)
}
