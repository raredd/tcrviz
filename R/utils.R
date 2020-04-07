### utils for tracking plots
# %||%, %inside%, cbindx, int, oig, rescaler, rm_alpha, tcol, write_htmlTable
###


`%||%` <- function(x, y) {
  if (is.null(x))
    y else x
}

`%inside%` <- function(x, interval) {
  interval <- sort(interval)
  x >= interval[1L] & x <= interval[2L]
}

#' cbindx
#' 
#' \code{cbind} objects with varying numbers of rows.
#' 
#' @seealso
#' \code{rawr::cbindx}
#' 
#' @examples
#' cbindx(1, 1:2, 1, 1:5)
#' 
#' @noRd

cbindx <- function(..., fill = NA) {
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

#' int
#' 
#' Insert \code{inter} between non \code{value}-sandwiched sequences.
#' 
#' @examples
#' x <- c(0, 1, 0, 0, 1, 0, 1)
#' int(x, inter = 0.1)
#' 
#' @noRd

int <- function(x, value = 0, inter = 0.1, dups = FALSE) {
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

#' Order-in-group
#' 
#' Get index of firsts or lasts by a grouping variable.
#' 
#' @examples
#' g <- c(1, 2, 2, 2, 3, 3, 3)
#' sapply(seq_along(g), function(x) oig(x, g, 'first'))
#' sapply(seq_along(g), function(x) oig(x, g, 'last')) 
#' 
#' @noRd

oig <- function(i, g, which = c('first', 'last')) {
  which <- match.arg(which)
  stopifnot(
    length(i) == 1L,
    length(g) > 0L,
    i %in% seq_along(g)
  )
  any(i == seq_along(g) & !duplicated(g, fromLast = which == 'last'))
}

#' Rescale a numeric vector
#' 
#' @seealso
#' \code{rawr::rescaler}
#' 
#' @examples
#' rescaler(1:5)
#' 
#' @noRd

rescaler <- function (x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
  (x - from[1L]) / diff(from) * diff(to) + to[1L]
}

#' Remove alpha from hex color
#' 
#' @examples
#' col <- tcol(1:3, c(NA, 0.5, 1))
#' rm_alpha(col)
#' rm_alpha(palette())
#' 
#' @noRd

rm_alpha <- function(x) {
  gsub('(^#[A-F0-9]{6})[A-F0-9]{2}$', '\\1', x, ignore.case = TRUE)
}

#' Add alpha to colors
#' 
#' @seealso
#' \code{rawr::tcol}
#' 
#' @examples
#' tcol(1:3)
#' tcol('red')
#' 
#' @noRd

tcol <- function(colors, alpha = 1) {
  v <- Vectorize(adjustcolor)
  res <- v(colors, replace(alpha, is.na(alpha) | !alpha %inside% 0:1, 1))
  unname(replace(res, is.na(colors), NA))
}

#' Write an htmlTable object (or similar) to file (or console)
#' 
#' @seealso
#' \code{rawr::write_htmlTable}
#' 
#' @examples
#' ht <- htmlTable::htmlTable(head(mtcars))
#' write_htmlTable(ht)
#' 
#' @noRd

write_htmlTable <- function(x, file = '', attributes = TRUE) {
  if (attributes) {
    x <- paste(
      '<!DOCTYPE html>\n<html>\n<body>',
      gsub('gr[ea]y\\s*(?=;)', '#bebebe', x, perl = TRUE),
      '</body>\n</html>', sep = '\n'
    )
    
    x <- structure(x, class = 'htmlTable', html = TRUE)
  }
  
  if (!is.character(file))
    x else cat(x, file = file)
}
