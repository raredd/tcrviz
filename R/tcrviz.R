### main tcrviz fns
# tcrviz, tcrdata, diversity
###


#' tcrviz
#'
#' TCR tracking visualizations.
#'
#' @param data an \code{n x t} count matrix with \code{n} unique sequences
#' for \code{t} time points
#' @param inter logical; use intermediate tracking for sequences; otherwise,
#' tracking only shown between columns with > 0 counts
#' @param n.group,group grouping of time points: n per group and group label
#' @param col.group,col.gmain colors for each group and main group color
#' @param labels individual time point labels
#' @param show_bar,show_pie,show_text logicals; determine which sections will
#' be drawn: bar plots of counts, pie charts of proportions, and text
#' summaries, respectively
#' @param denom a vector (or list of vectors for each column of \code{data})
#' of column indices to get the denominator for number of unique sequences;
#' see examples
#' @param show_counts logical; if \code{TRUE}, counts are shown beside bars
#' @param n_min optional y-axis tick mark added at \code{n_min}, e.g., for
#' cutting singletons from figure 
#' @param inner_pie logical; if \code{TRUE}, an inner pie chart added inside
#' main; the main pie is the proportion of each originating time point, and
#' the inner pie also includes undetected
#' @param track_type type of tracking arrows to draw between bars; one of
#' \code{"all"} (default, detected and undetected), \code{"detected"} (for
#' only detected), or \code{"none"} to suppress all arrows
#' @param which integer index of the column to highlight (i.e., all other
#' columns will be semi-transparent); use \code{0} to highlight all
#' @param xlim,ylim x- and y-axis limits
#' @param ylim_type extent of ylim shown (\code{"all"} - all unique sequences,
#' \code{"detected"} - max of number detected, \code{"post"} - max of number
#' detected after post) -- \code{all > detected > post}; ignored if
#' \code{ylim} is given
#' @param col.undetected color of undetected sequences (used in pie charts)
#' @param cex.main size of main labels
#' @param asterisk logical; if \code{TRUE}, new sequences are marked with
#' \code{*} at the originating time point
#' @param ids,col.ids sequences to track and highlight with \code{col.ids}
#' @param heights,widths see \code{\link{layout}}
#' @param space number of lines between top of bars from rest of plot
#' @param file output file name to write table; if none, no table is written
#' @param top for table output, only report \code{top} unique sequences by
#' frequency (use \code{Inf} for all)
#' @param highlight index of timepoint(s) to highlight (use 0 for all); e.g.,
#' if \code{highlight = 2}, the second column color is fully opaque with
#' transparency applied to other groups; see examples
#' @param reset_par logical; if \code{TRUE} (default), graphical parameters
#' are reset to state before plot call; use \code{FALSE} to continue adding
#' to figure
#' @param indices diversity indices to add to table output (all selected are
#' calculated for each column of \code{data}); one or more of \code{"shannon"},
#' \code{"simpson"}, or \code{"invsimpson"}; use \code{"none"} to exclude
#' 
#' @seealso
#' \code{\link{diversity}}
#'
#' @examples
#' ## generate example data (matrix must have class "tcr" and dimnames)
#' set.seed(1)
#' l <- 5:1 / 10
#' dat <- structure(
#'   sapply(l, function(x)
#'     sample(sort(rnbinom(50, 1, x)), prob = seq.int(50))),
#'   .Dimnames = list(1:50, letters[1:5])
#' )
#' 
#' ## intermediate tracking of undetected sequences (default)
#' tcrviz(dat, inter = TRUE)
#' ## no tracking undetected sequences, ie, lines only between >0 counts
#' tcrviz(dat, inter = FALSE)
#' 
#' 
#' ## highlight individual sequences
#' ids <- rownames(dat)[rowSums(!!dat) == 2]
#' tcrviz(dat, ids = ids, col.ids = 'purple', track_type = 'none')
#' tcrviz(dat, ids = ids, col.ids = 'purple',
#'        col.group = adjustcolor(palette(), 0.25))
#' 
#' 
#' ## highlighting time points
#' tcrviz(dat, highlight = 0) ## same as above
#' tcrviz(dat, highlight = 2) ## only second column is fully opaque
#' 
#' \dontrun{
#' ## all opaque then highlight each
#' par(ask = TRUE)
#' tcrviz(dat, highlight = 0:ncol(dat))
#' par(ask = FALSE)
#' 
#' tcrviz(
#'   dat, n.group = c(1, 2, 2), inner_pie = TRUE, n_min = 1,
#'   group = c('first', 'second', 'third'), track_type = 'det',
#'   # denom = list(1:5, 2:5, 2:5, 2:5, 2:5), ## cols 2-5 independent of 1
#'   col.group = c('grey', 'tomato', 'orange'),
#'   ## use different xlim for each panel
#'   xlim = lapply(1:5, function(x) c(0, max(dat[, x])))
#' )
#' }
#'
#' @export

tcrviz <- function(data, inter = TRUE, n.group = rep_len(1L, ncol(data)),
                   group = rep_len('', length(n.group)),
                   col.group = seq_along(group), col.gmain = col.group,
                   labels = rep_len('', length(group)),
                   show_bar = TRUE, show_pie = TRUE, show_text = TRUE,
                   denom = seq.int(ncol(data)), show_counts = FALSE,
                   n_min = NULL, inner_pie = TRUE,
                   track_type = c('all', 'detected', 'none'),
                   which = NULL, xlim = NULL, ylim = NULL,
                   ylim_type = c('all', 'detected', 'post'),
                   col.undetected = 'grey90', cex.main = 2,
                   asterisk = FALSE, ids = NULL, col.ids = 'coral3',
                   heights = NULL, widths = NULL, space = 5,
                   file = '', top = Inf, highlight = NULL, reset_par = TRUE,
                   indices = c('shannon', 'simpson', 'invsimpson', 'none')) {
  odata <- data
  force(n.group); force(denom)
  if (isTRUE(inter))
    data[] <- t(apply(data, 1L, int))
  
  data <- list(
    track = data,
    first = max.col(+!!odata, ties.method = 'first'),
    last  = max.col(+!!odata, ties.method = 'last')
  )
  
  ## if highlighting, loop over desired figures and exit (no table)
  if (!is.null(highlight)) {
    if (is.numeric(col.group))
      col.group <- palette()[as.integer(col.group)]
    col.group <- rep_len(col.group, ncol(odata))
    
    palette(col.group)
    for (color in highlight) {
      if (color > ncol(odata))
        palette(tcol(palette(), 0.25))
      Recall(
        odata, inter, n.group, group, col.group, col.gmain, labels, show_bar,
        show_pie, show_text, denom, show_counts, n_min, inner_pie, track_type,
        color, xlim, ylim, ylim_type, col.undetected, cex.main, asterisk, ids,
        col.ids, heights, widths, space, file, top, NULL, TRUE, indices
      )
    }
    palette('default')
    return(invisible(NULL))
  }
  
  write_table <- nzchar(file)
  
  op <- par(no.readonly = TRUE)
  oc <- palette()
  on.exit({
    palette(oc)
    par(op)
  })
  
  ## matrix for tracking, vector with first appearance of sequence
  ## rounded version of track for some figs (pie, text) -- useful with inter > 0
  track_type <- match.arg(track_type)
  
  track  <- data$track
  # if (!track_undet)
  #   track <- round(track)
  first  <- data$first
  rtrack <- round(track, 4L)
  
  ## number of time points should equal group n total
  nc <- ncol(track)
  stopifnot(nc == sum(n.group))
  sin <- seq.int(nc)
  
  ## denominator for pie charts (ie, detected + undetected for each time point)
  ## assumes first time point separate from others (eg, pre-vax)
  denom <- rep_len(if (inherits(denom, 'list')) denom else list(denom), nc)
  names(denom) <- colnames(track)
  
  ## vectors of group labels, group id, colors for looping
  groups <- rep(seq_along(group), n.group)
  group  <- rep(group, n.group)
  
  ## if one color given for each group, apply alpha
  col.group <- if (length(col.group) == length(n.group)) {
    if (is.numeric(col.group))
      col.group <- palette()[as.integer(col.group)]
    sapply(seq_along(n.group), function(ii) {
      if (n.group[ii] == 1L)
        col.group[ii]
      else {
        a <- rescaler(sequence(n.group[ii]), c(0.5, 1))
        tcol(col.group[ii], a)
      }
    })
  } else rep_len(col.group, nc)
  
  col.group <- unlist(col.group)
  col.gmain <- rm_alpha(col.group)
  
  
  ## initiate res for table output
  res <- NULL
  
  ## calculate ylim depending on ylim_type -- ignored if ylim given
  ylim_type <- match.arg(ylim_type)
  yl <- switch(
    ylim_type,
    ## use all unique sequences to set ylim
    all = NULL,
    ## use only max detected over all time points to set ylim
    detected = pmax(1, c(min(colSums(!rtrack)), nrow(track)) + c(-5, 1)),
    ## dont include first column of track to get ylim (includes break)
    post = pmax(1, c(min(colSums(!rtrack[, -1L])), nrow(track)) + c(-1, 1))
  )
  
  
  ## change palette for ease of coloring
  palette(col.group)
  
  ## add alpha to all other time points _except_ "which" to highlight
  if (!is.null(which)) {
    col_alpha <- ifelse(seq.int(nc) != which, tcol(oc, 0.25), oc)
    if (!which %in% seq.int(nc))
      col_alpha <- oc
    palette(col_alpha)
    write_table <- FALSE
  }
  
  ## set up layout -- add row(s) for pie, text figures
  nr  <- show_pie + show_text + show_bar
  mat <- matrix(seq.int(nc * nr), nr)
  hts <- c(show_pie * 1, show_text * 0.5, show_bar * 7)
  layout(mat, heights = heights %||% hts[hts > 0], widths = widths %||% 1)
  
  ## x/y axis styles, padding
  par(xaxs = 'i', yaxs = 'i', oma = c(1, 3, 3, 1))
  
  
  ## loop over each time point
  ## bar plot for each col, optional pie and/or text rows
  for (ii in sin) {
    ## current time point's placement in grouping (first or last)
    fig_1 <- oig(ii, groups, 'first') & ii > 1L
    fig_n <- oig(ii, groups, 'last')  & ii < nc
    
    ## add factor for better ordering of bars, rounded out when needed
    fudge <- first * 1e-4
    
    ## this time point, ordering of bars
    this <- track[, ii] - fudge
    oo_this <- order(this, decreasing = FALSE)
    
    ## the next time point and ordering, boolean if each sequence persists
    if (ii < nc) {
      that <- track[, ii + 1L] - fudge
      oo_that <- order(that, decreasing = FALSE)
      persist <- rowSums(!!rtrack[, c(ii, ii + 1L)]) == 2L
    } else {
      that <- oo_that <- NULL
      persist <- rep_len(0, nrow(track))
    }
    
    
    ## used to scale lwd/alpha trans in connections
    ## longer persistence == thicker, darker lines
    n_links <- rowSums(!!round(track))
    
    
    ## optional pie chart -- top row of each time point
    if (show_pie) {
      # browser()
      ## pie chart of not detected + timepoints
      tbl <- round(this[oo_this]) >= 1
      tbl <- tbl * first[oo_this]
      ## only use certain cols of track to get denominator
      tbl <- tbl[rowSums(!!round(track)[, denom[[ii]], drop = FALSE]) >= 1]
      ## if none/not detected, use palette()[nc + 1] -- col.undetected
      tbl <- replace(tbl, tbl < 1, max(nc) + 1L)
      col <- c(palette()[sin], col.undetected)
      tbl <- table(factor(col[tbl], unique(col)))
      
      ## exclude not detected
      tbl2 <- tbl[-length(tbl)]
      
      ## outer pie ("halo") with only detected
      par(mar = c(1, 1, 1, 1))
      pie(tbl2, col = names(tbl), border = 'white', labels = NA,
          clockwise = TRUE, radius = 1)
      
      ## inner pie chart with also undetected
      if (inner_pie) {
        ## spacer ring
        par(new = TRUE)
        pie(1, col = 'white', border = NA, labels = NA, radius = 0.75)
        
        par(new = TRUE)
        pie(tbl, col = names(tbl), border = NA, labels = NA,
            clockwise = TRUE, radius = 0.725)
      }
      
      ## short label above each pie - sample name
      mtext(labels[ii], line = 1, font = 2L)
      
      ## add lines separating groups of samples
      if (fig_1)
        abline(v = grconvertX(0, 'nfc'), col = 'black', lwd = 2, xpd = TRUE)
      if (fig_n)
        abline(v = grconvertX(1, 'nfc'), col = 'black', lwd = 2, xpd = TRUE)
    }
    
    
    ## optional text summary -- middle row of each time point
    if (show_text) {
      ## text box with counts of total cells, unique sequences, number
      ## undetected for each time point
      
      ## remove fudging to get correct counts
      rth <- round(this)
      cex <- 1.5
      leg <- function(...)
        legend(..., bty = 'n', xpd = NA, cex = cex)
      
      l <- list(
        'No. total cells' = sum(rth[rth >= 1]),
        'No. sequences'   = length(rth[rth >= 1]),
        'No. undetected'  = length(rth[rth < 1])
      )
      
      if (!inner_pie)
        l <- c(' '  = NA, l[1:2])
      
      nl <- format(unlist(l), big.mark = ',')
      nl[grepl('NA', nl)] <- ''
      
      par(mar = c(0, 1, 0, 1))
      plot.new()
      m <- max(strwidth(paste(names(l), '  '), cex = cex))
      leg(0, 0.9, names(l))
      leg(0 + m, 0.9, nl)
    }
    
    ## calculate x/y limits
    xl <- rep(list(c(0, max(c(track)))), nc)
    xlim <- if (inherits(xlim, 'list'))
      rep_len(xlim, nc) else xlim
    yl <- yl %||% c(0, nrow(track) + 1L)
    na <- sequence(yl[1L])
    if (!ylim_type %in% 'all')
      na <- c(na, max(na) + if (ylim_type %in% 'post') 3:4 else 3:4)
    
    if (!show_bar)
      next
    
    ## main figure -- ordered barplots for each time point
    ##   bars colored by time point each sequence first appears
    ##   connected by tracking lines with same color
    par(mar = c(3, 0.25, space, 0.25))
    col <- replace(first[oo_this], na, NA)
    col[names(this[oo_this]) %in% ids] <- col.ids
    
    bp <- barplot(
      this[oo_this], col = col, space = 0, border = 'white',
      xlim = xlim[[ii]] %||% xl[[ii]] * c(1, 1.15), ylim = ylim %||% yl,
      axes = FALSE, horiz = TRUE, names.arg = FALSE
    )
    if (show_counts) {
      lbl <- track[oo_this, ii]
      lbl[duplicated(lbl, fromLast = TRUE)] <- ''
      text(this[oo_this], bp, lbl, pos = 4L)
    }
    if (!is.null(n_min)) {
      at <- max(bp[which(track[oo_this, ii] < n_min)])
      axis(2L, at, sprintf('< %s', n_min), las = 1L)
      # abline(h = at, lty = 'dashed')
    }
    if (asterisk && ii > 1L) {
      rn <- diff(par('usr')[1:2]) / 50
      col <- replace(first[oo_this], first[oo_this] != ii, NA)
      col[names(this[oo_this]) %in% ids & !is.na(col)] <- col.ids
      
      points(this[oo_this] + rn, bp, pch = '*', xpd = NA, col = col)
      legend('bottomright', legend = '* new', bty = 'n', cex = 1.5, text.col = ii)
    }
    
    ## axes/labels for bar plots
    if (is.null(xlim[[ii]]))
      axis(1L, range(pretty(xl[[ii]])), if (ii == 1L)
        NULL else c(NA, max(pretty(xl[[ii]]))), lwd = 0.25, cex.axis = 1.5)
    else axis(1L, lwd = 0.25, cex.axis = 1.5, xaxp = c(0, max(pretty(xlim[[ii]])), 2))
    # box(bty = 'l', lwd = 0.5)
    abline(v = 0, lwd = 0.5)
    
    seg <- function(...) {
      segments(..., xpd = NA, lend = 3L)
    }
    
    ## add break in y-axis if some sequences are not being shown in ylim
    if (!ylim_type %in% 'all' && is.null(ylim)) {
      if (ii == 1L) {
        xx <- par('usr')[1L] + diff(par('usr')[1:2]) / c(-50, 50)
        df <- diff(par('usr')[3:4]) / 50
        yy <- max(na) + df * c(-0.5, 0.5)
        yf <- df / 2 - 0.25

        rect(xx[1L], yy[1L], xx[2L], yy[2L],
             col = 'white', xpd = NA, border = NA)
        seg(xx[1L], yy[1L] - yf, xx[2L], yy[1L] + yf, lwd = 0.25)
        seg(xx[1L], yy[2L] - yf, xx[2L], yy[2L] + yf, lwd = 0.25)
      }
    }
    
    ## this is a pain -- aesthetic lines for samples and groups
    xx <- c(0.15, 0.85)
    yy <- grconvertY(0.975, 'nfc')
    
    ## first in group adds each group label; otherwise, just line
    if (fig_1 | ii == 1L) {
      fr <- grconvertX(xx[1L], 'nfc')
      # to <- grconvertX(if (ii == 1L & (sum(groups[ii] == groups) != 2L) |
      #                      (length(groups) == 2L & ii == 2L))
      #   xx[2L] else 1, 'nfc')
      to <- grconvertX(xx[2L], 'nfc')
      seg(fr, yy, to, yy, col = ii, lwd = 10)
      text(fr, grconvertY(0.95, 'nfc'), group[ii], col = col.gmain[ii],
           adj = 0, xpd = NA, cex = 2.5, font = 2L)
    } else if (oig(ii, groups, 'last') & ii > 1L) {
      fr <- grconvertX(0, 'nfc')
      to <- grconvertX(xx[2L], 'nfc')
      seg(fr, yy, to, yy, col = ii, lwd = 10)
    } else {
      fr <- grconvertX(0, 'nfc')
      to <- grconvertX(1, 'nfc')
      fr <- grconvertX(xx[1L], 'nfc')
      to <- grconvertX(xx[2L], 'nfc')
      seg(fr, yy, to, yy, col = ii, lwd = 10)
    }
    
    ## main label for each time point
    mtext(colnames(track)[ii], col = ii, cex = cex.main,
          font = 2L, at = grconvertX(0.5, 'nfc'), adj = 0.5)
    
    ## tracking lines for persisting sequences
    if (!is.null(that)) {
      ## get location of sequences at following time point
      fr <- names(this[oo_this])
      to <- names(that[oo_that])
      to <- match(fr, to)
      
      ## make sure from and to sequences are identical
      stopifnot(
        all(names(this[oo_this]) == names(that[oo_that][to]))
      )
      
      ## extra padding at top of bar, end of arrow
      pad <- rep_len(diff(par('usr')[1:2]) / 50, 2L) * c(asterisk + 1L, 1L)
      
      ## scale lwd/alpha by number of total links (not used)
      nla <- rescaler(n_links, c(0.25, 1))
      nll <- rescaler(n_links, c(0.25, 2))
      
      ## identify undetected sequences (0 at current and/or next time point)
      undet <- if (track_type %in% 'all')
        round(this)[oo_this] < 1 | round(that)[oo_that][to] < 1
      else FALSE
      
      if (!track_type %in% 'none') {
        col <- ifelse(persist[oo_this], first[oo_this], NA)
        col[names(this[oo_this]) %in% ids] <- col.ids
        
        arrows(
          this[oo_this] + pad[1L], c(bp), par('usr')[2L] + pad[2L], c(bp)[to],
          col = col, lwd = ifelse(undet, 0.5, 1), lty = ifelse(undet, 4L, 1L),
          xpd = NA, code = 2L, length = 0.025
        )
      } else {
        col <- rep_len(col.ids[1L], length(oo_this))
        col[(!(names(this[oo_this]) %in% ids)) |
              (data$last[oo_this] <= ii) |
              data$first[oo_this] > ii] <- NA
        
        arrows(
          this[oo_this] + pad[2L], c(bp), par('usr')[2L] + pad[2L], c(bp)[to],
          col = col, xpd = NA, code = 2L, length = 0.025
        )
      }
    }
    
    ## create matrix of data for each time point:
    ## sequence, number, persist to next time (y/n), originating time point
    tbl <- cbind(
      Sequence = names(this),
      N = round(this),
      Persist = c('No', 'Yes')[persist + 1L],
      Origin = colnames(track)[first],
      Origin_idx = first
    )[rev(oo_this), ]
    
    if (ii == nc)
      tbl[, 'Persist'] <- '<font color=black>N/A</font>'
    tbl <- tbl[tbl[, 'N'] > 0, ]
    tbl <- tbl[seq.int(if (is.finite(top)) top else nrow(tbl)), ]
    
    tbl[, 'Origin'] <- sprintf(
      '<font color=%s>%s</font>',
      rm_alpha(tcol(palette()))[as.integer(tbl[, 'Origin_idx'])],
      tbl[, 'Origin']
    )
    
    tbl[, 'Persist'] <- sprintf(
      '<font color=%s>%s</font>',
      c('black', 'green')[(tbl[, 'Persist'] %in% 'Yes') + 1L],
      tbl[, 'Persist']
    )
    
    rownames(tbl) <- NULL
    ## drop origin_idx
    tbl <- tbl[, -ncol(tbl)]
    
    ## add diversity
    if (add_diversity <-
        !(identical('none', indices) | identical(indices, FALSE))) {
      counts <- as.integer(tbl[, 'N'])
      div <- diversity(counts, setdiff(indices, c('none')))
      
      tbl <- cbindx(
        tbl,
        Frequency = round(counts / sum(counts), 3L),
        Index = names(div),
        Measure = sprintf('%.4f', div)
      )
    }
    
    res <- if (ii == 1L)
      tbl else cbindx(res, tbl)
  }
  
  if (write_table) {
    ht <- htmlTable::htmlTable(
      res, align = if (add_diversity)
        strrep('lcccccc', nc) else strrep('lccc', nc),
      css.cell = 'padding: 0px 5px 0px; white-space: nowrap;',
      cgroup = sprintf('<font color=%s>%s</font>',
                       rm_alpha(tcol(palette()[sin])), colnames(track)),
      n.cgroup = rep_len(ncol(tbl), nc),
      caption = sprintf(
        'Summary of tcrviz: Sequence, count at current time, persisting to
        next time point, time point of origin%s.', if (add_diversity)
          'diversity indices' else ''
      )
    )
    
    write_htmlTable(ht, file = file)
  }
  
  invisible(res)
}

#' tcr data
#' 
#' Convenience function to format data for \code{tcrviz}.
#' 
#' @param data a data frame
#' @param seq,cnt column names or indices of \code{data} of sequences  and
#' counts, respectively
#' @param labels optional time point labels
#' 
#' @return
#' A matrix in the proper format for \code{\link{tcrviz}}.
#' 
#' @examples
#' ## alternative data format
#' dat <- tcrviz:::cbindx(
#'   time1 = letters[1:12], n1 = 1:12,
#'   space = NA,
#'   time2 = letters[2:10], n2 = 2:10,
#'   space = NA,
#'   time3 = letters[1:10], n3 = 1:10,
#'   space = NA,
#'   time4 = letters[6:15], n4 = 6:15
#' )
#' dat <- data.frame(dat)
#' dat[] <- lapply(dat, type.convert)
#' 
#' ## data in tcrviz format
#' tcrdata(dat)
#' tcrdata(dat, labels = paste0('week', 1:4))
#' 
#' @export

tcrdata <- function(data, seq = NULL, cnt = NULL, labels = NULL) {
  data <- data[, colSums(is.na(data)) < nrow(data)]
  
  cn <- colnames(data)
  cl <- sapply(data, class)
  
  seq <- data[, if (is.null(seq))
    cl %in% c('character', 'factor') else seq, drop = FALSE]
  cnt <- data[, if (is.null(cnt))
    cl %in% c('integer', 'numeric') else cnt, drop = FALSE]
  
  f <- function(x, y) {
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    rep(x, y)
  }
  
  ll <- Map(f, seq, cnt)
  ss <- levels(as.factor(unlist(ll)))
  dd <- data.frame(
    sequence = unlist(ll),
    timepoint = rep(names(ll), lengths(ll)),
    stringsAsFactors = FALSE
  )
  dd <- within(dd, {
    timepoint <- factor(timepoint, names(ll))
    sequence <- factor(sequence, ss)
  })
  
  res <- as.matrix(unclass(table(dd$sequence, dd$timepoint)))
  colnames(res) <- labels %||% colnames(seq)
  
  res
}

#' Diversity
#' 
#' Calculate diversity indices.
#' 
#' @param x a vector of counts
#' @param which the diversity index to calculate, one or more of
#' \code{"shannon"}, \code{"simpson"}, or \code{"invsimpson"}
#' 
#' @examples
#' x <- rpois(50, 5)
#' diversity(x, 'shannon')
#' diversity(x)
#' 
#' @export

diversity <- function(x, which = c('shannon', 'simpson', 'invsimpson')) {
  x <- x[is.finite(x)]
  x <- x / sum(x)
  
  ## shannon, simpson
  sh <- -x * log(x)
  si <- x * x
  
  div <- c(
    shannon = sum(sh, na.rm = TRUE),
    simpson = 1 - sum(si),
    invsimpson = 1 / sum(si)
  )
  
  div[match.arg(which, several.ok = TRUE)]
}
