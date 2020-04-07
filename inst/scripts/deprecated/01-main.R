### main
# S3 method:
# plot.tcr
#
# ignored:
# g
###


#' tcr viz
#'
#' tcr viz
#'
#' @param data an object of class \code{"tcr"}, usually returned by
#'   \code{get_track}
#'
#'   alternatively, an \code{n x t} integer matrix of class \code{"tcr"} with
#'   \code{n} unique sequences, \code{t} time points, and counts of each
#'   sequence at each timepoint; see examples
#' @param n.group,group grouping of time points: n per group and group label
#' @param col.group,col.gmain colors for each group and main group color
#' @param labels individual time point labels
#' @param show_bar,show_pie,show_text logicals; determines which sections
#'   will be drawn: bar plots of counts, pie charts of proportions, and text
#'   summaries, respectively
#' @param separate_first logical; if \code{TRUE}, first time point is treated
#'   as separate from later time points (only affects pie charts)
#' @param show_counts logical; if \code{TRUE}, counts are shown beside bars
#' @param inner_pie logical; if \code{TRUE}, a smaller pie chart is drawn
#'   inside the main
#'   
#'   the main pie is the proportion of each originating time point, and
#'   the inner pie is the same as the main but also with the undetected
#' @param track logical; if \code{TRUE}, arrows are drawn between bars
#' @param track_type type of tracking arrows to draw between bars; one of
#' \code{"all"} (default, detected and undetected), \code{"detected"} (for
#' only detected), or \code{"none"} to suppress all arrows
#' @param which integer index of the column to highlight (i.e., all other
#'   columns will be semi-transparent)
#' @param xlim,ylim x- and y-axis limits
#' @param ylim_type extent of ylim shown (\code{"all"} - all unique sequences,
#'   \code{"detected"} - max of number detected, \code{"post"} - max of number
#'   detected after post) -- all > detected > post; ignored if \code{ylim} is
#'   given
#' @param col.undetected color of undetected sequences (used in pie charts)
#' @param cex.main size of main labels
#' @param heights,widths see \code{\link{layout}}
#' @param file output file name to write table; if none, no table is written
#' @param top for table output, only report \code{top} clonotypes by
#'   frequency (use \code{Inf} for all)
#' @param highlight index of timepoint(s) to highlight (use 0 for all); e.g.,
#'   if \code{highlight = 2}, the second column color is fully opaque with
#'   transparency applied to other groups; see examples
#'
#' @examples
#' ## generate example data (matrix must have class "tcr")
#' set.seed(1)
#' dat <- structure(
#'   matrix(rpois(100, 0.5), 20, dimnames = list(LETTERS[1:20], letters[1:5])),
#'   class = 'tcr'
#' )
#'
#' plot(dat)
#' plot(dat, highlight = 0) ## same as above
#' plot(dat, highlight = 2) ## only second column is fully opaque
#' plot(dat, highlight = 0:ncol(dat)) ## all opaque then highlight each
#' 
#' plot(
#'   dat, n.group = c(1, 2, 2),
#'   group = c('first', 'second', 'third'),
#'   col.group = c('black', 'red', 'orange'),
#'   ylim = c(min(colSums(dat)), nrow(dat))
#' )
#'
#' @export

plot.tcr <- function(data,
                     n.group = rep_len(1L, ncol(data$track)),
                     group = rep_len('', length(n.group)),
                     col.group = seq_along(group), col.gmain = col.group,
                     labels = rep_len('', length(group)),
                     show_bar = TRUE, show_pie = TRUE, show_text = TRUE,
                     separate_first = TRUE, show_counts = FALSE, n_min = 2L,
                     inner_pie = FALSE, track_type = c('all', 'detected', 'none'),
                     which = NULL, xlim = NULL, ylim = NULL,
                     ylim_type = c('all', 'detected', 'post'),
                     col.undetected = 'grey90', cex.main = 2,
                     heights = NULL, widths = NULL, space = 10,
                     file = '', top = Inf, highlight = NULL,
                     add_diversity = TRUE,
                     indices = c('shannon', 'simpson', 'invsimpson')) {
  ## set up proper format if only a matrix is given
  if (is.matrix(data))
    data <- structure(
      list(track = data, first = max.col(+!!data, ties.method = 'first')),
      class = 'tcr'
    )
  
  if (is.numeric(col.group))
    col.group <- palette()[as.integer(col.group)]
  
  ## if highlighting, loop over desired figures and exit (no table)
  if (!is.null(highlight)) {
    palette(col.group)
    for (color in highlight) {
      Recall(data, n.group, group, col.group, col.gmain, labels, show_bar,
             show_pie, show_text, inner_pie, track_type, color,
             xlim, ylim, ylim_type, col.undetected, file, top, NULL)
    }
    palette('default')
    return(invisible(NULL))
  }
  
  write_table <- nzchar(file)
  
  oc <- palette()
  on.exit(palette(oc))
  
  ## matrix for tracking, vector with first appearance of sequence
  ## rounded version of track for some figs (pie, text) -- useful with inter > 0
  track_type <- match.arg(track_type)
  
  track  <- data$track
  # if (!track_undet)
  #   track <- round(track)
  first  <- data$first
  rtrack <- round(track, 4L)
  
  ## number of time points should equal group n total
  n <- ncol(track)
  stopifnot(n == sum(n.group))
  sin <- seq.int(n)
  
  ## denominator for pie charts (ie, detected + undetected for each time point)
  ## assumes first time point separate from others (eg, pre-vax)
  denom <- if (separate_first)
    c(list(sin), rep(list(sin[-1L]), n - 1)) else rep(list(sin), n)
  names(denom) <- colnames(track)
  
  ## vectors of group labels, group id, colors for looping
  force(col.gmain)
  groups <- rep(seq_along(group), n.group)
  group  <- rep(group, n.group)
  
  ## if one color given for each group, apply alpha
  col.group <- if (length(col.group) == length(n.group)) {
    sapply(seq_along(n.group), function(ii) {
      if (n.group[ii] == 1L)
        col.group[ii]
      else {
        a <- rescaler(sequence(n.group[ii]), c(0.5, 1))
        tcol(col.group[ii], a)
      }
    })
  } else rep_len(col.group, n)
  
  col.group <- unlist(col.group)
  col.gmain <- tcol(rep(col.gmain, n.group), 1)
  
  
  ## initiate res for table output
  res <- NULL
  
  ## calculate ylim depending on ylim_type -- ignored if ylim given
  ylim_type <- match.arg(ylim_type)
  yl <- switch(
    ylim_type,
    ## use all unique sequences to set ylim
    all = NULL,
    ## use only max detected over all time points to set ylim
    detected = c(min(colSums(!rtrack)), nrow(track)) + c(-5, 1),
    ## dont include first column of track to get ylim (includes break)
    post = c(min(colSums(!rtrack[, -1L])), nrow(track)) + c(-1, 1)
  )
  
  ## change palette for ease of coloring
  palette(col.group)
  
  ## add alpha to all other time points _except_ "which" to highlight
  if (!is.null(which)) {
    col_alpha <- ifelse(seq.int(n) != which, tcol(oc, 0.25), oc)
    if (!which %in% seq.int(n))
      col_alpha <- oc
    palette(col_alpha)
    write_table <- FALSE
  }
  
  ## set up layout -- add row(s) for pie, text figures
  nr  <- show_pie + show_text + show_bar
  mat <- matrix(seq.int(n * nr), nr)
  hts <- c(show_pie * 1, show_text * 0.5, show_bar * 5)
  layout(mat, heights = heights %||% hts[hts > 0], widths = widths %||% 1)
  
  ## x/y axis styles, padding
  par(xaxs = 'i', yaxs = 'i', oma = c(1, 3, 3, 1))
  
  
  ## loop over each time point
  ## bar plot for each col, optional pie and/or text rows
  for (ii in sin) {
    ## current time point's placement in grouping (first or last)
    fig_1 <- oig(ii, groups, 'first') & ii > 1L
    fig_n <- oig(ii, groups, 'last')  & ii < n
    
    ## add factor for better ordering of bars, rounded out when needed
    fudge <- first * 1e-4
    
    ## this time point, ordering of bars
    this <- track[, ii] - fudge
    oo_this <- order(this, decreasing = FALSE)
    
    ## the next time point and ordering, boolean if each sequence persists
    if (ii < n) {
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
      tbl <- tbl[rowSums(!!round(track)[, denom[[ii]]]) >= 1]
      ## if none/not detected, use n + 1 col.undetected
      tbl <- replace(tbl, tbl < 1, max(n) + 1L)
      col <- c(palette()[sin], col.undetected)
      tbl <- table(factor(col[tbl], col))
      
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
      ## text box with counts of total cells, unique clonotypes,
      ## number undetected for each time point
      
      ## remove fudging to get correct counts
      rth <- round(this)
      cex <- 1.5
      leg <- function(...)
        legend(..., bty = 'n', xpd = NA, cex = cex)
      
      l <- list(
        'No. total cells' = sum(rth[rth >= 1]),
        'No. clonotypes'  = length(rth[rth >= 1]),
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
    xl <- rep(list(c(0, max(c(track)))), n)
    xlim <- if (inherits(xlim, 'list'))
      rep_len(xlim, n) else xlim
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
    bp <- barplot(
      this[oo_this], col = replace(first[oo_this], na, NA),
      space = 0, border = 'white',
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
    if (ii > 1L) {
      rn <- diff(par('usr')[1:2]) / 50
      points(this[oo_this] + rn, bp, pch = '*', cex = 1, xpd = NA,
             col = replace(first[oo_this], first[oo_this] != ii, NA))
      legend('bottomright', legend = '* new', bty = 'n', cex = 1.5, text.col = ii)
    }
    
    ## axes/labels for bar plots
    if (is.null(xlim[[ii]]))
      axis(1L, range(pretty(xl[[ii]])), if (ii == 1L)
        NULL else c(NA, max(pretty(xl[[ii]]))), lwd = 0.25, cex.axis = 1.5)
    else axis(1L, lwd = 0.25, cex.axis = 1.5, xaxp = c(0, max(pretty(xlim[[ii]])), 2))
    # box(bty = 'l', lwd = 0.5)
    abline(v = 0, lwd = 0.5)
    
    # axis(2L, bp, names(this[oo_this]), cex.axis = 0.2, lwd = 0, las = 1L)
    
    ## add break in y-axis if some sequences are not being shown in ylim
    if (!ylim_type %in% 'all' && is.null(ylim)) {
      if (ii == 1L) {
        co <- c(-1, max(na) - 0.5 - 0.5, 1, max(na) - 0.5)
        sl <- (co[4L] - co[2L]) / (co[3L] - co[1L])
        seg <- function(...)
          segments(..., xpd = NA, col = 1L, lend = 3L)
        
        rect(co[1L], co[2L] - sl, co[2L], co[4L] - sl,
             border = NA, col = 'white', xpd = NA)
        seg(co[1L], co[2L], co[3L], co[4L])
        seg(co[1L], co[2L] - 0.5, co[3L], co[4L] - 0.5)
      }
    }
    
    ## this is a pain -- aesthetic lines for samples and groups
    xx <- c(0.15, 0.85)
    yy <- grconvertY(0.975, 'nfc')
    
    seg <- function(...)
      segments(..., xpd = NA, lend = 3L, lwd = 10)
    
    ## first in group adds each group label; otherwise, just line
    if (fig_1 | ii == 1L) {
    # if (0) {
      fr <- grconvertX(xx[1L], 'nfc')
      to <- grconvertX(if (ii == 1L & (sum(groups[ii] == groups) != 2L) |
                           (length(groups) == 2L & ii == 2L))
        xx[2L] else 1, 'nfc')
      seg(fr, yy, to, yy, col = col.group[ii])
      text(fr, grconvertY(0.95, 'nfc'), group[ii], col = col.gmain[ii],
           adj = 0, xpd = NA, cex = 2.5, font = 2L)
    } else if (oig(ii, groups, 'last') & ii > 1L) {
      fr <- grconvertX(0, 'nfc')
      to <- grconvertX(xx[2L], 'nfc')
      # seg(fr, yy, to, yy, col = col.gmain[ii])
      seg(fr, yy, to, yy, col = col.group[ii])
    } else {
      fr <- grconvertX(0, 'nfc')
      to <- grconvertX(1, 'nfc')
      # seg(fr, yy, to, yy, col = col.gmain[ii])
      seg(fr, yy, to, yy, col = col.group[ii])
    }
    
    ## main label for each time point
    mtext(colnames(track)[ii], col = ii, cex = cex.main, font = 2L,
          at = grconvertX(0.5, 'nfc'), adj = 0.5)
    
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
      pad <- c(0.25, 0)
      pad <- rep_len(diff(par('usr')[1:2]) / 75, 2L)
      
      ## scale lwd/alpha by number of total links
      nla <- rescaler(n_links, c(0.25, 1))
      nll <- rescaler(n_links, c(0.25, 2))
      
      ## identify undetected sequences (0 at current and/or next time point)
      undet <- if (track_type %in% 'all')
        round(this)[oo_this] < 1 | round(that)[oo_that][to] < 1
      else FALSE
      
      if (!track_type %in% 'none')
        arrows(
          this[oo_this] + pad[1L], c(bp), par('usr')[2L] + pad[2L], c(bp)[to],
          # col = tcol(ifelse(persist[oo_this], first[oo_this], 0), nla[oo_this]),
          # lwd = ifelse(undet, 0.5, nll[oo_this]),
          col = ifelse(persist[oo_this], first[oo_this], 0),
          lwd = ifelse(undet, 0.5, 1), lty = ifelse(undet, 4L, 1L),
          xpd = NA, code = 2L, length = 0.025
        )
    }
    
    ## create matrix of data for each time point:
    ## sequence, number, persist to next time (y/n), originating time point
    tbl <- cbind(
      Sequence = names(this), N = round(this), Persist = c('No', 'Yes')[persist + 1L],
      Origin = colnames(track)[first], Origin_idx = first
    )[rev(oo_this), ]
    
    if (ii == n)
      tbl[, 'Persist'] <- '<font color=black>N/A</font>'
    tbl <- head(tbl[tbl[, 'N'] > 0, ], top)
    tbl[, 'Origin'] <- sprintf(
      '<font color=%s>%s</font>',
      rm_alpha(tcol(palette()))[as.integer(tbl[, 'Origin_idx'])], tbl[, 'Origin']
    )
    tbl[, 'Persist'] <- sprintf(
      '<font color=%s>%s</font>',
      c('black', 'green')[(tbl[, 'Persist'] %in% 'Yes') + 1L], tbl[, 'Persist']
    )
    
    rownames(tbl) <- NULL
    ## drop origin_idx
    tbl <- tbl[, -ncol(tbl)]
    
    ## add diversity
    if (add_diversity) {
      counts <- as.integer(tbl[, 'N'])
      div <- diversity(counts, indices)
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
        strrep('lcccccc', n) else strrep('lccc', n),
      css.cell = 'padding: 0px 5px 0px; white-space: nowrap;',
      cgroup = sprintf('<font color=%s>%s</font>',
                       rm_alpha(tcol(palette()[sin])), colnames(track)),
      n.cgroup = rep_len(ncol(tbl), n),
      caption = 'Time point summary: Sequence, N at current time, Persisting
                 to next time point, time point of origin.'
    )
    structure(ht, class = 'htmlTable')
    
    write_htmlTable(ht, file = sprintf('~/desktop/%s.html', file))
  }
  
  invisible(res)
}

# g <- function(..., pal = palette(), order = seq_along(pal)) {
#   palette(pal)
#   for (color in order) {
#     plot(..., which = color)
#   }
#   palette('default')
#   invisible(NULL)
# }
