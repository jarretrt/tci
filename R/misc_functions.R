# --------------------------------------------------------------------------------------------------------------------------------
# - Misc. Functions --------------------------------------------------------------------------------------------------------------
# - Helper functions that aren't central to use in the package -------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------


#' Color palate for tci plotting functions
#'
#' @examples pal["navy"]
#'
#' @export
pal  <- c(black = "#020201",
          navy  = "#6f859b",
          brown = "#7c4c4d",
          teal  = "#60ccd9",
          pink  = "#e3bab3",
          darkgrey = "#97898a",
          grey  = "#c2c5cb")


#' Extract last element or column
#'
#' Function to extract the last element from a vector or the last column from a matrix
#'
#' @param x Vector or matrix
#'
#' @export
tail_vec <- function(x){

  args <- list(x)[[1]]

  if(is.null(dim(args))){
    return(args[length(args)])
  } else{
    return(args[,ncol(args)])
  }
}
#' @examples
#' tail_vec(1:8)
#' tail_vec(matrix(1:8,2,4))



#' Modified 'cut' function
#'
#' @name cut2
#'
#' Modified cut function to create equal-length intervals. See documentation
#' for base::cut for description of function arguments. Not currently used.
#'
#' @export
cut2 <- function(x, breaks, labels = NULL, include.lowest = FALSE, right = TRUE,
                  dig.lab = 3L, ordered_result = FALSE, ...)
{
  if (!is.numeric(x))
    stop("'x' must be numeric")
  if (length(breaks) == 1L) {
    if (is.na(breaks) || breaks < 2L)
      stop("invalid number of intervals")
    nb <- as.integer(breaks + 1)
    dx <- diff(rx <- range(x, na.rm = TRUE))
    if (dx == 0) {
      dx <- if (rx[1L] != 0)
        abs(rx[1L])
      else 1
      breaks <- seq.int(rx[1L] - dx/1000, rx[2L] + dx/1000,
                        length.out = nb)
    }
    else {
      breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
      breaks[c(1L, nb)] <- c(rx[1L] - dx/1000, rx[2L] +
                               dx/1000)
    }
  }
  else nb <- length(breaks <- sort.int(as.double(breaks)))
  if (anyDuplicated(breaks))
    stop("'breaks' are not unique")
  codes.only <- FALSE
  if (is.null(labels)) {
    for (dig in dig.lab:max(12L, dig.lab)) {
      # ch.br <- format(0 + breaks, digits = dig, width = 1L)
      ch.br <- sprintf(0+breaks, fmt = "%#.3f")
      if (ok <- all(ch.br[-1L] != ch.br[-nb]))
        break
    }
    labels <- if (ok)
      paste0(if (right)
        "("
        else "[", ch.br[-nb], ",", ch.br[-1L], if (right)
          "]"
        else ")")
    else paste0("Range_", seq_len(nb - 1L))
    if (ok && include.lowest) {
      if (right)
        substr(labels[1L], 1L, 1L) <- "["
      else substring(labels[nb - 1L], nchar(labels[nb -
                                                     1L], "c")) <- "]"
    }
  }
  else if (is.logical(labels) && !labels)
    codes.only <- TRUE
  else if (length(labels) != nb - 1L)
    stop("lengths of 'breaks' and 'labels' differ")
  code <- .bincode(x, breaks, right, include.lowest)
  if (codes.only)
    code
  else factor(code, seq_along(labels), labels, ordered = ordered_result)
}




#' Modified seq function to include bounds
#'
#' @param from sequence starting value
#' @param to sequence end value
#' @param by increment of the sequence
#'
#' @export
seqby <- function(from, to, by)
  sort(union(seq(from, to, by), c(from, to)))
#' @examples
#' seq(0,0.767,1/6)
#' seqby(0,0.767,1/6)



