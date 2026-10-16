#' Limnotrack QC colour palette — curated from pals::alphabet()
#'
#' 21 distinct colours for multi-lake plots. Darkest shades and most
#' similar pairs removed from the original 26-colour alphabet palette.
#'
#' @param n integer, number of colours needed (max 21)
#' @return character vector of hex colours
limno_colours <- function(n = 21) {
  pal <- c('#F0A0FF', '#0075DC', '#993F00', '#005C31', '#2BCE48',
           '#808080', '#94FFB5', '#8F7C00', '#9DCC00', '#C20088',
           '#003380', '#FFA405', '#426600', '#FF0010', '#5EF1F2',
           '#00998F', '#740AFF', '#990000', '#FFFF80', '#FFE100',
           '#FF5005')
  if (n > length(pal)) warning("Only 21 colours available, recycling")
  rep_len(pal, n)
}

#' ggplot2 colour scale using limnotrack palette
#'
#' @param ... passed to scale_colour_manual
#' @return ggplot2 scale
scale_colour_limno <- function(...) {
  ggplot2::scale_colour_manual(values = limno_colours(), ...)
}

#' ggplot2 fill scale using limnotrack palette
#'
#' @param ... passed to scale_fill_manual
#' @return ggplot2 scale
scale_fill_limno <- function(...) {
  ggplot2::scale_fill_manual(values = limno_colours(), ...)
}
