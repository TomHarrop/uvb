facet_adjust <- function(x, pos = c("up", "down"),
                         newpage = is.null(vp), vp = NULL) {
  # part of print.ggplot
  ggplot2:::set_last_plot(x)
  if(newpage)
    grid::grid.newpage()
  pos <- match.arg(pos)
  p <- ggplot2::ggplot_build(x)
  gtable <- ggplot2::ggplot_gtable(p)
  # finding dimensions
  dims <- apply(p$panel$layout[2:3], 2, max)
  nrow <- dims[1]
  ncol <- dims[2]
  # number of panels in the plot
  panels <- sum(grepl("panel", names(gtable$grobs)))
  space <- ncol * nrow
  # missing panels
  n <- space - panels
  # checking whether modifications are needed
  if(panels != space){
    # indices of panels to fix
    idx <- (space - ncol - n + 1):(space - ncol)
    # copying x-axis of the last existing panel to the chosen panels
    # in the row above
    gtable$grobs[paste0("axis_b",idx)] <-
      list(gtable$grobs[[paste0("axis_b",panels)]])
    if(pos == "down"){
      # if pos == down then shifting labels down to the same level as
      # the x-axis of last panel
      rows <- grep(paste0("axis_b\\-[", idx[1], "-", idx[n], "]"),
                   gtable$layout$name)
      lastAxis <- grep(paste0("axis_b\\-", panels), gtable$layout$name)
      gtable$layout[rows, c("t","b")] <- gtable$layout[lastAxis, c("t")]
    }
  }
  # again part of print.ggplot, plotting adjusted version
  if(is.null(vp)){
    grid::grid.draw(gtable)
  }
  else{
    if (is.character(vp))
      grid::seekViewport(vp)
    else grid::pushViewport(vp)
    grid::grid.draw(gtable)
    grid::upViewport()
  }
  invisible(p)
}