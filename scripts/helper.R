# #' Create a subchunk calling DT::datatable(dataframe) or plotly::ggplotly(plot) 
#' 
#' @param name string. subchunk name
#' @param input string. 'df' for data.frame or 'p' for ggplot object
#' @param width integer. plot width
#' @param height integer. plot height
#' @param ggplotly logical. TRUE for interactive plotting using plotly::ggplot,
#'                          FALSE for non-interactive plotting
#' @return string to print a subchunk
subchunkify <- function(name, input='df', width=7, height=5, ggplotly=TRUE) {
    if (input == 'df') {
        t_deparsed <- paste0("DT::datatable(df)")
        more <- ""
    } else if (input == 'plot') {
        t_deparsed <- paste0("plotly::ggplotly(p)")
        if (!ggplotly) {
            t_deparsed <- paste0("print(p)")
        }
        more <- paste0(", fig.width=", width, ", fig.height=", height)
    } else {
        stop('Incorrect argument: input')
    }
    sub_chunk <- paste0("```{r sub_chunk_", name, ", results='asis', echo=FALSE", more, "}",
        "\n",
        t_deparsed,
        "\n\n```\n\n\n")

    cat(knitr::knit(text = sub_chunk, quiet=TRUE))
}


#' Print a Markdown link for a plot
#'
#' @param fn filepath
#' @return String ready to be printed
link.plot <- function(fn){
    cat(paste('\n\n- Download Plot: [', fn, '](', fn, ')\n\n'))
}

#' Print a Markdown link for a table
#'
#' @param fn filepath
#' @return String ready to be printed
link.table <- function(fn){
    cat(paste('\n\n- Download Table: [', fn, '](', fn, ')\n\n'))
}


#' Create a list storing column names for plotting
#'
#' @param modality rna or protein
#' @return list
plot_labels <- function(modality) {
    list(
        MA=c(
            xmetric=paste0('baseMean_', modality),
            ymetric=paste0('log2FoldChange_', modality),
            xlabel='Average Expression',
            ylabel='log2FoldChange'),
        Volcano=c(
            xmetric=paste0('log2FoldChange_', modality),
            ymetric='log_odds',
            xlabel='log2FoldChange',
            ylabel='-log10(FDR)')
    )
}

#' Create a list determining max values for given columns
#'
#' @param df input data frame
#' @param xcol column name for x-axis
#' @param ycol column name for y-axis
#' @param plottype type of plot, either "MA" or "Volcano"
#' @return list
xyend <- function(df, xcol, ycol, plottype="MA") {

    # Create a vector storing +/- max values for x
    xvec <- df[[xcol]][!is.na(df[[xcol]])]
    xend <- max(abs(xvec))
    xend <- c(-xend, xend) %>% sort()

    # Create a vector storing +/- max values for y
    yvec <- df[[ycol]][!is.na(df[[ycol]])]
    yend <- max(abs(yvec))
    yend <- c(-yend, yend) %>% sort()
    if (plottype == "Volcano") {
        yend[1] <- 0
    }
    return(list(xend=xend, yend=yend))
}

