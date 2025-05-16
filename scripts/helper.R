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


#' Get the OrgDb for the specified organism, using the cached AnnotationHub.
#'
#' @param species Case-sensitive genus and species
#' @param cache Directory in which the AnnotationHub cache is stored
#' @param annotation_key_override If not NA, forces the hub to use this
#'        accession. Use this when you know exactly which OrgDb to use.
#'
#' @return OrgDb object
get_orgdb <- function(species, cache, annotation_key_override=NA){

    # Workaround to allow AnnotationHub to use proxy. See
    # https://github.com/Bioconductor/AnnotationHub/issues/4, and thanks
    # Wolfgang!
    proxy <- Sys.getenv('http_proxy')
    if (proxy == ""){
        proxy <- NULL
    }

    if (!dir.exists(cache)){
        dir.create(cache, recursive=TRUE)
    }

    ah <- AnnotationHub(hub=getAnnotationHubOption('URL'),
             cache=cache,
             proxy=proxy,
             localHub=FALSE)

    find.annotationhub.name <- function(species.name, override.code) { #autodetect ah names based on loaded database
        if (is.na(override.code)) {
        ah.query <- query(ah, "OrgDb")
        ah.query.speciesmatch <- grepl(paste("^", species.name, "$", sep=""), ah.query$species)
        ah.query.which <- which(ah.query.speciesmatch)
        stopifnot(length(ah.query.which) > 0) #require at least one match
        if (length(ah.query.which) > 1) { #warn of duplicate matches
           print("WARNING: found multiple candidate species in AnnotationHub: ");
           print(ah.query.speciesmatch)
        }
        names(ah.query)[ah.query.which[1]]
        } else {
        override.code
        }
    }
    annotation_key <- find.annotationhub.name(species, annotation_key_override)
    orgdb <- ah[[annotation_key]]
    return(orgdb)
}


#' Writes out original and split clusterprofiler results
#'
#' @param res clusterProfiler results
#' @param cprof.folder Directory in which to save files. Will be created if needed.
#' @param Label to use for the results. Will generate a filename
#' cprof.folder/label.txt and cprof.folder/label_split.txt
#'
#' @return List of the two files that are created on disk.
write.clusterprofiler.results <- function(res, cprof.folder, label){
    dir.create(cprof.folder, showWarnings=FALSE, recursive=TRUE)
    filename.orig <- file.path(cprof.folder, paste0(label, '.tsv'))
    write.table(res, file=filename.orig, sep='\t', quote=FALSE, row.names=FALSE)
    filename.split <- file.path(cprof.folder, paste0(label, '_split.tsv'))
    res.split <- split.clusterProfiler.results(res)
    write.table(res.split, file=filename.split, sep='\t', quote=FALSE, row.names=FALSE)
    return(list(orig=filename.orig, split=filename.split))
}

#' Split clusterProfiler output into one line per gene
#'
#' @param x Results from clusterProfiler. It is expected that the
#' clusterProfiler enrichment function was called with "readable=TRUE"
#'
#' @return data.frame with genes one per line instead of "/" separated in one
#' line. The rest of the original line is repeated for each gene.
#'
split.clusterProfiler.results <- function(x){
    df <- x@result
    # loop over all rows
    df.parse <- NULL
    for(k in 1:dim(df)[1]){
        g <- strsplit(as.character(df$geneID[k]), "/")[[1]]
        gg <- df[rep(k, length(g)),]
        gg$geneParse <- g
        if(is.null(df.parse)){
            df.parse <- gg
        } else {
            df.parse <- rbind(df.parse, gg)
        }
    }
    return(df.parse)
}

# ------------------------------------------------------------------------------
# NOTE: Following functions are customized from the package `enrichplot`:
# Source versions are 1.10.2 (installed in env-r) and 1.16.2 (the newest).
# References:
# - https://rdrr.io/bioc/enrichplot/src/R/emapplot.R (for v1.10.2)
# - https://github.com/YuLab-SMU/enrichplot/tree/master/R (for 1.16.2)
update_n <- function(x, showCategory) {
    if (!is.numeric(showCategory)) {
        if (inherits(x, 'list')) {
            showCategory <- showCategory[showCategory %in% names(x)]
        } else {
            showCategory <- intersect(showCategory, x$Description)
        }
        return(showCategory)
    }

    ## geneSets <- geneInCategory(x) ## use core gene for gsea result
    n <- showCategory
    if (inherits(x, 'list')) {
        nn <- length(x)
    } else {
        nn <- nrow(x)
    }
    if (nn < n) {
        n <- nn
    }
    return(n)
}

pairwise_termsim <- function(x, method = "JC", semData = NULL, showCategory = 200) {
    y <- as.data.frame(x)
    geneSets <- geneInCategory(x)
    n <- update_n(x, showCategory)
    if (n == 0) stop("no enriched term found...")
    if (is.numeric(n)) {
        y <- y[1:n, ]
    } else {
        y <- y[match(n, y$Description),]
        n <- length(n)
    }

    x@termsim <- get_similarity_matrix(y = y, geneSets = geneSets, method = method,
                semData = semData)
    x@method <- method
    return(x)
}

overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}


# ------------------------------------------------------------------------------


##' Get the similarity matrix
##'
##' @param y a data.frame of enrichment result
##' @param geneSets a list, the names of geneSets are term ids,
##' and every object is a vertor of genes
##' @param method method of calculating the similarity between nodes,
##' one of "Resnik", "Lin", "Rel", "Jiang" , "Wang"  and
##' "JC" (Jaccard similarity coefficient) methods
##' @param semData GOSemSimDATA object
##' @noRd
get_similarity_matrix <- function(y, geneSets, method, semData = NULL) {
    id <- y[, "ID"]
    geneSets <- geneSets[id]
    n <- nrow(y)
    y_id <- unlist(strsplit(y$ID[1], ":"))[1]
    ## Choose the method to calculate the similarity
    if (method == "JC") {
        w <- matrix(NA, nrow=n, ncol=n)
        colnames(w) <- rownames(w) <- y$Description
        for (i in seq_len(n-1)) {
            for (j in (i+1):n) {
                w[i,j] <- overlap_ratio(geneSets[id[i]], geneSets[id[j]])
            }
        }
        return(w)
    }

    if (y_id == "GO") {
        if(is.null(semData)) {
            stop("The semData parameter is missing,
                and it can be obtained through godata function in GOSemSim package.")
        }
        w <- GOSemSim::mgoSim(id, id, semData=semData, measure=method,
                              combine=NULL)
    }

    if (y_id == "DOID") w <- DOSE::doSim(id, id, measure=method)
    return(w)
}


#' Truncate the names of an enrichment results object
#'
#' @param obj DOSE::enrichResult object
#' @param truncate_to Max number of characters in a label
#'
#' @return enrichResult object with labels truncated
truncate <- function(obj, truncate_to){
  f <- function(x){
    if (nchar(x) > truncate_to){
      x <- paste0(substr(x, 1, truncate_to), '...')
    }
    return(x)
  }
  obj@result$Description <- sapply(obj@result$Description, f)
  return(obj)
}

