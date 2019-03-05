library(tidyverse)
library(plotly)
library(RColorBrewer)

vline <- function(x = 0, color = "red") {
    list(
        type = "line", 
        y0 = 0, 
        y1 = 1, 
        yref = "paper",
        x0 = x, 
        x1 = x, 
        line = list(color = color)
    )
}

build_graph <- function(
    x, 
    fontsize, 
    colorpalette, 
    pvaluecolumn,
    termnamecolumn,
    geneoverlapcolumn,
    vcut = NULL
) {
    if (vcut == Inf) {
        vcut = NULL
    }
    x[["p.value"]] <- x[[pvaluecolumn]]
    x[["term.name"]] <- x[[termnamecolumn]]
    x[["intersection"]] <- x[[geneoverlapcolumn]]
    p <- plot_ly(
        data = x,
        y = ~ reorder(term.name, -p.value), x = ~ -log10(p.value), type = "bar",
        hoverinfo = "text",
        text = ~ paste(term.name, "( P-Value =", p.value, ")"),
        textfont = list(color = "black"),
        color = ~ domain, colors = colorpalette, showlegend = TRUE
    ) %>%
        add_trace(
            data = x,
            y = ~ reorder(term.name, -p.value), x = 0, type = "scatter",
            mode = "text", 
            text = x$term.name, textposition = "right", showlegend = FALSE
        ) %>%
        layout(
            font = list(size = fontsize),
            xaxis = list(
                title = "-Log10(FDR)",
                titlefont = list(
                    size = fontsize
                ),
                # tickfont = list(
                #     size = fontsize
                # ),
                automargin = TRUE
            ),
            yaxis = list(
                title = "",
                # tickfont = list(
                #     size = fontsize
                # ),
                showticklabels = FALSE,
                automargin = TRUE
            )
        )
    if (!is.null(vcut)) {
        p <- p %>%
            layout(
                shapes = list(vline(vcut))
            )
    }
    return(p)
}


build_graph_gg <- function(
    x, 
    fontsize_labels,
    fontsize_axes, 
    colorpalette, 
    pvaluecolumn,
    termnamecolumn,
    geneoverlapcolumn,
    pncolumn = NULL, 
    pflag = NULL,
    nflag = NULL,
    vcut = NULL
) {
    if (vcut == Inf) {
        vcut = NULL
    }
    if (!is.null(pncolumn)) {
        if (sum(is.null(c(pflag, nflag))) == 1 | sum(c(pflag, nflag) == "") == 1) {
            stop("Both +/- flags must be specified")
        }
    }
    x[["p.value"]] <- x[[pvaluecolumn]]
    x[["term.name"]] <- x[[termnamecolumn]]
    x[["intersection"]] <- x[[geneoverlapcolumn]]
    if (!is.null(pncolumn)) {
        x[["pnflag"]] <- sapply(x[[pncolumn]], function(p) {
            if (p == pflag) { return(1) }
            else if (p == nflag) { return(-1) }
            else { return(NA) }
        })
    } else {
        x[["pnflag"]] <- 1
    }
    x[["-Log10(P-Value)"]] <- (-log10(x[[pvaluecolumn]]) * x[["pnflag"]])
    p <- ggplot(
        data = x,
        aes(
            x = reorder(term.name, `-Log10(P-Value)`), 
            y = `-Log10(P-Value)`, 
            fill = domain, label = term.name,
            hjust = ifelse(pnflag == 1, 0, 1)
        )
    ) +
        coord_flip() +
        geom_bar(stat = "identity") +
        geom_text(y = 1 * x[["pnflag"]], size = (fontsize_labels), nudge_y = 1) +
        scale_fill_manual(values = colorpalette) +
        theme_bw() +
        theme(
            text = element_text(size = fontsize_axes),
            axis.title.y = element_blank(),
            axis.text.y = element_blank()
        )
    
    if (!is.null(vcut)) {
        p <- p + geom_hline(yintercept = vcut, colour = "red")
        if (!is.null(pncolumn)) {
            p <- p + geom_hline(yintercept = -vcut, colour = "red")
        }
    }
    return(p)
}