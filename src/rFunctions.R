### Various R Functions ###
convertDF <- function(df, meta) {
    dfDat <- list(x=df, metaS=as.data.frame(lapply(meta[rownames(df), ], as.factor)))
    rownames(dfDat$metaS) <- rownames(meta)
    return(dfDat)
}

# This function calculates the bray-curtis dissimilarity matrix for a dataset and runs a PERMANOVA for between- and within-subject factors.
# The PERMANOVA does 10000 permutations. Here are good sources for learning more:
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/
# https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/restricting-permutations/
# https://stats.stackexchange.com/questions/590510/repeated-measures-permanova-nowhere-to-find
runPermStress <- function(df) {
    library('vegan')
    library('labdsv')
    set.seed(13)
    D <- dsvdis(df$x, index="bray/curtis")
    within_perm <- how(nperm=9999, plots = Plots(strata=df$metaS$Individual, type = "none"), within=Within(type='free'), observed=TRUE)
    bet_perm <- how(nperm=9999, plots=Plots(strata=df$metaS$Individual, type='free'), within=Within(type='none'), observed=TRUE)
    withinRes <- adonis2(D ~ Individual + Condition, data=df$metaS, permutation=within_perm)
    betRes <- adonis2(D ~ Sex + Cage + Individual, data=df$metaS, permutation=bet_perm)
    set.seed(Sys.time())
    return(list(betRes, withinRes))
}

runPermNoStress <- function(df) {
    library('vegan')
    library('labdsv')
    set.seed(13)
    D <- dsvdis(df$x, index="bray/curtis")
    within_perm <- how(nperm=9999, plots = Plots(strata=df$metaS$Individual, type = "none"), within=Within(type='free'), observed=TRUE)
    bet_perm <- how(nperm=9999, plots=Plots(strata=df$metaS$Individual, type='free'), within=Within(type='none'), observed=TRUE)
    withinRes <- adonis2(D ~ Individual + Relative_Day, data=df$metaS, permutation=within_perm)
    betRes <- adonis2(D ~ Sex + Cage + Individual, data=df$metaS, permutation=bet_perm)
    set.seed(Sys.time())
    return(list(betRes, withinRes))
}


# Calculates bray-curtis distance (optional) and performs a PCoA
PCoA <- function(dat, D=NA, k=2) {
    library('vegan')
    library('labdsv')
    x <- dat$x
    k <- max(k, 2)
    D <- dsvdis(x, index="bray/curtis")
    pc <- pco(D, k)
    ord <- pc
    toteig <- sum(pc$eig[pc$eig>0])
    ord$ordnames <- sprintf("PCo %d (%0.1f%%)", 1:k, 100 * pc$eig[1:k] / toteig)
    return (ord)
}

ordplot <- function(dat, ord, pcos = 2, pointoutline = T,
        colour = NA, colour_title = NA, colour_names = NA, colour_override = NA,
        shape = NA,  shape_title = NA,  shape_names = NA,  shape_override = NA,
        size = NA, size_title = NA, size_names = NA, size_override = NA,
        size_abs = NA, outline_size = NA,
        surface = NA, surf_maj_lines = 8, surf_min_lines = 4, surf_extent = 1,
        surf_smoothness = 1, surf_quality = 100,
        centroid = NA, loading = NA, sigloadings = NA,
        enriched = NA, text_halo = T,
        colour_log = F, colour_log_pseudocount = 0, colour_log_base = 10,
        arrows_fixed = NULL, arrows = NULL, arrow_text = T, arrow_sqrtnorm = T,
        connect = NA, sequence = NA, connectwidth = NA, sortby = NA, decreasing = F) {

    library('ggplot2')
    library('ggpubr')
    # ggplot2-based function to visualize the results of an ordination
    dat$meta <- dat$metaS
    
    

    # Set up the plot
    ggp <- ggplot() + theme_classic() +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
    gptopt <- list()
    gptaes <- list(x="dim1", y="dim2")

    # Get the dimensions to plot, and set the axes appropriately
    if (length(pcos) == 1) {
        pcos <- max(pcos,2)
        pcos <- c(pcos-1, pcos)
    } else {
        stopifnot(length(pcos) == 2)
    }
    stopifnot(max(pcos) <= dim(ord$points)[2])
    pts <- as.data.frame(ord$points[,pcos])
    if (is.null(ord$ordnames)) {
        dims <- colnames(pts)
    } else {
        dims <- ord$ordnames[pcos]
    }
    colnames(pts) <- c("dim1", "dim2")
    ggp <- ggp + xlab(dims[1]) + ylab(dims[2])

    # Use this mapping to allow the ordination to be done on a subset of the
    # data without requiring that dat also be that subset
    metamatch <- match(rownames(pts), rownames(dat$x))
    if (any(is.na(metamatch))) {
        pts <- pts[!is.na(metamatch),]
        metamatch <- metamatch[!is.na(metamatch)]
    }

    getmeta <- function(name) {
        tgt <- name == colnames(dat$meta)
        if (any(tgt)) {
            return (dat$meta[metamatch,tgt])
        }
        tgt <- name == colnames(dat$x)
        if (any(tgt)) {
            return (dat$x[metamatch,min(which(tgt))])
        }

        stop(sprintf("%s is not a metadata or feature.", name))
    }

    # Apply sorting
    if (!is.na(sortby)) {
        s <- getmeta(sortby)
        pts <- pts[order(s, decreasing=decreasing),]
        metamatch <- match(rownames(pts), rownames(dat$x))
    }

    # Draw connectors first
    if (!is.na(connect) || !is.na(sequence)) {
        df <- data.frame(x=pts$dim1, y=pts$dim2)

        if (!is.na(connect)) {
            connectby <- getmeta(connect)
            df$group <- connectby
        }
        if (!is.na(sequence)) {
            seqby <- getmeta(sequence)
            df <- df[order(seqby),]
        }
        if (!is.na(connect)) {
            ggp <- ggp + geom_path(data=df, aes(x=x, y=y, group=group), size=ifelse(is.na(connectwidth), 1, connectwidth))
        } else {
            ggp <- ggp + geom_path(data=df, aes(x=x, y=y), size=ifelse(is.na(connectwidth), 1, connectwidth))
        }
    }

    # Draw the surface
    if (!is.na(surface)) {
        surfby <- getmeta(surface)

        # pick an appropriate bandwidth
        bw1 <- surf_smoothness * 1.06 * sd(pts$dim1) * length(pts$dim1)^-0.2
        bw2 <- surf_smoothness * 1.06 * sd(pts$dim2) * length(pts$dim2)^-0.2

        library(signal) # for unwrap

        # only make the surface from finite values
        keep <- !is.na(surfby) & is.finite(surfby)
        surfby <- surfby[keep]
        xp <- pts$dim1[keep]
        yp <- pts$dim2[keep]

        # distance from the points to draw
        extent <- 0.5 * surf_extent^2 / surf_smoothness

        # surface from a Gaussian weighted mean
        gkmeansurf <- function(x, y) {
            dx <- x - xp
            dy <- y - yp
            chi2 <- (dx / bw1)^2 + (dy / bw2)^2
            kern <- exp(-chi2)

            theta <- unwrap(atan2(dy, dx))
            mnchi2 <- min(chi2)
            if (mnchi2 > extent && (max(theta) - min(theta)) < pi) { return (NA) }

            return (-sum(surfby * kern) / sum(kern))
        }

        # drop a grid on the ordination
        n <- surf_quality + 1
        border <- surf_extent
        xx <- matrix(seq(from=min(pts$dim1) - border * bw1,
                         to  =max(pts$dim1) + border * bw1,
                         length=n), n, n, byrow=F)
        yy <- matrix(seq(from=min(pts$dim2) - border * bw2,
                         to  =max(pts$dim2) + border * bw2,
                         length=n), n, n, byrow=T)

        # measure the surface at all points
        z <- mapply(gkmeansurf, as.vector(xx), as.vector(yy))

        # only keep points in range
        keep <- !is.na(z)

        # plot contours
        zrange <- max(z[keep]) - min(z[keep])
        cdata <- data.frame(x=as.vector(xx)[keep], y=as.vector(yy)[keep], z=z[keep])
        ggp <- ggp + stat_contour(aes(x=x, y=y, z=z), data=cdata,
               size=0.5, colour="grey50", binwidth=zrange/(surf_maj_lines*(surf_min_lines+1)))
        ggp <- ggp + stat_contour(aes(x=x, y=y, z=z), data=cdata,
               size=1, colour="black", binwidth=zrange/surf_maj_lines)
    }

    # Set up the plot to shape by some metadata
    if (!is.na(shape)) {
        shapeby <- getmeta(shape)
        if (length(shape_names) > 1 || !is.na(shape_names)) {
            # reorder the factors to match the order in shape_names
            # so that the legend is displayed in the correct order
            shapeby <- factor(levels(shapeby)[as.numeric(shapeby)],
                              levels=names(shape_names))
        }
        pts$shape <- shapeby
        gptaes$shape <- "shape"
        if (is.na(shape_title)) {
            shape_title <- shape
        }

        ggp <- ggp + guides(shape = guide_legend(title = shape_title,
                                                 override.aes = list(size=4, fill="gray")))

        if (length(shape_override) > 1 || !is.na(shape_override)) {
            # User-defined shapes
            if (length(shape_names) == 1 && is.na(shape_names)) {
                shape_names <- shape_override
                shape_names[names(shape_override)] <- names(shape_override)
            }
            ggp <- ggp + scale_shape_manual(values=unlist(shape_override),
                                            breaks=names(shape_names),
                                            labels=unlist(shape_names))
        } else if (pointoutline) {
            # Automatic shapes - select the "fillable" ones
            shps <- list()
            shapeby <- factor(shapeby)
            for (i in 1:length(levels(shapeby))) {
                shps[[levels(shapeby)[i]]] <- 20+i
            }
            ggp <- ggp + scale_shape_manual(values=unlist(shps))
        } else {
            # Automatic shapes - no outlines
            shps <- list()
            shapeby <- factor(shapeby)
            for (i in 1:length(levels(shapeby))) {
                shps[[levels(shapeby)[i]]] <- 15+i
            }
            ggp <- ggp + scale_shape_manual(values=unlist(shps))
        }

    } else if (length(shape_override) == 1 && !is.na(shape_override)) {
        gptopt$shape <- shape_override
    } else {
        gptopt$shape <- if (pointoutline) {21} else {19}
    }

    # Set up the plot to size by some metadata
    if (!is.na(size)) {
        sizeby <- getmeta(size)

        gptaes$size <- "size"
        if (is.na(size_title)) {
            size_title <- size
        }

        if (is.numeric(sizeby)) {
            if (!is.na(size_abs)) {
                sizeby <- sizeby * size_abs
            }
        } else if (length(size_override) > 1 || !is.na(size_override)) {
            # User-defined shapes
            if (length(size_names) == 1 && is.na(size_names)) {
                size_names <- size_override
                size_names[names(size_override)] <- names(size_override)
            }
            ggp <- ggp + scale_size_manual(values=unlist(size_override),
                                            breaks=names(size_names),
                                            labels=unlist(size_names))
        }

        pts$size <- sizeby
        ggp <- ggp + guides(size = guide_legend(title = size_title))
    } else {
        if (is.na(size_abs)) {
            gptopt$size <- min(6, 50 / sqrt(length(metamatch)))
        } else {
            gptopt$size <- size_abs
        }
    }

    # Set up the plot to colour based on metadata
    if (pointoutline) {
        ggscale_manual <- scale_fill_manual
        ggscale_gradientn <- scale_fill_gradientn
        ggscale_brewer <- scale_fill_brewer
        ggcolorfield <- "fill"
        if (!is.na(outline_size)) {
            gptopt[["stroke"]] <- outline_size
        }
    } else {
        ggscale_manual <- scale_colour_manual
        ggscale_gradientn <- scale_colour_gradientn
        ggscale_brewer <- scale_colour_brewer
        ggcolorfield <- "colour"
    }
    if (!is.na(colour)) {
        colby <- getmeta(colour)
        if (length(colour_names) > 1 || !is.na(colour_names)) {
            # reorder the factors to match the order in colour_names
            # so that the legend is displayed in the correct order
            colby <- factor(as.character(colby),
                            levels=names(colour_names))
        }
        if (is.numeric(colby) && colour_log) {
            colby <- log(colby + colour_log_pseudocount, base=colour_log_base)
        }
        pts$colour <- colby
        gptaes[[ggcolorfield]] <- "colour"

        if (is.na(colour_title)) {
            colour_title <- colour
        }

        if (length(colour_override) > 1 || !is.na(colour_override)) {
            # Manual colouring
            if (length(colour_names) == 1 && is.na(colour_names)) {
                colour_names <- colour_override
                colour_names[names(colour_override)] <- names(colour_override)
            }
            if (is.numeric(colby)) {
                ggp <- ggp + ggscale_gradientn(colours = colour_override, na.value="grey80")
            } else {
                ggp <- ggp + ggscale_manual(values=unlist(colour_override),
                                            breaks=names(colour_names),
                                            labels=unlist(colour_names))
            }
            col_guide <- guide_legend(title = colour_title,
                                      override.aes = list(size=3, shape=21))
        } else if (is.numeric(colby)) {
            # Continuous data uses the RdYlGn palette
            library(viridis)
            ggp <- ggp + ggscale_gradientn(colours = viridis(100), na.value="grey80")
            col_guide <- guide_colourbar(title = colour_title)
        } else if (length(levels(colby)) > 9) {
            library(RColorBrewer)
            ggp <- ggp + ggscale_manual(values = colorRampPalette(
                rev(brewer.pal(n = 7, name = "RdYlBu")))(length(levels(colby))))
            col_guide <- guide_legend(title = colour_title,
                                    override.aes = list(size=3, shape=21))
        } else {
            # Discrete data uses the Set1 palette
            ggp <- ggp + scale_fill_brewer(palette = "Set1")
            col_guide <- guide_legend(title = colour_title,
                                override.aes = list(size=3, shape=21))
        }
        if (pointoutline) {
            ggp <- ggp + guides(fill = col_guide)
        } else {
            ggp <- ggp + guides(colour = col_guide)
        }
    } else if (length(colour_override) == 1 && !is.na(colour_override)) {
        gptopt[[ggcolorfield]] <- colour_override
    } else {
        gptopt[[ggcolorfield]] <- "darkslategray4"
    }

    # Draw the points
    gptopt$data <- pts
    gptopt <- c(list(do.call(aes_string, gptaes)), gptopt)
    ggp <- ggp + do.call(geom_point, gptopt)

    # Helper to add text annotations
    xr <- 0.0018
    dx <- xr * (max(pts[,1]) - min(pts[,1]))
    dy <- 1.1 * xr * (max(pts[,2]) - min(pts[,2]))
    halo_quality <- 12
    put_text <- function(x, y, text, ...) {
        if (text_halo) {
            phi <- 2 * pi / halo_quality
            for (i in 1:halo_quality) {
                ggp <<- ggp + annotate("text", label=text, color="white",
                                       x = x + dx*cos((i - 0.5) * phi),
                                       y = y + dy*sin((i - 0.5) * phi), ...)
            }
        }
        ggp <<- ggp + annotate("text", label=text, x=x, y=y, ..., color="black")
    }

    # Add metadata centroids
    if (length(centroid) > 1 || !is.na(centroid)) {
        for (i in seq_along(centroid)) {
            meta <- getmeta(centroid[i])

            if (is.numeric(meta)) {
                d1 <- sum(pts$dim1 * meta) / sum(meta)
                d2 <- sum(pts$dim2 * meta) / sum(meta)
                put_text(d1, d2, centroid[i])
            } else {
                meta <- factor(meta) # remove unused
                for (j in seq_along(levels(meta))) {
                    d1 <- mean(pts$dim1[meta == levels(meta)[j]])
                    d2 <- mean(pts$dim2[meta == levels(meta)[j]])
                    put_text(d1, d2, levels(meta)[j])
                }
            }
        }
    }

    # Add metadata enrichment modes
    if (length(enriched) > 1 || !is.na(enriched)) {
        library("neldermead")
        for (i in seq_along(enriched)) {
            meta <- getmeta(enriched[i])

            if (!is.numeric(meta)) {
                stop(sprintf("%s must be numeric to find its enrichment mode"))
            }

            # select appropriate bandwidths
            bw1 <- 1.06 * sd(pts$dim1) * length(pts$dim1)^-0.2
            bw2 <- 1.06 * sd(pts$dim2) * length(pts$dim2)^-0.2

            # only use finite values
            keep <- !is.na(meta) & is.finite(meta)
            meta <- meta[keep]
            xp <- pts$dim1[keep]
            yp <- pts$dim2[keep]

            enr <- function(xy) {
                chi2 <- ((xy[1] - xp) / bw1)^2 + ((xy[2] - yp) / bw2)^2
                kern <- exp(-chi2)
                nsamps <- sum(kern)
                modulator <- nsamps / (nsamps + 1)
                return (modulator * -sum(meta * kern) / sum(kern))
            }

            bestxy <- 0
            bestnenr <- Inf
            for (x0 in seq(from=min(pts$dim1), to=max(pts$dim1), length=5)) {
                for (y0 in seq(from=min(pts$dim2), to=max(pts$dim2), length=5)) {
                    res <- fminsearch(enr, c(x0, y0))
                    if (neldermead.get(res, "fopt") < bestnenr) {
                        bestxy <- neldermead.get(res, "xopt")
                        bestnenr <- neldermead.get(res, "fopt")
                    }
                }
            }

            put_text(bestxy[1], bestxy[2], enriched[i])
        }
    }

    # Add arrows
    if (!is.null(arrows) || !is.null(arrows_fixed)) {
        if (is.null(arrows_fixed)) {
            arr_pcos <- matrix(0, 0, 2, dimnames=list(c(), dims))
        } else {
            arr_pcos <- arrows_fixed[,pcos]
        }
        if (!is.null(arrows)) {
            if (is.numeric(arrows)) {
                arrSet <- c()
                feats <- c(colnames(dat$x))
                for (i in seq_along(feats)) {
                    meta <- getmeta(feats[i])

                    if (is.numeric(meta)) {
                        d1 <- sum(pts$dim1 * meta) / sum(meta)
                        d2 <- sum(pts$dim2 * meta) / sum(meta)
                        arr <- matrix(c(d1, d2), 1, 2)
                        rownames(arr) <- feats[i]
                        arrSet <- rbind(arrSet, arr)
                    } else {
                        meta <- factor(meta) # remove unused
                        for (j in seq_along(levels(meta))) {
                            d1 <- mean(pts$dim1[meta == levels(meta)[j]])
                            d2 <- mean(pts$dim2[meta == levels(meta)[j]])
                            arr <- matrix(c(d1, d2), 1, 2)
                            rownames(arr) <- levels(meta)[j]
                            arrSet <- rbind(arrSet, arr)
                        }
                    }
                }
                arrSz <- arrSet[,1]^2 + arrSet[,2]^2
                I <- order(-arrSz)
                arr_pcos <- rbind(arr_pcos, arrSet[I[1:arrows],])
            } else {
                for (i in seq_along(arrows)) {
                    meta <- getmeta(arrows[i])

                    if (is.numeric(meta)) {
                        d1 <- sum(pts$dim1 * meta) / sum(meta)
                        d2 <- sum(pts$dim2 * meta) / sum(meta)
                        arr <- matrix(c(d1, d2), 1, 2)
                        rownames(arr) <- arrows[i]
                        arr_pcos <- rbind(arr_pcos, arr)
                    } else {
                        meta <- factor(meta) # remove unused
                        for (j in seq_along(levels(meta))) {
                            d1 <- mean(pts$dim1[meta == levels(meta)[j]])
                            d2 <- mean(pts$dim2[meta == levels(meta)[j]])
                            arr <- matrix(c(d1, d2), 1, 2)
                            rownames(arr) <- levels(meta)[j]
                            arr_pcos <- rbind(arr_pcos, arr)
                        }
                    }
                }
            }
        }

        ptsRS <- rowSums(pts[,1:2]^2)
        maxRadius <- sqrt(max(ptsRS[is.finite(ptsRS)]))
        maxArrRadius <- sqrt(max(rowSums(arr_pcos^2)))
        arr_pcos <- arr_pcos * (0.85 * maxRadius / maxArrRadius)
        if (arrow_sqrtnorm) {
            arr_len <- sqrt(rowSums(arr_pcos^2))
            maxarr_len <- max(arr_len)
            narr_len <- maxarr_len * sqrt(arr_len / maxarr_len)
            arr_pcos <- arr_pcos * (narr_len / arr_len)
        }
        colnames(arr_pcos) <- paste("dim", 1:length(pcos), sep="")

        ggp <- ggp + geom_segment(data=as.data.frame(arr_pcos),
                                  aes(x=0, y=0, xend=dim1, yend=dim2),
                                  colour="red", arrow=arrow())

        if (arrow_text) {
            for (i in seq_len(nrow(arr_pcos))) {
                rr <- arr_pcos[i,]
                hjust <- 0.5
                vjust <- 0.5
                if (abs(rr[1]) > abs(rr[2])) {
                    hjust <- if (rr[1] > 0) {0} else {1}
                } else {
                    vjust <- if (rr[2] > 0) {0} else {1}
                }
                tc <- rr * (1 + 0.015 * maxRadius / sqrt(sum(rr^2)))
                put_text(tc[1], tc[2], rownames(arr_pcos)[i],
                         hjust=hjust, vjust=vjust, color='red')
            }
        }
    }

    return (ggp)
}

# Calculates the log2 fold change given a row of data and a row of metadata
# by subtracting the mean of the first factor (assumed to be a reference or "control")
# from the means of the other factors. Required data to be log2 transformed beforehand.
calculateLog2FC <- function(y, factor) {
    yMeans <- tapply(y, factor, mean)
    return(yMeans[2:dim(yMeans)] - yMeans[1])
}

# This function fits a linear mixed model to the data and extracts the pearson residuals for downstream analysis.
runMixedModel <- function(df) {
    library('lme4')
    library('lmerTest')
    
    # Residual matrix is created
    resid <- matrix(NA, nrow(df$x), ncol(df$x))
    rownames(resid) <- rownames(df$x)
    colnames(resid) <- colnames(df$x)
    
    # A dataframe of p-value is created
    sig <- data.frame('Feature_Name'=colnames(df$x))
    sig$pvaluePre <- NA
    sig$pvalueStress <- NA
    sig$pvaluePost <- NA
    sig$log2FCPre <- NA
    sig$log2FCStress <- NA
    sig$log2FCPost <- NA

    for (i in seq_along(colnames(df$x))) {
        y <- df$x[rownames(df$metaS),i]
        groupFC <- calculateLog2FC(y, df$metaS$Condition)
        mdl <- suppressMessages(lmer((y ~ Condition + Sex + (1 | Individual)), data=df$metaS, REML=TRUE))

        sig$pvaluePre[i] <- summary(mdl)$coefficients['ConditionPre', 'Pr(>|t|)']
        sig$pvalueStress[i] <- summary(mdl)$coefficients['ConditionStress', 'Pr(>|t|)']
        sig$pvaluePost[i] <- summary(mdl)$coefficients['ConditionPost', 'Pr(>|t|)']
        sig$log2FCPre[i] <- groupFC['Pre']
        sig$log2FCStress[i] <- groupFC['Stress']
        sig$log2FCPost[i] <- groupFC['Post']
        
        # Save full residuals for downstream analysis
        r <- residuals(mdl, type="pearson")
        resid[names(r),i] <- r

    }
    sig$qvaluePre <- p.adjust(sig$pvaluePre, 'fdr')
    sig$qvalueStress <- p.adjust(sig$pvalueStress, 'fdr')
    sig$qvaluePost <- p.adjust(sig$pvaluePost, 'fdr')
    rownames(sig) <- sig$Feature_Name
    sig <- sig[,-1]
    resid <- as.data.frame(resid)
    result <- list(sig=sig, res=resid)
    return(result)
}

# This is essentially the same function as the runMixedModel function, but it only uses a normal fixed-effects model
# (because the number of replicates per individual is too small to run a mixed model)
# and it extracts residuals without differential feature testing
runFixedModel <- function(df) {
    library('lme4')
    # Allocate output matrix
    resid <- matrix(NA, nrow(df$x), ncol(df$x))
    rownames(resid) <- rownames(df$x)
    colnames(resid) <- colnames(df$x)
    
    for (i in seq_along(colnames(df$x))) {
        y <- df$x[rownames(df$metaS),i]

        mdl <- lm((y ~ Individual + Sex + Cage), data=df$metaS)

        # Save full residuals for downstream analysis
        r <- residuals(mdl, type="pearson")
        resid[names(r),i] <- r

    }
    resid <- as.data.frame(resid)
    return(resid)
}
