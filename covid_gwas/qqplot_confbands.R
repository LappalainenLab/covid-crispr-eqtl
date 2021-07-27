#------------------------------------------------
# Function for QQ plot with confidence intervals
#------------------------------------------------

# Q-Q plot with confidence bands ----------------

# library(ggplot2)
# library(ggrastr)

qqplot_confbands <- function(pval,
                             confidence = 0.95, confbands = TRUE,
                             main = "Q-Q plot", color_group_name = "Trait 1",
                             add_p = NULL, add_group = NULL,
                             point_size = 2.5,
                             raster_plot = FALSE,
                             show_legend = FALSE) {
  #-----------------
  # add_p- NULL or list() of p-value vectors to add to the plot, add_p = list(x = sel_crispr_hits$all_inv_var_meta_p, y = random$all_inv_var_meta_p)
  # add_group - NULL or list() of group names, add_group = list(x = "CRISPR", y = "Random")
  # raster_plot - uses ggrastr::geom_point_rast to create rasterized plots and keep text as vectors
  #----------------

  # Exclude missing values
  if (sum(is.na(pval)) > 0) {
    cat("Exclude missing P-values - ", sum(is.na(pval)), fill = TRUE)
    pval <- pval[!is.na(pval)]
  }

  alpha <- 1 - confidence
  n <- length(pval)
  pval <- sort(pval)
  k <- c(1:n)

  lower <- qbeta(alpha/2, k, n + 1 - k)
  upper <- qbeta((1 - alpha/2), k, n + 1 - k)
  expected <- (1:n - 0.5)/n

  max_p <- ceiling(max(-log10(pval), -log10(expected), -log10(lower), -log10(upper)))
  if (!is.null(add_p)) {
    max_p <- max(max_p, ceiling(sapply(add_p, function(x) max(-log10(x), na.rm = T))))
  }

  # Prepare dataframe for plotting
  df <- data.frame(expected = -log10(expected),
                   observed = -log10(pval),
                   lower = -log10(lower),
                   upper = -log10(upper))

  # Plot Q-Q plot
  g <- ggplot(data = df, aes(x = expected, y = observed)) +
      labs(title = main,
           x = bquote(-log[10] ~ "(expected" ~ italic(P) ~ "-value)"),
           y = bquote(-log[10] ~ "(observed" ~ italic(P) ~ "-value)")) +
      geom_abline(slope = 1, intercept = 0) +
      scale_x_continuous(limits = c(0, ceiling(max(-log10(df$expected))))) +
      scale_y_continuous(limits = c(0, max_p))

  if (raster_plot) {
    g <- g +
      geom_point_rast(aes(color = color_group_name), pch = 21, size = point_size)
  } else {
    g <- g +
      geom_point(aes(color = color_group_name), pch = 21, size = point_size)
  }

  if (confbands) {
    g <- g +
      geom_line(data = df, aes(x = expected, y = lower), col = "grey65", lty = 2) +
      geom_line(data = df, aes(x = expected, y = upper), col = "grey65", lty = 2)
  }

  if (!is.null(add_p)) {
    show_legend <- TRUE # show legend, otherwise keep `show_legend = FALSE`
    l <- length(add_p)
    for (i in 1:l) {
      obs <- add_p[[i]]
      if (sum(is.na(obs)) > 0) {
        cat("Exclude missing P-values in `add_p` element", i, "-", sum(is.na(obs)), fill = TRUE)
        obs <- obs[!is.na(obs)]
      }

      # Prepare dataframe for plotting
      n <- length(obs)
      expected <- (1:n - 0.5)/n
      df_add <- data.frame(expected = -log10(expected),
                           observed = -log10(sort(obs)))

      if (i == 1) {
        add_group_name_1 <- ifelse(is.null(add_group[[i]]), paste0("Trait", i + 1), add_group[[i]])
        if (raster_plot) {
          g <- g +
            geom_point_rast(data = df_add, aes(x = expected, y = observed, color = add_group_name_1), pch = 21, size = point_size)
        } else {
          g <- g +
            geom_point(data = df_add, aes(x = expected, y = observed, color = add_group_name_1), pch = 21, size = point_size)
        }
      }
      if (i == 2) {
        add_group_name_2 <- ifelse(is.null(add_group[[i]]), paste0("Trait", i + 1), add_group[[i]])
        if (raster_plot) {
          g <- g +
            geom_point_rast(data = df_add, aes(x = expected, y = observed, color = add_group_name_2), pch = 21, size = point_size)
        } else {
          g <- g +
            geom_point(data = df_add, aes(x = expected, y = observed, color = add_group_name_2), pch = 21, size = point_size)
        }
      }
      if (i == 3) {
        add_group_name_3 <- ifelse(is.null(add_group[[i]]), paste0("Trait", i + 1), add_group[[i]])
        if (raster_plot) {
          g <- g +
            geom_point_rast(data = df_add, aes(x = expected, y = observed, color = add_group_name_3), pch = 21, size = point_size)
        } else {
          g <- g +
            geom_point(data = df_add, aes(x = expected, y = observed, color = add_group_name_3), pch = 21, size = point_size)
        }
      }
      if (i == 4) {
        add_group_name_4 <- ifelse(is.null(add_group[[i]]), paste0("Trait", i + 1), add_group[[i]])
        if (raster_plot) {
          g <- g +
            geom_point_rast(data = df_add, aes(x = expected, y = observed, color = add_group_name_4), pch = 21, size = point_size)
        } else {
          g <- g +
            geom_point(data = df_add, aes(x = expected, y = observed, color = add_group_name_4), pch = 21, size = point_size)
        }
      }
      if (i == 5) {
        add_group_name_5 <- ifelse(is.null(add_group[[i]]), paste0("Trait", i + 1), add_group[[i]])
        if (raster_plot) {
          g <- g +
            geom_point_rast(data = df_add, aes(x = expected, y = observed, color = add_group_name_5), pch = 21, size = point_size)
        } else {
          g <- g +
            geom_point(data = df_add, aes(x = expected, y = observed, color = add_group_name_5), pch = 21, size = point_size)
        }
      }
      if (i > 5) {
        warning("Up to 5 different groups allowed! Modify the function to add more!")
      }
    }
  }
  if (!show_legend) {
    g <- g +
      scale_color_manual(values = "royalblue") +
      theme(legend.position = "none")
  }
  return(g)
}

# Calculate lambda for different chi-square distribution quantiles -------------

lambda_quantile <- function(pval, q = 0.5) {
  # Exclude missing values
  if (sum(is.na(pval)) > 0) {
    cat("Exclude missing P-values - ", sum(is.na(pval)), fill = TRUE)
    pval <- pval[!is.na(pval)]
  }

  l <- quantile(qchisq(pval, df = 1, lower.tail = FALSE), probs = q) / qchisq(q, df = 1)
  return(l)
}
