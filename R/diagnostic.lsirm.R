#' Diagnostic the result of LSIRM
#'
#' @description \code{diagnostic} is used to diagnostic the result of LSIRM.
#'
#' @param object object of class \code{lsirm}.
#' @param plot If \code{TRUE}, MCMC diagnostic plots are returned
#' @param draw.item Select items for diagnosis. Parameter "alpha", however, uses "second" as a default value because the first alpha has an estimation issue. Positions and item names are supported. For instance, if the name of item is "i1", "i2", "i3" and its positions is in order, the result of \code{beta = c("i1","i2","i3")} and \code{beta = c(1,2,3)} are equivalent.
#'
#' @return \code{diagnostic} returns plots for checking MCMC convergence for selected parameters.
#'
#' @examples
#' \donttest{
#' # generate example item response matrix
#' data     <- matrix(rbinom(500, size = 1, prob = 0.5), ncol=10, nrow=50)
#' lsirm_result <- lsirm(data ~ lsirm1pl(spikenslab = FALSE, fixed_gamma = FALSE))
#'
#' # 1PL model
#' diagnostic(lsirm_result, plot=TRUE)
#'
#' # 1PL model, multiple items
#' diagnostic(lsirm_result, plot=TRUE, draw.item=list(beta = c(1,2,3)))
#'
#' lsirm_result <- lsirm(data ~ lsirm2pl(spikenslab = FALSE, fixed_gamma = FALSE))
#'
#' # 2PL model
#' diagnostic(lsirm_result, plot=TRUE)
#' }

#' @export diagnostic
diagnostic <- function(object,
                       draw.item=list(beta="first",
                                      theta="first",
                                      alpha="second",
                                      zw.dist = matrix(c(1, 1), ncol = 2)),
                       gelman.diag = FALSE){
  UseMethod("diagnostic")
}

#' @export
diagnostic.lsirm <- function(object,
                             draw.item = list(beta="first",
                                              theta="first",
                                              alpha="second",
                                              zw.dist = matrix(c(1, 1), ncol = 2)),
                             gelman.diag = FALSE)
{

  orders = data.frame(idx = c(1:7),
                      param = c("beta", "theta", "gamma", "alpha", "sigma", "sigma_sd", "zw.dist"))
  which.draw = names(draw.item)
  porder = orders[orders$param %in% which.draw,]

  if(object$chains > 1){
    multi_chain = TRUE
    nmcmc <- dim(object[[1]]$w)[1]
    which.draw.temp = which.draw
  }else{
    multi_chain = FALSE
    nmcmc <- dim(object$w)[1]
    which.draw.temp = which.draw
  }
  if(multi_chain){

    if((is.null(object[[1]]$beta)) & ("beta" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object[[1]]$gamma)) & ("gamma" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object[[1]]$theta)) & ("theta" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object[[1]]$z) | is.null(object[[1]]$w)) & ("zw.dist" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object[[1]]$alpha)) & ("alpha" %in% which.draw)) stop("MCMC sample was not stored. The parameter alpha is only stored in method lsirm2pl")

  }else{
    if((is.null(object$beta)) & ("beta" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object$gamma)) & ("gamma" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object$theta)) & ("theta" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object$z) | is.null(object$w)) & ("zw.dist" %in% which.draw)) stop("MCMC sample was not stored")
    if((is.null(object$alpha)) & ("alpha" %in% which.draw)) stop("MCMC sample was not stored. The parameter alpha is only stored in method lsirm2pl")

  }

  if(is.null(draw.item$beta)&("beta" %in% which.draw)) draw.item$beta = "first"
  if(is.null(draw.item$theta)&("theta" %in% which.draw)) draw.item$theta = "first"
  if(is.null(draw.item$alpha)&("alpha" %in% which.draw)) draw.item$alpha = "second"
  if(is.null(draw.item$zw.dist)&("zw.dist" %in% which.draw)) draw.item$zw.dist = matrix(c(1, 1), ncol = 2)
  chain_list_all <- vector("list", length = object$chains)

  #### Beta  ------------
  if(((!is.null(object[[1]]$beta))|(!is.null(object$beta))) & ("beta" %in% which.draw)){

    if(multi_chain){
      if(is.null(colnames(object[[1]]$data))){
        colnames(object[[1]]$beta) = 1:ncol(object[[1]]$beta)
      }else{colnames(object[[1]]$beta) <- colnames(object[[1]]$data)}



      if((length(draw.item$beta) == 1) & (draw.item$beta[1] == "first")){
        draw.item.temp = colnames(object[[1]]$beta)[1]
        draw.item.num <- 1
      }else if(is.numeric(draw.item$beta)){
        draw.item.temp = colnames(object[[1]]$beta)[draw.item$beta]
        draw.item.num <- draw.item$beta
      }else{
        draw.item.temp = draw.item$beta
        draw.item.num = which(colnames(object[[1]]$beta) %in% draw.item.temp)
      }


      pnames <- c(paste('beta [', draw.item.temp, ']', sep = ''))
      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$beta[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{
      if(is.null(colnames(object$data))){
        colnames(object$beta) = 1:ncol(object$beta)
      }else{colnames(object$beta) <- colnames(object$data)}



      if((length(draw.item$beta) == 1) & (draw.item$beta[1] == "first")){
        draw.item.temp = colnames(object$beta)[1]
        draw.item.num <- 1
      }else if(is.numeric(draw.item$beta)){
        draw.item.temp = colnames(object$beta)[draw.item$beta]
        draw.item.num <- draw.item$beta
      }else{
        draw.item.temp = draw.item$beta
        draw.item.num = which(colnames(object$beta) %in% draw.item.temp)
      }


      pnames <- c(paste('beta [', draw.item.temp, ']', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$beta[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))
    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        data.frame(iteration = seq_along(chain_data[[j]]), value = chain_data[[j]], Chain = j)
      }))


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        # autocorrelation_plot <- ggplot(acf_combined,
        #                                aes(x = Lag, y = ACF,
        #                                    color = factor(Chain), group = Chain)) +
        #   geom_hline(aes(yintercept = 0)) +
        #   geom_segment(mapping = aes(xend = Lag, yend = 0)) +
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")

        pdf(file = NULL)
        # Call gelman.plot
        gelman_plot <- gelman.plot(chain_data)
        # Close the null device
        dev.off()



        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter

        # Create a data frame
        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )

        # Reshape for plotting
        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))
        # Check the data

        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +  # Specify line type and color in aes to merge legends
          # geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", size = 1) +
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position = c(1, 1),  # Place legend inside the plot at top-right
                legend.justification = c("right", "top")  # Anchor the legend at its top-right corner

          )


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "beta",]$idx == max(porder$idx)){
          if(((length(draw.item.num) > 1)&(i < length(draw.item.num))) | !is.null(which.draw.temp) ){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })



        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        # autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
        #   geom_line(color = "#268bd2") +  # Changed to line plot
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "beta",]$idx == max(porder$idx)){
          if((length(draw.item.num) > 1)&(i < length(draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }
    which.draw.temp <- setdiff(which.draw.temp, c("beta"))
    # coda::gelman.diag(mcmclist)
  } ## Beta

  ## Theta ----------
  if(((!is.null(object[[1]]$theta))|(!is.null(object$theta))) & ("theta" %in% which.draw)){

    if(multi_chain){


      if(is.null(rownames(object[[1]]$data))){
        colnames(object[[1]]$theta) = 1:ncol(object[[1]]$theta)
      }else{colnames(object[[1]]$theta) <- rownames(object[[1]]$data)}



      if((length(draw.item$theta) == 1) & (draw.item$theta[1] == "first")){
        draw.item.temp = colnames(object[[1]]$theta)[1]
        draw.item.num <- 1
      }else if(is.numeric(draw.item$theta)){
        draw.item.temp = colnames(object[[1]]$theta)[draw.item$theta]
        draw.item.num <- draw.item$theta
      }else{
        draw.item.temp = draw.item$theta
        draw.item.num = which(colnames(object[[1]]$theta) %in% draw.item.temp)
      }


      pnames <- c(paste('theta [', draw.item.temp, ']', sep = ''))
      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$theta[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{
      if(is.null(rownames(object$data))){
        colnames(object$theta) = 1:ncol(object$theta)
      }else{colnames(object$theta) <- rownames(object$data)}



      if((length(draw.item$theta) == 1) & (draw.item$theta[1] == "first")){
        draw.item.temp = colnames(object$theta)[1]
        draw.item.num <- 1
      }else if(is.numeric(draw.item$theta)){
        draw.item.temp = colnames(object$theta)[draw.item$theta]
        draw.item.num <- draw.item$theta
      }else{
        draw.item.temp = draw.item$theta
        draw.item.num = which(colnames(object$theta) %in% draw.item.temp)
      }


      pnames <- c(paste('theta [', draw.item.temp, ']', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$theta[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))

    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        data.frame(iteration = seq_along(chain_data[[j]]), value = chain_data[[j]], Chain = j)
      }))


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")
        # autocorrelation_plot <- ggplot(acf_combined,
        #                                aes(x = Lag, y = ACF,
        #                                    color = factor(Chain), group = Chain)) +
        #   geom_hline(aes(yintercept = 0)) +
        #   geom_segment(mapping = aes(xend = Lag, yend = 0)) +
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")

        pdf(file = NULL)
        # Call gelman.plot
        gelman_plot <- gelman.plot(chain_data)
        # Close the null device
        dev.off()



        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter

        # Create a data frame
        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )

        # Reshape for plotting
        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))
        # Check the data

        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +  # Specify line type and color in aes to merge legends
          # geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", size = 1) +
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position = c(1, 1),  # Place legend inside the plot at top-right
                legend.justification = c("right", "top")  # Anchor the legend at its top-right corner

          )


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "theta",]$idx == max(porder$idx)){
          if((length(draw.item.num) > 1)&(i < length(draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        # autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
        #   geom_line(color = "#268bd2") +  # Changed to line plot
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "theta",]$idx == max(porder$idx)){
          if((length(draw.item.num) > 1)&(i < length(draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }
    which.draw.temp <- setdiff(which.draw.temp, c("theta"))
  } ## Theta

  ## Gamma =====
  if(((!is.null(object[[1]]$gamma))|(!is.null(object$gamma))) & ("gamma" %in% which.draw)){

    if(multi_chain){

      draw.item.num <- 1
      pnames <- c(paste('gamma', sep = ''))

      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$gamma[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{


      draw.item.num <- 1
      pnames <- c(paste('gamma', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$gamma[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))
    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        data.frame(iteration = seq_along(chain_data[[j]]), value = chain_data[[j]], Chain = j)
      }))


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")
        # autocorrelation_plot <- ggplot(acf_combined,
        #                                aes(x = Lag, y = ACF,
        #                                    color = factor(Chain), group = Chain)) +
        #   geom_hline(aes(yintercept = 0)) +
        #   geom_segment(mapping = aes(xend = Lag, yend = 0)) +
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")

        pdf(file = NULL)
        # Call gelman.plot
        gelman_plot <- gelman.plot(chain_data)
        # Close the null device
        dev.off()



        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter

        # Create a data frame
        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )

        # Reshape for plotting
        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))
        # Check the data

        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +  # Specify line type and color in aes to merge legends
          # geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", size = 1) +
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position = c(1, 1),  # Place legend inside the plot at top-right
                legend.justification = c("right", "top")  # Anchor the legend at its top-right corner

          )


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "gamma",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        # autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
        #   geom_line(color = "#268bd2") +  # Changed to line plot
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "gamma",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }
    which.draw.temp <- setdiff(which.draw.temp, c("gamma"))
  } ## Gamma


  ## alpha -----
  if(((!is.null(object[[1]]$alpha))|(!is.null(object$alpha))) & ("alpha" %in% which.draw)){

    if(multi_chain){

      draw.item.num <- 1
      pnames <- c(paste('alpha', sep = ''))
      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$alpha[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{


      draw.item.num <- 1
      pnames <- c(paste('alpha', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$alpha[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))
    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        data.frame(iteration = seq_along(chain_data[[j]]), value = chain_data[[j]], Chain = j)
      }))


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")
        # autocorrelation_plot <- ggplot(acf_combined,
        #                                aes(x = Lag, y = ACF,
        #                                    color = factor(Chain), group = Chain)) +
        #   geom_hline(aes(yintercept = 0)) +
        #   geom_segment(mapping = aes(xend = Lag, yend = 0)) +
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")

        pdf(file = NULL)
        # Call gelman.plot
        gelman_plot <- gelman.plot(chain_data)
        # Close the null device
        dev.off()



        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter

        # Create a data frame
        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )

        # Reshape for plotting
        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))
        # Check the data

        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +  # Specify line type and color in aes to merge legends
          # geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", size = 1) +
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position = c(1, 1),  # Place legend inside the plot at top-right
                legend.justification = c("right", "top")  # Anchor the legend at its top-right corner

          )


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "alpha",]$idx == max(porder$idx)){
          if((length(draw.item.num) > 1)&(i < length(draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        # autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
        #   geom_line(color = "#268bd2") +  # Changed to line plot
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))
        if(porder[porder$param == "alpha",]$idx == max(porder$idx)){
          if((length(draw.item.num) > 1)&(i < length(draw.item.num))){
            readline(prompt = "Hit <Return> to see next plot")
          }
        }else{
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }

  } ## Alpha


  ## sigma --------
  if(((!is.null(object[[1]]$sigma))|(!is.null(object$sigma))) & ("sigma" %in% which.draw)){

    if(multi_chain){

      draw.item.num <- 1
      pnames <- c(paste('sigma', sep = ''))
      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$sigma[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{


      draw.item.num <- 1
      pnames <- c(paste('sigma', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$sigma[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))
    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        data.frame(iteration = seq_along(chain_data[[j]]), value = chain_data[[j]], Chain = j)
      }))


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")
        # autocorrelation_plot <- ggplot(acf_combined,
        #                                aes(x = Lag, y = ACF,
        #                                    color = factor(Chain), group = Chain)) +
        #   geom_hline(aes(yintercept = 0)) +
        #   geom_segment(mapping = aes(xend = Lag, yend = 0)) +
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")

        pdf(file = NULL)
        # Call gelman.plot
        gelman_plot <- gelman.plot(chain_data)
        # Close the null device
        dev.off()



        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter

        # Create a data frame
        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )

        # Reshape for plotting
        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))
        # Check the data

        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +  # Specify line type and color in aes to merge legends
          # geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", size = 1) +
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position = c(1, 1),  # Place legend inside the plot at top-right
                legend.justification = c("right", "top")  # Anchor the legend at its top-right corner

          )


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "sigma",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        # autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
        #   geom_line(color = "#268bd2") +  # Changed to line plot
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "sigma",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }
    which.draw.temp <- setdiff(which.draw.temp, c("sigma"))
  } ## Sigma

  ## Sigma Theta ---------
  if(((!is.null(object[[1]]$theta_sd))|(!is.null(object$theta_sd))) & ("theta_sd" %in% which.draw)){

    if(multi_chain){

      draw.item.num <- 1
      pnames <- c(paste('theta_sd', sep = ''))
      chains <- object$chains

      chain_list <- list()
      for(i in 1:chains){
        chain_list[[i]] <- matrix(c(object[[i]]$theta_sd[,draw.item.num]), ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{


      draw.item.num <- 1
      pnames <- c(paste('theta_sd', sep = ''))
      chains <- 1

      chain_list <- list()
      chain_list[[1]] <- matrix(c(object$theta_sd[,draw.item.num]), ncol = length(pnames),
                                dimnames= list(NULL,pnames))
    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:length(draw.item.num)){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        data.frame(iteration = seq_along(chain_data[[j]]), value = chain_data[[j]], Chain = j)
      }))


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")
        # autocorrelation_plot <- ggplot(acf_combined,
        #                                aes(x = Lag, y = ACF,
        #                                    color = factor(Chain), group = Chain)) +
        #   geom_hline(aes(yintercept = 0)) +
        #   geom_segment(mapping = aes(xend = Lag, yend = 0)) +
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")

        pdf(file = NULL)
        # Call gelman.plot
        gelman_plot <- gelman.plot(chain_data)
        # Close the null device
        dev.off()



        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter

        # Create a data frame
        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )

        # Reshape for plotting
        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))
        # Check the data

        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +  # Specify line type and color in aes to merge legends
          # geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", size = 1) +
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position = c(1, 1),  # Place legend inside the plot at top-right
                legend.justification = c("right", "top")  # Anchor the legend at its top-right corner

          )


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "theta_sd",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        # autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
        #   geom_line(color = "#268bd2") +  # Changed to line plot
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if(porder[porder$param == "theta_sd",]$idx != max(porder$idx)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }
    which.draw.temp <- setdiff(which.draw.temp, c("theta_sd"))
  } ## Sigma theta

  ## zw.dist ---------
  if(((!is.null(object[[1]]$z) & !is.null(object[[1]]$w))|
      ((!is.null(object$w)) & (!is.null(object$z)))) &
     ("zw.dist" %in% which.draw)){

    if(multi_chain){


      draw.item.temp <- apply(draw.item$zw.dist, 1, function(row) paste(row, collapse = "-"))
      draw.item.num <- length(draw.item.temp)

      pnames <- c(paste('zw.dist [', draw.item.temp, ']', sep = ''))
      chains <- object$chains

      chain_list <- list()
      nmcmc <- dim(object[[1]]$w)[1]

      for(i in 1:chains){
        dist_temp <- matrix(0, nrow = nmcmc, ncol = nrow(draw.item$zw.dist))
        for(j in 1:nrow(draw.item$zw.dist)){

          iz = draw.item$zw.dist[j,1]
          iw = draw.item$zw.dist[j,2]

          if(is.vector(object[[i]]$z[,j,])){
            dist_temp[,j] <- sum((object[[i]]$z[,iz,] - object[[i]]$w[,iw,])^2)
          }else{
            dist_temp[,j] <- apply(apply(object[[i]]$z[,iz,] - object[[i]]$w[,iw,], 1, function(x) x^2), 2, sum)
          }
        }
        chain_list[[i]] <- matrix(dist_temp, ncol = length(pnames),
                                  dimnames= list(NULL,pnames))
        if(length(chain_list_all[[i]]) == 0){
          chain_list_all[[i]] <- chain_list[[i]]
        }else{
          chain_list_all[[i]] <- cbind(chain_list_all[[i]], chain_list[[i]])
        }
      }
    }else{

      draw.item.temp <- apply(draw.item$zw.dist, 1, function(row) paste(row, collapse = "-"))
      draw.item.num <- length(draw.item.temp)

      pnames <- c(paste('zw.dist [', draw.item.temp, ']', sep = ''))

      chain_list <- list()
      nmcmc <- dim(object$w)[1]

      dist_temp <- matrix(0, nrow = nmcmc, ncol = nrow(draw.item$zw.dist))

      for(j in 1:nrow(draw.item$zw.dist)){

        iz = draw.item$zw.dist[j,1]
        iw = draw.item$zw.dist[j,2]

        if(is.vector(object$z[,j,])){
          dist_temp[,j] <- sum((object$z[,iz,] - object$w[,iw,])^2)
        }else{
          dist_temp[,j] <- apply(apply(object$z[,iz,] - object$w[,iw,], 1, function(x) x^2), 2, sum)
        }

      }
      chain_list[[1]] <- matrix(dist_temp, ncol = length(pnames),
                                dimnames= list(NULL,pnames))


    }

    # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
    mcmc_chains <- lapply(chain_list, function(x) coda::mcmc(x, thin = 1))
    # Convert the list of mcmc objects to an mcmc.list object
    mcmclist <- coda::as.mcmc.list(mcmc_chains)
    params <- dimnames(mcmclist[[1]])[[2]]

    plots_list <- list()

    for(i in 1:draw.item.num){
      # Extract data for each parameter from each chain
      param <- colnames(mcmclist[[1]])[i]
      chain_data <- lapply(mcmclist, function(chain) chain[, param])

      # Create data frame
      combined_data <- do.call(rbind, lapply(seq_along(chain_data), function(j) {
        data.frame(iteration = seq_along(chain_data[[j]]), value = chain_data[[j]], Chain = j)
      }))


      if(multi_chain){

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line() +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1, fill = factor(Chain), color = factor(Chain))) +
          geom_density(alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
          geom_line() +  # Changed to line plot
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")

        # autocorrelation_plot <- ggplot(acf_combined,
        #                                aes(x = Lag, y = ACF,
        #                                    color = factor(Chain), group = Chain)) +
        #   geom_hline(aes(yintercept = 0)) +
        #   geom_segment(mapping = aes(xend = Lag, yend = 0)) +
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")

        pdf(file = NULL)
        # Call gelman.plot
        gelman_plot <- gelman.plot(chain_data)
        # Close the null device
        dev.off()



        median_values <- gelman_plot$shrink[,, "median"]
        upper_values <- gelman_plot$shrink[,, "97.5%"]
        iterations <- gelman_plot$last.iter

        # Create a data frame
        psrf_data <- data.frame(
          Iteration = iterations,
          MedianPSRF = as.vector(median_values),
          UpperPSRF = as.vector(upper_values)
        )

        # Reshape for plotting
        psrf_long <- pivot_longer(psrf_data, cols = c("MedianPSRF", "UpperPSRF"), names_to = "Type", values_to = "PSRF") %>%
          mutate(Type = recode(Type, "UpperPSRF" = "97.5%", "MedianPSRF" = "Median"))
        # Check the data

        psrf_plot <- ggplot(psrf_long, aes(x = Iteration, y = PSRF, group = Type)) +
          geom_line(aes(linetype = Type, color = Type)) +  # Specify line type and color in aes to merge legends
          # geom_hline(yintercept = 1.1, linetype = "dashed", color = "black", size = 1) +
          scale_linetype_manual(values = c("97.5%" = "dashed", "Median" = "solid")) +
          scale_color_manual(values = c("97.5%" = "#377EB8", "Median" = "#E41A1C")) +
          labs(title = "", x = "Last Iteration in chain", y = "Shrink Factor") +
          theme_minimal() +
          theme(legend.title = element_blank(),
                legend.position = c(1, 1),  # Place legend inside the plot at top-right
                legend.justification = c("right", "top")  # Anchor the legend at its top-right corner

          )


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  psrf_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))

        if((draw.item.num > 1)&(i < draw.item.num)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }else{

        # Trace Plot
        trace_plot <- ggplot(combined_data, aes(x = iteration, y = var1, group = Chain, color = factor(Chain))) +
          geom_line(color = "#268bd2") +
          xlab("Iterations") +
          ylab("Value") + theme(legend.position = "none")

        # Density Plot
        density_plot <- ggplot(combined_data, aes(x = var1)) +
          geom_density(color = "#268bd2", fill = "#268bd2", alpha = 0.1) +
          xlab("Value") +
          ylab("Density") +
          theme(legend.position = "none")

        # Autocorrelation Plot
        autocorrelation_plots <- lapply(seq_along(chain_data), function(j) {
          acf_data <- acf(chain_data[[j]], plot = FALSE, lag.max = 60)
          data.frame(Lag = acf_data$lag, ACF = acf_data$acf, Chain = j)
        })

        # Merge autocorrelation data and plot
        acf_combined <- do.call(rbind, autocorrelation_plots)
        # autocorrelation_plot <- ggplot(acf_combined, aes(x = Lag, y = ACF, color = factor(Chain), group = Chain)) +
        #   geom_line(color = "#268bd2") +  # Changed to line plot
        #   xlab("Lag") +
        #   ylab("Autocorrelation") +
        #   theme(legend.position = "none")
        autocorrelation_plot <- ggplot(acf_combined,
                                       aes(x = Lag, y = ACF,
                                           color = factor(Chain), group = Chain)) +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = Lag, yend = 0)) +
          xlab("Lag") +
          geom_hline(aes(yintercept = (exp(2*1.96/sqrt(nmcmc-3)-1)/(exp(2*1.96/sqrt(nmcmc-3)+1)))), linetype='dotted', col = 'blue') +
          ylab("Autocorrelation") +
          theme(legend.position = "none")


        # Combine plots using grid.arrange
        combined_plot <- grid.arrange(arrangeGrob(trace_plot,
                                                  density_plot,
                                                  autocorrelation_plot,
                                                  top = textGrob(param, vjust = 1, gp = gpar(fontface = "bold", cex = 1)),
                                                  ncol = 2))


        if((draw.item.num > 1)&(i < draw.item.num)){
          readline(prompt = "Hit <Return> to see next plot")
        }
      }
    }

    # coda::gelman.diag(mcmclist)
    which.draw.temp <- setdiff(which.draw.temp, c("zw.dist"))
  } ## Dist

  if(gelman.diag == TRUE){
    if(object$chains > 1){
      # mcmc_chains <- lapply(chain_list, coda::mcmc, thin = 1)
      mcmc_chains <- lapply(chain_list_all, function(x) coda::mcmc(x, thin = 1))
      # Convert the list of mcmc objects to an mcmc.list object
      mcmclist <- coda::as.mcmc.list(mcmc_chains)
      print(coda::gelman.diag(mcmclist))
    }else{
      stop("Gelman and Rubin's convergence diagnostic requires at least two chains for computation.")
    }
  }


}


