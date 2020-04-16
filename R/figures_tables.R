library(fishflux)
library(purrr)
library(png)
library(grid)
library(fishualize)
library(ggplot2)
library(egg)
library(dplyr)
library(tidybayes)
library(cowplot)

##### table ingestion rates #####
make_table_ingestion <- function(out){

  ##data
  
  param <- read.csv("data/param_modelsp.csv")
  param_bu <- as.list(param[param$species=="Balistapus undulatus",1:37])
  param_zs <- as.list(param[param$species=="Zebrasoma scopas",1:37])
  param_em <- as.list(param[param$species=="Epinephelus merra",1:37])
  
  ## run model
  mod_bu <- fishflux::cnp_model_mcmc(2:25, param = param_bu, iter = 5000,  
                                     cor = list(ro_Qc_Qn = 0.5, ro_Qc_Qp = -0.3, ro_Qn_Qp = -0.2,
                                                ro_Dc_Dn = 0.8, ro_Dc_Dp = -0.3, ro_Dn_Dp = -0.3,
                                                ro_lwa_lwb = 0.9, ro_alpha_f0 = 0.9))
  mod_zs <- fishflux::cnp_model_mcmc(2:25, param = param_zs, iter = 5000)
  mod_em <- fishflux::cnp_model_mcmc(2:25, param = param_em, iter = 5000,  
                                     cor = list(ro_Qc_Qn = 0.5, ro_Qc_Qp = -0.3, ro_Qn_Qp = -0.2,
                                                ro_Dc_Dn = 0.8, ro_Dc_Dp = -0.3, ro_Dn_Dp = -0.3,
                                                ro_lwa_lwb = 0.9, ro_alpha_f0 = 0.9))
  extract_ing <- function(sp, mod){
    e <- fishflux::extract(mod, par = c("l1","w1","IN"))
    e$inn <- e$IN_median/e$w1_median
    e$inn1 <- e$`IN_2.5%`/e$w1_median
    e$inn2 <- e$`IN_97.5%`/e$w1_median
    
    return(data.frame(
      species = sp,
      TL = e$TL,
      weight = e$w1_median,
      w1 = e$`w1_2.5%`,
      w2 = e$`w1_97.5%`,
      ing_m = e$inn * 1000,
      ing_1 = e$inn1 * 1000,
      ing_2 = e$inn2 * 1000
    ))
  }
  
  e_zs <- extract_ing("Z. scopas", mod_zs)
  e_bu <- extract_ing("B. undulatus", mod_bu)
  e_em <- extract_ing("E. merra", mod_em)
  e <- bind_rows(e_zs, e_bu, e_em)
  
  e[,3:8] <- apply(e[,3:8], 1:2, round, digits = 1)
  
  result <- 
    data.frame(
      Species = e$species,
      TL = e$TL,
      Biomass = paste(e$weight, " (", e$w1, "-", e$w2, ")", sep = ""),
      Ingestion = paste(e$ing_m, " (", e$ing_1, "-", e$ing_2, ")", sep = ""))
  
  write.csv(result, out, row.names = FALSE)
}


#### Figure 2 ####

makeplot2 <- function(out){

  ## fish 
  img_bu <- readPNG("figures/Bundulatus_grey.png")
  g_bu <- rasterGrob(img_bu, interpolate=TRUE)
  img_em <- readPNG("figures/Emerra_grey.png")
  g_em <- rasterGrob(img_em, interpolate=TRUE)
  img_zs <- readPNG("figures/Zscopas_grey.png")
  g_zs <- rasterGrob(img_zs, interpolate=TRUE)
  
  
  ## data
  
  param <- read.csv("data/param_modelsp.csv")
  param_bu <- as.list(param[param$species == "Balistapus undulatus",1:37])
  param_zs <- as.list(param[param$species == "Zebrasoma scopas",1:37])
  param_em <- as.list(param[param$species == "Epinephelus merra",1:37])
  
  ## run model
  mod_bu <- fishflux::cnp_model_mcmc(2:25, param = param_bu, iter = 5000,  cor = list(ro_Qc_Qn = 0.5, ro_Qc_Qp = -0.3, ro_Qn_Qp = -0.2,
                                                                                      ro_Dc_Dn = 0.8, ro_Dc_Dp = -0.3, ro_Dn_Dp = -0.3,
                                                                                      ro_lwa_lwb = 0.9, ro_alpha_f0 = 0.9))
  mod_zs <- fishflux::cnp_model_mcmc(2:25, param = param_zs, iter = 5000)
  mod_em <- fishflux::cnp_model_mcmc(2:25, param = param_em, iter = 5000,  cor = list(ro_Qc_Qn = 0.5, ro_Qc_Qp = -0.3, ro_Qn_Qp = -0.2,
                                                                                      ro_Dc_Dn = 0.8, ro_Dc_Dp = -0.3, ro_Dn_Dp = -0.3,
                                                                                      ro_lwa_lwb = 0.9, ro_alpha_f0 = 0.9))
  
  plot_lim <- function(mod){
    
    lim <- lapply(mod$stanfit, function(x){
      ee <- rstan::extract(x,"lim")[[1]]
      c <- length(which(ee==1))/length(ee)
      n <- length(which(ee==2))/length(ee)
      p <- length(which(ee==3))/length(ee)
      return(data.frame(c = c,
                        n = n,
                        p = p))
    }) %>%
      
      dplyr::bind_rows() %>%
      dplyr::mutate(tl = unique(mod$summary$TL)) %>%
      
      tidyr::gather("nutrient", "prop_lim", - tl)
    
    plot <- ggplot(lim) +
      geom_point(aes(x = tl, y = prop_lim, color = nutrient)) +
      geom_line(aes(x = tl, y = prop_lim, color = nutrient)) +
      labs(x = "Total length (cm)", y = "Proportion of iterations", color = "Limiting element") +
      theme_bw() +
      scale_color_fish_d(option = "Trimma_lantana", end = 0.8)+
      theme(legend.position = "bottom")
    print(plot)
    
    return(plot)
    
  }
  
  p1 <- plot_lim(mod_zs) + labs(title = "Z. scopas") +
    theme(legend.position = "none", plot.title = element_text(face = "italic")) +
    annotation_custom(g_zs, xmin=21, xmax=25, ymin=-Inf, ymax=Inf)
  p2 <- plot_lim(mod_bu) + labs(title = "B. undulatus") +
    theme(legend.position = "none", plot.title = element_text(face = "italic"))+
    annotation_custom(g_bu, xmin=21, xmax=25, ymin=-Inf, ymax=Inf)
  p3 <- plot_lim(mod_em) + labs(title = "E. merra")+
    theme(legend.position = "none", plot.title = element_text(face = "italic"))+
    annotation_custom(g_em, xmin=21, xmax=25, ymin=-Inf, ymax=Inf)
  g_legend <- function(plot){ 
    tmp <- ggplot_gtable(ggplot_build(plot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    return(legend)} 
  
  legend <- g_legend(plot_lim(mod_bu)) 
  
  plot <- 
    cowplot::plot_grid(p1,p2,p3, legend, ncol = 1, rel_heights = c(5,5,5,1))
  plot
  ggsave(out, width = 8, height = 10)
}

#### Figure 3 ####

##### new plot full model vs c model
makeplot3 <- function(out){
  
  ##data
  param <- read.csv("data/param_modelsp.csv")
  param_bu <- as.list(param[param$species=="Balistapus undulatus",1:37])
  param_zs <- as.list(param[param$species=="Zebrasoma scopas",1:37])
  param_em <- as.list(param[param$species=="Epinephelus merra",1:37])
  
  ## fish 
  img_bu <- readPNG("figures/Bundulatus_grey.png")
  g_bu <- rasterGrob(img_bu, interpolate=TRUE)
  img_em <- readPNG("figures/Emerra_grey.png")
  g_em <- rasterGrob(img_em, interpolate=TRUE)
  img_zs <- readPNG("figures/Zscopas_grey.png")
  g_zs <- rasterGrob(img_zs, interpolate=TRUE)
  
  
  mod_zs <- fishflux::cnp_model_mcmc(TL = seq(2, 16, by = 2), param = param_zs, iter = 6000)
  
  res <- 
    lapply(1:8, function(x){
      ee <- rstan::extract(mod_zs$stanfit[[x]]) %>% as.data.frame() %>% 
        mutate(Fp_2 = (0.7 * Sc * Dp/Dc) - Gp,
               Fn_2 = (0.7 * Sc * Dn/Dc) - Gn) %>% mutate(lt = round(lt)) %>%
        group_by(lt) %>%
        summarise(lt = unique(lt),
                  Fp_m = median(Fp), Fp_q1 = quantile(Fp, 0.025), Fp_q3 = quantile(Fp, 0.975), 
                  Fp_q51 = quantile(Fp, 0.25), Fp_q53 = quantile(Fp, 0.75),
                  Fp_q81 = quantile(Fp, 0.1), Fp_q83 = quantile(Fp, 0.9),
                  Fn_m = median(Fn), Fn_q1 = quantile(Fn, 0.025), Fn_q3 = quantile(Fn, 0.975),
                  Fn_q51 = quantile(Fn, 0.25), Fn_q53 = quantile(Fn, 0.75),
                  Fn_q81 = quantile(Fn, 0.1), Fn_q83 = quantile(Fn, 0.9),
                  Ic_m = median(Ic), Ic_q1 = quantile(Ic, 0.025), Ic_q3 = quantile(Ic, 0.975),
                  Ic_q51 = quantile(Ic, 0.25), Ic_q53 = quantile(Ic, 0.75),
                  Ic_q81 = quantile(Ic, 0.1), Ic_q83 = quantile(Ic, 0.9),
                  Sc_m = median(Sc), Sc_q1 = quantile(Sc, 0.025), Sc_q3 = quantile(Sc, 0.975),
                  Sc_q51 = quantile(Sc, 0.25), Sc_q53 = quantile(Sc, 0.75),
                  Sc_q81 = quantile(Sc, 0.1), Sc_q83 = quantile(Sc, 0.9),
                  Fp_2_m = median(Fp_2), Fp_2_q1 = quantile(Fp_2, 0.025), Fp_2_q3 = quantile(Fp_2, 0.975),
                  Fp_2_q51 = quantile(Fp_2, 0.25), Fp_2_q53 = quantile(Fp_2, 0.75),
                  Fp_2_q81 = quantile(Fp_2, 0.1), Fp_2_q83 = quantile(Fp_2, 0.9),
                  Fn_2_m = median(Fn_2), Fn_2_q1 = quantile(Fn_2, 0.025), Fn_2_q3 = quantile(Fn_2, 0.975),
                  Fn_2_q51 = quantile(Fn_2, 0.25), Fn_2_q53 = quantile(Fn_2, 0.75),
                  Fn_2_q81 = quantile(Fn_2, 0.1), Fn_2_q83 = quantile(Fn_2, 0.9)
        )
      return(ee)
    }) %>% plyr::ldply()
  
  
  col <- fishualize::fish(option = "Trimma_lantana", n = 2, end = 0.7)
  
  zs_p1 <- 
    ggplot(res) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Ic_q1, yend = Ic_q3), size = 1, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Ic_q81, yend = Ic_q83), size = 2, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Ic_q51, yend = Ic_q53), size = 3, color = col[1]) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Sc_q1, yend = Sc_q3), size = 1, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Sc_q81, yend = Sc_q83), size = 2, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Sc_q51, yend = Sc_q53), size = 3, color = col[2]) +
    geom_line(aes(x = lt, y = Ic_m), color = col[1], linetype = 3, size = 1) +
    geom_point(aes(x = lt + 0.2 , y = Ic_m), shape = "-", size = 6) +
    geom_line(aes(x = lt, y = Sc_m), color = col[2], linetype = 3, size = 1) +
    geom_point(aes(x = lt - 0.2 , y = Sc_m), shape = "-", size = 6) +
    labs(x = "Total length (cm)", y = "Ingestion rate (g C/day)") +
    annotation_custom(g_zs, xmin = 3 , xmax = 8, ymax = max(res$Ic_q3), ymin =  max(res$Ic_q3)/2) +
    theme_bw()
  
  zs_p3 <-
    ggplot(res) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fp_q1, yend = Fp_q3), size = 1, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fp_q81, yend = Fp_q83), size = 2, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fp_q51, yend = Fp_q53), size = 3, color = col[1]) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fp_2_q1, yend = Fp_2_q3), size = 1, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fp_2_q81, yend = Fp_2_q83), size = 2, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fp_2_q51, yend = Fp_2_q53), size = 3, color = col[2]) +
    geom_line(aes(x = lt, y = Fp_m), color = col[1], linetype = 3, size = 1) +
    geom_point(aes(x = lt + 0.2,  y = Fp_m), shape = "-", size = 6) +
    geom_line(aes(x = lt,  y = Fp_2_m), color = col[2], linetype = 3, size = 1) +
    geom_point(aes(x = lt - 0.2,  y = Fp_2_m), shape = "-", size = 6) +
    labs(x = "Total length (cm)", y = "Excretion rate (g P/day)") +
    annotation_custom(g_zs, xmin = 3 , xmax = 8, ymax = max(res$Fp_q3), ymin =  max(res$Fp_q3)/2) +
    theme_bw()
  
  zs_p2 <-
    ggplot(res) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fn_q1, yend = Fn_q3), size = 1, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fn_q81, yend = Fn_q83), size = 2, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fn_q51, yend = Fn_q53), size = 3, color = col[1]) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fn_2_q1, yend = Fn_2_q3), size = 1, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fn_2_q81, yend = Fn_2_q83), size = 2, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fn_2_q51, yend = Fn_2_q53), size = 3, color = col[2]) +
    geom_line(aes(x = lt, y = Fn_m), color = col[1], linetype = 3, size = 1) +
    geom_point(aes(x = lt + 0.2, y = Fn_m), shape = "-", size = 6) +
    geom_line(aes(x = lt, y = Fn_2_m), color = col[2], linetype = 3, size = 1) +
    geom_point(aes(x = lt - 0.2, y = Fn_2_m), shape = "-", size = 6) +
    labs(x = "Total length (cm)", y = "Excretion rate (g N/day)") +
    annotation_custom(g_zs, xmin = 3 , xmax = 8, ymax = max(res$Fn_q3), ymin =  max(res$Fn_q3)/2) +
    theme_bw()
  
  
  mod_bu <- fishflux::cnp_model_mcmc(TL = seq(2, 16, by = 2), param = param_bu, iter = 6000)
  
  res <- 
    lapply(1:8, function(x){
      ee <- rstan::extract(mod_bu$stanfit[[x]]) %>% as.data.frame() %>% 
        mutate(Fp_2 = (0.7 * Sc * Dp/Dc) - Gp,
               Fn_2 = (0.7 * Sc * Dn/Dc) - Gn) %>% mutate(lt = round(lt)) %>%
        group_by(lt) %>%
        summarise(lt = unique(lt),
                  Fp_m = median(Fp), Fp_q1 = quantile(Fp, 0.025), Fp_q3 = quantile(Fp, 0.975), 
                  Fp_q51 = quantile(Fp, 0.25), Fp_q53 = quantile(Fp, 0.75),
                  Fp_q81 = quantile(Fp, 0.1), Fp_q83 = quantile(Fp, 0.9),
                  Fn_m = median(Fn), Fn_q1 = quantile(Fn, 0.025), Fn_q3 = quantile(Fn, 0.975),
                  Fn_q51 = quantile(Fn, 0.25), Fn_q53 = quantile(Fn, 0.75),
                  Fn_q81 = quantile(Fn, 0.1), Fn_q83 = quantile(Fn, 0.9),
                  Ic_m = median(Ic), Ic_q1 = quantile(Ic, 0.025), Ic_q3 = quantile(Ic, 0.975),
                  Ic_q51 = quantile(Ic, 0.25), Ic_q53 = quantile(Ic, 0.75),
                  Ic_q81 = quantile(Ic, 0.1), Ic_q83 = quantile(Ic, 0.9),
                  Sc_m = median(Sc), Sc_q1 = quantile(Sc, 0.025), Sc_q3 = quantile(Sc, 0.975),
                  Sc_q51 = quantile(Sc, 0.25), Sc_q53 = quantile(Sc, 0.75),
                  Sc_q81 = quantile(Sc, 0.1), Sc_q83 = quantile(Sc, 0.9),
                  Fp_2_m = median(Fp_2), Fp_2_q1 = quantile(Fp_2, 0.025), Fp_2_q3 = quantile(Fp_2, 0.975),
                  Fp_2_q51 = quantile(Fp_2, 0.25), Fp_2_q53 = quantile(Fp_2, 0.75),
                  Fp_2_q81 = quantile(Fp_2, 0.1), Fp_2_q83 = quantile(Fp_2, 0.9),
                  Fn_2_m = median(Fn_2), Fn_2_q1 = quantile(Fn_2, 0.025), Fn_2_q3 = quantile(Fn_2, 0.975),
                  Fn_2_q51 = quantile(Fn_2, 0.25), Fn_2_q53 = quantile(Fn_2, 0.75),
                  Fn_2_q81 = quantile(Fn_2, 0.1), Fn_2_q83 = quantile(Fn_2, 0.9)
        )
      return(ee)
    }) %>% plyr::ldply()
  
  
  bu_p1 <- 
    ggplot(res) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(aes(x = lt, y = Ic_m), color = col[1], linetype = 3, size = 1) +
    geom_line(aes(x = lt, y = Sc_m), color = col[2], linetype = 3, size = 1) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Ic_q1, yend = Ic_q3), size = 1, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Ic_q81, yend = Ic_q83), size = 2, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Ic_q51, yend = Ic_q53), size = 3, color = col[1]) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Sc_q1, yend = Sc_q3), size = 1, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Sc_q81, yend = Sc_q83), size = 2, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Sc_q51, yend = Sc_q53), size = 3, color = col[2]) +
    geom_point(aes(x = lt + 0.2, y = Ic_m), shape = "-", size = 6) +
    geom_point(aes(x = lt - 0.2, y = Sc_m), shape = "-", size = 6) +
    labs(x = "Total length (cm)", y = "Ingestion rate (g C/day)") +
    annotation_custom(g_bu, xmin = 3 , xmax = 8, ymax = max(res$Ic_q3), ymin =  max(res$Ic_q3)/2) +
    theme_bw()
  
  bu_p3 <-
    ggplot(res) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(aes(x = lt, y = Fp_m), color = col[1], linetype = 3, size = 1) +
    geom_line(aes(x = lt, y = Fp_2_m), color = col[2], linetype = 3, size = 1) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fp_q1, yend = Fp_q3), size = 1, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fp_q81, yend = Fp_q83), size = 2, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fp_q51, yend = Fp_q53), size = 3, color = col[1]) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fp_2_q1, yend = Fp_2_q3), size = 1, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fp_2_q81, yend = Fp_2_q83), size = 2, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fp_2_q51, yend = Fp_2_q53), size = 3, color = col[2]) +
    geom_point(aes(x = lt + 0.2, y = Fp_m), shape = "-", size = 6) +
    geom_point(aes(x = lt - 0.2, y = Fp_2_m), shape = "-", size = 6) +
    labs(x = "Total length (cm)", y = "Excretion rate (g P/day)") +
    annotation_custom(g_bu, xmin = 3 , xmax = 8, ymax = max(res$Fp_q3), ymin =  max(res$Fp_q3)/2) +
    theme_bw()
  
  bu_p2 <-
    ggplot(res) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(aes(x = lt, y = Fn_m), color = col[1], linetype = 3, size = 1) +
    geom_line(aes(x = lt, y = Fn_2_m), color = col[2], linetype = 3, size = 1) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fn_q1, yend = Fn_q3), size = 1, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fn_q81, yend = Fn_q83), size = 2, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fn_q51, yend = Fn_q53), size = 3, color = col[1]) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fn_2_q1, yend = Fn_2_q3), size = 1, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fn_2_q81, yend = Fn_2_q83), size = 2, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fn_2_q51, yend = Fn_2_q53), size = 3, color = col[2]) +
    geom_point(aes(x = lt + 0.2, y = Fn_m), shape = "-", size = 6) +
    geom_point(aes(x = lt - 0.2,  y = Fn_2_m), shape = "-", size = 6) +
    labs(x = "Total length (cm)", y = "Excretion rate (g N/day)") +
    annotation_custom(g_bu, xmin = 3 , xmax = 8, ymax = max(res$Fn_q3), ymin =  max(res$Fn_q3)/2) +
    theme_bw()
  
  mod_em <- fishflux::cnp_model_mcmc(TL = seq(2, 16, by = 2), param = param_em, iter = 6000)
  
  res <- 
    lapply(1:8, function(x){
      ee <- rstan::extract(mod_em$stanfit[[x]]) %>% as.data.frame() %>% 
        mutate(Fp_2 = (0.7 * Sc * Dp/Dc) - Gp,
               Fn_2 = (0.7 * Sc * Dn/Dc) - Gn) %>% mutate(lt = round(lt)) %>%
        group_by(lt) %>%
        summarise(lt = unique(lt),
                  Fp_m = median(Fp), Fp_q1 = quantile(Fp, 0.025), Fp_q3 = quantile(Fp, 0.975), 
                  Fp_q51 = quantile(Fp, 0.25), Fp_q53 = quantile(Fp, 0.75),
                  Fp_q81 = quantile(Fp, 0.1), Fp_q83 = quantile(Fp, 0.9),
                  Fn_m = median(Fn), Fn_q1 = quantile(Fn, 0.025), Fn_q3 = quantile(Fn, 0.975),
                  Fn_q51 = quantile(Fn, 0.25), Fn_q53 = quantile(Fn, 0.75),
                  Fn_q81 = quantile(Fn, 0.1), Fn_q83 = quantile(Fn, 0.9),
                  Ic_m = median(Ic), Ic_q1 = quantile(Ic, 0.025), Ic_q3 = quantile(Ic, 0.975),
                  Ic_q51 = quantile(Ic, 0.25), Ic_q53 = quantile(Ic, 0.75),
                  Ic_q81 = quantile(Ic, 0.1), Ic_q83 = quantile(Ic, 0.9),
                  Sc_m = median(Sc), Sc_q1 = quantile(Sc, 0.025), Sc_q3 = quantile(Sc, 0.975),
                  Sc_q51 = quantile(Sc, 0.25), Sc_q53 = quantile(Sc, 0.75),
                  Sc_q81 = quantile(Sc, 0.1), Sc_q83 = quantile(Sc, 0.9),
                  Fp_2_m = median(Fp_2), Fp_2_q1 = quantile(Fp_2, 0.025), Fp_2_q3 = quantile(Fp_2, 0.975),
                  Fp_2_q51 = quantile(Fp_2, 0.25), Fp_2_q53 = quantile(Fp_2, 0.75),
                  Fp_2_q81 = quantile(Fp_2, 0.1), Fp_2_q83 = quantile(Fp_2, 0.9),
                  Fn_2_m = median(Fn_2), Fn_2_q1 = quantile(Fn_2, 0.025), Fn_2_q3 = quantile(Fn_2, 0.975),
                  Fn_2_q51 = quantile(Fn_2, 0.25), Fn_2_q53 = quantile(Fn_2, 0.75),
                  Fn_2_q81 = quantile(Fn_2, 0.1), Fn_2_q83 = quantile(Fn_2, 0.9)
        )
      return(ee)
    }) %>% plyr::ldply()
  
  
  em_p1 <- 
    ggplot(res) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(aes(x = lt, y = Ic_m), color = col[1], linetype = 3, size = 1) +
    geom_line(aes(x = lt, y = Sc_m), color = col[2], linetype = 3, size = 1) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Ic_q1, yend = Ic_q3), size = 1, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Ic_q81, yend = Ic_q83), size = 2, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Ic_q51, yend = Ic_q53), size = 3, color = col[1]) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Sc_q1, yend = Sc_q3), size = 1, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Sc_q81, yend = Sc_q83), size = 2, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Sc_q51, yend = Sc_q53), size = 3, color = col[2]) +
    geom_point(aes(x = lt + 0.2, y = Ic_m), shape = "-", size = 6) +
    geom_point(aes(x = lt - 0.2, y = Sc_m), shape = "-", size = 6) +
    labs(x = "Total length (cm)", y = "Ingestion rate (g C/day)") +
    annotation_custom(g_em, xmin = 3 , xmax = 8, ymax = max(res$Ic_q3), ymin =  max(res$Ic_q3)/2) +
    theme_bw()
  
  em_p3 <-
    ggplot(res) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(aes(x = lt, y = Fp_m), color = col[1], linetype = 3, size = 1) +
    geom_line(aes(x = lt, y = Fp_2_m), color = col[2], linetype = 3, size = 1) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fp_q1, yend = Fp_q3), size = 1, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fp_q81, yend = Fp_q83), size = 2, color = col[1], alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fp_q51, yend = Fp_q53), size = 3, color = col[1]) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fp_2_q1, yend = Fp_2_q3), size = 1, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fp_2_q81, yend = Fp_2_q83), size = 2, color = col[2], alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fp_2_q51, yend = Fp_2_q53), size = 3, color = col[2]) +
    geom_point(aes(x = lt + 0.2, y = Fp_m), shape = "-", size = 6) +
    geom_point(aes(x = lt - 0.2, y = Fp_2_m), shape = "-", size = 6) +
    labs(x = "Total length (cm)", y = "Excretion rate (g P/day)") +
    annotation_custom(g_em, xmin = 3 , xmax = 8, ymax = max(res$Fp_q3), ymin =  max(res$Fp_q3)/2) +
    theme_bw()
  
  em_p2 <-
    ggplot(res) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(aes(x = lt, y = Fn_m), color = col[1], linetype = 3, size = 1) +
    geom_line(aes(x = lt, y = Fn_2_m), color = col[2], linetype = 3, size = 1) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fn_q1, yend = Fn_q3, color = "col1"), size = 1, alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fn_q81, yend = Fn_q83, color = "col1"), size = 2, alpha = 0.5) +
    geom_segment(aes(x = lt + 0.2, xend = lt + 0.2, y = Fn_q51, yend = Fn_q53, color = "col1"), size = 3) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fn_2_q1, yend = Fn_2_q3, color = "col2"), size = 1, alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fn_2_q81, yend = Fn_2_q83, color = "col2"), size = 2, alpha = 0.5) +
    geom_segment(aes(x = lt - 0.2, xend = lt - 0.2, y = Fn_2_q51, yend = Fn_2_q53, color = "col2"), size = 3) +
    geom_point(aes(x = lt + 0.2, y = Fn_m), shape = "-", size = 6) +
    geom_point(aes(x = lt - 0.2, y = Fn_2_m), shape = "-", size = 6) +
    labs(x = "Total length (cm)", y = "Excretion rate (g N/day)") +
    scale_color_manual(name = "", values = c(col1 = col[1], col2 = col[2]), labels = c("Full model", "C-only model")) +
    annotation_custom(g_em, xmin = 3 , xmax = 8, ymax = max(res$Fn_q3), ymin =  max(res$Fn_q3)/2) +
    theme_bw() +
    theme(legend.position = "bottom", legend.text = element_text(size = 12))
  
  
  p_all <- egg::ggarrange(zs_p1, zs_p2, zs_p3, bu_p1, bu_p2, bu_p3, em_p1, em_p2, em_p3, nrow = 3, 
                          labels = c("A", "B", "C", "D", "E","F", "G", "H", "I"),  
                          label.args = list(gp = grid::gpar(font = 2, cex = 1.2)))
  
  ggsave(out, p_all, width = 14, height = 12 )
  
}

#### Figure 4 ####
makeplot4 <- function(out){

  ## fish 
  img_bu <- readPNG("figures/Bundulatus_grey.png")
  g_bu <- rasterGrob(img_bu, interpolate=TRUE)
  img_em <- readPNG("figures/Emerra_grey.png")
  g_em <- rasterGrob(img_em, interpolate=TRUE)
  img_zs <- readPNG("figures/Zscopas_grey.png")
  g_zs <- rasterGrob(img_zs, interpolate=TRUE)
  
  
  ##data
  param <- read.csv("data/param_modelsp.csv")
  param_bu <- as.list(param[param$species=="Balistapus undulatus",1:37])
  param_zs <- as.list(param[param$species=="Zebrasoma scopas",1:37])
  param_em <- as.list(param[param$species=="Epinephelus merra",1:37])
  
  exp <- read.csv("data/excretion_data.csv")
  exp_bu <- exp[exp$sp == "Balistapus undulatus",]
  exp_zs <- exp[exp$sp == "Zebrasoma scopas",]
  exp_em <- exp[exp$sp == "Epinephelus merra",]
  
  ## tl
  tl_bu <- exp((log(exp_bu$biomass/param_bu$lwa_m)/param_bu$lwb_m))
  tl_zs <- unique(exp((log(exp_zs$biomass/param_zs$lwa_m)/param_zs$lwb_m)))
  tl_em <- exp((log(exp_em$biomass/param_em$lwa_m)/param_em$lwb_m))
  
  ## run model
  mod_bu <- fishflux::cnp_model_mcmc(tl_bu, param = param_bu, iter = 5000,  cor = list(ro_Qc_Qn = 0.5, ro_Qc_Qp = -0.3, ro_Qn_Qp = -0.2,
                                                                                       ro_Dc_Dn = 0.8, ro_Dc_Dp = -0.3, ro_Dn_Dp = -0.3,
                                                                                       ro_lwa_lwb = 0.9, ro_alpha_f0 = 0.9))
  mod_zs <- fishflux::cnp_model_mcmc(tl_zs, param = param_zs, iter = 5000)
  mod_em <- fishflux::cnp_model_mcmc(tl_em, param = param_em, iter = 5000,  cor = list(ro_Qc_Qn = 0.5, ro_Qc_Qp = -0.3, ro_Qn_Qp = -0.2,
                                                                                       ro_Dc_Dn = 0.8, ro_Dc_Dp = -0.3, ro_Dn_Dp = -0.3,
                                                                                       ro_lwa_lwb = 0.9, ro_alpha_f0 = 0.9))
  
  
  get_iter <- function(x){
    get <- t(plyr::ldply(x))
    colnames(get) <- get[1,]
    get <- data.frame(apply(get[-1,],2,as.numeric))
    get$iter <- 1:nrow(get)
    return(get)
  }
  
  iter_bu <- (lapply(mod_bu$stanfit, FUN = function(x){
    rstan::extract(x, c("Fn", "Fp", "w1", "Ic", "IN"))})) %>%
    lapply( FUN = get_iter) %>%
    dplyr::bind_rows()
  iter_bu$tl <- rep(tl_bu, each = 2500)
  iter_bu$biomass <- rep(exp_bu$biomass, each = 2500)
  
  
  plot_bu_N <- 
    ggplot(group_by(iter_bu, iter), aes(x = biomass, y = Fn))+
    stat_lineribbon(alpha = 0.8, show.legend = FALSE) +
    scale_fill_brewer()+
    theme_bw()+
    geom_point(aes(x = biomass, y = ex_N), data = exp_bu, size = 1)+
    labs(x = "biomass (g)", y = "N excretion (g/day)") +
    annotation_custom(g_bu, xmin = 0, xmax = 150, 
                      ymin = quantile(filter(iter_bu, tl == max(iter_bu$tl))$Fn, 0.975)/2, 
                      ymax = quantile(filter(iter_bu, tl == max(iter_bu$tl))$Fn, 0.975))
  
  plot_bu_P <- 
    ggplot(group_by(iter_bu, iter), aes(x = biomass, y = Fp)) +
    stat_lineribbon(alpha = 0.8, show.legend = FALSE) +
    scale_fill_brewer() +
    theme_bw() +
    geom_point(aes(x=biomass,y=ex_P),data=exp_bu, size = 1) +
    labs(x = "biomass (g)", y = "P excretion (g/day)") +
    annotation_custom(g_bu, xmin = 0, xmax = 150, 
                      ymin = quantile(filter(iter_bu, tl == max(iter_bu$tl))$Fp, 0.975)/2, 
                      ymax = quantile(filter(iter_bu, tl == max(iter_bu$tl))$Fp, 0.975))
  
  
  iter_zs <- (lapply(mod_zs$stanfit, FUN = function(x){
    rstan::extract(x, c("Fn", "Fp", "w1", "Ic", "IN"))})) %>%
    lapply( FUN = get_iter) %>%
    dplyr::bind_rows()
  iter_zs$tl <- rep(tl_zs, each = 2500)
  iter_zs$biomass <- rep(unique(exp_zs$biomass), each = 2500)
  
  plot_zs_N <- 
    ggplot(group_by(iter_zs, iter), aes(x = biomass, y = Fn)) +
    stat_lineribbon(alpha = 0.8, show.legend = FALSE) +
    scale_fill_brewer() +
    theme_bw() +
    geom_point(aes(x=biomass,y=ex_N),data=exp_zs, size = 1) +
    labs(x = "biomass (g)", y = "N excretion (g/day)") +
    annotation_custom(g_zs, xmin = 0, xmax = 50, 
                      ymin = quantile(filter(iter_zs, tl == max(iter_zs$tl))$Fn, 0.975)/2, 
                      ymax = quantile(filter(iter_zs, tl == max(iter_zs$tl))$Fn, 0.975))
  
  
  plot_zs_P <- 
    ggplot(group_by(iter_zs, iter), aes(x = biomass, y = Fp)) +
    stat_lineribbon(alpha = 0.8, show.legend = FALSE) +
    scale_fill_brewer() +
    theme_bw() +
    geom_point(aes(x=biomass,y=ex_P),data=exp_zs, size = 1) +
    labs(x = "biomass (g)", y = "P excretion (g/day)") +
    annotation_custom(g_zs, xmin = 0, xmax = 50, 
                      ymin = max(exp_zs$ex_P, na.rm = TRUE)/2, 
                      ymax = max(exp_zs$ex_P, na.rm = TRUE))
  
  
  iter_em <- (lapply(mod_em$stanfit, FUN = function(x){
    rstan::extract(x, c("Fn", "Fp", "w1", "Ic", "IN"))})) %>%
    lapply( FUN = get_iter) %>%
    dplyr::bind_rows()
  iter_em$tl <- rep(tl_em, each = 2500)
  iter_em$biomass <- rep(exp_em$biomass, each = 2500)
  
  plot_em_N <- 
    ggplot(group_by(iter_em, iter), aes(x = biomass, y = Fn))+
    stat_lineribbon(alpha = 0.8, show.legend = FALSE) +
    scale_fill_brewer()+
    theme_bw()+
    geom_point(aes(x=biomass,y=ex_N),data=exp_em, size = 1)+
    labs(x = "biomass (g)", y = "N excretion (g/day)") +
    annotation_custom(g_em, xmin = 0, xmax = 90, 
                      ymin = quantile(filter(iter_em, tl == max(iter_em$tl))$Fn, 0.975)/2, 
                      ymax = quantile(filter(iter_em, tl == max(iter_em$tl))$Fn, 0.975))
  
  
  plot_em_P <- 
    ggplot(group_by(iter_em, iter), aes(x = biomass, y = Fp))+
    stat_lineribbon(alpha = 0.8, show.legend = FALSE) +
    scale_fill_brewer()+
    theme_bw()+
    geom_point(aes(x=biomass,y=ex_P),data=exp_em, size = 1)+
    labs(x = "biomass (g)", y = "P excretion (g/day)")+
    annotation_custom(g_em, xmin = 0, xmax = 90, 
                      ymin = quantile(filter(iter_em, tl == max(iter_em$tl))$Fp, 0.975)/2, 
                      ymax = quantile(filter(iter_em, tl == max(iter_em$tl))$Fp, 0.975))
  
  
  plot <- cowplot::plot_grid(plot_zs_N,plot_zs_P,plot_bu_N,plot_bu_P,plot_em_N,plot_em_P, nrow = 3, 
                             labels = c("A", "B", "C", "D", "E","F"))
  plot
  ggsave(out, width = 8, height = 10)
}

##### Figure 5 #####
makeplot5 <- function(out){
  
  param <- read.csv("data/param_modelsp.csv")
  params <- as.list(param[param$species=="Zebrasoma scopas",1:37])
  #only averages for this simulation
  params <- params[grep("_m" ,names(params))]  
  
  mod <- cnp_model_mcmc(10:11, params)
  
  ex <- fishflux::extract(mod, c("Ic", "In", "Ip", "Fn", "Fp", 
                                 "Sc", "Sn", "Sp",
                                 "st_cn", "st_np")) %>% dplyr::filter(TL < 10.5)
  
  # TERs
  Scn <- ex$st_cn_median
  Snp <- ex$st_np_median
  
  # fix Dc at 20
  params$Dc_m <- 20
  Dcm <- params$Dc_m
  Dnm <- Dcm/Scn
  Dpm <- Dnm / Snp
 
  
  DN <- seq(from = 1.2, to = 3.2, by = 0.1)
  DP <- DN/Snp
  DNP <- expand.grid(DN, DP)
  colnames(DNP) <- c("Dn", "Dp")
  
  # Run model for multiple Dn and Dp values
  sim <- 
    lapply(1:nrow(DNP), function(i){
      params$Dn_m <- DNP$Dn[i]
      params$Dp_m <- DNP$Dp[i]
      mod <- fishflux::cnp_model_mcmc(TL = 10:11, param = params)
      ex <- fishflux::extract(mod, c("Ic", "In", "Ip", "Fn", "Fp", "Sc", "Sn", "Sp", 
                                     "lim", "st_cn", "st_np")) %>% 
        dplyr::filter(TL < 10.5)
      result <- cbind(DNP[i,], ex)
      print(i/nrow(DNP))
      return(result)
    })
  
  simd <- plyr::ldply(sim)
  
  p1 <- 
    ggplot()+
    geom_segment(aes(x = 1.2, y = 1.2/Snp, xend = Dnm, yend = Dpm), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 3.2/Snp), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 3.2, yend = Dpm), size = 2) +
    geom_text(aes(x = 1.7, y = 1, label = "N"), size = 15) +
    geom_text(aes(x = 2.7, y = 1.1, label = "C"), size = 15) +
    geom_text(aes(x = 2.4, y = 0.7, label = "P"), size = 15) +
    labs( x = "Dn (%)", y = "Dp (%)") +
    theme_bw() +
    theme(aspect.ratio=1, title = element_text(size = 16), 
          text = element_text(size = 14)) 
  
  p2 <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = Ic_median)) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    geom_segment(aes(x = 1.2, y = 1.2/Snp, xend = Dnm, yend = Dpm), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 3.2/Snp), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 3.2, yend = Dpm), size = 2) +
    theme_bw() +
    theme(aspect.ratio=1, title = element_text(size = 16), 
          text = element_text(size = 14)) +
    labs(x = "Dn (%)", y = "Dp (%)", fill = "Ic (g/day)") 
  
  p3 <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = Fn_median)) +
    geom_segment(aes(x = 1.2, y = 1.2/Snp, xend = Dnm, yend = Dpm), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 3.2/Snp), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 3.2, yend = Dpm), size = 2) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    theme_bw() +
    theme(aspect.ratio=1, title = element_text(size = 16), 
          text = element_text(size = 14))  +
    labs(x = "Dn (%)", y = "Dp (%)", fill = "Fn (g/day)") 
  
  p4 <-
    ggplot(simd) +
    geom_raster(aes(x = Dn, y = Dp, fill = Fp_median)) +
    scale_fill_fish(option = "Trimma_lantana", trans = "sqrt") +
    geom_segment(aes(x = 1.2, y = 1.2/Snp, xend = Dnm, yend = Dpm), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = Dnm, yend = 3.2/Snp), size = 2) +
    geom_segment(aes(x = Dnm, y = Dpm, xend = 3.2, yend = Dpm), size = 2) +
    theme_bw() +
    theme(aspect.ratio=1, title = element_text(size = 16), 
          text = element_text(size = 14)) +
    labs(x = "Dn (%)", y = "Dp (%)", fill = "Fp (g/day)") 
  
  
  all <- 
    egg::ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, labels = c("A", "B", "C", "D") ,
                   label.args = list(gp = grid::gpar(font = 2, cex =1.2)))
  
  ggsave(out, all, width = 8, height = 6)
}

