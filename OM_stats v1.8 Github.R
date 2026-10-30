##################################################################
######## Script for data processing, statistical        ##########
######## analysis, and plot generation of organic       ##########
######## matter fractions for burned sites on Chile     ##########
######## Developed by: Luis Daniel Olivares Martínez    ##########
######## Dec 2024 - Oct 2026                            ##########
######## Universidad Miguel Hernandez, Elche, Spain     ##########
######## Contact email: lolivares@umh.es                ##########
##################################################################

## STEP 0: Loading of project packages & working directory - RUN ALL THE STEP ####
library(readxl)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # establish project's working directory

## STEP 1: Data loading from working directory - RUN ALL THE STEP ####
sit <- read_excel("harmmonized_DB_Chile.xlsx", sheet="standard_DB"); sit <- as.data.frame(sit)
sit$sev_dpt_0 <- factor(sit$sev_dpt, levels=c("high","medium","low"))
sit$sev_dpt <- dplyr::recode(sit$sev_dpt, high = "medium-high", medium = "medium-high")
sit$sev_dpt <- factor(sit$sev_dpt, levels=c("low","medium-high"))
sit$sev_veg <- factor(sit$sev_veg, levels=c("low","high","medium"))
sit$Zonas <- factor(sit$Zonas)
sit$burned <- factor(sit$burned)
sit$clim <- factor(sit$clim)
# sit$age_veg <- factor(sit$age_veg, levels = c('stand initiation','stem exclusion','steady state'))
sit$veg_1 <- factor(sit$veg_1)
sit$veg_2 <- factor(sit$veg_2)
sit$asl <- cut(sit$masl, 
                breaks = seq(0,300, by=100), 
                labels = c("0-100", "100-200", "200-300"),
                right = F)  # incluye el límite inferior

cob <- read_excel("BD versiones/Datos 2012-2019.xlsx", sheet="Coberturas 2012-2019") %>% 
  rename(pto = `Puntos finales`,
         soil = `suelo desnudo %`,
         slope = `pendiente %`,
         yr = Año) %>% 
  mutate(pto = gsub("[^0-9]", "",pto) %>% as.numeric())
cob <- sit %>% select(key_1:pto,Sev.Cobertura:sev_dpt) %>% 
  distinct(pto, .keep_all = TRUE) %>% 
  right_join(., cob, by = "pto") #, relationship = "one-to-many"

nmr <- read_excel("harmmonized_DB_Chile.xlsx", sheet="Feuil1"); nmr <- as.data.frame(nmr)
cn3f <- read_excel("harmmonized_DB_Chile.xlsx", sheet="Datos_fracc_Patagonia", skip = 1) # , range = cell_cols("A:AA")
nmr_fg <- read_excel("harmmonized_DB_Chile.xlsx", sheet="4regions (3)"); nmr_fg <- as.data.frame(nmr_fg)
ftir <- read_excel("harmmonized_DB_Chile.xlsx", sheet="ftir") %>% as.data.frame() %>% 
  mutate(across(Amide:Aromatic, ~ replace_na(., 0)))

cn <- cn3f %>% 
  filter(is.na(olyr)) %>%
  group_by(key_1) %>%
  summarise(across(where(is.numeric), ~na_if(mean(., na.rm = TRUE), NaN)),
            across(where(is.character), unique)) %>% 
  select(-olyr) %>% 
  mutate(numeric_key = as.numeric(sub("R.*", "", key_1))) %>%  # Extraer el componente numérico
  arrange(numeric_key) %>% # Ordenar por el componente numérico
  select(-numeric_key) %>% 
  rename(TN = N, T_CN = CN) %>%
  pivot_longer(
    cols = starts_with("N_"):last_col(),
    names_to = c(".value", "fraction"),
    names_pattern = "(.*)_(fPOM|oPOM|MAOM)$"#,values_drop_na = T
  ) %>% 
  mutate(across(where(is.numeric), ~ round(., 2) )) %>% 
  filter(!is.na(ww))
cn$fraction <- factor(cn$fraction, levels = c('MAOM', 'oPOM', 'fPOM'))

# Function to replace 0 with NA in numeric columns using dplyr
rep_zero_na <- function(df) {
  df %>%
    mutate(across(where(is.numeric) & !all_of("fire"), ~ na_if(., 0)))
}

# Apply the function to each data frame
cn <- rep_zero_na(cn)

ftir <- ftir %>% mutate(across(Amide:Aromatic, ~ifelse(.==0,1e-6,.)))

full_db <- list(cn %>% select(-ZONA,-fire), 
                nmr %>% select(-ZONA,-fire, -name), 
                nmr_fg %>% select(-num,-pto, -name), 
                ftir %>% select(-ZONA,-fire, -name)) %>% 
  reduce(full_join, by = c("key_1", "fraction")) %>% 
  right_join(sit %>% select(key_1, Zonas, pto, burned, asl, sev_veg:pH), .,
             by = "key_1") %>% 
  mutate(fraction = factor(fraction, levels = c('fPOM', 'oPOM', 'MAOM')))

cn <- merge(cn,sit, by = "key_1", all.x = T)
nmr <- merge(nmr,sit, by = "key_1", all.x = T)
nmr_fg <- merge(nmr_fg,sit, by = "key_1", all.x = T)
ftir <- merge(ftir,sit, by = "key_1", all.x = T)

### Generating Variable-Distribution Table for step 3 ####
dis <- data.frame(var_ori = c(names(cn)[c(4:6,8:13)], 
                              names(nmr)[c(6:22)], names(ftir)[c(6:9)], 
                              names(nmr_fg)[c(6:10)]),
                  var = NA, dist = NA, bd = c(rep('cn',9),
                                              rep('nmr',17), rep('ftir',4),
                                              rep('nmr_fg',5)), neg = NA)
dis$dist <- c('normal', 'log-normal', 'normal', 'log-normal', 'log-normal', 'log-normal', 'gamma', 'exponential', 'gamma', 'normal', 'normal', 'normal', 'normal', 'gamma', 'gamma', 'normal', 'normal', 'uniform', 'normal', 'gamma', 'log-normal', 'normal', 'normal', 'gamma', 'normal', 'gamma', 'gamma', 'normal', 'normal', 'normal', 'normal', 'log-normal', 'normal', 'normal', 'normal')
normi <- which(dis$dist != 'log-normal')
dis[normi,2] <- dis[normi,1]
normi <- which(dis$dist == 'log-normal')
dis[normi,2] <- paste0('log(',dis[normi,1],')')
dis$dec <- c(1,2,1,rep(1,4),2,1, #tabla CN
             rep(3,6),5,2,0,4,2,3,3,4,4,3,4,rep(2,4), #tabla NMR y FTIR
             rep(1,4),3) #tabla NMR_fg
dis$neg[which(dis$var_ori == 'NOSC')] <- T #all the parameters with negative values (for Gamma distributions they need to be positive)

lom <- 9+21+5 # aumento en 4 para nmr y cn bajo de 21 a 9

### Generating Variable-Distribution Table for step 4 ####
full_db <- full_db %>% mutate(tca = tca/100 %>% round(5),
                              sa = sa/100 %>% round(5),
                              ang = ifelse(ang > 50, ang, 0), # Transision zone to hydrophillic phase at contact angle of 50° according to Shirtcliffe, N. J., McHale, G., Newton, M. I., Pyatt, F. B., & Doerr, S. H. (2006). Critical conditions for the wetting of soils. Applied Physics Letters, 89(9). https://doi.org/10.1063/1.2339072
                              wdpt = ifelse(wdpt < 5, 0,wdpt) )# Hydrophillic

full_db <- full_db %>% mutate(
  wdpt_0 = ifelse(wdpt != 0, 1, 0),
  ang_0 = ifelse(ang != 0, 1, 0),
  wdpt = ifelse(wdpt != 0, wdpt, NA),
  ang = ifelse(ang != 0, ang, NA)
)

fire_db <- subset(full_db, burned == 'yes')

sd_fire_db <- fire_db %>% 
  mutate(across(where(is.numeric) & -matches("_0$"), ~ as.vector(scale(.x, scale = T, center = 0))) )

sd_full_db <- full_db %>% 
  mutate(across(where(is.numeric) & -matches("_0$"), ~ as.vector(scale(.x, scale = T, center = 0))) )

response <- data.frame(var = NA,
                       var_ori = c("wdpt","wdpt_0","ang","ang_0","sa","tca","qCO2","Cmic","BSR"),
                       dist = c("log-normal",'binomial',"normal",'binomial',rep("beta",2),rep("gamma",2),"normal"),
                       dec = c(0,2,1,2,3,3,4,0,3) #,
                       # zero_mod = c(T,T,rep(F,5)), # activar si se eliminan _0
                       # loop = c(1,3,5:9) # activar si se eliminan _0
) %>% mutate(var = ifelse(dist == 'log-normal', 
                          paste0('log(',var_ori,')'), var_ori))

### Generate functions for post-processing data ####
post_process <- function(mods) {
  mods <- mods %>%
    filter(!is.na(REML)) %>%
    mutate(across(1:2, ~ ifelse(var %in% dis$var[fi.neg], . -0.4, .)))
}

p_vals <- function(textos) {
  vectores <- sapply(textos, function(x) eval(parse(text = x)))
  
  resultados <- data.frame()
  
  for (vector in vectores) {
    if (length(vector) == 4) {
      resultados <- rbind(resultados, as.data.frame(t(vector[-1])))
    } else if (length(vector) == 3) {
      resultados <- rbind(resultados, as.data.frame(t(c(vector[-1],NA))))
    } else if (length(vector) == 2) {
      resultados <- rbind(resultados, as.data.frame(t(c(NA,NA,vector[-1]))))
    }
  }
  
  colnames(resultados) <- paste0('p_',c("fPOM", "oPOM", "MAOM"))
  # resultados <- resultados[, !colnames(resultados) %in% "Intercept"]
  
  return(resultados)
}

p_vals2 <- function(textos) {
  vectores <- sapply(textos, function(x) eval(parse(text = x)))
  
  resultados <- data.frame()
  
  for (vector in vectores) {
    if (length(vector) == 4) {
      resultados <- rbind(resultados, as.data.frame(t(rev(vector[-1]))))
    } else if (length(vector) == 3) {
      resultados <- rbind(resultados, as.data.frame(t(c(vector[-1],NA))))
    } else if (length(vector) == 2) {
      resultados <- rbind(resultados, as.data.frame(t(c(NA,NA,vector[-1]))))
    }
  }
  
  colnames(resultados) <- paste0('p_',c("fPOM", "oPOM", "MAOM"))
  # resultados <- resultados[, !colnames(resultados) %in% "Intercept"]
  
  return(resultados)
}

post_process2 <- function(mods) {
  mods <- mods %>%
    mutate(AIC = as.numeric(AIC),BIC = as.numeric(BIC),
           fe_val = as.numeric(fe_val),
           fe_sde = as.numeric(fe_sde),
           fe_z = as.numeric(fe_z),
           fe_p = as.numeric(fe_p)
    ) %>% 
    filter(!is.na(AIC) & 
             !grepl("^Error", mod) & 
             !grepl("NaN", fe_p_vals))
  
  mods <- cbind(mods, p_vals(mods$fe_p_vals) )
  
  return(mods)
}

post_process3 <- function(mods) {
  mods <- mods %>%
    mutate(AIC_REML = as.numeric(AIC_REML),
    ) %>% 
    filter(!is.na(AIC_REML) & 
             !grepl("^Error", mod) & 
             !grepl("NaN", fe_p_vals))
  
  mods <- cbind(mods, p_vals2(mods$fe_p_vals) )
  
  return(mods)
}

min_mod <- function(mods) {
  mods <- mods %>%
    group_by(var) %>%
    # filter(if (all(fe1_p > 0.05)) TRUE else fe1_p < 0.05) %>%
    slice(which.min(AIC)) %>%
    arrange(match(var,dis$var_ori))
}

min_mod_BIC <- function(mods) {
  mods <- mods %>%
    group_by(var) %>%
    slice(which.min(BIC)) %>%
    arrange(match(var,dis$var_ori))
}

# min_mod <- function(mods) {
#   mods <- mods %>%
#     group_by(var) %>%
#     filter(if (all(fe1_p >= 0.05)) TRUE else
#       if (any(fe1_p < 0.05 &
#               REML - min(REML) < 10 & REML - min(REML) > 0) ) {
#         fe1_p < 0.05 & REML - min(REML) < 10
#       }  else fe1_p < 0.05
#     ) %>%
#     slice(which.min(REML)) %>%
#     arrange(match(var, dis$var_ori))
# }

check_p_values <- function(fe_p_vals, condition) {
  sapply(fe_p_vals, function(x) {
    p_values <- eval(parse(text = x))
    condition(p_values[-1])  # Discard the first element of the c()
  })
}

min_mod2 <- function(mods) {
  mods <- mods %>%
    # mutate(var = recode(var, "log(wdpt)" = "ln_wdpt")) %>% 
    group_by(var, x_var) %>%
    slice(which.min(AIC)) %>%
    mutate(
      var = factor(var, levels = response$var),
      x_var = factor(x_var, levels = dis$var_ori)
    ) %>%
    arrange(var, x_var)
  
  return(mods)
}

# min_mod2 <- function(mods) {
#   mods <- mods %>%
#     # mutate(var = recode(var, "log(wdpt)" = "ln_wdpt")) %>% 
#     group_by(var, x_var) %>%
#     filter(
#       if (all(check_p_values(fe_p_vals, function(p) all(p >= 0.05)))) {
#         TRUE
#       } else if (any(check_p_values(fe_p_vals, function(p) any(p < 0.05)) &
#                      (BIC - min(BIC) < 10 & BIC - min(BIC) > 0))) {
#         check_p_values(fe_p_vals, function(p) any(p < 0.05)) & BIC - min(BIC) < 10
#       } else {
#         check_p_values(fe_p_vals, function(p) any(p < 0.05))
#       }
#     ) %>%
#     slice(which.min(BIC)) %>%
#     # filter( # uncomment if you want just the significant ones
#     #   check_p_values(fe_p_vals, function(p) any(p < 0.05) )
#     # ) %>%
#     mutate(
#       var = factor(var, levels = response$var),
#       x_var = factor(x_var, levels = dis$var_ori)
#     ) %>%
#     arrange(var, x_var)
#   
#   return(mods)
# }


check_ps <- function(fe_p_vals, condition) {
  sapply(fe_p_vals, function(x) {
    p_values <- eval(parse(text = x))
    condition(p_values)  # All elements of the c()
  })
}

min_mod3 <- function(mods) {
  mods <- mods %>%
    group_by(var) %>%
    filter(
      if (all(check_ps(fe_p_vals, function(p) any(p >= 0.05)))) {
        TRUE
      } else if (any(check_ps(fe_p_vals, function(p) any(p < 0.05)) & 
                     (REML - min(REML) < 10 & REML - min(REML) > 0))) {
        check_ps(fe_p_vals, function(p) any(p < 0.05)) & REML - min(REML) < 10
      } else {
        check_ps(fe_p_vals, function(p) any(p < 0.05))
      }
    ) %>%
    slice(which.min(REML)) %>%
    arrange(match(var, dis$var_ori))
  
  return(mods)
}

min_mod4 <- function(mods) {
  mods <- mods %>%
    group_by(var) %>%
    filter(
      if (all(check_ps(fe_p_vals, function(p) any(p >= 0.05)))) {
        TRUE
      } else if (any(check_ps(fe_p_vals, function(p) any(p < 0.05)) & 
                     (AIC_REML - min(AIC_REML) < 10 & AIC_REML - min(AIC_REML) > 0))) {
        check_ps(fe_p_vals, function(p) any(p < 0.05)) & AIC_REML - min(AIC_REML) < 10
      } else {
        check_ps(fe_p_vals, function(p) any(p < 0.05))
      }
    ) %>%
    slice(which.min(AIC_REML)) %>%
    arrange(match(var, dis$var_ori))
  
  return(mods)
}

# Definir una función para transformar valores negativos a positivos
# neg_to_pos <- function(df, cols) {
#   df <- df %>%
#     mutate(across(all_of(cols), ~ abs(.)))
#   return(df)
# }
neg_pls4 <- function(df, cols) {
  df <- df %>%
    mutate(across(all_of(cols), ~ . + 0.4))
  return(df)
}
# Aplicar la función a cada data frame con valores negativos
fi.neg <- which(dis$neg)
lapply(dis$bd[fi.neg], function(df_name) {
  df <- get(df_name)
  cols_to_transform <- dis$var_ori[fi.neg][dis$bd[fi.neg] == df_name]
  df <- neg_pls4(df, cols_to_transform)
  assign(df_name, df, envir = .GlobalEnv)
}) %>% invisible()

## STEP 2. Generating statistical models ####
### environment ####
library(lme4)
library(lmerTest)
library(car)
library(nlme)
library(MuMIn)
library(emmeans)
# library(geepack)
# library(ggeffects)

### data.frame models to evaluate effects for presence - absence of fire p1 ####
np1 <- 12 #12 essays
mods <- array(NA,c(lom*np1,12))
mods <- as.data.frame(mods)
names(mods) <- c('var','dis.y','REML','AIC','BIC', 'mod','fe_p_vals','fe_p_names','re_names','re_p_vals', 'R2m', 'R2c')

mods$dis.y <- rep(dis$dist, each=np1)
mods$var <- rep(dis$var_ori, each=np1)

#### Function generate presence - absence models with random effects ####
run_model <- function(i, n) {
  rowfi <- np1*i-(np1-n)
  
  # Select random effects based on n
  raneff <- switch(as.character(n),
                   "1" = "(1|Zonas)",
                   "2" = "(1|clim:Zonas)",
                   "3" = "(1|veg_1:Zonas)",
                   "4" = "(1|veg_2:Zonas)",
                   "5" = "(1|veg_1)",
                   "6" = "(1|veg_2)",
                   "7" = "(1|clim:asl)",
                   "8" = "(1|asl:clim:Zonas)",
                   "9" = "(1|veg_2:veg_1)",
                   "10" = "(1|asl:clim:veg_1:Zonas)",
                   "11" = "(1|asl:clim:veg_2:Zonas)")
  
  # Add fixed effects and set fraction when applicable
  if (dis$bd[i] == "ftir") {
    fixeff <- "burned"
  } else {
    fixeff <- "fraction * burned"
    # raneff <- gsub("1\\|","fraction|", raneff) # removed due to loss of stat. power
    }
  
  formula_str <- paste(dis$var[i], "~", fixeff)
  
  if (n == 12) {
  
  mod.0 <- gls(as.formula(formula_str), data = get(dis$bd[i]))
  mods[rowfi, "AIC"] <- AIC(mod.0) %>% round(1)
  mods[rowfi, "BIC"] <- BIC(mod.0) %>% round(1)
  mods[rowfi, "mod"] <- formula_str
  
  mods[rowfi, "fe_p_names"] <- mod.0 %>% anova() %>% 
    row.names() %>% 
    paste(collapse = ', ')
  mods[rowfi, "fe_p_vals"] <- mod.0 %>% anova() %>% 
    pull('p-value') %>% round(5) %>% 
    paste(collapse = ', ')

  return(mods[rowfi,])  }
  
  formula_str <- paste(formula_str, "+", raneff)
  mods[rowfi, "mod"] <- formula_str
  
  # Generate model... omit when defined the if cases for different distributions
  # mod.cn <- lmerTest::lmer(as.formula(formula_str), data = get(dis$bd[i]))
  # s.mod <- summary(mod.cn)
  
  # Generate model according to variable distribution in dis$dist[i]
  if (dis$dist[i] %in% c("normal","log-normal") ) {
    mod.cn <- lmer(as.formula(formula_str), data = get(dis$bd[i]))
    s.mod <- summary(mod.cn)
    mods[rowfi, "REML"] <- s.mod[["AICtab"]] %>% round(1)
    
    reff <- mod.cn %>% ranova()
    mods[rowfi, "re_names"] <- reff %>% 
      row.names() %>% 
      .[-1] %>% paste(collapse = ', ')
    mods[rowfi, "re_p_vals"] <- reff %>% 
      pull('Pr(>Chisq)') %>% 
      .[-1] %>% round(5) %>% paste(collapse = ', ')
  } else if (dis$dist[i] %in% c("gamma", "exponential") ) {
    mod.cn <- glmer(as.formula(formula_str), 
                    data = get(dis$bd[i]), family = Gamma(link='log'))
  } 
  
  mods[rowfi, "AIC"] <- AIC(mod.cn) %>% round(1)
  mods[rowfi, "BIC"] <- BIC(mod.cn) %>% round(1)
  
  a.mod <- mod.cn %>% Anova()
  mods[rowfi, "fe_p_names"] <- a.mod %>% row.names() %>% 
    paste(collapse = ', ')
  mods[rowfi, "fe_p_vals"] <- a.mod %>% pull('Pr(>Chisq)') %>% 
    round(5) %>% paste(collapse = ', ')
  
  mods[rowfi, c('R2m', 'R2c')] <- mod.cn %>% r.squaredGLMM() %>% round(4)
  
  return(mods[rowfi,])
}

#### Loop through presence - absence models ####
for (i in 1:lom) {
  for (n in 1:np1) {
    if (dis$dist[i] == "uniform") { # 
      next  # Omitir y saltar al siguiente valor del bucle
    }
    mods[np1*i-(np1-n),] <- try({run_model(i, n)}, silent = TRUE)  # Usar silent = TRUE para suprimir mensajes de error
  }
}

# Generate and export summary tables
mods <- subset(mods, !is.na(AIC))
min.mods <- min_mod(mods)

write.csv(mods, file = "mod_abs_p1_reff.csv", row.names = F)
write.csv(min.mods, file = "mod_abs_p1_bst.csv", row.names = F)

### fixed effects function to evaluate presence - absence of fire p2 ####
run_mod2 <- function(i) {
  ll <- dim(mods2)[2]
  fi <- which(dis$var_ori == min.mods$var[i])
  df <- setNames(data.frame(matrix(ncol = ll, nrow = 1)), names(mods2))
  df$var <- min.mods$var[i]
  df$dis.y <- min.mods$dis.y[i]
  
  # Extract the fixed and random effects
  formula_str <- min.mods$mod[i]
  
  # all_eff <- gsub(".*(~.*)\\s*(.*)", "\\1", formula_str)
  # fix_eff <- gsub(".*~\\s*(.*?)\\s*\\+\\s*\\(1\\|.*", "\\1", formula_str)
  
  ran_eff <- ifelse(grepl("\\(1\\|", formula_str), TRUE, FALSE)
  
  # Check for the random effect, if any doing gls(method = 'ML')
  if(!ran_eff){
    mod <- gls(as.formula(formula_str), data = get(dis$bd[fi]),
               method = 'ML')
    mod$call$data <- dis$bd[fi] %>% as.name()
    # mod$call$model <- formula_str %>% as.formula() %>% as.call()
    
    a_mod <- mod %>% dredge(options(na.action = 'na.fail'))
    b_mod <- a_mod %>% get.models(1)
    b_mod <- b_mod[[1]]
    f_mod <- b_mod %>% formula()
    
    mod_out <- gls(f_mod, data = get(dis$bd[fi]))
    fix_eff <- ifelse(grepl("~\\s1", f_mod %>% deparse()), FALSE, TRUE)
    # df$fe_p_names <- mod_out %>% anova() %>%
    #   row.names() %>% .[-1] %>%
    #   paste(collapse = ', ')
    # df$fe_p_vals <- mod_out %>% anova() %>%
    #   pull('p-value') %>% .[-1] %>% round(5) %>%
    #   paste(collapse = ', ')
    
  }else{
    ran_eff <- gsub(".*?(\\(1\\|.*)", "\\1", formula_str)
    
    if (dis$dist[fi] %in% c("normal","log-normal") ) {
      mod <- lmer(as.formula(formula_str), data = get(dis$bd[fi]),
                  REML = FALSE)
    } else if (dis$dist[fi] %in% c("gamma", "exponential") ) {
      mod <- glmer(as.formula(formula_str), 
                      data = get(dis$bd[fi]), family = Gamma(link='log'))
    } 
    mod@call$data <- dis$bd[fi] %>% as.name()
    
    a_mod <- mod %>% dredge(options(na.action = 'na.fail'))
    b_mod <- a_mod %>% get.models(1)
    b_mod <- b_mod[[1]]
    f_mod <- b_mod %>% formula()
    
    fix_eff <- ifelse(grepl("~\\s\\s*[^()]+\\s*\\(1", f_mod %>% deparse()), TRUE, FALSE)
    
    if (dis$dist[fi] %in% c("normal","log-normal") ){
      mod_out <- lmer(f_mod, data = get(dis$bd[fi]))
      s.mod <- summary(mod_out)
      df$REML <- s.mod[["AICtab"]] %>% round(1)
      reff <- mod_out %>% ranova()
      df$re_names <- reff %>% 
        row.names() %>% 
        .[-1] %>% paste(collapse = ', ')
      df$re_p_vals <- reff %>% 
        pull('Pr(>Chisq)') %>% 
        .[-1] %>% round(5) %>% paste(collapse = ', ')
    } else if (dis$dist[fi] %in% c("gamma", "exponential") ) {
      mod_out <- b_mod
    }
    df[1, c('R2m', 'R2c')] <- mod_out %>% r.squaredGLMM() %>% round(4)
  }
  if (fix_eff){
    anova_mod <- mod_out %>% Anova()
    df$fe_p_names <- anova_mod %>% row.names() %>%
      paste(collapse = ', ')
    df$fe_p_vals <- anova_mod %>% pull('Pr(>Chisq)') %>%
      round(5) %>% paste(collapse = ', ')
  }

  df$AIC <- AIC(mod_out) %>% round(1)
  df$BIC <- BIC(mod_out) %>% round(1)
  df$mod <- f_mod %>% deparse()

  print(df$mod)
  print(a_mod)
  return(df[1,])
}

#### Loop through presence - absence models p2 ####
# min.mods <- read.csv("mod_abs_p1_bst.csv")
# mods <- read.csv("mod_abs_p1_reff.csv")
mods2 <- min.mods
mods2[1:dim(mods2)[1],3:dim(mods2)[2]] <- NA

for (i in 1:dim(mods2)[1]) {
  mods2[i,] <- try({run_mod2(i)}, silent = TRUE)  # Usar silent = TRUE para suprimir mensajes de error
}

# Generate and export summary tables
write.csv(mods2, file = "mod_abs_p2_bst.csv", row.names = F)

### data.frame models to evaluate effects for two grades of severity of fire p1 ####
ns1 <- 12 # essays
mods_sev_s <- array(NA,c(lom*ns1,12))
mods_sev_s <- as.data.frame(mods_sev_s)
names(mods_sev_s) <- c('var','dis.y','REML','AIC','BIC', 'mod','fe_p_names','fe_p_vals','re_names','re_p_vals', 'R2m', 'R2c')

mods_sev_s$dis.y <- rep(dis$dist, each=ns1)
mods_sev_s$var <- rep(dis$var_ori, each=ns1)

#### Function to run the model, soil severity 2 levels p1  ####
sev_model_2 <- function(fixed_effect = "sev_dpt", i, n) {
  ll <- 12 # number of columns
  df <- setNames(data.frame(matrix(ncol = ll, nrow = 1)), c('var','dis.y','REML','AIC','BIC', 'mod','fe_p_names','fe_p_vals','re_names','re_p_vals', 'R2m', 'R2c'))
  df$var <- dis$var_ori[i]
  df$dis.y <- dis$dist[i]
  
  df_mod <- get(dis$bd[fi]) %>% # podría ser i en lugar de fi
    filter(!is.na(!!sym(fixed_effect)))

  # Select random effects based on n
  raneff <- switch(as.character(n),
                   "1" = "(1|Zonas)",
                   "2" = "(1|clim:Zonas)",
                   "3" = "(1|veg_1:Zonas)",
                   "4" = "(1|veg_2:Zonas)",
                   "5" = "(1|veg_1)",
                   "6" = "(1|veg_2)",
                   "7" = "(1|clim:asl)",
                   "8" = "(1|asl:clim:Zonas)",
                   "9" = "(1|veg_2:veg_1)",
                   "10" = "(1|asl:clim:veg_1:Zonas)",
                   "11" = "(1|asl:clim:veg_2:Zonas)")

  # Add fixed effects and set fraction when applicable
  if (dis$bd[i] == "ftir") {
    fixeff <- fixed_effect
  } else {
    fixeff <- paste("fraction *",fixed_effect)
  }

  formula_str <- paste(dis$var[i], "~", fixeff)

  if (n == 12) {
    mod.0 <- gls(as.formula(formula_str), data = df_mod)
    df$AIC <- AIC(mod.0) %>% round(1)
    df$BIC <- BIC(mod.0) %>% round(1)
    df$mod <- formula_str

    df$fe_p_names <- mod.0 %>% anova() %>%
      row.names() %>% .[-1] %>%
      paste(collapse = ', ')
    df$fe_p_vals <- mod.0 %>% anova() %>%
      pull('p-value') %>% .[-1] %>% round(5) %>%
      paste(collapse = ', ')

    return(df[1,])  }

  formula_str <- paste(formula_str, "+", raneff)
  df$mod <- formula_str

  # Generate model according to variable distribution in dis$dist[i]
  if (dis$dist[i] %in% c("normal","log-normal") ) {
    mod.cn <- lmer(as.formula(formula_str), data = df_mod)
    s.mod <- summary(mod.cn)
    df$REML <- s.mod[["AICtab"]] %>% round(1)

    reff <- mod.cn %>% ranova()
    df$re_names <- reff %>%
      row.names() %>%
      .[-1] %>% paste(collapse = ', ')
    df$re_p_vals <- reff %>%
      pull('Pr(>Chisq)') %>%
      .[-1] %>% round(5) %>% paste(collapse = ', ')
  } else if (dis$dist[i] %in% c("gamma", "exponential") ) {
    mod.cn <- glmer(as.formula(formula_str),
                    data = df_mod, family = Gamma(link='log'))
  }

  df$AIC <- AIC(mod.cn) %>% round(1)
  df$BIC <- BIC(mod.cn) %>% round(1)

  a.mod <- mod.cn %>% Anova()
  df$fe_p_names <- a.mod %>% row.names() %>%
    paste(collapse = ', ')
  df$fe_p_vals <- a.mod %>% pull('Pr(>Chisq)') %>%
    round(5) %>% paste(collapse = ', ')

  df[1, c('R2m', 'R2c')] <- mod.cn %>% r.squaredGLMM() %>% round(4)
  
  return(df[1,])
}

#### Loop trough models with sev_dpt & sev_veg as fixed variables ####
for (i in 1:lom) {
  for (n in 1:ns1) {
    if (dis$dist[i] == "uniform") { # 
      next  # Omitir y saltar al siguiente valor del bucle
    }
    mods_sev_s[ns1*i-(ns1-n),] <- try({sev_model_2("sev_dpt", i, n)}, silent = T)
  }
}

### Generate and export summary tables sev_dpt p1 ####
mods_sev_s <- subset(mods_sev_s, !is.na(AIC))
min.sev_s <- min_mod(mods_sev_s)

write.csv(mods_sev_s, file = "mod_sev_p1_reff.csv", row.names = F)
write.csv(min.sev_s, file = "mod_sev_p1_bst.csv", row.names = F)

### Function to run the model, soil severity 2 levels p2  ####
# min.sev_s <- read.csv("mod_sev_p1_bst.csv")
# mods_sev_s <- read.csv("mod_sev_p1_reff.csv")
sev_mod_p2 <- function(i) { #fixed_effect no definido en la función
  ll <- dim(min.sev_s)[2]
  fi <- which(dis$var_ori == min.sev_s$var[i])
  df <- setNames(data.frame(matrix(ncol = ll, nrow = 1)), names(min.sev_s))
  df$var <- min.sev_s$var[i]
  df$dis.y <- min.sev_s$dis.y[i]
  
  df_mod <- get(dis$bd[fi]) %>% 
    filter(!is.na(!!sym(fixed_effect)))
  
  # Extract the fixed and random effects
  formula_str <- min.sev_s$mod[i]
  
  ran_eff <- ifelse(grepl("\\(1\\|", formula_str), TRUE, FALSE)
  
  # Check for the random effect, if any doing gls(method = 'ML')
  if(!ran_eff){
    mod <- gls(as.formula(formula_str), data = df_mod,
               method = 'ML')
    # mod$call$data <- dis$bd[fi] %>% as.name()
    # mod$call$model <- formula_str %>% as.formula() %>% as.call()
    
    a_mod <- mod %>% dredge(options(na.action = 'na.fail'))
    b_mod <- a_mod %>% get.models(1)
    b_mod <- b_mod[[1]]
    f_mod <- b_mod %>% formula()
    
    mod_out <- gls(f_mod, data = df_mod)
    fix_eff <- ifelse(grepl("~\\s1", f_mod %>% deparse()), FALSE, TRUE)
  }else{
    ran_eff <- gsub(".*?(\\(1\\|.*)", "\\1", formula_str)
    
    if (dis$dist[fi] %in% c("normal","log-normal") ) {
      mod <- lmer(as.formula(formula_str), data = df_mod,
                  REML = FALSE)
    } else if (dis$dist[fi] %in% c("gamma", "exponential") ) {
      mod <- glmer(as.formula(formula_str), 
                   data = df_mod, family = Gamma(link='log'))
    } 
    # mod@call$data <- dis$bd[fi] %>% as.name()
    
    a_mod <- mod %>% dredge(options(na.action = 'na.fail'))
    b_mod <- a_mod %>% get.models(1)
    b_mod <- b_mod[[1]]
    f_mod <- b_mod %>% formula()
    
    fix_eff <- ifelse(grepl("~\\s\\s*[^()]+\\s*\\(1", f_mod %>% deparse()), TRUE, FALSE)
    
    if (dis$dist[fi] %in% c("normal","log-normal") ){
      mod_out <- lmer(f_mod, data = df_mod)
      s.mod <- summary(mod_out)
      df$REML <- s.mod[["AICtab"]] %>% round(1)
      reff <- mod_out %>% ranova()
      df$re_names <- reff %>% 
        row.names() %>% 
        .[-1] %>% paste(collapse = ', ')
      df$re_p_vals <- reff %>% 
        pull('Pr(>Chisq)') %>% 
        .[-1] %>% round(5) %>% paste(collapse = ', ')
    } else if (dis$dist[fi] %in% c("gamma", "exponential") ) {
      mod_out <- b_mod
    }
    df[1, c('R2m', 'R2c')] <- mod_out %>% r.squaredGLMM() %>% round(4)
  }
  if (fix_eff){
    anova_mod <- mod_out %>% Anova()
    df$fe_p_names <- anova_mod %>% row.names() %>%
      paste(collapse = ', ')
    df$fe_p_vals <- anova_mod %>% pull('Pr(>Chisq)') %>%
      round(5) %>% paste(collapse = ', ')
  }
  
  df$AIC <- AIC(mod_out) %>% round(1)
  df$BIC <- BIC(mod_out) %>% round(1)
  df$mod <- f_mod %>% deparse()
  
  print(df$mod)
  print(a_mod)
  return(df[1,])
}

#### Table and loop for sev_dpt p2 ####
mods_sev_s2 <- min.sev_s
mods_sev_s2[1:dim(mods_sev_s2)[1],3:dim(mods_sev_s2)[2]] <- NA

for (i in 32:dim(mods_sev_s2)[1]) { 
  # df_mod <- get(dis$bd[fi]) %>%  # si da problemas descomentar
  #   filter(!is.na(!!sym(fixed_effect)))
  mods_sev_s2[i,] <- try({sev_mod_p2(i)}, silent = TRUE)  # Usar silent = TRUE para suprimir mensajes de error
}

# Generate and export summary tables #
write.csv(mods_sev_s2, file = "mod_sev_p2_bst.csv", row.names = F)

## STEP 3. Correlation plots and analysis for repellency data ####
### environment ####
library(lme4)
library(lmerTest)
library(car)
library(nlme)
library(MuMIn)
library(emmeans)
library(vegan)
library(glmmTMB)

### Generation of data.frame for model ####
nessays <- 12
mods_rel <- array(NA,c((dim(response)[1])*lom*nessays,11)) # the last number must correspond with length(names(mods_rel)) ; dim(response)[1]+2
mods_rel <- as.data.frame(mods_rel)

names(mods_rel) <- c('var','dis.y','AIC','BIC', 'mod','fe_p_names','fe_p_vals','re_names','re_p_vals', 'R2m', 'R2c')

mods_rel$dis.y <- rep(response$dist, each=nessays*lom)
mods_rel$var <- rep(response$var, each=nessays*lom)
mods_rel$x_var <- rep(dis$var_ori, by = nessays*9)
  
### Functions to extract effects on repellency+aggregation+microbiology ####
rel_model <- function(i, n, j) { #i-response, n-essays/raneff, j=dis/xvar, mods_rel = mods_rel
  row_fi <- lom*(nessays*(i-1)+n-1)+j
  df <- setNames(data.frame(matrix(ncol = 9, nrow = 1)), c('AIC','BIC', 'mod','fe_p_names','fe_p_vals','re_names','re_p_vals', 'R2m', 'R2c'))

    df_mod <- sd_fire_db %>% 
      filter(!is.na(!!sym(dis$var_ori[j])) & !is.na(!!sym(response$var_ori[i])))
  
  # Select random effects based on n
  raneff <- switch(as.character(n),
                   "1" = "(1|Zonas)",
                   "2" = "(1|clim:Zonas)",
                   "3" = "(1|veg_1:Zonas)",
                   "4" = "(1|veg_2:Zonas)",
                   "5" = "(1|veg_1)",
                   "6" = "(1|veg_2)",
                   "7" = "(1|clim:asl)",
                   "8" = "(1|asl:clim:Zonas)",
                   "9" = "(1|veg_2:veg_1)",
                   "10" = "(1|asl:clim:veg_1:Zonas)",
                   "11" = "(1|asl:clim:veg_2:Zonas)")
  
  # Add fixed effects and set fraction when applicable
  if (dis$bd[j] == "ftir") {
    fixeff <- dis$var_ori[j]
  } else {
    fixeff <- paste(dis$var_ori[j],": fraction")
  }
  
  formula_str <- paste(response$var[i], "~", fixeff)
  
  if (n == 12) {
    
    mod.0 <- gls(as.formula(formula_str), data = df_mod)
    df$AIC <- AIC(mod.0) %>% round(1)
    df$BIC <- BIC(mod.0) %>% round(1)
    df$mod <- formula_str
    
    df$fe_p_names <- mod.0 %>% anova() %>%
      row.names() %>% .[-1] %>%
      paste(collapse = ', ')
    df$fe_p_vals <- mod.0 %>% anova() %>%
      pull('p-value') %>% .[-1] %>% round(5) %>%
      paste(collapse = ', ')
    
    return(df[1,])  }
  
  formula_str <- paste(formula_str, "+", raneff)
  df$mod <- formula_str
  
  # Generate model according to variable distribution in dis$dist[i]
  if (response$dist[i] != "beta"){ # & response$dist[i] != "gamma"
    if (response$dist[i] %in% c("normal","log-normal") ) {
      mod.cn <- lmer(as.formula(formula_str), data = df_mod)
      reff <- mod.cn %>% ranova()
      df$re_names <- reff %>%
        row.names() %>%
        .[-1] %>% paste(collapse = ', ')
      df$re_p_vals <- reff %>%
        pull('Pr(>Chisq)') %>%
        .[-1] %>% round(5) %>% paste(collapse = ', ')
    } else {
      if (response$dist[i] == 'binomial') {
        mod.cn <- glmer(as.formula(formula_str),
                        data = df_mod, family = 'binomial')
      } else if (response$dist[i] %in% c("gamma", "exponential") ) {
        mod.cn <- glmer(as.formula(formula_str),
                        data = df_mod, family = Gamma(link='log'))
        }
    }
    df[1, c('R2m', 'R2c')] <- mod.cn %>% r.squaredGLMM() %>% round(4)
    
  } else {
      mod.cn <- glmmTMB(as.formula(formula_str),
                        data = df_mod, family = beta_family())
  }
  a.mod <- mod.cn %>% Anova()
  df$fe_p_names <- a.mod %>% row.names() %>%
    paste(collapse = ', ')
  df$fe_p_vals <- a.mod %>% pull('Pr(>Chisq)') %>%
    round(5) %>% paste(collapse = ', ')
  
  
  df$AIC <- AIC(mod.cn) %>% round(1)
  df$BIC <- BIC(mod.cn) %>% round(1)
  
  return(df[1,])
  }

rel_model_2 <- function(i) {
  var_y <- min.rel$var[i]
  var_x <- min.rel$x_var[i]
  fi <- which(dis$var_ori == var_x)
  fj <- which(response$var == var_y)
  
  df_mod <- full_db %>% 
    filter(!is.na(!!sym(var_x)) & 
             !is.na(!!sym(response$var_ori[fj])))
  
  # Extract the fixed and random effects
  formula_str <- min.rel$mod[i]
  
  ran_eff <- ifelse(grepl("\\(1\\|", formula_str), TRUE, FALSE)
  
  # Check for the random effect, if any doing gls()
  if(!ran_eff){
    mod <- gls(as.formula(formula_str), data = df_mod)
  }else{
    if (response$dist[fj] != "beta"){
      if (response$dist[fj] %in% c("normal","log-normal") ) {
        mod <- lmer(as.formula(formula_str), data = df_mod)
      } else {
        if (response$dist[fj] == 'binomial') {
          mod <- glmer(as.formula(formula_str),
                       data = df_mod, family = 'binomial')
        } else if (response$dist[fj] %in% c("gamma", "exponential") ) {
          mod <- glmer(as.formula(formula_str), 
                       data = df_mod, family = Gamma(link='log'))
        }
      }
    } else {
      mod <- glmmTMB(as.formula(formula_str),
                     data = df_mod, family = beta_family())
    }
  }
  
  if(dis$bd[fi] == "ftir"){
    a.mod <- mod %>% 
      emtrends(as.formula(paste('~', var_x)), 
               var = var_x, adjust = "sidak", data = df_mod) %>%
                 test()
    a.mod[1] <- "MAOM"
    names(a.mod)[1] <- "fraction"
  } else {
    a.mod <- mod %>% 
      emtrends(~ fraction, var = var_x, 
               adjust = "sidak", data = df_mod) %>%
      test()
  }
  
  a.mod$var <- response$var_ori[fj]
  a.mod$x_var <- var_x
  a.mod$dist <- response$dist[fj]
  names(a.mod)[2] <- 'trend'
  return(a.mod)
}

### Loop trough models for repellency+aggregation+microbiology ####
system.time({ # if you uncomment this also do it at +14 lines

for (i in 1:dim(response)[1]) { # 5:6
  for (n in 1:nessays) { #
    for(j in 1:lom){ # 2
      if (dis$dist[j] == "uniform") { # 
        next  # Omitir y saltar al siguiente valor del bucle
        }
      row_fil <- lom*(nessays*(i-1)+n-1)+j
      mods_rel[row_fil,3:11] <- try({rel_model(i, n, j)}, silent = TRUE) #, mods_rel
    }
  }
} #; View(mods_rel)

})

### ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
# i = 1:9 with nessays and lom
# (3780 models in 8 min)
#   user  system elapsed 
# 461.91    2.51  467.60
#
# 0.12 sec by model
# Ran on 01/aug/2025
### ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

### Generate and export summary tables for repellency+aggregation+microbiology ####
write.csv(mods_rel, file = "mod_rel_all_na.csv", row.names = FALSE)
# mods_rel <- read.csv("mod_rel_all_na.csv")
mods_rel <- subset(mods_rel, !is.na(AIC))
min.rel <- min_mod2(mods_rel)

write.csv(mods_rel, file = "mod_rel_all.csv", row.names = FALSE)
write.csv(min.rel, file = "mod_rel_bst_p1.csv", row.names = FALSE)

### Generate and export comparison analysis checking slopes p2 ####
min.rel <- read.csv("mod_rel_bst_p1.csv")
ls_min.rel <- lapply(1:dim(min.rel)[1], rel_model_2)
min.relp <- do.call(rbind, ls_min.rel)
min.relp <- min.relp %>% 
  dplyr::select(-Amide:-Aromatic)

write.csv(min.relp, file = "mod_rel_bst_p2.csv", row.names = FALSE)
## STEP 4 - Post-hoc analysis ####
### environment ####
library(nlme) # a esta función no le gustan los NA's gls lmer
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(multcomp)

### emmeans ####
mods2 <- read.csv("mod_abs_p2_bst.csv")
mods_sev_s2 <- read.csv("mod_sev_p2_bst.csv")

lista_burned <- as.list(mods2$var); names(lista_burned) <- mods2$var
lista_sev_s <- as.list(mods_sev_s2$var); names(lista_sev_s) <- mods_sev_s2$var

sum_emm <- function(i, fixed_effect = "sev_dpt", mods2 = mods_sev_s2) {
  fi <- which(dis$var_ori == mods2$var[i])
  formula_str <- mods2$mod[i]
  
  ran_eff <- grepl("\\(1", formula_str) #true or false
  
  if (mods2$var[i] %in% names(ftir)) {
    fixef <- fixed_effect
    formu <- paste("~", fixed_effect)
  } else {
    fixef <- paste(fixed_effect, "* fraction")
    formu <- paste("~", fixed_effect, "| fraction")
  }
  formula_str <- paste(mods2$var[i], "~", fixef, "+", mods2$re_names[i])
  
  if (fixed_effect == "sev_dpt") {
    df_mod <- get(dis$bd[fi]) %>% 
      filter(burned == "yes")
  } else {
    df_mod <- get(dis$bd[fi])
  }
  
  if (!ran_eff) {
    mod <- gls(as.formula(formula_str), data = df_mod)
  } else {
    dist_mod <- dis$dist[fi]
    if (dist_mod %in% c("normal", "log-normal")) {
      mod <- lmer(as.formula(formula_str), data = df_mod)
    } else if (dist_mod %in% c("gamma", "exponential")) {
      mod <- glmer(as.formula(formula_str), data = df_mod, family = Gamma(link='log'))
    } else {
      stop("Distribución no soportada para modelo con efectos aleatorios")
    }
  }
  
  emmeans_res <- emmeans(mod, specs = as.formula(formu), 
                         adjust = "sidak") %>% cld()
  
  emmeans_res <- data.frame(emmeans_res)
  emmeans_res$var <- mods2$var[i]
  if(!("fraction" %in% names(emmeans_res))) 
    {emmeans_res$fraction <- 'MAOM'}
  emmeans_res$dis.y <- mods2$dis.y[i]
  emmeans_res$formula <- formula_str
  
  print(formula_str)
  return(emmeans_res)
}
list_sum_br <- list()
list_sum_sv <- list()
for(i in 1:dim(mods2)[1]){
  if(dis$var_ori[i] == "ww") next
  print(i)
  list_sum_br[[i]] <- sum_emm(i,fixed_effect = "burned", mods2 = mods2)
  list_sum_sv[[i]] <- sum_emm(i,fixed_effect = "sev_dpt", mods2 = mods_sev_s2)
}  

names(list_sum_br) <- names(list_sum_sv) <- mods2$var
list_sum_br$ww <- list_sum_sv$ww <- NULL

df_sum_br <- list_sum_br %>%
  imap(~ mutate(.x, var = .y)) %>%
  bind_rows()
df_sum_sv <- df_sv <- list_sum_sv %>%
  imap(~ mutate(.x, var = .y)) %>%
  bind_rows()
trans_emm <- function(df) {
  df %>%
    mutate(
      emm = if_else(dis.y %in% c("normal", "log-normal"), emmean, exp(emmean)),
      se = if_else(dis.y %in% c("normal", "log-normal"), SE, SE * emm),
      LCL = if_else(dis.y %in% c("normal", "log-normal"), lower.CL, exp(asymp.LCL)),
      UCL = if_else(dis.y %in% c("normal", "log-normal"), upper.CL, exp(asymp.UCL)),
      emm = if_else(var %in% names(nmr_fg), emm / 100, emm),
      se = if_else(var %in% names(nmr_fg), se / 100, se),
      UCL = if_else(var %in% names(nmr_fg), UCL / 100, UCL),
      LCL = if_else(var %in% names(nmr_fg), LCL / 100, LCL),
      fraction = factor(fraction, levels = c("fPOM", "oPOM", "MAOM")),
      var = factor(var, levels = mods2$var)
      )
}
df_sum_br <- trans_emm(df_sum_br) %>% 
  mutate(
    burned = factor(burned, levels = c("no","yes"),
                    labels = c("Unburned", "Burned"))) %>% 
  arrange(var, fraction, burned)
df_sum_sv <- trans_emm(df_sum_sv) %>% 
  mutate(
    sev_dpt = factor(sev_dpt, levels = c("low", "medium-high"))) %>% 
  arrange(var, fraction, sev_dpt)

### Saving the emmeans datasets ##
write.csv(df_sum_br, "emm_abs.csv", row.names = F)
write.csv(df_sum_sv, "emm_sev.csv", row.names = F)

### post-hocs emmeans ####
post_hocs <- function(i, fixed_effect = "sev_dpt", 
                      mods2 = mods_sev_s2) {
  fi <- which(dis$var_ori == mods2$var[i])
  formula_str <- mods2$mod[i]
  
  ran_eff <- ifelse(grepl("\\(1", formula_str), TRUE, FALSE)
  fix_eff <- !(grepl("~\\s1", formula_str) | grepl("~\\s\\(1", formula_str))

  fixef <- gsub(".*~\\s*(.*?)\\s*\\+\\s*\\(?1.*", "\\1", formula_str)
  
  # Check for the fixed effect, if any return NA, if not continue
  if(!fix_eff){
    return(NA)
  }else{
    if(fixed_effect == "sev_dpt"){
      df_mod <- get(dis$bd[fi]) %>%
        filter(burned == "yes")
    }else{
      df_mod <- get(dis$bd[fi])
    }
    # df_mod <- get(dis$bd[fi]) %>%
    #   filter(!is.na(!!sym(fixed_effect)))
    
    # Check for the random effect
    if(!ran_eff){
      mod <- gls(as.formula(formula_str), data = df_mod)
    }else{
      if (dis$dist[fi] %in% c("normal","log-normal") ) {
        mod <- lmer(as.formula(formula_str), data = df_mod)
      } else if (dis$dist[fi] %in% c("gamma", "exponential") ) {
        mod <- glmer(as.formula(formula_str), 
                     data = df_mod, family = Gamma(link='log'))
      } 
    }
    
    group_mod <- emmeans(mod, as.formula(paste("~", fixef)), adjust = "sidak") %>%
      cld()
    group_mod$var <- mods2$var[i]
    print(formula_str)
    return(group_mod)
  }
  
}

# lista_burned <- lapply(seq_along(lista_burned), function(i)post_hocs(i,fixed_effect = "burned", mods2 = mods2))
# lista_sev_s <- lapply(seq_along(lista_sev_s), post_hocs(fixed_effect = "sev_dpt", mods2 = mods_sev_s2))

for(i in c(1:6,8:31,33:dim(mods2)[1])){ # 7 y 32 dan mucha lata
  print(i)
  lista_burned[[i]] <- post_hocs(i,fixed_effect = "burned", mods2 = mods2)
  lista_sev_s[[i]] <- post_hocs(i,fixed_effect = "sev_dpt", mods2 = mods_sev_s2)
}

for(i in c(7,32)){
  fixef <- gsub(".*~\\s*(.*?)\\s*\\+\\s*\\(?1.*", "\\1", mods2$mod[i])
  lista_burned[[i]] <- gls(as.formula(mods2$mod[i]), data = subset(get(dis$bd[i]),!is.na("burned"))) %>% 
    emmeans(as.formula(paste("~", fixef)), adjust = "sidak") %>%
    cld()
  lista_burned[[i]]$var <- mods2$var[i]
  
  fixefs <- gsub(".*~\\s*(.*?)\\s*\\+\\s*\\(?1.*", "\\1", mods_sev_s2$mod[i])
  lista_sev_s[[i]] <- gls(as.formula(mods_sev_s2$mod[i]), data = get(dis$bd[i]) %>% filter(!is.na(!!sym("sev_dpt"))) ) %>% 
    emmeans(as.formula(paste("~", fixefs)), adjust = "sidak") %>%
    cld()
  lista_sev_s[[i]]$var <- mods_sev_s2$var[i]
}

non_b <- lista_burned[!sapply(lista_burned, is.data.frame)] %>% names()
non_s <- lista_sev_s[!sapply(lista_sev_s, is.data.frame)] %>% names()
lista_burned <- lista_burned[sapply(lista_burned, is.data.frame)]
lista_sev_s <- lista_sev_s[sapply(lista_sev_s, is.data.frame)]

post_abs <- bind_rows(lista_burned)
post_sev <- bind_rows(lista_sev_s)

post_abs <- post_abs %>% relocate(var, burned, fraction, .group, .before = 1)
post_sev <- post_sev %>% relocate(var, sev_dpt, fraction, .group, .before = 1)

### Save the post-hoc datasets ##
write.csv(post_abs, "post_abs.csv", row.names = F)
write.csv(post_sev, "post_sev.csv", row.names = F)

### contrasts ####
dif_em <- function(i, fixed_effect = "sev_dpt", 
                   mods2 = mods_sev_s2){
  fi <- which(dis$var_ori == mods2$var[i])
  formula_str <- mods2$mod[i] # formula original
  
  ran_eff <- ifelse(grepl("\\(1", formula_str), TRUE, FALSE)
  
  if(mods2$var[i] %in% names(ftir)){
    fixef <- fixed_effect
    formu <- paste("trt.vs.ctrl1 ~", fixed_effect)
  } else {
    fixef <- paste(fixed_effect, "* fraction")
    formu <- paste("trt.vs.ctrl1 ~", fixed_effect, "| fraction")
  }
  formula_str <- paste(mods2$var[i], '~', fixef, '+', mods2$re_names[i])
  
  if(fixed_effect == "sev_dpt"){
    df_mod <- get(dis$bd[fi]) %>%
      filter(burned == "yes")
  }else{
    df_mod <- get(dis$bd[fi])
  }
  
  if(!ran_eff){
    mod <- gls(as.formula(formula_str), data = df_mod)
  }else{
    if (dis$dist[fi] %in% c("normal","log-normal") ) {
      mod <- lmer(as.formula(formula_str), data = df_mod)
    } else if (dis$dist[fi] %in% c("gamma", "exponential") ) {
      mod <- glmer(as.formula(formula_str),
                   data = df_mod, family = Gamma(link='log'))
    }
  }
  
  con_mod1 <- emmeans(mod, specs = as.formula(formu), adjust = "sidak")
  contr_mod <- contrast(con_mod1, ref = "no", method = "revpairwise") %>% summary(type = "response")
  contr_mod <- contr_mod %>%
    mutate(
      estimate = if ("estimate" %in% colnames(.)) estimate else ratio ,
      upper = estimate + SE,
      lower = estimate - SE
    )
  
  print(formula_str)
  return(contr_mod)
}

list_d_br <- list()
list_d_sv <- list()

for(i in 1:dim(mods2)[1]){
  if(dis$var_ori[i] == "ww") next
  print(i)
  list_d_br[[i]] <- dif_em(i,fixed_effect = "burned", mods2 = mods2)
  list_d_sv[[i]] <- dif_em(i,fixed_effect = "sev_dpt", mods2 = mods_sev_s2)
}  

names(list_d_br) <- names(list_d_sv) <- mods2$var
list_d_br$ww <- list_d_sv$ww <- NULL

df_d_br <- list_d_br %>%
  imap(~ mutate(.x, var = .y)) %>%
  bind_rows()
df_d_sv <- df_sv <- list_d_sv %>%
  imap(~ mutate(.x, var = .y)) %>%
  bind_rows()

df_d_br <- df_d_br %>%
  mutate(
    estimate = if_else(var %in% names(nmr_fg), estimate / 100, estimate),
    SE = if_else(var %in% names(nmr_fg), SE / 100, SE),
    upper = if_else(var %in% names(nmr_fg), upper / 100, upper),
    lower = if_else(var %in% names(nmr_fg), lower / 100, lower),
    fraction = replace_na(fraction, "MAOM")
  )

df_d_sv <- df_d_sv %>%
  mutate(
    estimate = if_else(var %in% names(nmr_fg), estimate / 100, estimate),
    SE = if_else(var %in% names(nmr_fg), SE / 100, SE),
    upper = if_else(var %in% names(nmr_fg), upper / 100, upper),
    lower = if_else(var %in% names(nmr_fg), lower / 100, lower),
    fraction = replace_na(fraction, "MAOM")
  )


