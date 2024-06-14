if (!("pacman" %in% rownames(installed.packages()))) {install.packages("pacman")}
pacman::p_load(pacman,
               tidyverse,
               readxl,
               lubridate,
               mgcv,
               gratia,
               furrr,
               tictoc,
               magrittr,
               beepr,
               RcppRoll,
               scales)

source("Gam_functions.R")

N03N_filtered <- readRDS("~/LAWA/LAWAgit/N03N_filtered.rds")

testdf <- N03N_filtered[[62]]

#testdf$yearmon

#testdf$ymon <-  tsibble::yearmonth(testdf$myDate)

test_gam <- screeningmodeling(.data = testdf , datevar = myDate,values = RawValue )

plot_individual_trend(test_gam)

plot_individual_trend_log(test_gam)
###


N03N_filtered_df <- map_df(N03N_filtered,~.)

N03N_filtered_df %>%
  group_by(lawa_site_id, measurement) %>% 
  #filter(measurement == "bd") %>% 
  #drop_na( value) %>% #TRY WITH NAs too
  mutate(n=length(yearmon )) %>% 
  screeningmodeling(values=  RawValue ,
                    datevar = myDate , 
                    link = "identity", 
                    conf.type = "conf",
                    conf.level=0.95,
                    beep = TRUE, 
                    tdist = F,
                    autocor = TRUE,
                    lawa_site_id, measurement) -> 
  LAWA_screenmodel_out




LAWA_screenmodel_out %>%
  #arrange(desc(SWEREFID))%>%
  plot_screeningtrends(sorting =lawa_site_id,
                       y_id = lawa_site_id,
                       wrappingvar = measurement)


LAWA_screenmodel_out$fderiv[[1]] %>% tail()

LAWA_screenmodel_out$fderiv[[1]] %>%
  mutate(deriv = .derivative, 
         deriv_se = .se, 
         lower=.lower_ci, 
         upper=.upper_ci) %>%
  as_tibble %>%
  rowwise %>%
  mutate(signif = !between(0, lower, upper),
         sign=sign(deriv),
         signif_sign = signif*sign) %>% 
  view()
