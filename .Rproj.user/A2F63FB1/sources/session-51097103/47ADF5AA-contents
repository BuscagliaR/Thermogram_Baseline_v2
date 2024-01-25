Same_Inner_Endpoint_Graphs <- read.csv("~/GitHub/Thermogram_Baseline_v2/generated_output/Same_Inner_Endpoint_Graphs.csv")
Same_Inner_Endpoint_Graphs <- Same_Inner_Endpoint_Graphs[-1]

devtools::install_github("BuscagliaR/tlbparam")
library(tlbparam)
library(dplyr)
library(stringr)
?clean_thermograms()

untamptered_data <-read.csv("~/GitHub/Thermogram_Baseline_v2/generated_output/innermost_data.csv")

inner <- read.csv("~/GitHub/Thermogram_Baseline_v2/generated_output/innermost_data.csv")
inner <- inner[-1]

clean_thermograms(inner)
head(t(inner))
inner.t <- t(inner)


inner <- inner %>% mutate(Temperature = str_c("T", Temperature, sep=""))

colnames(inner.t) <- inner.t[1,]
inner.t <- inner.t[-1,]
inner.t <- as.data.frame(inner.t)
inner.t <- inner.t %>% mutate(SampleCode = row.names(inner.t))
thing <- inner.t$SampleCode

urine_char <- clean_thermograms(inner.t, column = "SampleCode",summary = c("tarea", "fwhm", "max", "tpeakf", "peakf", "tpeak1", "peak1", "tpeak2",
                                                                           "peak2", "tpeak3", "peak3", "P1P2.trough", "tP1P2.trough", "tmax", "tfm",
                                                                           "peak12ratio", "min", "tmin", "median"))
urine_char <- load_thermograms(inner.t)

inner.t <- inner.t %>% as.numeric(inner.t)

inner.t <- data.frame(lapply(inner.t[,1:451], function(x) as.numeric(as.character(x))))
inner.t$SampleCode <- thing
str(inner.t)

data <- urine_char[,c(3,456:474)]


library(ggplot2)
first_obs <- untamptered_data %>% select(c("Temperature", "X1a"))
ggplot(first_obs, aes(x=Temperature,y=X1a))+geom_point()
