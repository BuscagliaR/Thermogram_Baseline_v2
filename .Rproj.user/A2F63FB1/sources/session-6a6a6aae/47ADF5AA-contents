Same_Inner_Endpoint_Graphs <- read.csv("~/GitHub/Thermogram_Baseline_v2/generated_output/Same_Inner_Endpoint_Graphs.csv")
Same_Inner_Endpoint_Graphs <- Same_Inner_Endpoint_Graphs[-1]

devtools::install_github("BuscagliaR/tlbparam")
library(tlbparam)
library(dplyr)
library(stringr)


untamptered_data <-read.csv("~/GitHub/Thermogram_Baseline_v2/generated_output/innermost_data.csv")

inner <- read.csv("~/GitHub/Thermogram_Baseline_v2/generated_output/innermost_data.csv")
inner <- inner[-1]

clean_thermograms(inner)
head(t(inner))
inner <- inner %>% mutate(Temperature = str_c("T", Temperature, sep=""))

inner.t <- t(inner)



colnames(inner.t) <- inner.t[1,]
inner.t <- inner.t[-1,]
inner.t <- as.data.frame(inner.t)
inner.t <- inner.t %>% mutate(SampleCode = row.names(inner.t))
thing <- inner.t$SampleCode
inner.t$SampleCode <- thing
inner.t <- data.frame(lapply(inner.t[,1:451], function(x) as.numeric(as.character(x))))
inner.t$SampleCode <- thing

inner.t <- inner.t %>% select("SampleCode",everything())

urine_char <- clean_thermograms(inner.t, column = "SampleCode",summary = c('Width','Area','Max','TMax','TFM','Peak 1','Peak 2','Peak 3','TPeak 1','TPeak 2',
                                                                           'TPeak 3','Peak 1 / Peak 2', 'Peak 1 / Peak 3', 'Peak 2 / Peak 3',
                                                                           'Median', 'V1.2', 'TV1.2', 'V1.2 / Peak 1', 'V1.2 / Peak 2', 'V1.2 / Peak 3',
                                                                           'Min', 'TMin', 'Peak F', 'TPeak F'))

data <- urine_char[,c(3,456:479)]

urine_classes <- readxl::read_excel("~/GitHub/Thermogram_Baseline_v2/data_raw/Lups urine classification 2.1.2024.xlsx", sheet = "Sheet1")
urine_classes <- urine_classes[,c(1,11)]
urine_classes$SampleNumber <- as.numeric(urine_classes$Sample)

data <- data %>% mutate(SampleNumber = as.numeric(str_extract(SampleCode, '\\d+')))

final_table <- left_join(urine_classes, data, by = "SampleNumber")
final_table <- final_table[!is.na(final_table$Disease),]

###ANOVA###
for(i in 5:(ncol(final_table)+5))
  test <- aov(final_table[i] ~ Disease, data = final_table)
  summary(test)

###Graphing###

library(ggplot2)
first_obs <- untamptered_data %>% select(c("Temperature", "X1a"))
ggplot(first_obs, aes(x=Temperature,y=X1a))+geom_point()

Sample_Names <- colnames(untamptered_data)[-1:-2]

pdf('generated_output/graphs.pdf')
{
  for(j in 1:length(Sample_Names))
  {
    g1 <- untamptered_data %>% select(Temperature, Sample_Names[j]) %>% 
      ggplot(aes(x = Temperature, y = .[,2])) + 
      geom_line() + 
      labs(title = paste0('Final Automated Sample for ', str_sub(Sample_Names[j], 2)))
    print(g1)
  }
}
dev.off()

write.csv(final_table, file = "generated_output/characteristics_urine.csv")
save.image("generated_output/characteristics_urine.R")
