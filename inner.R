grid.temp <- seq(45, 90, 0.1)

### Setup frame to store all data generated
Plasma.Final.Data.Inner <- data.frame(Temperature = grid.temp)

### loop over all samples and store the dCp with chosen temperature grid
for(j in 1:n.samples)
{
  cat('Working on Sample ', All.IDs[j], 'element ', j,' of', n.samples,' \n')
  ### select a sample
  working.sample <- final.plasma %>% 
    filter(SampleID == All.IDs[j]) %>% 
    select(Temperature, dCp)
  ### get a baseline-subtracted and interpolated final result!
  auto.output <- auto.baseline(x = working.sample, grid.temp = grid.temp, point = "innermost")
  Plasma.Final.Data.Inner <- Plasma.Final.Data.Inner %>% cbind(out = auto.output$dCp)
  cat("\014")
}

colnames(Plasma.Final.Data.Inner)[-1] <- All.IDs

thing.inner <- Plasma.Final.Data.Inner %>% pivot_longer(names_to="SampleID", values_to = "dCp", 2:73)

thing.inner <- thing.inner %>% mutate(SampleNumber = str_extract(SampleID, '\\d+'), SampleIteration = str_extract(SampleID, '\\D$'))

numbers <- thing.inner %>% pull(SampleNumber) %>% unique()

other.thing <- Plasma.Working.Corrected %>% filter(between(Temperature, 45, 90))

pdf('generated_output/Plasma_Final_Auto_Graphics_Eight_Inner.pdf')
{
  g<-NULL
  for(i in 1:length(numbers)){
    sample1 <- thing.inner %>% filter(SampleNumber == numbers[i]) %>% mutate(method = str_c("Auto",SampleIteration,sep=" ")) %>% mutate(type = "Auto")
    sample2 <- other.thing %>% filter(SampleNumber == numbers[i]) %>% mutate(method = str_c("Tech",SampleIteration,sep=" ")) %>% mutate(type = "Not Auto")
    
    combinded <- sample1 %>% rbind(sample2)
    
    g[[i]]<- ggplot(combinded, aes(Temperature, dCp, color = method, linetype = type))+
      geom_line()+
      scale_color_manual(breaks = c("Auto a", "Auto b", "Tech a", "Tech b"),
                         values=c("green1", "red1", "green4", "red4"))+
      labs(title = paste0('Final Automated Sample for ', str_sub(numbers[i], 1)))+
      theme(legend.position = "top",
            legend.key.height = unit(0.10, 'cm'),
            legend.key.width = unit(0.30, 'cm'),
            legend.box.spacing = unit(0.20, 'cm'),
            legend.margin = margin(1, 1, 1, 1),
            legend.title = element_text(size=5),
            legend.text=element_text(size=5),
            title = element_text(size = 8))
  }
  marrangeGrob(g,nrow=4,ncol=2)
}
dev.off()

endpoints.inner <- data.frame()
for(j in 1:n.samples)
{
  ### select a sample
  working.sample <- final.plasma %>% 
    filter(SampleID == All.IDs[j]) %>% 
    select(Temperature, dCp)
  ### get a baseline-subtracted and interpolated final result!
  auto.output <- moving.window(x = working.sample,point.selection = "innermost")
  endpoints.inner <- endpoints.inner %>% rbind(data.frame(sample = All.IDs[j],lower = auto.output[1], upper=auto.output[2]))
}

write.csv(x=endpoints.inner, file = 'generated_output/Plasma_Endpoints_Inner.csv')

new_thing.inner <- thing.inner

new_thing.inner <- new_thing.inner %>% mutate(new_title = str_c(SampleID,"Auto",sep="_"))

new_thing_wide.inner <- new_thing.inner[c("Temperature", "new_title", "dCp")] %>% pivot_wider(names_from = new_title, values_from = dCp)

final.plasma.inner <- inner_join(new_thing_wide.inner,data.set,by="Temperature")

write.csv(x=final.plasma.inner, file = 'generated_output/Plasma_Comparison_Inner.csv')
