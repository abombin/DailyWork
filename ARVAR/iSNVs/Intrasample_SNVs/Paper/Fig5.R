library(ggplot2)





combShannon = read.csv("IntraSnv_results/AllSamplesCombShannon.csv")
combShannon$NormShan = combShannon$Shannon
combShannon$NormShan[combShannon$NormShan == 0] = 0.005900276
combShannon$NormShan = scale(log(combShannon$NormShan, base = 10))


combFilt = combShannon[combShannon$days_post_symptom_onset < 100 & combShannon$days_post_symptom_onset  >=0,]
# scatter_quadratic <- ggplot(combFilt, aes(x = days_post_symptom_onset, y = NormShan)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula  = y ~ x + I(x^2), se = FALSE, color = "red")


symptomnms = ggplot(combFilt , aes(x = days_post_symptom_onset, y = NormShan)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + theme_classic() + 
  theme(text = element_text(size = 24))

ggsave(file = "Paper/Figs/Shannon_Symptoms.png", plot = symptomnms, height = 12, width = 18, units = "in", dpi = 300)

cor.test(combFilt$days_post_symptom_onset, combFilt$Shannon, method = "spearman")
modle1 = lm(NormShan~days_post_symptom_onset, data = combFilt)
summary(modle1)

# vaccination
daysVac = ggplot(combShannon , aes(x = days_since_last_vax, y = NormShan)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + theme_classic() + 
  theme(text = element_text(size = 24))

daysVacSq = ggplot(combFilt, aes(x = days_since_last_vax, y = NormShan)) +
geom_point() +
geom_smooth(method = "lm", formula  = y ~ x + I(x^2), se = FALSE, color = "red")+  theme_classic() + 
  theme(text = element_text(size = 24))

ggsave(file = "Paper/Figs/Shannon_DaysVax.png", plot = daysVac, height = 12, width = 18, units = "in", dpi = 300)
ggsave(file = "Paper/Figs/Shannon_DaysVaxSq.png", plot = daysVacSq, height = 12, width = 18, units = "in", dpi = 300)

modle1 = lm(NormShan~days_since_last_vax + I(days_since_last_vax^2), data = combFilt)
summary(modle1)

# 

combFilt = combShannon[combShannon$vax_doses_received < 5 ,]
combShannon$vax_doses_received = as.numeric(as.character(combShannon$vax_doses_received))
 VaxDoses <- ggplot(combShannon, aes(x = vax_doses_received, y = NormShan)) +
    geom_point() +
   geom_smooth(method = "lm", formula  = y ~ x + I(x^2), se = FALSE, color = "red")
# 
 VaxDosesLin <- ggplot(combShannon, aes(x = vax_doses_received, y = NormShan)) +
   geom_point() +
  geom_smooth(method = "lm",  se = FALSE, color = "blue") +
   theme_classic() + 
   theme(text = element_text(size = 24))

combFilt$vax_doses_received = as.character(combFilt$vax_doses_received)
combFilt= combFilt[!is.na(combFilt$vax_doses_received),]
VaxDoses <- ggplot(combFilt, aes(x=vax_doses_received, y=NormShan, fill = vax_doses_received)) + 
  geom_violin(trim=FALSE)+
  geom_jitter() + 
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 18,        # Use the shape code for a diamond
    size = 8,          # Adjust the size of the diamond
    color = "white"  # Adjust dodge width if necessary
  ) + theme_classic() + theme(text = element_text(size = 26))

ggsave(file = "Paper/Figs/Shannon_VaxDoses.png", plot = VaxDoses, height = 10, width = 14, units = "in", dpi = 300)
ggsave(file = "Paper/Figs/Shannon_VaxDosesLin.png", plot =  VaxDosesLin, height = 12, width = 18, units = "in", dpi = 300)

# disease severity

disSevLin = ggplot(combShannon , aes(x = disease_severity, y = NormShan)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") + theme_classic() + 
  theme(text = element_text(size = 24))

combShannon$disease_severity = as.character(combShannon$disease_severity)
disSev = ggplot(combShannon, aes(x=disease_severity, y=NormShan, fill = disease_severity)) + 
  geom_violin(trim=FALSE)+
  geom_jitter() + 
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 18,        # Use the shape code for a diamond
    size = 8,          # Adjust the size of the diamond
    color = "white"  # Adjust dodge width if necessary
  ) + theme_classic() + theme(text = element_text(size = 26))

ggsave(file = "Paper/Figs/Shannon_Disease_severity.png", plot = disSev, height = 10, width = 14, units = "in", dpi = 300)
ggsave(file = "Paper/Figs/Shannon_Disease_severityLin.png", plot =  disSevLin, height = 12, width = 18, units = "in", dpi = 300)

# Vax Bin
combFilt = combShannon[!is.na(combShannon$Vaccinated) ,]

vaxBin = ggplot(combFilt, aes(x=Vaccinated, y=NormShan, fill = Vaccinated)) + 
  geom_violin(trim=FALSE)+
  geom_jitter() + 
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 18,        # Use the shape code for a diamond
    size = 8,          # Adjust the size of the diamond
    color = "white"  # Adjust dodge width if necessary
  ) + theme_classic() + theme(text = element_text(size = 26))

 ggplot(combFilt, aes(x = Vaccinated, y = Shannon)) +
   geom_boxplot() 

ggsave(file = "Paper/Figs/Shannon_VaxBin.png", plot = vaxBin, height = 10, width = 14, units = "in", dpi = 300)

# who variants
colnames(combShannon)
table(combShannon$WHO_variant)

combFilt = combShannon[combShannon$WHO_variant != "Gamma" & combShannon$WHO_variant != "Mu" ,]

whoVar = ggplot(combFilt, aes(x=WHO_variant, y=NormShan, fill = WHO_variant)) + 
  geom_violin(trim=FALSE)+
  geom_jitter() + 
  stat_summary(
    fun = "median",
    geom = "point",
    shape = 18,        # Use the shape code for a diamond
    size = 8,          # Adjust the size of the diamond
    color = "white"  # Adjust dodge width if necessary
  ) + theme_classic() + theme(text = element_text(size = 26))

ggplot(combFilt, aes(x = WHO_variant, y = NormShan)) +
  geom_boxplot() 

ggsave(file = "Paper/Figs/Shannon_WHO.png", plot = whoVar, height = 10, width = 14, units = "in", dpi = 300)