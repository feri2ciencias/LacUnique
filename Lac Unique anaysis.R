install.packages("readr")
install.packages("paleoMAS")
install.packages('paleoMAS')
install.packages("remotes")
remotes::install_version("paleoMAS", version = "2.0-1")
install.packages('MASS')
install.packages('lattice')



library(readr)
library(tidyverse)     
library(tidypaleo)     
library(vegan)       
library(ggfortify)
library(broom)        
library(cowplot)      
library(ggrepel)       
library(ggthemes)     
library(ggforce)       
library(MASS)
library(lattice)
library(vegan)
library(paleoMAS)

data_set <- read_csv("https://raw.githubusercontent.com/feri2ciencias/LacUnique/refs/heads/main/Full_dataset.csv")


# Environmental variables
env <- data_set[, which(names(data_set) == "MnFe"):which(names(data_set) == "pheophytin a")]

# Cladocera data set
clad <- data_set[, which(names(data_set) == "Bosmina sp."):which(names(data_set) == "Sida crystallina")]

# Quironomids
qui <- data_set[, which(names(data_set) == "Procladius"):which(names(data_set) == "Tanytarsussp.")]


#Figure 2

env.pca <- env %>% 
  decostand(., method = "standardize") %>%  
  scale() %>%                               
  prcomp()                                  

# Results PCA
env.pca

# Biplot  scores
env.pca %>%
  augment(env) %>%  # agrega el dataset original
  ggplot(aes(.fittedPC1, .fittedPC2)) + 
  geom_point(size = 1.5) +
  theme_half_open(12) + 
  background_grid()

# Vectors
env.pca %>%
  tidy(matrix = "rotation")

arrow_style <- arrow(
  angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt"))

# Biplot de vectores
env.pca %>%
  tidy(matrix = "rotation") %>%
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>% 
  ggplot(aes(PC1, PC2)) +
  geom_segment(xend = 0, yend = 0, arrow = arrow_style) +
  xlim(-1, 0.5) + ylim(-0.5, 0.7) +
  coord_fixed() +
  theme_minimal_grid(12) +
  geom_text_repel(aes(label = column))

# Scores
scores <- env.pca$x
scores

PCA1 <- scores[, 1]

combined <- data.frame(
  Depth = data_set$`Midpoint Depth (cm)`,
  env,
  PCA1 = PCA1
)

combined

data.piv <- combined %>%
  pivot_longer(
    cols = -Depth,
    names_to = "element",
    values_to = "abun"
  )

unique_age<-age_depth_model(depth=data_set$`Midpoint Depth (cm)`,
                            age=data_set$`Year (AD)`)



theme_set(theme_paleo())
alta_plot <- ggplot(data.piv, aes(x=abun, y =Depth)) +
  geom_lineh() +
  geom_point() +
  labs(x = NULL, y = "Depth (cm)") +
  scale_y_reverse() +
  facet_geochem_gridh(vars(element)) +
  labs(y = "Depth (cm)", x = NULL)+
  theme(strip.text.x = element_text(size = 8),
        axis.text = element_text(size = 6))

alta_plot


alta_plot_Wage<-alta_plot+ 
  scale_y_depth_age(
    unique_age,
    age_name = "Age (Year AD)",
    age_breaks = waiver(),
    age_labels = waiver())

alta_plot_Wage


##Figure 3

ISOT <- data_set %>% 
  ggplot(aes(x = D13C, y = CN)) +
  
  geom_point(aes(shape = factor(Zone)), 
             size = 4, color = "black") +
  
  geom_path(aes(group = 1),   # mantiene orden de filas (importante)
            color = "black",  
            linetype = "dashed", 
            linewidth = 0.8) +
  
  xlab(bquote(delta^13 * C)) +
  ylab("C/N") +
  
  scale_shape_manual(
    name = "Zones",
    labels = c("Z3=1980-2014", "Z2=1930-1980","Z1=1880-1930"),
    values = c(16, 17, 15)
  ) +
  
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.20),
    legend.background = element_rect(fill = "white", color = NA),
    legend.box.background = element_rect(color = NA)
  )

ISOT


#FIg. 4 
combined.clad <- data.frame(
  Depth = data_set$`Midpoint Depth (cm)`,
  clad)

data.piv.clad <- combined.clad %>%
  pivot_longer(
    cols = -Depth,
    names_to = "element",
    values_to = "taxon.clad")


shannon.ind.clad <- data.piv.clad %>%
  mutate(
    abun = as.numeric(taxon.clad),       # asegurar numérico
    abun = ifelse(is.na(abun), 0, abun)  # reemplazar NA por 0
  ) %>%
  group_by(Depth) %>%
  mutate(
    total = sum(abun),
    p = abun / total
  ) %>%
  summarise(
    Shannon = -sum(ifelse(p > 0, p * log(p), 0)),
    .groups = "drop"
  )


theme_set(theme_paleo())
alta_plot.sha.clad <- ggplot(shannon.ind.clad, aes(x=Shannon, y =Depth)) +
  scale_y_reverse() +
  geom_point()+
  labs(x = "Shannon index", y = "Depth (cm)")+
  theme(strip.text.x = element_text(size = 5),
        axis.text = element_text(size = 10))

alta_plot.sha.clad     

alta_plot_Wage<-alta_plot.sha.clad+ 
  scale_y_depth_age(
    unique_age,
    age_name = "Age (Year AD)")


data.piv.clad$element<-factor(data.piv.clad$element, levels=c("Bosmina.sp.", "D..longispina", "D..pulex",  "Chydorus.brevilabris",      
                                                              "Alona.quadrangularis", "Alona.affinis",  "Rynchotalona.falcata",      
                                                              "Alona.guttata",  "Alona.circumfibriata", "Graptoleberis.testudinaria",
                                                              "Pleuroxus.sp.",  "Acroperus.harpae", "Alonella.excisa", "Alona.intermedia",          
                                                              "Eurycercus.sp.", "Chydorus.gibbus", "Leydigia.leydigi",  "Alonella.nana",             
                                                              "Camptocercus.sp.", "Chydorus.faviformus", "Monospilus.dispar",         
                                                              "Sida.crystallina" ),          
                       ordered=TRUE)
         
 



theme_set(theme_paleo())
alta_plot_clad <- ggplot(data.piv.clad, aes(x=taxon.clad, y =Depth)) +
  geom_colh(width = 0.5, position = "dodgev")+ 
  scale_y_reverse() +
  labs(x = NULL, y = "Depth (cm)")+
  facet_abundanceh(vars(element),
                   labeller = label_species,
                   rotate_facet_labels= 45,
                   scale = "free_x",
                   dont_italicize = c("\\(.*?\\)", "spp?\\.", "-complex", "[Oo]ther"))+
  theme(strip.text.x = element_text(size = 9),
        axis.text = element_text(size = 9))

alta_plot_clad


#Fig. 5

poll2 <- qui %>%
  rowwise() %>%  # Para operar por fila
  mutate(
    total = sum(c_across(everything()), na.rm = TRUE)  # Suma solo columnas de especies
  ) %>%
  mutate(across(-total, ~ (.x / total) * 100)) %>%  # Abundancia relativa (%)
  ungroup() %>%
  select(-total) 


combined.qui <- data.frame(
  Depth = data_set$`Midpoint Depth (cm)`,
  poll2)

data.piv.qui <- combined.qui %>%
  pivot_longer(
    cols = -Depth,
    names_to = "element",
    values_to = "taxon.qui")


shannon.ind.qui <- data.piv.qui %>%
  mutate(
    abun = as.numeric(taxon.qui),       # asegurar numérico
    abun = ifelse(is.na(abun), 0, abun)  # reemplazar NA por 0
  ) %>%
  group_by(Depth) %>%
  mutate(
    total = sum(abun),
    p = abun / total
  ) %>%
  summarise(
    Shannon = -sum(ifelse(p > 0, p * log(p), 0)),
    .groups = "drop"
  )


theme_set(theme_paleo())
alta_plot.sha.qui <- ggplot(shannon.ind.qui, aes(x=Shannon, y =Depth)) +
  scale_y_reverse() +
  geom_point()+
  labs(x = "Shannon index", y = "Depth (cm)")+
  theme(strip.text.x = element_text(size = 5),
        axis.text = element_text(size = 10))

alta_plot.sha.qui   

alta_plot_Wage. qui<-alta_plot.sha.qui+ 
  scale_y_depth_age(
    unique_age,
    age_name = "Age (Year AD)")


data.piv.qui$element<-factor(data.piv.qui$element, levels=c(
  "Procladius",  "Microtendipes", "Polypedilum", "T.mendax",  "Einfeldia",                
  "Endochironomus", "Corynoneura.Thienemaniella", "Cryptochironomus",
  "Micropsectra",  "C.plumosus", "Cladopelma", "Psectrocladius",
  "Glyptotendipespallens", "C.anthracinus", "Tanytarsussp."),          
                              ordered=TRUE)



theme_set(theme_paleo())
alta_plot_qui <- ggplot(data.piv.qui, aes(x=taxon.qui, y =Depth)) +
  geom_colh(width = 0.5, position = "dodgev")+ 
  scale_y_reverse() +
  labs(x = NULL, y = "Depth (cm)")+
  facet_abundanceh(vars(element),
                   labeller = label_species,
                   rotate_facet_labels= 45,
                   scale = "free_x",
                   dont_italicize = c("\\(.*?\\)", "spp?\\.", "-complex", "[Oo]ther"))+
  theme(strip.text.x = element_text(size = 9),
        axis.text = element_text(size = 9))

alta_plot_qui


#Fig. 6 RDA Cladocera
comm.cla <- decostand(combined.clad[,-1], method = "total")
comm.cla*100->porc.cla
sqrt(porc.cla)->raiz.cla

env.rda <- env %>% 
  decostand(., method = "standardize")

rda.uni.cla<-rda(downweight(raiz.cla)~., data=env.rda)

step.forward<- ordistep(rda(downweight(raiz.cla)~1, data=env.rda),
                        scope=formula(rda.uni.cla), 
                        direction = "forward", pstep=1000)

rda.uni.cla2 <- rda(downweight(raiz.cla) ~ D13C  + `pheophytin a` + CN+ MnFe, data = env.rda)

plot(rda.uni.cla, scaling=3)

vif.cca(rda.uni.cla2)

anova.cca(rda.uni.cla2, by="terms")

colores <- ifelse(data_set$Zone == "1", "blue",
                  ifelse(data_set$Zone == "2", "darkgoldenrod2",
                         ifelse(data_set$Zone == "3", "darkgreen", "gray")))
formas <- ifelse(data_set$Zone == "1", 19,
                 ifelse(data_set$Zone == "2", 17,
                        ifelse(data_set$Zone == "3", 15, NA)))


par(mfrow=c(1,2))
plot(rda.uni.cla2, type="n", xlim=c(-2,2), ylim=c(-1,1),
     xlab="RDA 1 (87.07%)", ylab="RDA 2 (8.30%)" )
text(rda.uni.cla2, display = "species",  scaling =3, 
     cex = 0.8, col = "black",pos=1)

labl <- c("13C", "D15N", "CN", "Ti")   
text(rda.uni.cla2, display="bp", col="Blue", scaling=3, labels=labl, pos=3)


plot(rda.uni.cla2, type="n", xlim=c(-2,2), ylim=c(-1,1),
     xlab="RDA 1 (87.07%)", ylab="RDA 2 (8.30%)" )
points(rda.uni.cla2, display="sites", scaling=1,
       col=colores, pch=as.numeric(formas), cex=2.5)

labl <- c("13C", "D15N", "CN", "Ti")  
text(rda.polo.qui2, display="bp", col="Blue", scaling=3, labels=labl, pos=3)


#Fig. 6 RDA Chironomids

comm.qui <- decostand(combined.qui[,-1], method = "total")
comm.qui*100->porc.qui
sqrt(porc.qui)->raiz.qui

env.rda <- env %>% 
  decostand(., method = "standardize")

rda.uni.qui<-rda(downweight(raiz.qui)~., data=env.rda)

step.forward.qui<- ordistep(rda(downweight(raiz.qui)~1, data=env.rda),
                        scope=formula(rda.uni.qui), 
                        direction = "forward", pstep=1000)

#Con esto hacemos el forward para escoger las variables

rda.uni.qui2 <- rda(downweight(raiz.qui) ~ `pheophytin a` + MnFe+ `K/Ti` +CN  +  D13C,  data = env.rda)

plot(rda.uni.qui2, scaling=3)
X11()

vif.cca(rda.uni.qui2)

anova.cca(rda.uni.qui2, by="terms")

plot(rda.uni.qui2, type="n", xlim=c(-2,2), ylim=c(-1,1),
     xlab="RDA 1 (59.61%)", ylab="RDA 2 (19.68%)" )
text(rda.uni.qui2, display = "species",  scaling =3, 
     cex = 0.8, col = "black",pos=1, font= 3)

labl <- c("pheophytin.a", "MnFe", "K.Ti", "CN",  "D13C" )   
text(rda.uni.qui2, display="bp", col="Blue", scaling=3, labels=labl, pos=3)


plot(rda.uni.qui2, xlim=c(-2,2), ylim=c(-1.5,1.5), type="n",
     xlab="RDA 1 (59.61%)", ylab="RDA 2 (19.68%)",pch=as.numeric(formas), col=colores)
points(rda.uni.qui2, display="sites", scaling=1,
       col=colores, pch=as.numeric(formas), cex=2.5)


#Fig.7 Ecological change

##Summer temperature

summer.temp <- read_csv("https://raw.githubusercontent.com/feri2ciencias/LacUnique/refs/heads/main/Weather.csv")



decorana(clad)->clados.dca
decorana(qui)->quiros.dca
qui[,-8]->qui2
data_set$`Year (AD)`->años

scores(clados.dca,display="species")->scores.spp.cla
scores(quiros.dca,display="species")->scores.spp.qui

scores(clados.dca,display="species")->scores.spp.cla
scores(quiros.dca,display="species")->scores.spp.qui

change(clad, años, dca=T, roc=T )

par(mfrow=c(1,2))
analog.sing(clad, base=1, años, dca=T, method= "morisita")->cla.dist
analog.sing(qui2, base=1, años, dca=T, method= "morisita")->qui.dist

as.data.frame(cla.dist)->cla.dist.df
as.data.frame(qui.dist)->qui.dist.df

plot(cla.dist.df$Distance, cla.dist.df$Age, col="blue", type= "o", lwd=2,
     xlim=c(0,4), ylim=c(1884,2014), xlab= ("Ecological distance to present (SD)"),
     ylab=("Age (years CE)"))
lines(qui.dist.df$Distance, cla.dist.df$Age, col="red", type= "o", lwd=2)
abline(v=2, lty=2, col="slategrey")
abline(h=1980, lty=3, col="black")
abline(h=1930, lty=3, col="black")

#Ecological change
change(clad,años,dca=T)->cla.ecoch
change(qui, años, dca=T)->qui.ecoch

as.data.frame(cla.ecoch)->cla.ecoch.df
as.data.frame(qui.ecoch)->qui.ecoch.df

par(mfrow=c(1,2))

plot(qui.ecoch.df$Distance, qui.ecoch.df$Age, col="gray58", type= "o", lwd=2,
     xlim=c(0,1.5), ylim=c(1898,2014),
     ylab=("Age (years CE)"),  axes=FALSE, ann=FALSE)
box()
axis(side=2, las=1, at=z1)
axis(1)
lines(qui.ecoch.df$Distance, cla.ecoch.df$Age, col="gray22", type= "o", lwd=2)
mtext(text="Ecological change (SD)", side=1, line=2.3, cex=0.6)
abline(h=1980, lty=3, col="black")
abline(h=1930, lty=3, col="black")

##Combinado con temp

par(mfrow=c(1,4))
par(oma=c(2,2,2,2))


#Ecological change
par(mar=c(4,4,4,0))
plot(cla.ecoch.df$Distance, cla.ecoch.df$Age, col="gray58", type= "o", lwd=2,
     xlim=c(0,1.5), ylim=c(1884,2014),
     ylab=("Age (years CE)"),  axes=FALSE, ann=FALSE)
axis(side=2, las=1, at=cla.ecoch.df$Age)  #las= gira las etiquetas
box()
axis(1)
lines(qui.ecoch.df$Distance, cla.ecoch.df$Age, col="gray22", type= "o", lwd=2)
mtext(text="Ecological change (SD)", side=1, line=2.3, cex=0.6)
abline(h=1980, lty=3, col="black")
abline(h=1930, lty=3, col="black")

#Cladoceros

par(mar=c(4,0,4,0))
plot(cla.dist.df$Distance, cla.dist.df$Age, col="gray58", type= "o", lwd=2,
     xlim=c(0,4), ylim=c(1884,2014), xlab= ("Ecological distance to present (SD)"),
     axes=FALSE, ann=FALSE)
axis(1)
box()
aux=expression("Cladocerans")
mtext(aux,3,3.3)
mtext(text="Ecological distance to present (SD)", side=1, line=2.3, cex=0.6)
abline(v=1, lty=2, col="slategrey")
abline(h=1980, lty=3, col="black")
abline(h=1930, lty=3, col="black")

# quironomids

par(mar=c(4,0,4,0))
plot(qui.dist.df$Distance, qui.dist.df$Age, col="gray22", type= "o", lwd=2,
     xlim=c(0,4), ylim=c(1884,2014), xlab= ("Ecological distance to present (SD)"),
     ylab=("Age (years CE)"),axes=FALSE, ann=FALSE)
axis(1)
box()
aux=expression("Quironomids")
mtext(aux,3,3.3)
mtext(text="Ecological distance to present (SD)", side=1, line=2.3, cex=0.6)
abline(v=1, lty=2, col="slategrey")
abline(h=1980, lty=3, col="black")
abline(h=1930, lty=3, col="black")

#Summer
cbind(summer.temp$Summer,summer.temp$Year)->Summer
par(mar=c(4,0,4,4))
plot(Summer, type="o", lwd=2,
     xlim=c(14.73,17.43), ylim=c(1880,2014),
     axes=FALSE, ann=FALSE)
axis(side=4, las=1, at=cla.ecoch.df$Age)  #las= gira las etiquetas
axis(1)
box()
aux=expression("Summer")
mtext(aux,3,3.3)
mtext(text="°C", side=3, line=2.3)
abline(h=1980, lty=3, col="black")
abline(h=1930, lty=3, col="black")


