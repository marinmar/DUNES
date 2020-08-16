library(readxl)
library(gdata)
library(vegclust)
library(BiodiversityR)
library(officer)
library(flextable)
library(htmltools)
library(webshot)
library(ggpmisc)
library(multcomp)
library(ggplot2)
library(tibble)
library(RColorBrewer)
library(ggpubr)
library(vegan)
library(reshape2)
library(corrplot)
library(car)
library(cowplot)
library("FactoMineR")
library("factoextra")
library(tidyverse)
library(Rmisc)
library(dplyr)

set.seed(1L)

std_border = fp_border(color="gray")

formula <- y ~ x
specify_decimal <-
  function(x, k)
    trimws(format(round(x, k), nsmall = k))

mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(24)

my_comparisons <-
  list(c("Megaride", "Dunes"),
       c("Banc", "Dunes"),
       c("Banc", "Megaride"),
       c("Dunes", "Reference"),
       c("Megaride", "Reference"),
       c("Banc", "Reference"))


#### Dunes.env ####

Dunes.env <-
  read_excel(
    "G:/cloud/stage 2a/figures et tables/tables abio bio/explanatory_variables_Dunes_compiled.xlsx",
    na = "NA"
  ) %>% as_tibble()

Dunes.env <- select(Dunes.env, !Commentaire)
Dunes.env <-
  Dunes.env %>%  mutate_at("Box_name", str_replace_all, "Binnen_Ratel", "Banc")
Dunes.env <-
  Dunes.env %>%  mutate_at("Box_name", str_replace_all, "Nord_Breedt", "Megaride")
Dunes.env <-
  Dunes.env %>%  mutate_at("Box_name", str_replace_all, "Chenal", "Dunes")
Dunes.env <-
  Dunes.env %>%  mutate_at("Season", str_replace_all, "Autumn", "Automne")
Dunes.env <-
  Dunes.env %>%  mutate_at("Season", str_replace_all, "Spring", "Printemps")
Dunes.env$Box_name <-
  factor(Dunes.env$Box_name,
         levels = c("Banc", 'Megaride', "Dunes", "Reference"))
Dunes.env$Position <-
  factor(Dunes.env$Position,
         levels = c("Creux", "Pente", "Sommet", "Pente raide", "Reference"))
Dunes.env$Transect <- factor(Dunes.env$Transect, levels = 1:7)
Dunes.env$Position <-
  factor(Dunes.env$Position,
         levels = c("Creux", 'Pente', "Sommet", "Reference"))


Treated_Stations = c("S13",
                     "S14",
                     "S15",
                     "S9",
                     "S8",
                     "S10",
                     "S21",
                     "S23",
                     "S1",
                     "S2",
                     "S4")
Dunes.env <- Dunes.env %>% filter(Station %in% Treated_Stations)
Dunes.env$Station <-
  factor(Dunes.env$Station,
         levels = c("S15",
                    "S14",
                    "S13",
                    "S8",
                    "S9",
                    "S10",
                    "S2",
                    "S4",
                    "S1",
                    "S21",
                    "S23"))
Dunes.env.comp <- Dunes.env
Dunes.env <- na.omit(Dunes.env)

#### Dunes.bio ####

Dunes.bio <-
  Matrice_station_espece_Mai_2020_Matched <-
  read_excel(
    "G:/cloud/stage 2a/figures et tables/tables abio bio/Matrice_station_espece_compiled.xlsx"
  ) %>% as_tibble()
Dunes.bio <- Dunes.bio %>% filter(Station %in% Treated_Stations)

Species_name <-
  unique(pull(Dunes.bio, "ScientificName_accepted")) #liste de toutes les espèces présentes
StationxSpecies <-
  distinct(Dunes.bio[, 1]) #liste avec stations uniques

for (i in 1:length(Species_name)) {
  StationxSpecies = add_column(StationxSpecies, !!Species_name[i] := 0)
}
StationxSpecies_abundance <- StationxSpecies

for (i in 1:nrow(StationxSpecies)) {
  M = filter(Dunes.bio, Code == StationxSpecies[[i, "Code"]])
  for (j in 1:length(Species_name)) {
    for (k in 1:nrow(M)) {
      if (Species_name[j] == M[k, "ScientificName_accepted"]) {
        StationxSpecies_abundance[i, Species_name[j]] = M[k, "Number"]
      }
    }
  }
}

#Calcul des indices de diversité (Shannon, Richesse spécifique, indice d'équitabilité de Piélou)

Shannon <-
  diversity(subset(StationxSpecies_abundance, select = -c(Code)))
Rich_spe <-
  specnumber(subset(StationxSpecies_abundance, select = -c(Code)))
Pielou <-
  Shannon / log(specnumber(subset(
    StationxSpecies_abundance, select = -c(Code)
  )))

StationxSpecies_abundance <-
  mutate(
    StationxSpecies_abundance,
    Total_abundance = select_if(StationxSpecies_abundance, is.numeric) %>% rowSums(na.rm = TRUE)
  )
StationxSpecies_abundance = add_column(
  StationxSpecies_abundance,
  Shannon = !!Shannon,
  Rich_spe = !!Rich_spe,
  Pielou = !!Pielou
)
StationxSpecies_abundance <-
  left_join(StationxSpecies_abundance, distinct(subset(
    Dunes.bio, select = c(Code, Station, Code_prelevement)
  )))
StationxSpecies_abundance <-
  left_join(
    StationxSpecies_abundance,
    subset(
      Dunes.env,
      select = c(
        index,
        Box_name,
        Season,
        Depth,
        Longitude,
        Latitude,
        Position,
        Median,
        Sorting,
        Gravel,
        Sand,
        Silt,
        Clay,
        Pebble,
        Granules,
        Very.Coarse.Sand,
        Coarse.Sand,
        Medium.Sand,
        Fine.Sand,
        Very.Fine.Sand,
        Temperature,
        Salinity,
        OxygenPercent,
        Temperature.at.bottom.Mean,
        Salinity.Mean,
        Current.speed.Mean,
        Chlorophyll.a.Mean,
        MeanTrawlingSurface,
        MeanTrawlingSubsurfSAR,
        Organic_matter_content
      )
    ),
    by = c("Code_prelevement" = "index")
  )

StationxSpecies_abundance <- na.omit(StationxSpecies_abundance)

#### CTD ####

ctd <- read_excel("ctd.xlsx") %>% as_tibble()
ctd_t <- filter(ctd,ctd$"Salinité (g/l)" >= 1)
ctd_t <- aggregate(ctd_t, by = list(ctd_t$Box_name), FUN = sd)
colMeans(select_if(ctd_t,is.numeric))

ggplot(ctd) +
 aes(x = `Température (°C)`, y = `Profondeur (m)`, colour = Box_name) +
 geom_path(size = 1L) +
 scale_color_hue() +
 scale_y_continuous(trans = "reverse") +
 theme_minimal() +
 facet_wrap(vars(Box_name), scales = "free_x")

ggplot(ctd) +
 aes(x = `Salinité (g/l)`, y = `Profondeur (m)`, colour = Box_name) +
 geom_path(size = 1L) +
 scale_color_hue() +
 scale_y_continuous(trans = "reverse") +
 theme_minimal() +
 facet_wrap(vars(Box_name), scales = "free_x")

ggplot(ctd) +
 aes(x = `Oxygène (%)`, y = `Profondeur (m)`, colour = Box_name) +
 geom_path(size = 1L) +
 scale_color_hue() +
 scale_y_continuous(trans = "reverse") +
 theme_minimal() +
 facet_wrap(vars(Box_name), scales = "free_x")
ctd

#### Courbes granulométriques ####



granulometrie <- data.frame(index = "S1",
                            fraction = "X1",
                            proportion = 0)

G = cbind(select(Dunes.env, index), select(Dunes.env, starts_with("X")))
G <- na.omit(G)
granulometrie <- melt(G, value.name = "value")

granulometrie <-
  left_join(granulometrie, subset(
    Dunes.env,
    select = c(index, Station, Season, Transect, Box_name, Position, Depth)
  ))

granulometrie$fraction <-
  factor(
    granulometrie$variable,
    levels = c(
      "X0",
      "X63",
      "X80",
      "X100",
      "X125",
      "X160",
      "X200",
      "X250",
      "X315",
      "X400",
      "X500",
      "X630",
      "X800",
      "X1000",
      "X1250",
      "X1600",
      "X2000",
      "X2500",
      "X3150",
      "X4000",
      "X5000",
      "X6300",
      "X8000",
      "X10000",
      "X12500"
    )
  )
granulometrie$Station <-
  factor(
    granulometrie$Station,
    levels = c(
      "S1",
      "S2",
      "S3",
      "S4",
      "S5",
      "S6",
      "S7",
      "S8",
      "S9",
      "S10",
      "S11",
      "S12",
      "S13",
      "S14",
      "S15",
      "S16",
      "S17",
      "S18",
      "S19",
      "S20",
      "S21",
      "S22",
      "S23"
    )
  )
granulometrie$Box_name <-
  factor(granulometrie$Box_name,
         levels = c("Banc", 'Megaride', "Dunes", "Reference"))
granulometrie$Transect <-
  factor(granulometrie$Transect, levels = 1:7)


image = ggplot(granulometrie, aes(x = variable, y = value, group = index, lines)) +
  geom_line(aes(color = Transect, linetype = Position), size = 0.7) +
  facet_grid(rows = vars(Box_name), cols = vars(Season)) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  scale_fill_brewer(palette = "Dark2") +
  panel_border() +
  labs(x = "Fraction granulométrique", y = "Proportion")
image
#ggsave(file="granulo_unified.png", plot=image, width=8, height=6)

#### Caract. points ####

image = ggplot(Dunes.env.comp, aes(x = Station, y = Depth, fill = Position)) +
  geom_col() +
  facet_grid(
    cols = vars(Box_name),
    rows = vars(Season),
    scales = "free",
    space = "free_x") +
  theme(
    axis.text.x =  element_text(angle = 90, hjust = 0),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")) +
  scale_fill_brewer(palette = "Dark2",name = "Position") +
  labs(x = "Station", y = "Profondeur (m)")
image
#ggsave(file="StationDepth.png", plot=image, width=8, height=6)

#### Caract. sedimentaires ####

# modèle linéaire
lm.Median <- lm(data = Dunes.env, Median ~ Box_name + Season + Depth)
summary(lm.Median)
lm.MO <-
  lm(data = Dunes.env, Organic_matter_content ~ Box_name + Season + Depth)
summary(lm.MO)
lm.Clay <- lm(data = Dunes.env, Clay ~ Box_name + Season + Depth)
summary(lm.Clay)
lm.Sorting <-
  lm(data = Dunes.env, Sorting ~ Box_name + Season + Depth)
summary(lm.Sorting)

#ACP

res.pca <-
  PCA(
    select(
      Dunes.env,
      Median,
      Sorting,
      Clays,
      Organic_matter_content,
    ),
    scale.unit = TRUE,
    graph = TRUE
  )
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
res.hcpc <- HCPC(res.pca, graph = FALSE)
Plot.Var <- fviz_pca_var(
  res.pca,
  col.var = "cos2",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
)
Plot.Season <-
  fviz_pca_ind(
    res.pca,
    col.ind = Dunes.env$Season,
    geom.ind = "point",
    legend.title = "Saison",
    palette = brewer.pal(8, "Set2")
  )

Plot.Box_name <-
  fviz_pca_ind(
    res.pca,
    col.ind = Dunes.env$Box_name,
    geom.ind = "point",
    legend.title = "Boite d'échantillonage",
    palette = brewer.pal(8, "Set2")
  )

image = ggdraw() +
  draw_plot(Plot.Var, 0, .5, .5, .5) +
  draw_plot(Plot.Box_name, .5, .5, .5, .5) +
  draw_plot(Plot.Season, 0, 0, .5, .5) +
  draw_plot_label(c("A", "C", "B"), c(0, 0, 0.5), c(1, 0.5, 1), size = 15)
image
#ggsave(file="ACP.env.png", plot=image, width=10, height=8)



#Permanova
Dune.env.sub <-
  select(Dunes.env, Median, Organic_matter_content, Clay, Sorting)
adonis(data = Dunes.env, Dune.env.sub ~ Box_name + Season + Depth)

#Tests de significativité divers

compare_means(Median ~ Box_name, data = Dunes.env)
compare_means(Organic_matter_content ~ Season,method = "t.test", data = Dunes.env,group.by = "Box_name")
compare_means(Clay ~ Box_name, data = Dunes.env)

Dunes.env.sub <- select(Dunes.env, Season,Box_name,Clay,Organic_matter_content)
agg <- aggregate(select_if(Dunes.env.sub,is.numeric), list(Dunes.env.sub$Season,Dunes.env.sub$Box_name), sd)
agg

#### Median ####

Median_1 = ggplot(Dunes.env, aes(x = Box_name, y = Median, fill = Box_name)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Boite d'échantillonnage", y = "Taille médiane des grains (µm)") +
  stat_compare_means(comparisons = my_comparisons  , label =  "p.signif")
Median_1

Median_2 = ggplot(Dunes.env, aes(x = Season, y = Median, fill = Box_name)) +
  geom_boxplot() +
  facet_grid(cols = vars(Box_name),
             scales = "free",
             space = "free_x") +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Boite d'échantillonnage", y = "Taille médiane des grains (µm)") +
  stat_compare_means(aes(group = Season) , label =  "p.signif")
Median_2

Median_3 = ggplot(Dunes.env, aes(x = Depth, y = Median, group = Season)) +
  geom_point(aes(color = Season)) +
  facet_grid(cols = vars(Box_name), scales = "free") +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  stat_smooth(method = "lm", se = FALSE, aes(color = Season)) +
  labs(x = "Profondeur (m)", y = "Taille médiane des grains (µm)") +
  scale_color_brewer(palette = "Dark2", name = "Saison")
Median_3
#ggsave(file="medianXDepth.png", plot=Median_3, width=6, height=3)

#### Clay ####

Clay_1 = ggplot(Dunes.env, aes(x = Box_name, y = Clay, fill = Box_name)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Boite d'échantillonnage", y = "Taux d'argiles (% du poids total)") +
  stat_compare_means(comparisons = my_comparisons  , label =  "p.signif")
Clay_1


Clay_2 = ggplot(Dunes.env, aes(x = Season, y = Clay, fill = Box_name)) +
  geom_boxplot() +
  facet_grid(cols = vars(Box_name),
             scales = "free",
             space = "free_x") +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Saison") +
  labs(x = "Boite d'échantillonnage", y = "Taux d'argiles (% du poids total)") +
  stat_compare_means(aes(group = Season) , label =  "p.signif", method = "t.test")
Clay_2

Clay_3 = ggplot(Dunes.env, aes(x = Depth, y = Clay, group = Season)) +
  geom_point(aes(color = Season)) +
  facet_grid(cols = vars(Box_name)) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  labs(x = "Boite d'échantillonnage", y = "Taux d'argiles (% du poids total)") +
  stat_smooth(method = "lm", se = FALSE, aes(color = Season)) +
  scale_color_brewer(palette = "Dark2", name = "Saison")
Clay_3

#### MO ####

MO_1 = ggplot(Dunes.env,
              aes(x = Box_name, y = Organic_matter_content, fill = Box_name)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Boite d'échantillonnage", y = "Taux de matière organique (% du poids total)") +
  stat_compare_means(comparisons = my_comparisons  , label =  "p.signif")
MO_1


MO_2 = ggplot(Dunes.env,
              aes(x = Season, y = Organic_matter_content, fill = Box_name)) +
  geom_boxplot() +
  facet_grid(cols = vars(Box_name),
             scales = "free",
             space = "free_x") +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Saison") +
  labs(x = "Boite d'échantillonnage", y = "Taux de matière organique (% du poids total)") +
  stat_compare_means(aes(group = Season) , label =  "p.signif", method = "t.test")
MO_2

MO_3 = ggplot(Dunes.env,
              aes(x = Depth, y = Organic_matter_content, group = Season)) +
  geom_point(aes(color = Season)) +
  facet_grid(cols = vars(Box_name)) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  labs(x = "Boite d'échantillonnage", y = "Taux de matière organique (% du poids total)") +
  stat_smooth(method = "lm", se = FALSE, aes(color = Season)) +
  scale_color_brewer(palette = "Dark2", name = "Saison")
MO_3

#### Sorting ####

Sorting_1 = ggplot(Dunes.env, aes(x = Box_name, y = Sorting, fill = Box_name)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Boite d'échantillonnage", y = "Indice de tri") +
  stat_compare_means(comparisons = my_comparisons  , label =  "p.signif")
Sorting_1

Sorting_2 = ggplot(Dunes.env, aes(x = Season, y = Sorting, fill = Box_name)) +
  geom_boxplot() +
  facet_grid(cols = vars(Box_name),
             scales = "free",
             space = "free_x") +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Saison") +
  labs(x = "Boite d'échantillonnage", y = "Indice de tri") +
  stat_compare_means(aes(group = Season) , label =  "p.signif", method = "t.test")
Sorting_2

Sorting_3 = ggplot(Dunes.env, aes(x = Depth, y = Sorting, group = Season)) +
  geom_point(aes(color = Season)) +
  facet_grid(cols = vars(Box_name)) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  labs(x = "Boite d'échantillonnage", y = "Indice de tri") +
  stat_smooth(method = "lm", se = FALSE, aes(color = Season)) +
  scale_color_brewer(palette = "Dark2", name = "Saison")
Sorting_3

image = ggdraw() +
  draw_plot(Median_1, 0, .5, .5, .5) +
  draw_plot(Clay_1, .5, .5, .5, .5) +
  draw_plot(MO_1, 0, 0, .5, .5) +
  draw_plot(Sorting_1, 0.5, 0, .5, .5) +
  draw_plot_label(c("A", "C", "B", "D"), c(0, 0, 0.5, 0.5), c(1, 0.5, 1, 0.5), size = 15)
image
#ggsave(file="IndiceXBoites.png", plot=image, width=6, height=6)

image = ggdraw() +
  draw_plot(Median_3, 0, 0, 1, 0.5) +
  draw_plot(Clay_2, 0, 0.5, .5, .5) +
  draw_plot(MO_2, 0.5, 0.5, .5, .5) +
  draw_plot_label(c("A", "C", "B"), c(0, 0, 0.5), c(1, 0.5, 1), size = 15)
image
#ggsave(file="IndiceXBoites2.png", plot=image, width=6, height=6)

#### Tableau résumé sediments ####


tab <- tibble(name = "",
              unit ="",
              Moyenne = "",
              Ecart_type = "")

tab <-
  tab %>% add_row(
    name = "Taille médiane des grains",
    unit = "µm",
    Moyenne = mean(pull(Dunes.env, Median)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, Median)) %>% specify_decimal(3)
  )
tab <-
  tab %>% add_row(
    name = "Taux d'argiles",
    unit = "%",
    Moyenne = mean(pull(Dunes.env, Clay)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, Clay)) %>% specify_decimal(3)
  )
tab <-
  tab %>% add_row(
    name = "Taux de matière organique",
    unit = "%",
    Moyenne = mean(pull(Dunes.env, Organic_matter_content)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, Organic_matter_content)) %>% specify_decimal(3)
  )
tab <-
  tab %>% add_row(
    name = "Indice de tri",
    unit= "",
    Moyenne = mean(pull(Dunes.env, Sorting)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, Sorting)) %>% specify_decimal(3)
  )
tab <- tab[-c(1), ]
tab

tab <- flextable(data = tab, col_keys = c("name","col1","unit","col2","Moyenne","Ecart_type"))
tab <-  autofit(tab)
tab <- set_header_labels(tab, name = "",unit="Unité",Ecart_type="Ecart-type")
tab <- hline(tab, part="body", border = std_border )
tab <- empty_blanks(tab)
tab
#print(tab, preview = "docx")

tab <- tibble(origine="",
              name = "",
              unit = "",
              Moyenne = "",
              Ecart_type = "")

tab <-
  tab %>% add_row(
    origine = "Sonde CTD",
    name = "Température",
    unit = "°C",
    Moyenne = mean(pull(Dunes.env, Temperature)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, Temperature)) %>% specify_decimal(3)
  )
tab <-
  tab %>% add_row(
    origine = "Sonde CTD",
    name = "Salinité",
    unit = "g.l^-1",
    Moyenne = mean(pull(Dunes.env, Salinity)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, Salinity)) %>% specify_decimal(3)
  )
tab <-
  tab %>% add_row(
    origine = "Sonde CTD",
    name = "Taux d'oxygène",
    unit = "%",
    Moyenne = mean(pull(Dunes.env, OxygenPercent)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, OxygenPercent)) %>% specify_decimal(3)
  )
tab <-
  tab %>% add_row(
    origine = "Copernicus",
    name = "Vitesse des courants",
    unit = "m.s^-1",
    Moyenne = mean(pull(Dunes.env, Current.speed.Mean)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, Current.speed.Mean)) %>% specify_decimal(3)
  )
tab <-
  tab %>% add_row(
    origine = "Copernicus",
    name = "Température au fond",
    unit = "°C",
    Moyenne = mean(pull(Dunes.env, Temperature.at.bottom.Mean)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, Temperature.at.bottom.Mean)) %>% specify_decimal(3)
  )
tab <-
  tab %>% add_row(
    origine = "Copernicus",
    name = "Salinité",
    unit = "mg.l^-1",
    Moyenne = mean(pull(Dunes.env, Salinity.Mean)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, Salinity.Mean)) %>% specify_decimal(3)
  )
tab <-
  tab %>% add_row(
    origine = "Copernicus",
    name = "Chlorophylle a",
    unit = "mg.m^3",
    Moyenne = mean(pull(Dunes.env, Chlorophyll.a.Mean)) %>% specify_decimal(3),
    Ecart_type = sd(pull(Dunes.env, Chlorophyll.a.Mean)) %>% specify_decimal(3)
  )
 
tab <- tab[-c(1), ] 
tab <- flextable(data = tab, col_keys = c("origine","name","col2","unit","col1","Moyenne","Ecart_type"))
tab <-  autofit(tab)
tab <- set_header_labels(tab, name = "",unit="Unité",origine="",Ecart_type="Ecart-type")
tab <- merge_v(tab, j = ~ origine )
tab <- hline(tab, part="body", border = std_border )
tab <- empty_blanks(tab)
tab
#print(tab, preview = "docx")

#indices biodiv

tab <- tibble(name = "",
              unit = "",
              Moyenne = "",
              Ecart_type = "")

tab <-
  tab %>% add_row(
    name = "Richesse spécifique",
    unit = "espèces/0.1m²",
    Moyenne = mean(pull(StationxSpecies_abundance, Rich_spe)) %>% specify_decimal(1),
    Ecart_type = sd(pull(StationxSpecies_abundance, Rich_spe)) %>% specify_decimal(1)
  )
tab <-
  tab %>% add_row(
    name = "Abondance totale",
    unit = "ind/0.1m²",
    Moyenne = mean(pull(StationxSpecies_abundance, Total_abundance)) %>% specify_decimal(1),
    Ecart_type = sd(pull(StationxSpecies_abundance, Total_abundance)) %>% specify_decimal(1)
  )
tab <-
  tab %>% add_row(
    name = "Indice de Shannon",
    Moyenne = mean(pull(StationxSpecies_abundance, Shannon)) %>% specify_decimal(2),
    Ecart_type = sd(pull(StationxSpecies_abundance, Shannon)) %>% specify_decimal(2)
  )
tab <-
  tab %>% add_row(
    name = "Indice d'équitabilité de Piélou",
    Moyenne = mean(pull(StationxSpecies_abundance, Pielou)) %>% specify_decimal(2),
    Ecart_type = sd(pull(StationxSpecies_abundance, Pielou)) %>% specify_decimal(2)
  )

tab <- tab[-c(1), ] 
tab <- flextable(data = tab, col_keys = c("name","col1","unit","col2","Moyenne","Ecart_type"))
tab <-  autofit(tab)
tab <- set_header_labels(tab, name = "",Ecart_type="Ecart-type",unit="Unité")
tab <- hline(tab, part="body", border = std_border )
tab <- empty_blanks(tab)
tab
#print(tab, preview = "docx")




#### Analyses biodiv ####

MSE <- StationxSpecies_abundance #Matrice Station Espèces
MSE

#On retire l'outlier S2R3A1
#MSE <- MSE %>% filter(Code != "S2R3A1")

MSE.sub <-
  MSE %>% select(c(Total_abundance, Shannon, Rich_spe, Pielou))

#Permanova

adonis(data = MSE, MSE.sub ~ Box_name + Season + Depth)

# analyses moyennes diverses

compare_means(Pielou ~ Box_name, data = MSE, group.by = "Season")
compare_means(Clay ~ Box_name, data = Dunes.env)
MSE.sub <- select(MSE, Season,Box_name,Rich_spe,Total_abundance,Shannon,Pielou)
agg <- aggregate(select_if(MSE.sub,is.numeric), list(MSE.sub$Season,MSE.sub$Box_name), sd)
agg

# diagrammes rang fréquence
MSE.sub <- MSE %>% select(1, all_of(Species_name))
MSE.sub <- MSE.sub %>% column_to_rownames(var = "Code")
ra <- rankabundance(MSE.sub) 

ra <- ra  %>% data.frame() %>% rownames_to_column(var = "Species")
ra <- select(head(ra,10),Species,proportion)
ra

tab <- flextable(data = ra, col_keys = c("Species","col1","proportion"))
tab <-  autofit(tab)
tab <- set_header_labels(tab, Species = "Espèce",proportion="Fréquence")
tab <- hline(tab, part="body", border = std_border )
tab <- empty_blanks(tab)
tab
#print(tab, preview = "docx")

# analyse univariées indices

lm.Median <- lm(data = MSE, Rich_spe ~ Median + Organic_matter_content )
summary(lm.Median)

#### Nmds all ####

MSE.sub <- MSE %>% select(1, all_of(Species_name))
MSE.sub <- MSE.sub %>% column_to_rownames(var = "Code")
MSE.sub <- MSE.sub %>% log1p()

Nmds = metaMDS(MSE.sub, trymax = 100)

spp.scrs <- as.data.frame(na.omit(Nmds$species))

data.scores = as.data.frame(scores(Nmds))
data.scores$Box_name = MSE$Box_name
data.scores$Season = MSE$Season
data.scores$Position = MSE$Position
data.scores$Depth = MSE$Depth

simp <- simper(MSE.sub,MSE$Season,permutations=999)
simp

spp.scrs <- as.data.frame(na.omit(Nmds$species))

Automne <- tibble("Species",.name_repair = ~ c("Species"))
Printemps <- tibble("Species",.name_repair = ~ c("Species"))


simp_list <- list("Automne" = Automne,"Printemps" = Printemps)

for (i in 1:length(simp)){
  a <- str_split(names(simp)[i], "_")[[1]][1]
  b <- str_split(names(simp)[i], "_")[[1]][2]
  for (j in 1:length(simp_list)){
    if (a == names(simp_list)[j]){
      l <- cbind(rownames_to_column(data.frame(simp[[i]]$ava),var = "Species"),data.frame(simp[[i]]$p))
      names(l)=c("Species","sim","p")
      l <- filter(l,p<0.05)
      l <- filter(l,sim!=0)
      simp_list[[j]] <- full_join(simp_list[[j]],select(l,Species,sim),by = "Species")
    }
    if (b == names(simp_list)[j]){
      l <- cbind(rownames_to_column(data.frame(simp[[i]]$avb),var = "Species"),data.frame(simp[[i]]$p))
      names(l)=c("Species","sim","p")
      l <- filter(l,p<0.05)
      l <- filter(l,sim!=0)
      simp_list[[j]] <- full_join(simp_list[[j]],select(l,Species,sim),by = "Species")
    }
  }
}
spc_list <- tibble(Species="",Value=0,Code="",.name_repair = ~ c("Species","Value","Code"))
for (i in 1:length(simp_list)){
  simp_list[[i]]$avg <- rowMeans(simp_list[[i]][,c(2:length(simp_list[[i]]))],na.rm = TRUE)
  simp_list[[i]] <- select(simp_list[[i]],Species,avg) %>% na.omit()
  names(simp_list[[i]])=c("Species","Value")
  simp_list[[i]] <- simp_list[[i]] %>% arrange(desc(Value))
  simp_list[[i]]$Code = names(simp_list)[i]
  spc_list <- union_all(spc_list,simp_list[[i]])
}

spc_list <- spc_list %>% group_by(Species) %>% filter(Value == max(Value))
spc_list <- spc_list[-c(1), ] 

simp_list_res_all <- tibble(index= 1:30,.name_repair = ~ c("Index"))

for (i in 1:length(simp_list)){
  simp_list_res_all <- cbindX(simp_list_res_all,select(simp_list[[i]],Species,Value))
}

#write.csv(simp_list_res_all,file = "simp_list_res_all.csv")

spp.scrs <- rownames_to_column(spp.scrs,var="Species")
sig.spp.scrs <- spp.scrs %>% filter(Species %in% spc_list$Species)
sig.spp.scrs <- left_join(sig.spp.scrs,spc_list)

Nmds_graph = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    data = data.scores,
    aes(colour = Season),
    size = 3,
    alpha = 0.5
  ) +
  scale_color_brewer(palette = "Dark2", name = "Saison",breaks=c("Automne","Printemps")) +
  geom_text(
    data = sig.spp.scrs,
    aes(x = MDS1, y = MDS2,colour=Code),
    fontface = "bold",
    label = sig.spp.scrs$Species
  ) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 9, colour = "grey30")
  ) +
  stat_ellipse(aes(colour = Season), level = 0.95) +
  labs(colour = "Box_name")

Nmds_graph
#ggsave(file="Nmds_season.png", plot=Nmds_graph, width=14, height=8)

# Nmds par saisons ####

MSE.spring <- MSE %>% filter(Season == "Printemps")
MSE.sub <- MSE.spring %>% select(1, all_of(Species_name))
MSE.sub <- MSE.sub %>% column_to_rownames(var = "Code")
MSE.sub <- MSE.sub[, colSums(MSE.sub != 0) > 0]
MSE.sub <- MSE.sub %>% log1p()

Nmds = metaMDS(MSE.sub, trymax = 100)
Nmds

data.scores = as.data.frame(scores(Nmds))
data.scores$Box_name = MSE.spring$Box_name
data.scores$Position = MSE.spring$Position

simp <- simper(MSE.sub,MSE.spring$Box_name,permutations=999)
simp

spp.scrs <- as.data.frame(na.omit(Nmds$species))

Megaride <- tibble("Species",.name_repair = ~ c("Species"))
Dunes <- tibble("Species",.name_repair = ~ c("Species"))
Reference <- tibble("Species",.name_repair = ~ c("Species"))
Banc <- tibble("Species",.name_repair = ~ c("Species"))

simp_list <- list("Dunes" = Dunes,"Megaride" = Megaride,"Banc" = Banc,"Reference"=Reference)

for (i in 1:length(simp)){
  a <- str_split(names(simp)[i], "_")[[1]][1]
  b <- str_split(names(simp)[i], "_")[[1]][2]
  for (j in 1:length(simp_list)){
    if (a == names(simp_list)[j]){
      l <- cbind(rownames_to_column(data.frame(simp[[i]]$ava),var = "Species"),data.frame(simp[[i]]$p))
      names(l)=c("Species","sim","p")
      l <- filter(l,p<0.05)
      l <- filter(l,sim!=0)
      simp_list[[j]] <- full_join(simp_list[[j]],select(l,Species,sim),by = "Species")
    }
    if (b == names(simp_list)[j]){
      l <- cbind(rownames_to_column(data.frame(simp[[i]]$avb),var = "Species"),data.frame(simp[[i]]$p))
      names(l)=c("Species","sim","p")
      l <- filter(l,p<0.05)
      l <- filter(l,sim!=0)
      simp_list[[j]] <- full_join(simp_list[[j]],select(l,Species,sim),by = "Species")
    }
  }
}
spc_list <- tibble(Species="",Value=0,Code="",.name_repair = ~ c("Species","Value","Code"))
for (i in 1:length(simp_list)){
  simp_list[[i]]$avg <- rowMeans(simp_list[[i]][,c(2:4)],na.rm = TRUE)
  simp_list[[i]] <- select(simp_list[[i]],Species,avg) %>% na.omit()
  names(simp_list[[i]])=c("Species","Value")
  simp_list[[i]] <- simp_list[[i]] %>% arrange(desc(Value))
  simp_list[[i]]$Code = names(simp_list)[i]
  spc_list <- union_all(spc_list,simp_list[[i]])
}

spc_list <- spc_list %>% group_by(Species) %>% filter(Value == max(Value))
spc_list <- spc_list[-c(1), ] 

simp_list_res_spring <- tibble(index= 1:30,.name_repair = ~ c("Index"))

for (i in 1:length(simp_list)){
  simp_list_res_spring <- cbindX(simp_list_res_spring,select(simp_list[[i]],Species,Value))
}

spp.scrs <- rownames_to_column(spp.scrs,var="Species")
sig.spp.scrs <- spp.scrs %>% filter(Species %in% spc_list$Species)
sig.spp.scrs <- left_join(sig.spp.scrs,spc_list)

Nmds_graph_spring = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    data = data.scores,
    aes(colour = Box_name),
    size = 3,
    alpha = 0.5
  ) +
  scale_color_brewer(palette = "Dark2", name = "Boite d'échantillonnage",breaks=c("Automne","Printemps","Banc","Dunes","Megaride","Reference")) +
  geom_text(
    data = sig.spp.scrs,
    aes(x = MDS1, y = MDS2,colour=Code),
    fontface = "bold",
    label = sig.spp.scrs$Species
  ) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 9, colour = "grey30")
  ) +
  geom_path(
    data = data.scores %>% group_by(Box_name, Position) %>% summarise_at(vars(matches("Nmds")), mean),
    aes(colour = Box_name),
    size = 2
  ) +
  geom_point(
    data = data.scores %>% group_by(Box_name, Position) %>% summarise_at(vars(matches("Nmds")), mean),
    aes(shape = Position),
    size = 2
  ) +
  stat_ellipse(aes(colour = Box_name), level = 0.90) +
  labs(colour = "Boite")

Nmds_graph_spring

MSE.autumn <- MSE %>% filter(Season == "Automne")
MSE.sub <- MSE.autumn %>% select(1, all_of(Species_name))
MSE.sub <- MSE.sub %>% column_to_rownames(var = "Code")
MSE.sub <- MSE.sub[, colSums(MSE.sub != 0) > 0]
MSE.sub <- MSE.sub %>% log1p()

Nmds = metaMDS(MSE.sub, trymax = 100)
Nmds

data.scores = as.data.frame(scores(Nmds))
data.scores$Box_name = MSE.autumn$Box_name
data.scores$Position = MSE.autumn$Position

simp <- simper(MSE.sub,MSE.autumn$Box_name,permutations=999)
simp

spp.scrs <- as.data.frame(na.omit(Nmds$species))

Megaride <- tibble("Species",.name_repair = ~ c("Species"))
Dunes <- tibble("Species",.name_repair = ~ c("Species"))
Reference <- tibble("Species",.name_repair = ~ c("Species"))
Banc <- tibble("Species",.name_repair = ~ c("Species"))

simp_list <- list("Dunes" = Dunes,"Megaride" = Megaride,"Banc" = Banc,"Reference"=Reference)

for (i in 1:length(simp)){
  a <- str_split(names(simp)[i], "_")[[1]][1]
  b <- str_split(names(simp)[i], "_")[[1]][2]
  for (j in 1:length(simp_list)){
    if (a == names(simp_list)[j]){
      l <- cbind(rownames_to_column(data.frame(simp[[i]]$ava),var = "Species"),data.frame(simp[[i]]$p))
      names(l)=c("Species","sim","p")
      l <- filter(l,p<0.05)
      l <- filter(l,sim!=0)
      simp_list[[j]] <- full_join(simp_list[[j]],select(l,Species,sim),by = "Species")
    }
    if (b == names(simp_list)[j]){
      l <- cbind(rownames_to_column(data.frame(simp[[i]]$avb),var = "Species"),data.frame(simp[[i]]$p))
      names(l)=c("Species","sim","p")
      l <- filter(l,p<0.05)
      l <- filter(l,sim!=0)
      simp_list[[j]] <- full_join(simp_list[[j]],select(l,Species,sim),by = "Species")
    }
  }
}
spc_list <- tibble(Species="",Value=0,Code="",.name_repair = ~ c("Species","Value","Code"))
for (i in 1:length(simp_list)){
  simp_list[[i]]$avg <- rowMeans(simp_list[[i]][,c(2:4)],na.rm = TRUE)
  simp_list[[i]] <- select(simp_list[[i]],Species,avg) %>% na.omit()
  names(simp_list[[i]])=c("Species","Value")
  simp_list[[i]] <- simp_list[[i]] %>% arrange(desc(Value))
  simp_list[[i]]$Code = names(simp_list)[i]
  spc_list <- union_all(spc_list,simp_list[[i]])
}

spc_list <- spc_list %>% group_by(Species) %>% filter(Value == max(Value))
spc_list <- spc_list[-c(1), ] 

simp_list_res_autumn <- tibble(index= 1:30,.name_repair = ~ c("Index"))

for (i in 1:length(simp_list)){
  simp_list_res_autumn <- cbindX(simp_list_res_autumn,select(simp_list[[i]],Species,Value))
}

spp.scrs <- rownames_to_column(spp.scrs,var="Species")
sig.spp.scrs <- spp.scrs %>% filter(Species %in% spc_list$Species)
sig.spp.scrs <- left_join(sig.spp.scrs,spc_list)

Nmds_graph_autumn = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    data = data.scores,
    aes(colour = Box_name),
    size = 3,
    alpha = 0.5
  ) +
  scale_color_brewer(palette = "Dark2", name = "Boite d'échantillonnage",breaks=c("Automne","Printemps","Banc","Dunes","Megaride","Reference")) +
  geom_text(
    data = sig.spp.scrs,
    aes(x = MDS1, y = MDS2,color=Code),
    fontface = "bold",
    label = sig.spp.scrs$Species
  ) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 9, colour = "grey30"),
    legend.position = "none"
  ) +
  geom_path(
    data = data.scores %>% group_by(Box_name, Position) %>% summarise_at(vars(matches("Nmds")), mean),
    aes(colour = Box_name),
    size = 2
  ) +
  geom_point(
    data = data.scores %>% group_by(Box_name, Position) %>% summarise_at(vars(matches("Nmds")), mean),
    aes(shape = Position),
    size = 2
  ) +
  stat_ellipse(aes(colour = Box_name), level = 0.90) +
  labs(colour = "Box_name")


Nmds_graph_autumn

image <-
  plot_grid(
    Nmds_graph_autumn,
    Nmds_graph_spring,
    rel_widths = c(1, 1.4),
    labels = c('Automne', 'Printemps'),
    label_size = 12,
    label_y = 0.98,
    label_x = 0.35
  )
image
#ggsave(file="Nmds_synth.png", plot=image, width=14, height=8)

names(simp_list_res_all) = c("index","sp1","v1","sp2","v2")
tab <- simp_list_res_all %>% filter(!is.na(v2)) %>% as_tibble()

tab <- flextable(data = tab, col_keys = c("sp1","v1","col1","sp2","v2"))
tab <-  autofit(tab)
tab <- set_header_labels(tab, Species = "",V1="",Ecart_type="Ecart-type")
tab <- hline(tab, part="body", border = std_border )
tab <- colformat_num(tab,na_str = "",digits = 2)
tab <- empty_blanks(tab)
tab
#print(tab, preview = "docx")

#### dbRDA ####


# all  #

RDA.env <-
  MSE %>% select(
    Median,
    Clay,
    Temperature,
    Temperature.at.bottom.Mean,
    Salinity.Mean,
    Current.speed.Mean
  ) %>% scale() %>% data.frame()
RDA.env <- cbind(RDA.env, select(MSE))
MSE.sub <- MSE %>% select(1, all_of(Species_name))
MSE.sub <- MSE.sub %>% column_to_rownames(var = "Code")
RDA.bio <- MSE.sub %>% log1p()

mcor <- cor(select_if(RDA.env, is.numeric), use = "complete.obs")
corrplot(
  mcor,
  type = "upper",
  method = "number",
  tl.cex = 0.5,

)

dbrda <- dbrda(RDA.bio ~ ., data = RDA.env)
dbrda
RsquareAdj(dbrda)
plot(dbrda, )
anova(dbrda, by = "terms", permu = 200)

mod0 <- dbrda(RDA.bio ~ 1, RDA.env)
mod1 <- dbrda(RDA.bio ~ ., RDA.env)

step.res <- ordiR2step(mod0, mod1, perm.max = 200)
step.res$anova
capture.output(step.res$anova,file="test.doc")


en = envfit(step.res, RDA.bio , permutations = 999, na.rm = TRUE)

spp.dbrda <-
  as.data.frame(scores(en, display = "vectors")) #save species intrinsic values into dataframe
spp.dbrda <-
  cbind(spp.dbrda, Species = rownames(spp.dbrda)) #add species names to dataframe
spp.dbrda <-
  cbind(spp.dbrda, pval = en$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
sig.spp.dbrda <- subset(spp.dbrda, pval <= 0.05)
sig.spp.dbrda #espèces les plus significatives

sppscores(step.res) <- RDA.bio

smry <- summary(step.res)

df1  <- data.frame(smry$sites[, 1:2])
dist.mat <- vegdist(df1, method = "euclidean")
clust.res <- hclust(dist.mat)
clust.res <- cutree(clust.res, k = 5) %>% data.frame()
df1$clust = clust.res$.
df1$Box_name = MSE$Box_name
df1$Position = MSE$Position
df1$Season = MSE$Season
df1$clust <- factor(df1$clust, levels = 1:6)
df2  <- data.frame(smry$biplot[, 1:2])


image = ggplot(data = df1, aes(x = dbRDA1, y = dbRDA2)) +
  geom_point(
    data = df1,
    aes(),
    shape = 3,
    size = 3,
    alpha = 0.5
  ) +
  scale_color_brewer(palette = "Dark2", name = "Box_name") +
  geom_text(
    data = df2,
    aes(x = dbRDA1, y = dbRDA2),
    colour = "blue",
    fontface = "bold",
    label = row.names(df2)
  ) +
  geom_text(
    data = sig.spp.dbrda,
    aes(x = dbRDA1, y = dbRDA2),
    colour = "orange",
    fontface = "bold",
    label = row.names(sig.spp.dbrda)
  ) +
  theme(
    axis.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    axis.ticks = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 10,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 9, colour = "grey30")
  ) +
  labs(colour = "Box_name")
image
#ggsave(file="dbRDA.png", plot=image, width=10, height=8)



#### Shannon ####

Shannon_1 = ggplot(MSE, aes(x = Season, y = Shannon, fill = Season)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Saison", y = "Indice de Shannon") +
  stat_compare_means(aes(group = Season) , label =  "p.signif")
Shannon_1

Shannon_2 = ggplot(MSE, aes(x = Box_name, y = Shannon, fill = Box_name)) +
  geom_boxplot() +
  facet_grid(cols = vars(Season),
             scales = "free",
             space = "free_x") +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Boite d'échantillonnage", y = "Indice de Shannon") +
  stat_compare_means(comparisons = my_comparisons , label =  "p.signif",label.y = c(2.7,2.7,2.7,3,3,3,2.7,2.7,2.7,3,3,3))
Shannon_2

Shannon_3 = ggplot(MSE, aes(x = Depth, y = Shannon, group = Box_name)) +
  geom_point(aes(color = Box_name)) +
  facet_grid(cols = vars(Season)) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  stat_smooth(method = "lm", se = FALSE, aes(color = Box_name)) +
  labs(x = "Profondeur (m)", y = "Indice de Shannon") +
  scale_color_brewer(palette = "Dark2", name = "Saison")
Shannon_3


#### Rich_spe ####

Rich_spe_1 = ggplot(MSE, aes(x = Season, y = Rich_spe, fill = Season)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Saison", y = "Richesse spécifique") +
  stat_compare_means(aes(group = Season) , label =  "p.signif")
Rich_spe_1

Rich_spe_2 = ggplot(MSE, aes(x = Box_name, y = Rich_spe, fill = Box_name)) +
  geom_boxplot() +
  facet_grid(cols = vars(Season),
             scales = "free",
             space = "free_x") +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Boite d'échantillonnage", y = "Richesse spécifique") +
  stat_compare_means(comparisons = my_comparisons , label =  "p.signif",label.y = c(20,20,20,23,23,23,20,20,20,23,23,23))
Rich_spe_2

Rich_spe_3 = ggplot(MSE, aes(x = Depth, y = Rich_spe, group = Box_name)) +
  geom_point(aes(color = Box_name)) +
  facet_grid(cols = vars(Season)) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  stat_smooth(method = "lm", se = FALSE, aes(color = Box_name)) +
  labs(x = "Profondeur (m)", y = "Richesse spécifique") +
  scale_color_brewer(palette = "Dark2", name = "Saison")
Rich_spe_3


#### Abondance totale ####

Ab_tot_1 = ggplot(MSE, aes(x = Season, y = Total_abundance, fill = Season)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Saison", y = "Abondance totale") +
  stat_compare_means(aes(group = Season) , label =  "p.signif")
Ab_tot_1

Ab_tot_2 = ggplot(MSE, aes(x = Box_name, y = Total_abundance, fill = Box_name)) +
  geom_boxplot() +
  facet_grid(cols = vars(Season),
             scales = "free",
             space = "free_x") +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Boite d'échantillonnage", y = "Abondance totale") +
  stat_compare_means(comparisons = my_comparisons , label =  "p.signif",label.y = c(450,450,450,500,500,500,450,450,450,500,500,500))
Ab_tot_2

Ab_tot_3 = ggplot(MSE, aes(x = Depth, y = Total_abundance, group = Box_name)) +
  geom_point(aes(color = Box_name)) +
  facet_grid(cols = vars(Season)) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  stat_smooth(method = "lm", se = FALSE, aes(color = Box_name)) +
  labs(x = "Profondeur (m)", y = "Abondance totale") +
  scale_color_brewer(palette = "Dark2", name = "Saison")
Ab_tot_3


## Piélou ####

Pielou_1 = ggplot(MSE, aes(x = Season, y = Pielou, fill = Season)) +
  geom_boxplot() +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Saison", y = "Indice d'équitabilité de Piélou") +
  stat_compare_means(aes(group = Season) , label =  "p.signif")
Pielou_1

Pielou_2 = ggplot(MSE, aes(x = Box_name, y = Pielou, fill = Box_name)) +
  geom_boxplot() +
  facet_grid(cols = vars(Season),
             scales = "free",
             space = "free_x") +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Dark2", name = "Boite d'échantillonnage") +
  labs(x = "Boite d'échantillonnage", y = "Indice d'équitabilité de Piélou") +
  stat_compare_means(comparisons = my_comparisons , label =  "p.signif",label.y = c(1.1,1.1,1.1,1.2,1.2,1.2,1.1,1.1,1.1,1.2,1.2,1.2))
Pielou_2

Pielou_3 = ggplot(MSE, aes(x = Depth, y = Pielou, group = Box_name)) +
  geom_point(aes(color = Box_name)) +
  facet_grid(cols = vars(Season)) +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  ) +
  stat_smooth(method = "lm", se = FALSE, aes(color = Box_name)) +
  labs(x = "Profondeur (m)", y = "Indice d'équitabilité de Piélou") +
  scale_color_brewer(palette = "Dark2", name = "Saison")
Pielou_3


## Synth graph #####

image = ggdraw() +
  draw_plot(Rich_spe_1, 0, .5, .5, .5) +
  draw_plot(Ab_tot_1, .5, .5, .5, .5) +
  draw_plot(Shannon_1, 0, 0, .5, .5) +
  draw_plot(Pielou_1, 0.5, 0, .5, .5)
image
#ggsave(file="Indices_bioXSaisons.png", plot=image, width=6, height=6)
image = ggdraw() +
  draw_plot(Rich_spe_2, 0, .5, .5, .5) +
  draw_plot(Ab_tot_2, .5, .5, .5, .5) +
  draw_plot(Shannon_2, 0, 0, .5, .5) +
  draw_plot(Pielou_2, 0.5, 0, .5, .5)
image
##ggsave(file="Indices_bioXBox_name.png", plot=image, width=8, height=6)
image = ggdraw() +
  draw_plot(Rich_spe_3, 0, .5, .5, .5) +
  draw_plot(Ab_tot_3, .5, .5, .5, .5) +
  draw_plot(Shannon_3, 0, 0, .5, .5) +
  draw_plot(Pielou_3, 0.5, 0, .5, .5)
image
#ggsave(file="Indices_bioXDepth.png", plot=image, width=8, height=6)

#### nb species ####

spc <- select(Dunes.bio,12:20)
spc <- unique(spc)
spc
count(spc,"Phylum")
count(spc,"Species")


