#Quantify dental complexity using molaR
library(molaR); library(geomorph)

#Set path to import files
#Meshes can be found in the Trimmed_Mesh_Files folder in this repo
S10 = 'Trimmed_Mesh_Files/S10_CS.ply' #Sw.clarki
S12 = 'Trimmed_Mesh_Files/S12_CS.ply' #Sw.clarki
S45 = 'Trimmed_Mesh_Files/S45_CS.ply' #Sw.obliquidentatus
S1 = 'Trimmed_Mesh_Files/S1_CS.ply' #Sw.whitei
S9 = 'Trimmed_Mesh_Files/S9_CS.ply' #Sw.asymmetricus

#Import trimmed mesh files
Sweet_S10 =read.ply(S10, ShowSpecimen = FALSE)
Sweet_S12 =read.ply(S12, ShowSpecimen = FALSE)
Sweet_S45 =read.ply(S45, ShowSpecimen = FALSE)
Sweet_S1 =read.ply(S1, ShowSpecimen = FALSE)
Sweet_S9 =read.ply(S9, ShowSpecimen = FALSE)

#Calculate slope
slope_S10 = molaR::Slope(Sweet_S10, Guess = FALSE)
slope_S12 = molaR::Slope(Sweet_S12, Guess = FALSE)
slope_S45 = molaR::Slope(Sweet_S45, Guess = FALSE)
slope_S1 = molaR::Slope(Sweet_S1, Guess = FALSE)
slope_S9 = molaR::Slope(Sweet_S9, Guess = FALSE)


molaR::Slope3d(slope_S10)
molaR::Slope3d(slope_S12)
molaR::Slope3d(slope_S45)
molaR::Slope3d(slope_S1)
molaR::Slope3d(slope_S9)

histogram(slope_S9$plyFile$Face_Slopes, xlim=c(0,90)) # Sw.asymmetricus
histogram(slope_S10$plyFile$Face_Slopes, xlim=c(0,90)) # Sw.clarki
histogram(slope_S1$plyFile$Face_Slopes, xlim=c(0,90)) # Sw.whitei
histogram(slope_S45$plyFile$Face_Slopes, xlim=c(0,90)) # Sw.obliquidentatus

#Calculate DNE
DNE_S10 = molaR::DNE(Sweet_S10)
DNE_S45 = molaR::DNE(Sweet_S45)
DNE_S1 = molaR::DNE(Sweet_S1)
DNE_S9 = molaR::DNE(Sweet_S9)
DNE_S12 = molaR::DNE(Sweet_S12)

molaR::DNE3d(DNE_S1)
molaR::DNE3d(DNE_S9)
molaR::DNE3d(DNE_S10)
molaR::DNE3d(DNE_S45)
molaR::DNE3d(DNE_S12)

#Calcualte OPCR
OPCr_S10 = molaR::OPCr(Sweet_S10,Steps = 10)
OPCr_S45 = molaR::OPCr(Sweet_S45, Steps = 10)
OPCr_S1 = molaR::OPCr(Sweet_S1, Steps = 10)
OPCr_S9 = molaR::OPCr(Sweet_S9, Steps = 10)
OPCr_S12 = molaR::OPCr(Sweet_S12, Steps = 10)

#Calculate OPC to visualizatize OPCr results
OPC_S10 = molaR::OPC(Sweet_S10)
OPC_S45 = molaR::OPC(Sweet_S45)
OPC_S1 = molaR::OPC(Sweet_S1)
OPC_S9 = molaR::OPC(Sweet_S9)
OPC_S12 = molaR::OPC(Sweet_S12)

molaR::OPC3d(OPC_S10)
molaR::OPC3d(OPC_S45)
molaR::OPC3d(OPC_S1)
molaR::OPC3d(OPC_S9)
molaR::OPC3d(OPC_S12)

#Histograms of slope data
library(plotly)
library(ggplot2)
library(dplyr)
library(broom)

fig_slope_S1 <- plot_ly(alpha = 0.5)
fig_slope_S1 <- fig_slope_S1 %>%
  add_histogram(
    x = slope_S1$plyFile$Face_Slopes[data.frame(slope_S1$plyFile$Face_Slopes)[,1] <= 90])
fig_slope_S1 <- fig_slope_S1 %>% layout(yaxis = list(range = c(0,375)), xaxis = list(range = c(0,30)))
fig_slope_S1

fig_slope_S9 <- plot_ly(alpha = 0.5)
fig_slope_S9 <- fig_slope_S9 %>%
  add_histogram(
    x = slope_S9$plyFile$Face_Slopes[data.frame(slope_S9$plyFile$Face_Slopes)[,1] <= 90])
fig_slope_S9 <- fig_slope_S9 %>% layout(yaxis = list(range = c(0,375)), xaxis = list(range = c(0,30)))
fig_slope_S9

fig_slope_S10 <- plot_ly(alpha = 0.5)
fig_slope_S10 <- fig_slope_S10 %>%
  add_histogram(
    x = slope_S10$plyFile$Face_Slopes[data.frame(slope_S10$plyFile$Face_Slopes)[,1] <= 90])
fig_slope_S10 <- fig_slope_S10 %>% layout(yaxis = list(range = c(0,375)),xaxis = list(range = c(0,30)))
fig_slope_S10

fig_slope_S45 <- plot_ly(alpha = 0.5)
fig_slope_S45 <- fig_slope_S45 %>%
  add_histogram(
    x = slope_S45$plyFile$Face_Slopes[data.frame(slope_S45$plyFile$Face_Slopes)[,1] <= 90])
fig_slope_S45 <- fig_slope_S45 %>% layout(yaxis = list(range = c(0,375)), xaxis = list(range = c(0,30)))
fig_slope_S45

#Calculate slope of lines for histogram values
iter = list(slope_S1$plyFile$Face_Slopes
         ,slope_S9$plyFile$Face_Slopes
         ,slope_S10$plyFile$Face_Slopes
         ,slope_S45$plyFile$Face_Slopes)

breaks = seq(0,29,1)
y_vals <- lapply(iter, FUN = function(x) hist(x[data.frame(x)[,1] <= 30],breaks = seq(0,30,1))$counts)
fit <- lapply(y_vals[], FUN = function(x) lm(x~breaks))
loessFit <- fit <- lapply(y_vals[], FUN = function(x) loess(x~breaks))

#plot lm results
fig_lm <- plot_ly(x = ~breaks) %>%
  add_markers(y = y_vals[[1]], name = "Sw.whitei(LA)") %>%
  add_lines(y = ~fitted(loessFit[[1]]),  name = "Sw.whitei(LA)") %>%
  layout(xaxis = list(title='Slope (degrees)'),
         yaxis = list(title='Triangle Count')) %>%
  add_ribbons(data = augment(loessFit[[1]]),
              ymin = ~.fitted - 1.96*.se.fit,
              ymax = ~.fitted + 1.96*.se.fit,
              fillcolor = 'rgba(168,216,234,0.1)',
              line = list(color = 'rgba(168,216,234,0.1)'),
              showlegend = FALSE)

fig_lm <- fig_lm %>%
  add_markers(y = y_vals[[2]], name = "Sw.asymmetricus(S)") %>%
  add_lines(x = ~breaks, y = fitted(loessFit[[2]]), name = "Sw.asymmetricus(S)") %>%
  add_ribbons(data = augment(loessFit[[2]]),
              ymin = ~.fitted - 1.96*.se.fit,
              ymax = ~.fitted + 1.96*.se.fit,
              fillcolor = 'rgba(100,150,80,0.1)',
              line = list(color = 'rgba(100,150,80,0.1)'),
              showlegend = FALSE)

fig_lm <- fig_lm %>%
  add_markers(y = y_vals[[3]], name = "Sw.clarki(S)") %>%
 add_lines(x = ~breaks, y = fitted(loessFit[[3]]), name = "Sw.clarki(S)") %>%
  add_ribbons(data = augment(loessFit[[3]]),
              ymin = ~.fitted - 1.96*.se.fit,
              ymax = ~.fitted + 1.96*.se.fit,
              fillcolor = 'rgba(100,150,80,0.1)',
              line = list(color = 'rgba(100,150,80,0.1)'),
              showlegend = FALSE)

fig_lm <- fig_lm %>%
  add_markers(y = y_vals[[4]], name = "Sw.obliquidentatus(LA)") %>%
  add_lines(x = ~breaks, y = fitted(loessFit[[4]]), name = "Sw.obliquidentatus(LA)") %>%
  add_ribbons(data = augment(loessFit[[4]]),
              ymin = ~.fitted - 1.96*.se.fit,
              ymax = ~.fitted + 1.96*.se.fit,
              fillcolor = 'rgba(168,216,234,0.1)',
              line = list(color = 'rgba(168,216,234,0.1)'),
              showlegend = FALSE)
fig_lm
