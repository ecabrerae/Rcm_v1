#Numeracion eficiente de Aquimpe con PCA.

rm(list=ls()) #Eliminar todas las variables en memoria
#Entrada de datos-----------------------------------------
setwd('F:/Proj_R/Numbering_FEM_PCA/Rcm_v1') #Se selecciona un directorio de trabajo 

#--- Librerias -------
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(scales)
library(ggrepel)
library(igraph)
library(Matrix)
library(visNetwork)
#---------------------

#------ Aqui se cargan los datos -------------------------
Nodos_orig <-read.csv("data/Ejemplo_Tequesquitengo/Nodos.csv", header = T)           #Se cargan todos los Nodos (Principales y Secundarios) con la numeracion inicial.
Triangulos_orig <-read.csv("data/Ejemplo_Tequesquitengo/Triangulos.csv", header=T)   #Se cargan los triangulos con la numeracion inicial.

Nomb_proj<-c("Ejemplo Ariguanabo")
dx_ma<-1                                        #dx para graficos de matriz de adyacencia.
x_nudg=0.5; y_nudg=0.5                          #Cuanto se despega el label de nodos de los nodos en el mapa
nod_size=3; triang_size=4; n_size=2         #Tamano del label de los nodos y de los triangulos
epsg_proj=0                                  #Si se tiene se pone el numero sino se ubica un "0"
#xlim_izq= 330000; xlim_der= 365000; dx=5000 #Para mapas de triangulacion
#ylim_izq= 340000; xlim_der= 355000; dx=5000 #Para mapas de triangulacion
#---------------------------------------------------------

ne<-nrow(Triangulos_orig) #Cantidad de triangulos "elementos".
np<-nrow(Nodos_orig)      #Cantidad de nodos (Principales + Secundarios).


#-------------------------------------------------------------------------------#
#--------- Preparacion de dataframes para mapa de numeracion ORIGINAL ----------#
#-------------------------------------------------------------------------------#

Nodos_orig$Coord_X<-as.numeric(Nodos_orig$Coord_X)
Nodos_orig$Coord_Y<-as.numeric(Nodos_orig$Coord_Y)
Nodos_orig1<-Nodos_orig; Triangulos_orig1<-Triangulos_orig

#------ Columnas nuevas en Triangulos_orig y creación de data.frame Triangulos_av_r para mapas -----------#
Triangulos_orig1$Area=c(rep(0, times=ne))      #Aqui se le agrega al data.frame de Triangulos_orig1 la columna de area de los triangulos.
Triangulos_orig1$x_centr=c(rep(0, times=ne))   #Aqui se le agrega al data.frame de Triangulos_orig1 la columna de x centro.
Triangulos_orig1$y_centr=c(rep(0, times=ne))   #Aqui se le agrega al data.frame de Triangulos_orig1 la columna de y centro.

Nodprincip_orig1<-Nodos_orig1 %>% filter(Nodos_orig1$Tipo_Nodos =="NP")
nnp<-nrow(Nodprincip_orig1)
Nodprincip_orig1$Nro<-c(1:nnp)
Nodprincip_orig1<-Nodprincip_orig1 %>% select("Nro","Nodo", "Coord_X", "Coord_Y")


source("fun_Numbering_2DcFEMmesh.R")
Lista2<- fun_area_centroids_triang(Triangulos_orig1,ne,Nodprincip_orig1,nnp)
Triangulos_orig1=data.frame(Lista2[1])

Triangulos_orig1_r<-data.frame("Nro"=as.numeric(c(seq(1,6*ne,by=1))),
                            "Nodo"=c(rep(0,times=(6*ne))),
                            "Tipo_nodo"=as.character(rep("NP",times=(6*ne))))
source("fun_Numbering_2DcFEMmesh.R")
Lista3<- fun_triang_r(Triangulos_orig1,Triangulos_orig1_r, ne, Nodos_orig1)
Triangulos_orig1_r=data.frame(Lista3[1])

source("fun_Numbering_2DcFEMmesh.R")
Lista4<- fun_ancho_banda(Triangulos_orig1, Nodos_orig1, ne)
ancho_banda_uniq_gath_orig1=data.frame(Lista4[1])
ab_orig<-ancho_banda_uniq_gath_orig1[1,1]

Madj_orig<-matrix(0,nrow=np, ncol=np) #Matriz de adyacencia Adjacency Triangulacion original

for (i in 1:ne){
  Madj_orig[Triangulos_orig[i,2], c(Triangulos_orig[i,5], Triangulos_orig[i,7])]=1
  Madj_orig[Triangulos_orig[i,3], c(Triangulos_orig[i,5], Triangulos_orig[i,6])]=1
  Madj_orig[Triangulos_orig[i,4], c(Triangulos_orig[i,6], Triangulos_orig[i,7])]=1
  Madj_orig[Triangulos_orig[i,5], c(Triangulos_orig[i,2], Triangulos_orig[i,3])]=1
  Madj_orig[Triangulos_orig[i,6], c(Triangulos_orig[i,3], Triangulos_orig[i,4])]=1
  Madj_orig[Triangulos_orig[i,7], c(Triangulos_orig[i,2], Triangulos_orig[i,4])]=1
}

polig_ab_orig<- data.frame("Nro"=as.numeric(c(seq(1, 6, by=1))),
                           "Coord_x"=as.numeric(c(1, ab_orig+1, np, np, np-(ab_orig), 1)),
                           "Coord_y"=as.numeric(c(1, 1, np-(ab_orig), np, np, (ab_orig+1))))
li_ab_orig<- data.frame("Nro"= as.numeric(c(1, 1)),
                        "Coord_x"=as.numeric(c(1, np)),
                        "Coord_y"=as.numeric(c(1, np)))

# De esta manara se puede calcular el ancho de banda mas facil------------------#
bw_orig_vect<-as.data.frame(which(Madj_orig != 0, arr.ind = T))
bw_orig_vect$bw<-(bw_orig_vect$row-bw_orig_vect$col)
bw_orig<-max(bw_orig_vect$bw) #-------------------------------------------------#

#write.csv(Madj_orig, "F:\\Proj_R\\Numbering_FEM_PCA\\Rcm_v1\\data\\Madj_orig.csv", row.names=FALSE)

#-------------------------------------------------------------------------------#
#---------- FIN preparac de dataframes para mapa de numerac ORIGINAL -----------#
#-------------------------------------------------------------------------------#


#-------------------------------------------------------------------------------#
#------- Proceso de renumeracion de nodos segun metodo de AUTOVECTORES ---------#
#-------------------------------------------------------------------------------#
Nodos_av<-Nodos_orig; Triangulos_av<-Triangulos_orig

source("fun_Numbering_2DcFEMmesh.R")
Lista1<- fun_renumera_av(Triangulos_av,Nodos_av,np)

Triangulos_av<-data.frame(Lista1[1])
Nodos_av<-data.frame(Lista1[2])
Nodprincip_av<-data.frame(Lista1[3])
nnp<-nrow(Nodprincip_av)  #Cantidad de nodos principales.
#-------------------------------------------------------------------------------#

#------ Columnas nuevas en Triangulos_av y creación de data.frame Triangulos_av_r para mapas -----------#
Triangulos_av$Area=c(rep(0, times=ne))      #Aqui se le agrega al data.frame de Triangulos_av la columna de area de los triangulos.
Triangulos_av$x_centr=c(rep(0, times=ne))   #Aqui se le agrega al data.frame de Triangulos_av la columna de x centro.
Triangulos_av$y_centr=c(rep(0, times=ne))   #Aqui se le agrega al data.frame de Triangulos_av la columna de y centro.

source("fun_Numbering_2DcFEMmesh.R")
Lista2<- fun_area_centroids_triang(Triangulos_av,ne,Nodprincip_av,nnp)
Triangulos_av=data.frame(Lista2[1])

Triangulos_av_r<-data.frame("Nro"=as.numeric(c(seq(1,6*ne,by=1))),
                            "Nodo"=c(rep(0,times=(6*ne))),
                            "Tipo_nodo"=as.character(rep("NP",times=(6*ne))))
source("fun_Numbering_2DcFEMmesh.R")
Lista3<- fun_triang_r(Triangulos_av,Triangulos_av_r, ne, Nodos_av)
Triangulos_av_r=data.frame(Lista3[1])
#-------------------------------------------------------------------------------#

source("fun_Numbering_2DcFEMmesh.R")
Lista4<- fun_ancho_banda(Triangulos_av, Nodos_av, ne)
ancho_banda_uniq_gath_av=data.frame(Lista4[1])
ab_av<-ancho_banda_uniq_gath_av[1,1]


Madj_av<-matrix(0,nrow=np, ncol=np) #Matriz de adyacencia Adjacency Triangulacion metodo av

for (i in 1:ne){
  Madj_av[Triangulos_av[i,2], c(Triangulos_av[i,5], Triangulos_av[i,7])]=1
  Madj_av[Triangulos_av[i,3], c(Triangulos_av[i,5], Triangulos_av[i,6])]=1
  Madj_av[Triangulos_av[i,4], c(Triangulos_av[i,6], Triangulos_av[i,7])]=1
  Madj_av[Triangulos_av[i,5], c(Triangulos_av[i,2], Triangulos_av[i,3])]=1
  Madj_av[Triangulos_av[i,6], c(Triangulos_av[i,3], Triangulos_av[i,4])]=1
  Madj_av[Triangulos_av[i,7], c(Triangulos_av[i,2], Triangulos_av[i,4])]=1
}

polig_ab_av<- data.frame("Nro"=as.numeric(c(seq(1, 6, by=1))),
                         "Coord_x"=as.numeric(c(1, ab_av+1, np, np, np-(ab_av), 1)),
                         "Coord_y"=as.numeric(c(1, 1, np-(ab_av), np, np, (ab_av+1))))
li_ab_av<- data.frame("Nro"= as.numeric(c(1, 1)),
                      "Coord_x"=as.numeric(c(1, np)),
                      "Coord_y"=as.numeric(c(1, np)))

#-------------------------------------------------------------------------------#
#------ FIN Proceso de renumeracion de nodos segun metodo de AUTOVECTORES ------#
#-------------------------------------------------------------------------------#



#----------------------------------------------------------------------------------------#
#------- Proceso de renumeracion de nodos segun metodo de Reverse Cuthill-McKee ---------#
#----------------------------------------------------------------------------------------#

  grad<-data.frame(nodes=c(rowSums(Madj_orig)))    #Aquí está el grado de cada nodo
  vect_ord_grad_min<-which(grad==min(grad$nodes))  #Vector donde estan los nodos con el menor grado
  Vect_ancho_semiband<-vect_ord_grad_min
  
  for (i in 1:length(vect_ord_grad_min)){
    
    source("fun_Numbering_2DcFEMmesh.R")
    nstart=vect_ord_grad_min[i]
    Rcm<-fun_Rcm(Madj_orig, np, nstart) #Aqui este el metodo Reverse Cuthill Mc-Kee
    
    Triangulos_rcm<-Triangulos_orig1
    Nodos_Rcm<-Nodos_orig
    source("fun_Numbering_2DcFEMmesh.R")
    Lista5<- fun_num_Triang_Rcm(Triangulos_rcm, Rcm, Nodos_Rcm, np)
    Triangulos_rcm=data.frame(Lista5[1])
    Nodos_Rcm=data.frame(Lista5[2])
    Nodos_Rcm$Nodo<-Nodos_Rcm$Nro #Aqui ya se tiene el data frame de con los nodos corregido por RCM
    
    Nodprincip_rcm<-Nodos_Rcm %>% filter(Nodos_Rcm$Tipo_Nodos =="NP")
    nnp<-nrow(Nodprincip_rcm)
    Nodprincip_rcm$Nro<-c(1:nnp)
    Nodprincip_rcm<-Nodprincip_rcm %>% select("Nro","Nodo", "Coord_X", "Coord_Y")
    
    source("fun_Numbering_2DcFEMmesh.R")
    Lista2<- fun_area_centroids_triang(Triangulos_rcm,ne,Nodprincip_rcm,nnp)
    Triangulos_rcm=data.frame(Lista2[1])
    
    source("fun_Numbering_2DcFEMmesh.R")
    Lista4<- fun_ancho_banda(Triangulos_rcm, Nodos_Rcm, ne)
    ancho_banda_uniq_gath_rcm=data.frame(Lista4[1])
    ab_rcm<-ancho_banda_uniq_gath_rcm[1,1]
    
    Vect_ancho_semiband[i]<-ab_rcm
  }

  #Una vez conocido en nodo de inicio que garantiza el menor ancho de banda, se calcula todo para él
  nstart=vect_ord_grad_min[which(Vect_ancho_semiband==min(Vect_ancho_semiband))[1]]
  source("fun_Numbering_2DcFEMmesh.R")
  #nstart=vect_ord_grad_min[i]
  Rcm<-fun_Rcm(Madj_orig, np, nstart) #Aqui este el metodo Reverse Cuthill Mc-Kee
  
  Triangulos_rcm<-Triangulos_orig1
  Nodos_Rcm<-Nodos_orig

  source("fun_Numbering_2DcFEMmesh.R")
  Lista5<- fun_num_Triang_Rcm(Triangulos_rcm, Rcm, Nodos_Rcm, np)

  Triangulos_rcm=data.frame(Lista5[1])
  Nodos_Rcm=data.frame(Lista5[2])

  Nodos_Rcm$Nodo<-Nodos_Rcm$Nro #Aqui ya se tiene el data frame de con los nodos corregido por RCM

  Nodprincip_rcm<-Nodos_Rcm %>% filter(Nodos_Rcm$Tipo_Nodos =="NP")
  nnp<-nrow(Nodprincip_rcm)
  Nodprincip_rcm$Nro<-c(1:nnp)
  Nodprincip_rcm<-Nodprincip_rcm %>% select("Nro","Nodo", "Coord_X", "Coord_Y")

  source("fun_Numbering_2DcFEMmesh.R")
  Lista2<- fun_area_centroids_triang(Triangulos_rcm,ne,Nodprincip_rcm,nnp)
  Triangulos_rcm=data.frame(Lista2[1])

  Triangulos_rcm_r<-data.frame("Nro"=as.numeric(c(seq(1,6*ne,by=1))),
                             "Nodo"=c(rep(0,times=(6*ne))),
                             "Tipo_nodo"=as.character(rep("NP",times=(6*ne))))

  source("fun_Numbering_2DcFEMmesh.R")
  Lista3<- fun_triang_r(Triangulos_rcm, Triangulos_rcm_r, ne, Nodos_Rcm)
  Triangulos_rcm_r=data.frame(Lista3[1])

  source("fun_Numbering_2DcFEMmesh.R")
  Lista4<- fun_ancho_banda(Triangulos_rcm, Nodos_Rcm, ne)
  ancho_banda_uniq_gath_rcm=data.frame(Lista4[1])
  ab_rcm<-ancho_banda_uniq_gath_rcm[1,1]

  
  Madj_rcm<-matrix(0,nrow=np, ncol=np) #Matriz de adyacencia Adjacency Triangulacion metodo Reverse Cuthill-McKee algorithm
  
  for (i in 1:ne){
    Madj_rcm[Triangulos_rcm[i,2], c(Triangulos_rcm[i,5], Triangulos_rcm[i,7])]=1
    Madj_rcm[Triangulos_rcm[i,3], c(Triangulos_rcm[i,5], Triangulos_rcm[i,6])]=1
    Madj_rcm[Triangulos_rcm[i,4], c(Triangulos_rcm[i,6], Triangulos_rcm[i,7])]=1
    Madj_rcm[Triangulos_rcm[i,5], c(Triangulos_rcm[i,2], Triangulos_rcm[i,3])]=1
    Madj_rcm[Triangulos_rcm[i,6], c(Triangulos_rcm[i,3], Triangulos_rcm[i,4])]=1
    Madj_rcm[Triangulos_rcm[i,7], c(Triangulos_rcm[i,2], Triangulos_rcm[i,4])]=1
  }
  
  polig_ab_rcm<- data.frame("Nro"=as.numeric(c(seq(1, 6, by=1))),
                            "Coord_x"=as.numeric(c(1, ab_rcm+1, np, np, np-(ab_rcm), 1)),
                            "Coord_y"=as.numeric(c(1, 1, np-(ab_rcm), np, np, (ab_rcm+1))))
  li_ab_rcm<- data.frame("Nro"= as.numeric(c(1, 1)),
                         "Coord_x"=as.numeric(c(1, np)),
                         "Coord_y"=as.numeric(c(1, np)))
  
  # De esta manara se puede calcular el ancho de banda ----------------------
  bw_rcm_vect<-as.data.frame(which(Madj_rcm != 0, arr.ind = T))
  bw_rcm_vect$bw<-(bw_rcm_vect$row-bw_rcm_vect$col)
  bw_rcm<-max(bw_rcm_vect$bw)
  #--------------------------------------------------------------------------  
  
#  write.csv(Madj_rcm, "F:\\Proj_R\\Numbering_FEM_PCA\\Rcm_v1\\data\\Madj_rcm.csv", row.names=FALSE)
  
#----------------------------------------------------------------------------------------#
#------ FIN Proceso de renumeracion de nodos segun metodo de Reverse Cuthill-McKee ------#
#----------------------------------------------------------------------------------------#














# Grafica 1. Mapa de triagulacion con ancho de semibanda, Métodos de Numeración: Numeración original --------------#
source("fun_Numbering_2DcFEMmesh.R")
met<-c(1,0,0)
g1<- fun_map_triang(Nomb_proj, Triangulos_orig1_r, ancho_banda_uniq_gath_orig1, Nodos_orig1, met, 
                    x_nudg, y_nudg, epsg_proj, nod_size, triang_size, n_size)
#g1
# ggsave(plot = g1, filename = 'F:/Proj_R/Numbering_FEM_PCA/Imagenes/Triangulacion_original.png', 
#        units = 'in', width = 9, height = 6, dpi = 300)
#-------------------------------------------------------------------------------------------------------------------#

# Grafica 2. Mapa de triagulacion con ancho de semibanda, Método de Numeración: Autovectores ----------------------#
source("fun_Numbering_2DcFEMmesh.R")
met<-c(0,0,0); met<-c(0,1,0)
g2<- fun_map_triang(Nomb_proj, Triangulos_av_r, ancho_banda_uniq_gath_av, Nodos_av, met, 
                    x_nudg, y_nudg, epsg_proj, nod_size, triang_size, n_size)
#g2
#-------------------------------------------------------------------------------------------#

# Grafica 3. Mapa de triagulacion con ancho de semibanda, Método de Numeración: Reverse Cuthill-McKee algorithm ---#
source("fun_Numbering_2DcFEMmesh.R")
met<-c(0,0,0); met<-c(0,0,1)
g3<- fun_map_triang(Nomb_proj, Triangulos_rcm_r, ancho_banda_uniq_gath_rcm, Nodos_Rcm, met, 
                    x_nudg, y_nudg, epsg_proj, nod_size, triang_size, n_size)
#g3
#-------------------------------------------------------------------------------------------------------------------#

# Grafica 4. Matriz Sparce, matriz de adjacencia numeracion original ----------------------#
source("fun_Numbering_2DcFEMmesh.R")
caso<-c("triangulación con numeración original.")
g4<- fun_mat_sparce(Madj_orig, Nomb_proj, ab_orig, polig_ab_orig, li_ab_orig, caso, dx_ma)
#g4
#-------------------------------------------------------------------------------------------#

# Grafica 5. Matriz Sparce, matriz de adjacencia numeracion autovector --------------------#
source("fun_Numbering_2DcFEMmesh.R")
caso<-c("triangulación con numeración método de autovector.")
g5<- fun_mat_sparce(Madj_av, Nomb_proj, ab_av, polig_ab_av, li_ab_av, caso, dx_ma)
#g5
#-------------------------------------------------------------------------------------------#

# Grafica 6. Matriz Sparce, matriz de adjacencia Reverse Cuthill-McKee algorithm ----------#
source("fun_Numbering_2DcFEMmesh.R")
caso<-c("triangulación con numeración método Reverse Cuthill-McKee.")
g6<- fun_mat_sparce(Madj_rcm, Nomb_proj, ab_rcm, polig_ab_rcm, li_ab_rcm, caso, dx_ma)
#g6
#-------------------------------------------------------------------------------------------#


pdf('F:/Proj_R/Numbering_FEM_PCA/Rcm_v1/Imagenes/Ejemplo_1.pdf', height = 6, width = 9, onefile = TRUE)
g1; g2; g3; g4; g5; g6
dev.off()











  
  
#--------------------------------------------igraph (grafo NO dinamico) ----------------------------------------------#
Meth_graph<-c("rcm")             #Aqui debe definir "original", "av" o "rcm"
Madj_graph<-Madj_rcm                  #Aqui debe definir si entra Madj_orig, Madj_av, o Madj_rcm
Nodos_graph<-Nodos_Rcm                #Aqui debe definir si entra Nodos_orig, Nodos_av, o Nodos_Rcm
Triangulos_graph<-Triangulos_rcm      #Aqui debe definir si entra Triangulos_orig, Triangulos_av, o Triangulos_rcm


# Crear grafo con la matriz de adyacencia. -----------------------------------------------------------------------#
# Aqui se obtiene la matriz de adyacencia (en el formato de los datos de entrada pero con 0 y 1 en los valores)
M1<- Madj_graph %>% as.vector %>%
  tibble(value = ., row = rep(1:nrow(Madj_graph), times = ncol(Madj_graph)),
         col = rep(1: ncol(Madj_graph), each = nrow(Madj_graph)))

M1$value<-if_else (M1$value!= 0, 1, 0)

# Aqui se crea un grafo partiendo de la matriz de adjacencia qua hay que tranformarla al estilo del paquete Matrix

M1sp<-sparseMatrix(M1$row, M1$col, x=M1$value)

gr1<-graph_from_adjacency_matrix(
  M1sp,
  mode = c("undirected"),
  weighted = NULL,
  diag = TRUE,
  add.colnames = NULL,
  add.rownames = NA
)
l<-as.matrix(Nodos_graph[,c(3,4)]) #Aqui se pasan las coord x, y de los nodos para que el grafo quede OK

plot.igraph(gr1,
            vertex.size=2, vertex.label.dist=0.5, vertex.color="red", vertex.label.cex= 0.6,
            layout = l) #Esta funcion modifica las coordenadas en las escalas de -1 a 1 en x y y.

#layout.reingold.tilford(gr1) #Con esta instruccion se pueden ver las coordenadas x e y del grafo.

#------------------------------------------------------------------------------------ FIN Grafico estatico con igraph




#------------------------------------------ visNetwork (grafo dinamico)  ---------------------------------------------#
#Nota: en visNetwork el origen de coordenadas se encuentra en el extremo superior izquiero, por eso se hace un reajuste de y

lvn<-as.matrix(Nodos_graph[,c(3,4)]) #Aqui se pasan las coord x, y de los nodos para que el grafo quede OK
lvn[,2]<-max(l[,2])-l[,2]            #Esto se hace porque visNetwork tiene el origen de coordenadas en el extremo superior izquierdo

#-- Variante -1.

visIgraph(gr1, 
          idToLabel= FALSE,
#          layout = "layout_nicely", 
          physics=TRUE) %>% 
  visNodes(shadow=TRUE)
#----------------------------------

#-- Variante -2.

data<- toVisNetworkData(gr1)  #Pasa un grafo de igraph para ser usado por VisNetwork (es un metodo LENTO)

visNetwork(nodes = data$nodes, edges = data$edges, height = "500px", width = "100%") %>%
  visNodes(shadow = TRUE, label=data$nodes$label)
#---------------------------------  

#-- Variante -3  (en este caso si se respeta las coordenadas de los nodos)

opts <- . %>% visOptions(highlightNearest = TRUE) %>%
  visInteraction(navigationButtons = TRUE, dragNodes = FALSE, 
                 dragView = FALSE, zoomView = FALSE)

visIgraph(gr1, physics=TRUE) %>% 
  visIgraphLayout(layout = "layout.norm", layoutMatrix = lvn) %>% 
  opts
#--------------------------------------

#-- Variante -4.

edges_graph<- fun_edges_graph(ne, Triangulos_graph)

df.edges_graph<-data.frame(from=(edges_graph[,1]),
                          to=(edges_graph[,2]))
df.edges_graph<-df.edges_graph[!duplicated(df.edges_graph),]

df.edges_graph$dif<-abs(df.edges_graph$to - df.edges_graph$from)                 #anchos de banda en cada linea 
df.edges_graph$ord_dif<-if_else(df.edges_graph$dif==max(df.edges_graph$dif),6,1) #Columna para espesor de la linea en el grafo
df.edges_graph$lab_ord_dif<-as.character(if_else(df.edges_graph$dif==max(df.edges_graph$dif),as.character(max(df.edges_graph$dif)),NA))


nodes_graph<-Nodos_graph
nodes_graph<- nodes_graph %>% dplyr::rename("id"="Nodo")            #Tiene que ser "id" para que visNetwork lo entienda 
nodes_graph$Nro<-as.character(nodes_graph$Nro)                      #Tienen que ser character para que sean mostrados en el grafo
nodes_graph$Tipo_integ<-if_else(nodes_graph$Tipo_Nodos=="NP",1,2)   #Para poder asignarle un color a los nodos principales y otro a los secundarios.



# Opcion 1 Grafo basico: (sin nada practicamente)
visNetwork(nodes_graph, df.edges_graph, width="100%", height="400px")
#----------------------------------------------------

# Opcion 2 Grafo algo mas elaborado.
opts <- . %>% visOptions(highlightNearest = TRUE) %>%
  visInteraction(navigationButtons = TRUE, dragNodes = FALSE, 
                 dragView = FALSE, zoomView = FALSE)

visNetwork(nodes_graph, df.edges_graph, width = "100%") %>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = lvn) %>% 
  opts
#---------------------------------------------------------------


# Opcion 3 Grafo bastante bien logrado.

if (Meth_graph==c("original")){
  caso<-c("numeración original.")
  bw<-max(df.edges_graph$dif)
}else if(Meth_graph==c("av")){
  caso<-c("numeración de nodos según: Autovector.")
  bw<-max(df.edges_graph$dif)
}else {
  caso<-c("numeración de nodos según: Reverse Cuthill−McKee.")
  bw<-max(df.edges_graph$dif)
}


vis.nodes<-nodes_graph
vis.links<-df.edges_graph

vis.nodes$shape<-"dot"                             #Forma del nodo
vis.nodes$shadow<-TRUE                             #Nodes will drop shadow
vis.nodes$title<-vis.nodes$id                      #Text on click
vis.nodes$label<-vis.nodes$Nro                     #Node Label (tiene que ser character)
vis.nodes$size<-10                                 #Node size
vis.nodes$borderWidth <- 2                         #Node border width
vis.nodes$color.background <- c("tomato", "slategrey")[nodes_graph$Tipo_integ]
vis.nodes$color.border <- "black"
vis.nodes$color.highlight.background <- "orange"   #Cuando pincho un nodo se pone orange el relleno
vis.nodes$color.highlight.border <- "darkred"      #Cuando pincho un nodo se pone darkred el borde

vis.links$label<-vis.links$lab_ord_dif             #label de los arcos con mayor ancho de banda.
vis.links$color <- c("gray", "", "", "", "", "red")[df.edges_graph$ord_dif]                          # line color  
vis.links$width <- 1*(df.edges_graph$ord_dif)       # line width
vis.links$smooth <- FALSE                          # should the edges be curved?


visnet<-visNetwork(vis.nodes, vis.links,
                   main=iconv(paste0("Mapa de la triangulación con", " ", caso), "UTF-8", "latin1"), 
                   submain=iconv(paste0("Ancho de semibanda :", bw), "UTF-8", "latin1"),
                   footer= paste0("Fuente: Proyecto", " ", Nomb_proj))%>%
  visIgraphLayout(layout = "layout.norm", layoutMatrix = lvn)

visOptions(visnet, highlightNearest = TRUE, selectedBy = "Tipo_Nodos")

#---------------------------------------------------------------



#LEEME 
# Por hacer: 1-Crear grafos dinamicos con otros paquetes, ver: https://kateto.net/network-visualization
#            2- Aplicar ggraph para hacerlo en ggplot2
#-----------------------------------------------------------------------------------------------------------------#



    # pdf('F:/Proj_R/Numbering_FEM_PCA/Imagenes/Ejemplo_1.pdf', height = 6, width = 9, onefile = TRUE)
    # g1; g4
    # dev.off()
















# M1<-matrix(0,nrow=5,ncol=5)                             #Matriz de ceros
# M1[1,1]=1; M1[1,2]=2; M1[1,3]=3; M1[1,4]=4; M1[1,5]=5
# M1<-M1
# 
# 
# M2<- M1 %>% as.vector %>% 
#   tibble(value = ., row = rep(1:nrow(M1), times = ncol(M1)),
#          col = rep(1: ncol(M1), each = nrow(M1)))
# 
# M2$value<-if_else (M2$value!= 0, 1, 0)
# M2<-as.data.frame(M2)
# 
# ggplot(data= M2, 
#        aes(x = col, y = row, group = factor(value), colour = factor(value))) +
#   geom_point(size = 2) +
#   scale_color_manual(values = c('white','black'),
#                      name="Values",
#                      labels=c("False", "True"))+
#   scale_y_reverse()+
#   scale_x_continuous(position = "top")+
#   labs(title=iconv("Matriz de tipo SPARCE.", "UTF-8", "latin1"),
#      subtitle = iconv(paste0("Ancho de semibanda : 45"), "UTF-8", "latin1"),
#      caption = 'Fuente: Proyecto La Cana.',
#      x="Rows",
#      y="Columns")+
#   cowplot::theme_cowplot()+
#   theme(panel.grid.major=element_line(color='darkgrey',
#                                       linetype = 'dashed',
#                                       size = 0.1),
#         panel.grid.minor = element_blank(),
#         panel.ontop = TRUE,
#         panel.background = element_rect(fill=NA,
#                                         color='black'))+
#   theme(plot.title=element_text(size=13, face='bold'),
#         axis.title.x = element_text(size = 11),
#         axis.title.y = element_text(size = 11),
#         axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, color = "black"),
#         axis.text.y = element_text(size = 10, color = "black"),
#         plot.caption = element_text(size = 9, color = "black"))+
#   theme(aspect.ratio = 3/4)   
# 
# #theme_minimal()


# plot(bd2$xorig, bd2$yorig, pch=21, bg="black", cex=1, col="black", 
#      main="Regresión lineal", xlab = "Coordenadas x", ylab = "Coordenadas y")
# grid(NULL,NULL)
# abline(rect_1, col="gray", lwd=3)
# abline(a=b0, b=b1, col="red", lwd=3)
# mtext(ec1, 3, line=-2); 
# 
# 
# summary(rect_1)



# plot(bd2$xorig,bd2$yorig)
# abline(rect_1,col="red",lwd=2)
# 
# 
# par(mfrow=c(1, 2))
# 
# plot(bd$x, bd$y, pch=21, bg="black", cex=1, col="black", 
#      main="Regresión lineal con residuos", xlab = "x", ylab = "y")
# grid(NULL,NULL)
# abline(recta, col="gray", lwd=3)
# mtext(ec1, 3, line=-2); 
# mtext(bquote(R^2 == .(round(r2, digits = 4))), 3, line=-3, cex=0.8)
# mtext(bquote(Sr== .(round(Sr, digits = 2))), 3, line=-4, cex=0.8)
# for (i in 1:n){
#   lines(c(bd[i,1], bd[i,1]), c(bd[i,2], eval_recta[i]), col="red", lwd=2, lty=3)
# }
# points(bd$x, bd$y, pch=21, bg="black", cex=1, col="black")
# points(bd$x, bd$eval_rect, bg="red", cex=1, col="red", pch=21)
# 
# plot(bd$x, bd$y, pch=21, cex=1, col="red", main="Regresión lineal", xlab = "x", ylab = "y")
# abline(recta, col="red", lwd=2)
# mtext(ec2, 3, line=-2); 
# mtext(bquote(R^2 == .(round(r2, digits = 4))), 3, line=-3, cex=0.8)
# #------------------------------------------------------
# summary(recta)



#bd<-bd1[,c(2,3)]