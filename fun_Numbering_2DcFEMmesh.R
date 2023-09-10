

#--------------------------------------------Funcion ------ fun_mapa_ef_triang mapa de eficiencia de triangulos.
fun_renumera_av<-function(Triangulos_av,Nodos_av,np) {
  
  #------ Matriz de varianzas y covarianzas.
  Nodos_av$Coord_X2<-(Nodos_av$Coord_X)^2 #Cuadrados de cada valor de las xorig.
  Nodos_av$Coord_Y2<-(Nodos_av$Coord_Y)^2 #Cuadrados de cada valor de las yorig.
  Nodos_av$Coord_XY<-(Nodos_av$Coord_X)*(Nodos_av$Coord_Y) #multiplicacion de las xorig y yorig.
  
  vacov<-data.frame("c1"=c(0,0),"c2"=c(0,0))  #Matriz de varianzas y covarianzas.
  
  vacov[1,1]<-((np*sum(Nodos_av$Coord_X2))-(sum(Nodos_av$Coord_X)*sum(Nodos_av$Coord_X)))/(np*(np-1))
  vacov[1,2]<-((np*sum(Nodos_av$Coord_XY))-(sum(Nodos_av$Coord_X)*sum(Nodos_av$Coord_Y)))/(np*(np-1))
  vacov[2,1]<-((np*sum(Nodos_av$Coord_XY))-(sum(Nodos_av$Coord_Y)*sum(Nodos_av$Coord_X)))/(np*(np-1))
  vacov[2,2]<-((np*sum(Nodos_av$Coord_Y2))-(sum(Nodos_av$Coord_Y)*sum(Nodos_av$Coord_Y)))/(np*(np-1))
  #--------------------------------------------------------
  
  #------ Autovalores y autovectores 
  e<-eigen(vacov) #Funcion propia de r para calculo de autovalores y autovectores.
  #--------------------------------------------------------
  
  #Autovalores ---------
  a<-1; b<--(vacov[1,1]+vacov[2,2]); c<-(vacov[1,1]*vacov[2,2])-(vacov[1,2]*vacov[2,1])
  disc<-(b^2)-(4*a*c)  #Discriminante.
  av<-data.frame("a"=c(0,0))
  av[1,1]<-((-b)+sqrt(disc))/(2*a)
  av[2,1]<-((-b)-sqrt(disc))/(2*a)
  
  av1<-as.data.frame(e$values) #Obtenido con funcion eigen de r.
  names (av1) = c("a")
  #-----------------------------------------------------------------------
  
  #Autovectores --------
  avect<-as.data.frame(e$vectors)
  names (avect) = c("avu_x","avu_y") #avu: Autovector unitario.
  #-----------------------------------------------------------------------
  
  
  Nodos_av$xtransf<-rep(0,times=np)
  Nodos_av$ytransf<-rep(0,times=np)
  
  for (i in 1:np){
    Nodos_av[i,9]<-(Nodos_av[i,3]*avect[1,1])+(Nodos_av[i,4]*avect[2,1])
    Nodos_av[i,10]<-(Nodos_av[i,3]*avect[1,2])+(Nodos_av[i,4]*avect[2,2])
  }
  
  bd1a<-Nodos_av[order(Nodos_av$xtransf, decreasing = TRUE),]
  bd1a$Num_nodos<-c(1:np)
  bd1a$Nro<-c(1:np)
  
  
  Triangulos_av1<-Triangulos_av
  Triangulos_av1[,c(2:7)]=0
  
  bd1anp<-bd1a %>% filter(Tipo_Nodos=="NP")
  nnp<-nrow(bd1anp)#Cantidad de Nodos Principales.
  bd1anp$Nro<-c(1:nnp)
  
  for (i in 1:nnp){
    for (j in 1:3){
      Triangulos_av1[c(which(Triangulos_av[,j+1]==bd1anp[i,2])),j+1]= bd1anp[i,11]  
    }
  }
  
  bd1ans<-bd1a %>% filter(Tipo_Nodos=="NS")
  nns<-nrow(bd1ans)#Cantidad de Nodos Secundarios.
  bd1ans$Nro<-c(1:nns)
  
  for (i in 1:nns){
    for (j in 1:3){
      Triangulos_av1[c(which(Triangulos_av[,j+4]==bd1ans[i,2])),j+4]= bd1ans[i,11]  
    }
  }
  
  Triangulos_av<-Triangulos_av1
  
  bd1a$Nodo<-bd1a$Num_nodos
  Nodos_av<-bd1a %>% select("Nro","Nodo", "Coord_X", "Coord_Y", "Tipo_Nodos")
  
  Nodprincip_av<-Nodos_av %>% filter(Nodos_av$Tipo_Nodos =="NP")
  Nodprincip_av$Nro<-c(1:nnp)
  Nodprincip_av<-Nodprincip_av %>% select("Nro","Nodo", "Coord_X", "Coord_Y")
  

  Lista1<-list(Triangulos_av,Nodos_av, Nodprincip_av)
  return(Lista1)
}
#-----------------------------------------------------------------------------------


#--------------------------------------------Funcion ------ fun_area_centroids_triang

fun_area_centroids_triang<-function(Triangulos,ne,Nodprincip,nnp) {
  
  for (i in 1:ne){
    
    np1<-Triangulos[i,2] 
    j=which(Nodprincip$Nodo==np1)
    xnp1=Nodprincip[j,3]; ynp1=Nodprincip[j,4]  #coordenadas x,y del primer nodo principal del triangulo
    
    np2<-Triangulos[i,3] 
    j=which(Nodprincip$Nodo==np2)
    xnp2=Nodprincip[j,3];  ynp2=Nodprincip[j,4] #coordenadas x,y del segundo nodo principal del triangulo
    
    np3<-Triangulos[i,4]
    j=which(Nodprincip$Nodo==np3)
    xnp3=Nodprincip[j,3];  ynp3=Nodprincip[j,4] #coordenadas x,y del tercer nodo principal del triangulo
    
    Area_Triang<-abs((1/2)*(((xnp2*ynp3)-(xnp3*ynp2))+((xnp3*ynp1)-(xnp1*ynp3))+((xnp1*ynp2)-(xnp2*ynp1)))) #Aqui se calculan las Areas de los triangulos
    
    Triangulos[i,12]= Area_Triang
    Triangulos[i,13]=(xnp1+xnp2+xnp3)/3
    Triangulos[i,14]=(ynp1+ynp2+ynp3)/3
    
  }
  
  # xmin=min(Nodprincip[,3])
  # xmax=max(Nodprincip[,3])
  
  Lista2<-list(Triangulos) #,g4)
  return(Lista2)
  
}   
#-----------------------------------------------------------------------------------


#--------------------------------------------Funcion ------ Fun_triang_r Se llena el data.frame fun_triang_r

fun_triang_r<-function(Triangulos,Triangulos_r, ne,Todos_los_nodos) {
  
  for (i in 1:ne){
    for (j in 1:6){
      Triangulos_r[((i-1)*6)+j, 1]=i
      Triangulos_r[((i-1)*6)+j, 2]=Triangulos[i,j+1]
      if(j==4|j==5|j==6){Triangulos_r[((i-1)*6)+j, 3]="NS"}
      #    Triangulos_r[((i-1)*6)+j, 4]=Triangulos[i,8]
    }
  }
  
  Todos_los_nodos<-Todos_los_nodos %>% select("Nodo","Coord_X", "Coord_Y")
  
  Triangulos_r<-left_join(Triangulos_r, Triangulos, by = "Nro") %>% 
    select(-NP1, -NP2, -NP3,-NS1, -NS2, -NS3)
  Triangulos_r<-left_join(Triangulos_r, Todos_los_nodos, by = "Nodo") 
  #colnames(Propiedades)<-c("Nrog","Num_Mat","Kd","E")
  #Triangulos_r<-left_join(Triangulos_r, Propiedades, by = "Num_Mat") %>% 
  #  select(-Nrog)
  #Triangulos_r<-Triangulos_r %>% mutate(Td=Kd*Esp_Cota, .keep = "all")
  
  Lista3<-list(Triangulos_r)
  return(Lista3)
}
#--------------------------------------------------------------------------------#



#--------------------------------------------Funcion ------ fun_ancho_banda Se llena el data.frame ancho_banda_uniq_gath_av

fun_ancho_banda<-function(Triang, Nods, ne) {

  ancho_banda<-data.frame("N1"=as.integer(c(rep(0,times=(6*ne)))),
                          "N2"=as.integer(c(rep(0,times=(6*ne))))) 
  j=1
  for (i in 1:ne){
    
    if (Triang[i,2]< Triang[i,5]){
      ancho_banda[j,1]<-Triang[i,2]
      ancho_banda[j,2]<-Triang[i,5]
    }else{
      ancho_banda[j,2]<-Triang[i,2]
      ancho_banda[j,1]<-Triang[i,5]
    }
    j=j+1
    
    if (Triang[i,5]< Triang[i,3]){
      ancho_banda[j,1]<-Triang[i,5]
      ancho_banda[j,2]<-Triang[i,3]
    }else{
      ancho_banda[j,2]<-Triang[i,5]
      ancho_banda[j,1]<-Triang[i,3]
    }  
    j=j+1
    
    if(Triang[i,3]<Triang[i,6]){
      ancho_banda[j,1]<-Triang[i,3]
      ancho_banda[j,2]<-Triang[i,6]
    }else{
      ancho_banda[j,2]<-Triang[i,3]
      ancho_banda[j,1]<-Triang[i,6]
    }
    j=j+1
    
    if(Triang[i,6]<Triang[i,4]){
      ancho_banda[j,1]<-Triang[i,6]
      ancho_banda[j,2]<-Triang[i,4]
    }else{
      ancho_banda[j,2]<-Triang[i,6]
      ancho_banda[j,1]<-Triang[i,4] 
    }  
    j=j+1
    
    if(Triang[i,4]<Triang[i,7]){
      ancho_banda[j,1]<-Triang[i,4]
      ancho_banda[j,2]<-Triang[i,7] 
    }else{
      ancho_banda[j,2]<-Triang[i,4]
      ancho_banda[j,1]<-Triang[i,7]  
    }  
    j=j+1
    
    if(Triang[i,7]<Triang[i,2]){
      ancho_banda[j,1]<-Triang[i,7]
      ancho_banda[j,2]<-Triang[i,2]  
    }else{
      ancho_banda[j,2]<-Triang[i,7]
      ancho_banda[j,1]<-Triang[i,2] 
    }  
    j=j+1
  }
  
  ancho_banda$ab<-ancho_banda$N2-ancho_banda$N1
  ancho_banda_uniq <-ancho_banda[!duplicated(ancho_banda), ]
  ancho_banda_uniq<-ancho_banda_uniq [ancho_banda_uniq$ab == max(ancho_banda_uniq$ab), ]
  ancho_banda_uniq$line_r<-c(1:nrow(ancho_banda_uniq))
  
  ancho_banda_uniq_gath<-gather(ancho_banda_uniq,"Nodos","Nodos_num", 1:2)
  ancho_banda_uniq_gath$Coord_X<-c(rep(0, times=nrow(ancho_banda_uniq_gath)))
  ancho_banda_uniq_gath$Coord_Y<-c(rep(0, times=nrow(ancho_banda_uniq_gath)))
  
  for (i in 1:nrow(ancho_banda_uniq_gath)){
    ancho_banda_uniq_gath[i, 5]<-Nods[which(Nods$Nodo==c(ancho_banda_uniq_gath[i, 4])), 3]
    ancho_banda_uniq_gath[i, 6]<-Nods[which(Nods$Nodo==c(ancho_banda_uniq_gath[i, 4])), 4]
  }
  
  Lista4<-list(ancho_banda_uniq_gath)
  return(Lista4)  

}
#--------------------------------------------------------------------------------#



#--------------------------------------------Funcion ------ Funcion para el metodo Reverse Cuthill-Mckee

fun_Rcm<-function(Madj_orig, np, nstart) {

  #grad<-data.frame(nodes=c(rowSums(Madj_orig))) #Aquí está el grado de cada vector
  Rcm<-data.frame(nodes=c(rep(0, np)))
  #Rcm[1,1]=which(grad==min(grad$nodes))[1]
  Rcm[1,1]=nstart
  
  k=1; m=1
  
  while (Rcm[nrow(Rcm),1]==0){
    
    a=which(Madj_orig[Rcm[k,1],]==1, arr.ind=TRUE) #Nodos adyacentes al nodo ubicando en el indice 'k' en Rcm. 
    a=setdiff(a, Rcm[1:m,1]) #Idem al anterior pero se garantiza que no se repitan con los ya definidos (numerados)
    
    b=as.integer(grad[c(a),1]) #Se obtiene el grado de cada nodo de a
    a<-a[order(b)]
    b<-b[order(b)]
    
    if (length(a)!=0) {
      Rcm[m+1:length(a),1]<-a
    }
    
    m=m+length(a)
    k=k+1
  }
  #Rcm<-as.numeric(c(Rcm[[1]]))   #Cuthill-Mckee
  Rcm<-rev(Rcm[,1])               #Reverse Cuthill-Mckee
  
  return(Rcm)
 
}
#--------------------------------------------------------------------------------#



#--------------------------------------------Funcion ----------------------------#
fun_num_Triang_Rcm<-function(Triangulos_rcm, Rcm, Nodos_Rcm, np) {

  for (i in 1:np){
    Nodos_Rcm[which(Nodos_Rcm$Nodo==Rcm[i]), 1]<-i
  }
  Nodos_Rcm <- Nodos_Rcm %>% dplyr::arrange(Nro)
  
  
  for (i in 1:6){
    for (j in 1:ne){
      Triangulos_rcm[j,i+1]<-Nodos_Rcm[which(Nodos_Rcm$Nodo==Triangulos_rcm[j,i+1]), 1]     
    }
  }
  
Lista5<-list(Triangulos_rcm, Nodos_Rcm)
return(Lista5)  

}





#--------------------------------------------Funcion ------ fun_mapa_triang mapa de triangulacion con ancho de banda para todos los métodos.
fun_map_triang<-function(Nomb_proj, Triang, ab_uniq_gath, Nods, met, x_nudg, y_nudg, epsg_proj, 
                         nod_size, triang_size, n_size) {
  
  ab_uniq_gath$ab<-as.factor(ab_uniq_gath$ab)
  ab_uniq_gath$Nodos<-as.factor(ab_uniq_gath$Nodos)
 
  if (epsg_proj==0){
   
    if (which(met==1)==1){
      
      caso<-c("numeración original.")
      
      if (nrow(ab_uniq_gath)>2){
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, group=line_r, colour=ab), linewidth = 0.6)+
          geom_point(data = Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, 
                     shape = 21, fill="white", stroke = 0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x= x_nudg, y = y_nudg))+
          #coord_sf(datum = "EPSG:3795")+
          #scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          #scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv(paste0("Mapa de la triangulación con", " ", caso), "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        #+theme(aspect.ratio = 3/4)
        
      }else{
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, colour=ab), linewidth = 0.6)+
          geom_point(data=Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, 
                     shape = 21, fill="white", stroke = 0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x= x_nudg, y = y_nudg))+
          #coord_sf(datum = "EPSG:3795")+
          # scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          # scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv(paste0("Mapa de la triangulación con", " ", caso), "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        #+theme(aspect.ratio = 3/4)   
      } 
    }  
    
    if (which(met==1)==2){
      
      rect_1<-lm(formula=Coord_Y~Coord_X, data=Nods)
      eval_rect_1=predict(rect_1); Coeficientes=coef(rect_1)
      a0=round(Coeficientes[1], digits = 5)
      a1=round(Coeficientes[2], digits = 5)
      ec1=ifelse(a0<0,paste('y=',a1,'x',a0), paste('y=',a1,'x +',a0))
      
      xmed<-mean(Nods$Coord_X)
      ymed_mod<-(a1*xmed)+a0
      
      b1<-(-1/a1)
      b0<-ymed_mod-(b1*xmed)
      
      caso<-c("numeración de nodos según: Autovector.")
      
      if (nrow(ab_uniq_gath)>2){
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, group=line_r, colour=ab), linewidth = 0.6)+
          geom_point(data=Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, shape = 21, fill="white", stroke=0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+ 
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x= x_nudg, y = y_nudg))+
          #coord_sf(datum="EPSG:3795")+
          # scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          # scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+  
          geom_abline(intercept=a0, slope=a1, color="blue", lty=1, lwd=1)+
          geom_abline(intercept=b0, slope=b1, color="grey", lwd=1)+
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv(paste0("Mapa de la triangulación con"," ", caso), "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        #theme(aspect.ratio = 3/4)   
        
      }else{
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, colour=ab), linewidth = 0.6)+
          geom_point(data=Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, shape = 21, fill="white", stroke=0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+ 
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x= x_nudg, y = y_nudg))+
          geom_abline(intercept=a0, slope=a1, color="blue", lty=1, lwd=1)+
          geom_abline(intercept=b0, slope=b1, color="grey", lwd=1)+
          #coord_sf(datum="EPSG:3795")+
          # scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          # scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv("Mapa de la triangulación con método de numeración: Autovector.", "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        #theme(aspect.ratio = 3/4)   
      }    
    }
    
    if (which(met==1)==3){
      
      caso<-c("numeración de nodos según: Reverse Cuthill-McKee.")
      
      if (nrow(ab_uniq_gath)>2){
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, group=line_r, colour=ab), linewidth = 0.6)+
          geom_point(data = Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, shape = 21, fill="white", stroke=0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x=x_nudg, y = y_nudg))+
          #coord_sf(datum="EPSG:3795")+
          # scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          # scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+        
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv(paste0("Mapa de la triangulación con", " ", caso), "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        #theme(aspect.ratio = 3/4)   
        
      }else{
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, colour=ab), linewidth = 0.6)+
          geom_point(data=Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, shape = 21, fill="white", stroke=0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x=x_nudg, y = y_nudg))+
          #coord_sf(datum="EPSG:3795")+
          # scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          # scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv(paste0("Mapa de la triangulación con", " ", caso), "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        # theme(aspect.ratio = 3/4)   
      } 
      
    }  
    
    return(g) 
    
  }else{

    if (which(met==1)==1){
      
      caso<-c("numeración original.")
      
      if (nrow(ab_uniq_gath)>2){
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, group=line_r, colour=ab), linewidth = 0.6)+
          geom_point(data = Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, 
                     shape = 21, fill="white", stroke = 0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x= x_nudg, y = y_nudg))+
          coord_sf(datum = "EPSG:3795")+
          #scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          #scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv(paste0("Mapa de la triangulación con", " ", caso), "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        #+theme(aspect.ratio = 3/4)
        
      }else{
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, colour=ab), linewidth = 0.6)+
          geom_point(data=Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, 
                     shape = 21, fill="white", stroke = 0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x= x_nudg, y = y_nudg))+
          coord_sf(datum = "EPSG:3795")+
          # scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          # scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv(paste0("Mapa de la triangulación con", " ", caso), "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        #+theme(aspect.ratio = 3/4)   
      } 
    }  
    
    if (which(met==1)==2){
      
      rect_1<-lm(formula=Coord_Y~Coord_X, data=Nods)
      eval_rect_1=predict(rect_1); Coeficientes=coef(rect_1)
      a0=round(Coeficientes[1], digits = 5)
      a1=round(Coeficientes[2], digits = 5)
      ec1=ifelse(a0<0,paste('y=',a1,'x',a0), paste('y=',a1,'x +',a0))
      
      xmed<-mean(Nods$Coord_X)
      ymed_mod<-(a1*xmed)+a0
      
      b1<-(-1/a1)
      b0<-ymed_mod-(b1*xmed)
      
      caso<-c("numeración de nodos según: Autovector.")
      
      if (nrow(ab_uniq_gath)>2){
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, group=line_r, colour=ab), linewidth = 0.6)+
          geom_point(data=Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, shape = 21, fill="white", stroke=0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+ 
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x= x_nudg, y = y_nudg))+
          coord_sf(datum="EPSG:3795")+
          # scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          # scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+  
          geom_abline(intercept=a0, slope=a1, color="blue", lty=1, lwd=1)+
          geom_abline(intercept=b0, slope=b1, color="grey", lwd=1)+
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv(paste0("Mapa de la triangulación con"," ", caso), "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        #theme(aspect.ratio = 3/4)   
        
      }else{
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, colour=ab), linewidth = 0.6)+
          geom_point(data=Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, shape = 21, fill="white", stroke=0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+ 
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x= x_nudg, y = y_nudg))+
          geom_abline(intercept=a0, slope=a1, color="blue", lty=1, lwd=1)+
          geom_abline(intercept=b0, slope=b1, color="grey", lwd=1)+
          coord_sf(datum="EPSG:3795")+
          # scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          # scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv("Mapa de la triangulación con método de numeración: Autovector.", "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        #theme(aspect.ratio = 3/4)   
      }    
    }
    
    if (which(met==1)==3){
      
      caso<-c("numeración de nodos según: Reverse Cuthill-McKee.")
      
      if (nrow(ab_uniq_gath)>2){
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, group=line_r, colour=ab), linewidth = 0.6)+
          geom_point(data = Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, shape = 21, fill="white", stroke=0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x=x_nudg, y = y_nudg))+
          coord_sf(datum="EPSG:3795")+
          # scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          # scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+        
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv(paste0("Mapa de la triangulación con", " ", caso), "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        #theme(aspect.ratio = 3/4)   
        
      }else{
        
        g <- ggplot() +
          geom_polygon(data=Triang%>% filter(Tipo_nodo %in% c("NP")), 
                       aes(x = Coord_X, y = Coord_Y, group = Nro, fill = "black"), 
                       color="black", alpha=0.05, linetype=1, size=0.1)+
          geom_line(data=ab_uniq_gath, 
                    mapping = aes(x= Coord_X, y= Coord_Y, colour=ab), linewidth = 0.6)+
          geom_point(data=Triang, aes(x = Coord_X, y = Coord_Y), color="black", size=n_size, shape = 21, fill="white", stroke=0.1)+
          geom_text(data = Triang, aes(x = x_centr, y = y_centr, group = Nro, label = paste0(Nro)), size = triang_size, col = "gray70")+
          geom_text(data = Triang %>% distinct(Nodo, .keep_all=TRUE), aes(x = Coord_X, y = Coord_Y, label = paste0(Nodo)),
                    size = nod_size, col = "black", position = position_nudge(x=x_nudg, y = y_nudg))+
          coord_sf(datum="EPSG:3795")+
          # scale_x_continuous(limits = c(330000, 365000), breaks = seq(330000, 365000, 5000), labels = comma)+
          # scale_y_continuous(limits = c(340000, 355000), breaks = seq(340000, 355000, 5000), labels = comma)+
          scale_fill_discrete(name="Triángulos", labels=c(""))+
          scale_colour_manual(name="Ancho de banda", values=c("red"),
                              labels=c(paste0(ab_uniq_gath[1,1])))+
          labs(title=iconv(paste0("Mapa de la triangulación con", " ", caso), "UTF-8", "latin1"),
               subtitle = iconv(paste0("Ancho de semibanda :", ab_uniq_gath[1,1]), "UTF-8", "latin1"),
               caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
               x="Coordenas X",
               y="Coordenadas Y")+
          cowplot::theme_cowplot()+
          theme(panel.grid.major=element_line(color='darkgrey',
                                              linetype = 'dashed',
                                              size = 0.1),
                panel.grid.minor = element_blank(),
                panel.ontop = FALSE,
                panel.background = element_rect(fill=NA,
                                                color='black'))+
          theme(plot.title=element_text(size=11, face='bold'),
                plot.subtitle = element_text(size=10),
                legend.title = element_text(size=10),
                legend.text = element_text(size=9),
                axis.title.x = element_text(size = 9),
                axis.title.y = element_text(size = 9),
                axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5, color = "black"),
                axis.text.y = element_text(size = 9, angle=90, hjust = 0.5, color = "black"),
                plot.caption = element_text(size = 8, color = "black"))
        # theme(aspect.ratio = 3/4)   
      } 
      
    }  
    
    return(g)    
    
  }
  
  
}
#--------------------------------------------------------------------------------#



#--------------------------------------------Funcion ------ fun_mat_sparce.
fun_mat_sparce<-function(M1, Nomb_proj, ab, polig_ab, li_ab, caso, dx_ma) {

  # M1<-matrix(0,nrow=5,ncol=5)                             #Matriz de ceros
  # M1[1,1]=1; M1[1,2]=2; M1[1,3]=3; M1[1,4]=4; M1[1,5]=5
  # M1<-M1
  M1<-as.matrix(M1)

  M2<- M1 %>% as.vector %>%
    tibble(value = ., row = rep(1:nrow(M1), times = ncol(M1)),
           col = rep(1: ncol(M1), each = nrow(M1)))

  M2$value<-if_else (M2$value!= 0, 1, 0)
  #  M2<-as.data.frame(M2)

g<- ggplot() +
  geom_point(data=M2,
            mapping= aes(x = col, y = row, group = factor(value), colour = factor(value)), size=2, shape = 18)+
  geom_line(data=li_ab, mapping = aes(x= Coord_x, y= Coord_y), color="red", linewidth = 1)+
  geom_polygon(data=polig_ab, aes(x= Coord_x, y= Coord_y, fill="red"),
               color="black", alpha=0.03, linetype=1, size=0.4)+
  scale_fill_discrete(name="Area", labels=c(""))+
  scale_colour_manual(name="Values", values = c("white","black"), labels=c("False", "True"))+
  scale_y_reverse(breaks = seq(from = 1, to = nrow(M2), by = dx_ma))+
  scale_x_continuous(position = "top", breaks = seq(from = 1, to = nrow(M2), by = dx_ma))+
  labs(title=iconv(paste0("Matriz de tipo SPARCE,", " ", caso), "UTF-8", "latin1"),
         subtitle = iconv(paste0("Ancho de semibanda:", " ", ab), "UTF-8", "latin1"),
         caption = paste0("Fuente: Proyecto", " ", Nomb_proj),
         x="Rows",
         y="Columns")+
    cowplot::theme_cowplot()+
    theme(panel.grid.major=element_line(color='grey90',
                                        size = 0.1,
                                        linetype = 'dashed'),
          panel.grid.minor = element_blank(),
          panel.ontop = FALSE,
          panel.background = element_rect(fill=NA,
                                          color='black'))+
    theme(plot.title=element_text(size=13, face='bold'),
          axis.title.x = element_text(size = 11),
          axis.title.y = element_text(size = 11),
          axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          plot.caption = element_text(size = 9, color = "black"))+
    theme(aspect.ratio = 3/4)

  return(g)
}
#--------------------------------------------------------------------------------#





#--------------------------------------------Funcion ------ fun_edges_graph.
fun_edges_graph<-function(ne, Triangulos_graph) {
  
  #edges_orig<-matrix(0, nrow=6*ne, ncol=2)
  ed_graph_1<-matrix(0, nrow=6, ncol=2)
  for (i in 1:ne){
    
    if (Triangulos_graph[i,2]<Triangulos_graph[i,5]){
      ed_graph_1[1,1]<-Triangulos_graph[i,2]
      ed_graph_1[1,2]<-Triangulos_graph[i,5]
    }else{
      ed_graph_1[1,1]<-Triangulos_graph[i,5]
      ed_graph_1[1,2]<-Triangulos_graph[i,2]
    }
    
    if (Triangulos_graph[i,3]<Triangulos_graph[i,5]){
      ed_graph_1[2,1]<-Triangulos_graph[i,3]
      ed_graph_1[2,2]<-Triangulos_graph[i,5]
    }else{
      ed_graph_1[2,1]<-Triangulos_graph[i,5]
      ed_graph_1[2,2]<-Triangulos_graph[i,3]
    }  
    
    if (Triangulos_graph[i,3]<Triangulos_graph[i,6]){
      ed_graph_1[3,1]<-Triangulos_graph[i,3]
      ed_graph_1[3,2]<-Triangulos_graph[i,6]
    }else{
      ed_graph_1[3,1]<-Triangulos_graph[i,6]
      ed_graph_1[3,2]<-Triangulos_graph[i,3]
    }    
    
    if (Triangulos_graph[i,4]<Triangulos_graph[i,6]){
      ed_graph_1[4,1]<-Triangulos_graph[i,4]
      ed_graph_1[4,2]<-Triangulos_graph[i,6]
    }else{
      ed_graph_1[4,1]<-Triangulos_graph[i,6]
      ed_graph_1[4,2]<-Triangulos_graph[i,4]
    }    
    
    if (Triangulos_graph[i,4]<Triangulos_graph[i,7]){
      ed_graph_1[5,1]<-Triangulos_graph[i,4]
      ed_graph_1[5,2]<-Triangulos_graph[i,7]
    }else{
      ed_graph_1[5,1]<-Triangulos_graph[i,7]
      ed_graph_1[5,2]<-Triangulos_graph[i,4]
    }      
    
    if (Triangulos_graph[i,2]<Triangulos_graph[i,7]){
      ed_graph_1[6,1]<-Triangulos_graph[i,2]
      ed_graph_1[6,2]<-Triangulos_graph[i,7]
    }else{
      ed_graph_1[6,1]<-Triangulos_graph[i,7]
      ed_graph_1[6,2]<-Triangulos_graph[i,2]
    }        
    
    if (i==1){
      edges_graph=ed_graph_1
    }else{
      edges_graph=rbind(edges_graph, ed_graph_1)
    }
  }  
  
  
  return(edges_graph)
}



