
entera_experimento <- read_excel("datos/230224-LFQ_incl_isoforms_FILTRADAS.xlsx")
entera_experimento$Description <- gsub("\\[OS=Rattus norvegicus\\]", "", entera_experimento$Description)
sample_metadata <- read.csv("datos/sample-metadata.csv")

df_experimento<-entera_experimento

accesion_experimento<-dplyr::select(df_experimento, Accession)
gene_ID_experimento<-dplyr::select(df_experimento, "Gene Symbol")
indice_proteina<-cbind(accesion_experimento,gene_ID_experimento)
colnam<-colnames(entera_experimento)
# indices <- grep("Normalized)", colnam)
indices <- grep("Abundance: ", colnam)

abundancias_experimento<-df_experimento[, indices]
id<-indice_proteina[,1]

colnam<-colnames(abundancias_experimento)
colnam <- gsub(".*rat(\\d+).*", "\\1", colnam)
colnam <- paste("rat-", gsub("\\D+", "", colnam), sep = "")
colnames(abundancias_experimento) <- colnam

df_experimento<-cbind(accesion_experimento, gene_ID_experimento, abundancias_experimento)

df_experimento <- df_experimento %>% 
  arrange(Accession)

entera_experimento <- entera_experimento %>% 
  arrange(Accession)

#Grupos experimentales

#control

prot_control<- sample_metadata %>% 
  filter(alcohol == "abstemio" & stress == "relax") %>% 
  pull("sample.id")

tabla_control<- df_experimento %>% 
  select(1, 2, all_of(prot_control))



# alcohol tranqui
prot_abstemio_stress<- sample_metadata %>% 
  filter(alcohol == "abstemio" & stress == "stress") %>% 
  pull("sample.id")

tabla_abstemio_stress<- df_experimento %>% 
  select(1, 2, all_of(prot_abstemio_stress))


tabla_control<-tabla_control %>% select(starts_with("rat"))
tabla_abstemio_stress<-tabla_abstemio_stress %>% select(starts_with("rat"))
nprot<-nrow(tabla_abstemio_stress) 
nratas<-ncol(tabla_abstemio_stress) 


#Crear una lista para almacenar las líneas que no cumplen con el mínimo de 3 datos no vacíos
falta_control <- c()

# Identificar las líneas que no cumplen con el mínimo de 3 datos no vacíos en tabla_abstemio_stress
for (n in 1:nprot) {
  c_abstemio_relax <- as.numeric(tabla_control[n, ])
  
  # Verificar si hay menos de 3 datos no vacíos
  if (sum(!is.na(c_abstemio_relax)) < 3) {
    falta_control <- c(falta_control, n)
  }
}



# Crear una lista para almacenar las líneas que no cumplen con el mínimo de 3 datos no vacíos
falta_stress <- c()

# Identificar las líneas que no cumplen con el mínimo de 3 datos no vacíos en tabla_abstemio_stress
for (n in 1:nprot) {
  c_abstemio_stress <- as.numeric(tabla_abstemio_stress[n, ])
  
  # Verificar si hay menos de 3 datos no vacíos
  if (sum(!is.na(c_abstemio_stress)) < 3) {
    falta_stress <- c(falta_stress, n)
  }
}

valores_comunes <- union(falta_stress, falta_control)

######################################################################################


######################################################################################

entera_stress<-entera_experimento[-valores_comunes,]



accesion_stress<-dplyr::select(entera_stress, Accession)
gene_ID_stress<-dplyr::select(entera_stress, "Gene Symbol")
indice_proteina<-cbind(accesion_stress,gene_ID_stress)
colnam<-colnames(entera_stress)
# indices <- grep("Normalized)", colnam)
indices <- grep("Abundance: ", colnam)
# indices <- grep("Scaled)", colnam)

abundancias_experimento<-entera_stress[, indices]
id<-indice_proteina[,1]

colnam<-colnames(abundancias_experimento)
colnam <- gsub(".*rat(\\d+).*", "\\1", colnam)
colnam <- paste("rat-", gsub("\\D+", "", colnam), sep = "")
colnames(abundancias_experimento) <- colnam

# abundancias_experimento[is.na(abundancias_experimento)] <- 0

# df_stress<-cbind(accesion_stress, gene_ID_stress, abundancias_experimento)
rownames(abundancias_experimento)<-id
rownames(abundancias_experimento)<-id



prot_control<- sample_metadata %>% 
  filter(alcohol == "abstemio" & stress == "relax") %>% 
  pull("sample.id")


prot_abstemio_stress<- sample_metadata %>% 
  filter(alcohol == "abstemio" & stress == "stress") %>% 
  pull("sample.id")


tabla_control<- abundancias_experimento %>% 
  select(all_of(prot_control))

tabla_abstemio_stress<- abundancias_experimento %>% 
  select(all_of(prot_abstemio_stress))

stress<-cbind(tabla_control, tabla_abstemio_stress)



suma<-colSums(stress, na.rm = TRUE)


suma_total<-sum(suma)
media_total<-sum(suma)/ncol(stress)
factor_correlaccion<-media_total/suma



factor_correlaccion_rep <- lapply(factor_correlaccion, rep, times = nrow(stress))
factor_correlaccion_rep <- as.data.frame(do.call(cbind, factor_correlaccion_rep))
ajustado<-factor_correlaccion_rep*stress

suma_2<-colSums(ajustado, na.rm = TRUE)



prot_control<- sample_metadata %>% 
  filter(alcohol == "abstemio" & stress == "relax") %>% 
  pull("sample.id")


prot_abstemio_stress<- sample_metadata %>% 
  filter(alcohol == "abstemio" & stress == "stress") %>% 
  pull("sample.id")


tabla_control<- ajustado %>% 
  select(all_of(prot_control))

tabla_abstemio_stress<- ajustado %>% 
  select(all_of(prot_abstemio_stress))

rownames(tabla_control)<-id
rownames(tabla_abstemio_stress)<-id





n<-ncol(tabla_control)
tabla_control_datos<-tabla_control
tabla_control_datos$media<-(rowSums(tabla_control / n, na.rm = TRUE)) 

tabla_abstemio_stress_datos<-tabla_abstemio_stress
tabla_abstemio_stress_datos$media <- rowSums(tabla_abstemio_stress, na.rm = TRUE) / n

nprot<-nrow(tabla_abstemio_stress) 
nratas<-ncol(tabla_abstemio_stress) 


abundancias_stress_normalizada<-cbind(tabla_control,tabla_abstemio_stress)
rownames(abundancias_stress_normalizada)<-id

p_stress <- c()
normalidad_stress <- c()
normalidad_control <- c()
cohen_values <- c() 

# Realizar el test y comprobación de normalidad para cada proteína
for (n in 1:nprot) {
  c_abstemio_stress <- as.numeric(tabla_abstemio_stress[n, ])
  c_abstemio_relax <- as.numeric(tabla_control[n, ])
  
  # Calcular la d de Cohen
  resultado_cohen_stress <- cohen.d(c_abstemio_stress, c_abstemio_relax)
  cohen <- resultado_cohen_stress$estimate
  
  # Almacenar el valor de la d de Cohen
  cohen_values <- c(cohen_values, cohen)
  
  
  # Comprobar normalidad para el grupo "stress"
  if (length(unique(c_abstemio_stress)) == 1) {
    p_valor_stress <- 1
  } else {
    p_valor_stress <- shapiro.test(c_abstemio_stress)$p.value
  }
  normalidad_stress <- c(normalidad_stress, p_valor_stress)
  
  # Comprobar normalidad para el grupo control
  if (length(unique(c_abstemio_relax)) == 1) {
    p_valor_control <- 1
  } else {
    p_valor_control <- shapiro.test(c_abstemio_relax)$p.value
  }
  normalidad_control <- c(normalidad_control, p_valor_control)
  
  # Verificar si los valores de ambos grupos son idénticos
  if (all(c_abstemio_stress == c_abstemio_relax)) {
    p_valor <- 1
  } else {
    # Elegir el test estadístico en función de la normalidad de ambos grupos
    if (p_valor_stress > 0.05 && p_valor_control > 0.05) {
      # Ambos grupos son normales, realizar t-test
      resultado_stress <- t.test(c_abstemio_stress, c_abstemio_relax, var.equal = FALSE )
      
    } else {
      # Al menos uno de los grupos no es normal, realizar Wilcoxon test
      resultado_stress <- wilcox.test(c_abstemio_stress, c_abstemio_relax, var.equal = FALSE)
    }
    
    # Almacenar el p-value
    p_valor <- resultado_stress$p.value
    
  }
  
  
  # Almacenar el p-value
  p_stress <- c(p_stress, p_valor)
}

# Crear un data frame con los resultados
result_stress <- data.frame(Protein = rownames(tabla_abstemio_stress), 
                            P_Value = p_stress, 
                            Cohen_D = cohen_values, 
                            Normality_P_Value = normalidad_stress, 
                            row.names = NULL)



# Filtrar las de proteínas significativas
significativa_stress <- result_stress %>% filter(P_Value < 0.05)



write.xlsx(significativa_stress, "datos/estadística_stress_BR16_machos.xlsx", sheetName = "Desreguladas")





tabla_abstemio_stress_datos<-tabla_abstemio_stress
tabla_abstemio_stress_datos$media<- rowSums(tabla_abstemio_stress) / nratas

tabla_control_datos<-tabla_control
tabla_control_datos$media<- rowSums(tabla_control) / nratas



cociente_stress<-tabla_abstemio_stress_datos$media / tabla_control_datos$media
cociente_stress<-data.frame(cociente_stress)
log2cociente_stress<-apply(cociente_stress, 1, log2)
log2cociente_stress<-data_frame(log2cociente_stress)

p_valor_stress<-result_stress$P_Value
p_valor_stress<-data.frame(p_valor_stress)

p_valor_stress<-data.frame(p_valor_stress)
log10pvalor_stress<-apply(p_valor_stress, 1, log10)
log10pvalor_stress<-data.frame(log10pvalor_stress)
log10pvalor_stress<-apply(log10pvalor_stress, 1, function(x) -x)
log10pvalor_stress<-data.frame(log10pvalor_stress)

vulcanot<-data.frame(Accession = result_stress$Protein,
                     Gene.Symbol = entera_stress$`Gene Symbol`,
                     log2cociente_stress, 
                     log10pvalor_stress, 
                     P_Value = p_stress,
                     Cohen_D = cohen_values, 
                     Normality_P_Value = normalidad_stress, 
                     tabla_abstemio_stress_datos$media, 
                     tabla_control_datos$media,
                     entera_stress$Description)
rata_control<-ncol(tabla_control)
rata_stress<-ncol(tabla_abstemio_stress)
p_valor<-1.30103
LFC<-0.58
# LFC<-0.32192809488
a<-  rata_control + rata_stress
ref <- 1.35
vulcanot_etiquetas<- vulcanot %>% filter(log2cociente_stress < -LFC | log2cociente_stress > LFC)
vulcanot_etiquetas<- vulcanot_etiquetas %>% filter(log10pvalor_stress>p_valor)



# mas altos sobrexpresados


tops<- vulcanot_etiquetas %>% filter(log2cociente_stress>LFC)
tops <- tops %>% arrange(desc(log10pvalor_stress)) %>%  # Ordenar el DataFrame por 'log10pvalor_stress' en orden descendente
  head(3)                         # Seleccionar las primeras tres filas

top_rows<-tops


# Obtener los valores de la columna 'Accesion' de las tres filas
sobreexpresadas_evidentes_accesion <- top_rows$Accession

# mas altos subrexpresados


high_negative_rows <- vulcanot_etiquetas %>%
  filter(log2cociente_stress < -LFC) %>%   # Filtrar filas con 'log2cociente_stress' negativo
  arrange(desc(log10pvalor_stress)) %>%  # Ordenar el DataFrame resultante por 'log10pvalor_stress' en orden descendente
  head(3)                         # Seleccionar lasrownames(abundancias) <- id



# Obtener los valores de la columna 'Accesion' de las tres filas
subrrepresentadas_evidentes_accesion <- high_negative_rows$Accession

vulcanot_etiquetas_evidente<- vulcanot_etiquetas  %>% filter(log10pvalor_stress>p_valor)

far_rows <- vulcanot_etiquetas_evidente %>%
  arrange(desc(log2cociente_stress)) %>%  # Ordenar el DataFrame por 'log10pvalor_stress' en orden descendente
  head(3)                         # Seleccionar las primeras tres filas

# Obtener los valores de la columna 'Accesion' de las tres filas
mas_sobreexpresadas_accesion <- far_rows$Accession



# Obtener los tres valores máximos de 'log10pvalor_stress' para filas con 'log2cociente_stress' negativo
far_negative_rows <- vulcanot_etiquetas %>%
  filter(log2cociente_stress < (-LFC) & log10pvalor_stress > p_valor) %>%   # Filtrar filas con 'log2cociente_stress' negativo
  arrange(log2cociente_stress) %>%  # Ordenar el DataFrame resultante por 'log10pvalor_stress' en orden descendente
  head(3)                         # Seleccionar las primeras tres filas

mas_subrrepresentadas_accesion <- far_negative_rows$Accession

accession_values <- c(mas_sobreexpresadas_accesion, mas_subrrepresentadas_accesion, 
                      subrrepresentadas_evidentes_accesion, sobreexpresadas_evidentes_accesion)
accession_values <- data.frame("Accession" = accession_values)
subset_Desreguladas<- subset(vulcanot, Accession %in% accession_values$Accession)

########################################################################################################################
########################################################################################################################
#######################################################################################################################
########################################################################################################################
########################################################################################################################
#######################################################################################################################
############################################################
############################################################
########################################################################################################################
#######################################################################################################################
sobrexpresadas<-filter(vulcanot, log2cociente_stress>LFC & log10pvalor_stress>p_valor)
subexpresadas <- filter(vulcanot, log2cociente_stress < -LFC & log10pvalor_stress > p_valor )

resto<-anti_join(vulcanot, sobrexpresadas) %>% anti_join(subexpresadas)
sobrexpresadas$estado<-"Sobrexpresado"
subexpresadas$estado<-"Subexpresado"


df_Desreguladas<-rbind(subexpresadas, sobrexpresadas)
write.xlsx(df_Desreguladas, "datos/efecto_stress_BR16_machos.xlsx", sheetName = "Desreguladas")
write.xlsx(vulcanot, "datos/estadistica_entera_stress_BR16_machos.xlsx", sheetName = "todas")




# Calcula el número de puntos en tus datos
numero_de_puntos <- sum(
  nrow(sobrexpresadas),
  nrow(subexpresadas),
  nrow(resto)
)

# Crea el gráfico sin capas geométrimio
tu_grafico <- ggplot() +
  labs(
    title = "Volcano efecto de la proteómica estrés BR-16",
    x = bquote(log[2](Fold_Change)),
    y = bquote(-log[10](pvalue)),
    color = "Estado"
  ) +
  scale_color_manual(values = c("Sobrexpresado" = "lightsalmon", "Subexpresado" = "cadetblue3")) +
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")
  ) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 5))

# Agrega capas geométrimio
tu_grafico <- tu_grafico +
  geom_point(data = sobrexpresadas, aes(x = log2cociente_stress, y = log10pvalor_stress, color = estado), size = 0.5) +
  geom_point(data = subexpresadas, aes(x = log2cociente_stress, y = log10pvalor_stress, color = estado), size = 0.5) +
  geom_point(data = resto, aes(x = log2cociente_stress, y = log10pvalor_stress), color = "grey", size = 0.5) +
  geom_vline(xintercept = -LFC, linetype = "dashed", color = "black") +
  geom_vline(xintercept = LFC, linetype = "dashed", color = "black") +
  geom_hline(yintercept = p_valor, linetype = "dashed", color = "black") +
  geom_label_repel(data = subset_Desreguladas,
                   aes(x = log2cociente_stress, y = log10pvalor_stress, label = `Gene.Symbol`),
                   color = "black", box.padding = 0.2, point.padding = 0.1,
                   force = 10, nudge_x = 0.000000000000001, nudge_y = 0.0000000000001)

# Agrega el texto con el número de puntos
tu_grafico <- tu_grafico +
  geom_text(aes(x = 2, y = ref - 0.3, label = paste("p-valor(0.05)= 1.30103")), size = 3, hjust = 0)+
  geom_text(aes(x = 2, y = ref - 0.55, label = paste("N (de proteínas=", numero_de_puntos)), size = 3, hjust = 0) +
  geom_text(aes(x = 2, y = ref - 0.80, label = paste("N (ratas)  =", a)), size = 3, hjust = 0) +
  geom_text(aes(x = 2, y = ref - 1.05 , label = paste("subexpresadas =", nrow(subexpresadas))), size = 3, hjust = 0) +
  geom_text(aes(x = 2, y = ref - 1.30 , label = paste("sobrexpresadas =", nrow(sobrexpresadas))), size = 3, hjust = 0) +
  geom_text(aes(x = -5, y = 5, label = paste("Ratas consumo stress vs control")), size = 3, hjust = 0)

print(tu_grafico)  
ruta_archivo <- "datos/efecto_stress_BR16_machos.png"
ggsave(filename = ruta_archivo, plot = last_plot(), device = "png", width = 7.5, height = 7.5)



lista_red <- paste(df_Desreguladas$Accession, collapse = ", ")

lista_red <- strsplit(lista_red, ", ")[[1]]

# Agregar comillas a cada elemento
lista_red <- paste0('"', lista_red, '"')

# Cadena formateada
lista_red <- paste("[", paste(lista_red, collapse = ", "), "]", sep = "")

cat(lista_red)


br_16_stress <- read_excel("datos/RED_STRESS.xlsx")



br_16_interacciones <- br_16_stress[, (3:4)]

br_16_interacciones <- replace(br_16_interacciones, br_16_interacciones == "Rps10l1", "Rps10")

br_16_interacciones <- replace(br_16_interacciones, br_16_interacciones == "LOC102555453", "Rpl12")

br_16_interacciones <- replace(br_16_interacciones, br_16_interacciones == "Qk", "Qki")


# # Plot ----
# g <- graph_from_data_frame(br_16_interacciones, directed = FALSE) %>% simplify
# plot.igraph(g, vertex.label.dist = 3.5, layout = layout_nicely, edge.arrow.size = 0.5)

nodos <- unique(br_16_interacciones[, 1])

# Supongo que aquí estás utilizando el dataframe correcto (df_Desreguladas)
proteinas_all <-
  tibble(
    Gene.Symbol = df_Desreguladas$Gene.Symbol,
    P_Value = df_Desreguladas$log10pvalor_stress,
    Cociente = df_Desreguladas$log2cociente_stress  ,
  )

# Get unique genes from both columns
nodos <-
  unique(c(
    as.character(br_16_interacciones$preferredName_A),
    as.character(br_16_interacciones$preferredName_B)
  ))

g <- graph_from_data_frame(d = br_16_interacciones,
                           vertices = proteinas_all$Gene.Symbol,
                           directed = FALSE)

# Asignar datos a los nodos
V(g)$cociente <- proteinas_all$Cociente
V(g)$p_valor <- proteinas_all$P_Value

# Es mejor representar de manera binaria la expresion de los genes
categoria <- ifelse(V(g)$cociente > 0, "Sobreexpresado", "Subexpresado")
V(g)$categoria <- categoria

# Grafo ----
mi_paleta <- colorRampPalette(c("darkcyan", "white", "coral"))(n = 100)

ggraph(g, layout = "igraph", algorithm = "nicely") +
  geom_edge_link() +
  geom_node_point(aes(size = p_valor, fill = cociente, shape = categoria),
                  colour = "black") +
  geom_node_text(aes(label = name), nudge_x = 0, nudge_y = 0.5) + 
  theme_graph() +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_gradientn(colors = mi_paleta, limits = c(min(-5, 0), max(0,5))) +  scale_size(range = c(5, 12)) +  # Cambia el tamaño base de los nodos
  labs(fill = "Cociente",
       shape = "Categoría",
       size = "p-valor")+
  
  ggtitle("Red interacción de proteínas desreguladas por el consumo de estrés")


ruta_archivo <- "datos/red_stress.png"
ggsave(filename = ruta_archivo, plot = last_plot(), device = "png", width = 10, height = 10)









sample_metadata <- read.csv("datos/sample-metadata.3.csv")


prot_control<- sample_metadata %>% 
  filter(stress == "relax") %>% 
  pull("sample.id")

tabla_control_relleno<- abundancias_stress_normalizada %>% 
  select(all_of(prot_control))



# Calcular la media de cada fila excluyendo los valores vacíos
tabla_control_relleno$media <- rowMeans(tabla_control_relleno, na.rm = TRUE)

# Reemplazar los valores vacíos en cada fila con la media calculada
for (i in 1:nrow(tabla_control_relleno)) {
  tabla_control_relleno[i, is.na(tabla_control_relleno[i, ])] <- tabla_control_relleno$media[i]
}

# Eliminar la columna de medias si no la necesitas
tabla_control_relleno$media <- NULL


# stress tranqui
prot_abstemio_stress<- sample_metadata %>% 
  filter(stress == "stress") %>% 
  pull("sample.id")

tabla_abstemio_stress_relleno<- abundancias_stress_normalizada %>% 
  select(all_of(prot_abstemio_stress))


# Calcular la media de cada fila excluyendo los valores vacíos
tabla_abstemio_stress_relleno$media <- rowMeans(tabla_abstemio_stress_relleno, na.rm = TRUE)

# Reemplazar los valores vacíos en cada fila con la media calculada
for (i in 1:nrow(tabla_abstemio_stress_relleno)) {
  tabla_abstemio_stress_relleno[i, is.na(tabla_abstemio_stress_relleno[i, ])] <- tabla_abstemio_stress_relleno$media[i]
}

# Eliminar la columna de medias si no la necesitas
tabla_abstemio_stress_relleno$media <- NULL




abundancias_stress_relleno<-cbind(tabla_abstemio_stress_relleno, tabla_control_relleno)

rownames(abundancias_stress_relleno)<-id

suma<-colSums(abundancias_stress_relleno, na.rm = TRUE)
print(suma)
max<-which.max(suma)
figura<-barplot(suma, cex.names = 0.6, las = 2, xlab = "Muestras", ylab = "Expresión")
title("Diagrama de barras del total de lecturas por muestras normalizado con HK1" )
# abline(h=median(colSums(abundancias_stress_relleno)),col="blue")


n<-(abundancias_stress_relleno)

medias_gen<-colSums(abundancias_stress_relleno, na.rm = TRUE)

cpm<-cpm(abundancias_stress_relleno)
logcpm <- cpm(abundancias_stress_relleno, log=TRUE)
plotDensities(logcpm, legend = "topright")




boxplot(logcpm, xlab="", ylab="Recuento en Log2 por millón de lecturas",las=2)
abline(h=median(logcpm),col="blue")
title("Boxplots de los log(CPM)")



sample_metadata$stress <- factor(sample_metadata$stress, levels = c("relax", "stress"))
# asociamos un color a cada tipo de muestra
sample.color <- c("blue","red")[sample_metadata$stress]
# usamos ese código de colores para pintar el MDSplot
plotMDS(logcpm, col=sample.color)
# añadimos leyenda y título
legend("bottomright", fill=c("blue", "red"), legend=levels(sample_metadata$stress))
title("MDSplot con muestras coloreadas según estado")




varianza<-apply(abundancias_stress_relleno, 1, var )
logabundancias<-apply(abundancias_stress_relleno, 1, log10 )
logabundancias <- t(logabundancias)
# seleccionadas<-names(sort(varianza, decreasing = TRUE))[1:50]
Desreguladas<-df_Desreguladas$Accession

matriz_elegidas<-logabundancias[Desreguladas,]
Desreguladas_nombres <- indice_proteina[match(rownames(matriz_elegidas), indice_proteina$Accession),]
rownames(matriz_elegidas) <- Desreguladas_nombres$`Gene Symbol`

rownames(sample_metadata) <- sample_metadata$sample.id
pheatmap(matriz_elegidas,  annotation_col = sample_metadata)
