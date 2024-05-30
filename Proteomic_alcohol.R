

entera_alcohol <- read_excel("datos/230224-LFQ_incl_isoforms_FILTRADAS.xlsx")
entera_alcohol$Description <- gsub("\\[OS=Rattus norvegicus\\]", "", entera_alcohol$Description)
sample_metadata <- read.csv("datos/sample-metadata.csv")

df_alcohol<-entera_alcohol

accesion_alcohol<-dplyr::select(df_alcohol, Accession)
gene_ID_alcohol<-dplyr::select(df_alcohol, "Gene Symbol")
indice_proteina<-cbind(accesion_alcohol,gene_ID_alcohol)
colnam<-colnames(entera_alcohol)
# indices <- grep("Normalized)", colnam)
indices <- grep("Abundance: ", colnam)

abundancias_experimento<-df_alcohol[, indices]
id<-indice_proteina[,1]

colnam<-colnames(abundancias_experimento)
colnam <- gsub(".*rat(\\d+).*", "\\1", colnam)
colnam <- paste("rat-", gsub("\\D+", "", colnam), sep = "")
colnames(abundancias_experimento) <- colnam

df_alcohol<-cbind(accesion_alcohol, gene_ID_alcohol, abundancias_experimento)

df_alcohol <- df_alcohol %>% 
  arrange(Accession)

entera_alcohol <- entera_alcohol %>% 
  arrange(Accession)

#Grupos experimentales

#control

prot_control<- sample_metadata %>% 
  filter(alcohol == "abstemio" & stress == "relax") %>% 
  pull("sample.id")

tabla_control<- df_alcohol %>% 
  select(1, 2, all_of(prot_control))



# alcohol tranqui
prot_alcohol_relax<- sample_metadata %>% 
  filter(alcohol == "EtOH" & stress == "relax") %>% 
  pull("sample.id")

tabla_alcohol_relax<- df_alcohol %>% 
  select(1, 2, all_of(prot_alcohol_relax))


tabla_control<-tabla_control %>% select(starts_with("rat"))
tabla_alcohol_relax<-tabla_alcohol_relax %>% select(starts_with("rat"))
nprot<-nrow(tabla_alcohol_relax) 
nratas<-ncol(tabla_alcohol_relax) 


#Crear una lista para almacenar las líneas que no cumplen con el mínimo de 3 datos no vacíos
falta_control <- c()

# Identificar las líneas que no cumplen con el mínimo de 3 datos no vacíos en tabla_alcohol_relax
for (n in 1:nprot) {
  c_abstemio_relax <- as.numeric(tabla_control[n, ])
  
  # Verificar si hay menos de 3 datos no vacíos
  if (sum(!is.na(c_abstemio_relax)) < 3) {
    falta_control <- c(falta_control, n)
  }
}



# Crear una lista para almacenar las líneas que no cumplen con el mínimo de 3 datos no vacíos
falta_alcohol <- c()

# Identificar las líneas que no cumplen con el mínimo de 3 datos no vacíos en tabla_alcohol_relax
for (n in 1:nprot) {
  c_alcohol_relax <- as.numeric(tabla_alcohol_relax[n, ])
  
  # Verificar si hay menos de 3 datos no vacíos
  if (sum(!is.na(c_alcohol_relax)) < 3) {
    falta_alcohol <- c(falta_alcohol, n)
  }
}

valores_comunes <- union(falta_alcohol, falta_control)

######################################################################################

entera_alcohol<-entera_alcohol[-valores_comunes,]
df_alcohol<-entera_alcohol



accesion_alcohol<-dplyr::select(df_alcohol, Accession)
gene_ID_alcohol<-dplyr::select(df_alcohol, "Gene Symbol")
indice_proteina<-cbind(accesion_alcohol,gene_ID_alcohol)
colnam<-colnames(entera_alcohol)
# indices <- grep("Normalized)", colnam)
indices <- grep("Abundance: ", colnam)
# indices <- grep("Scaled)", colnam)

abundancias_experimento<-df_alcohol[, indices]
id<-indice_proteina[,1]

colnam<-colnames(abundancias_experimento)
colnam <- gsub(".*rat(\\d+).*", "\\1", colnam)
colnam <- paste("rat-", gsub("\\D+", "", colnam), sep = "")
colnames(abundancias_experimento) <- colnam

# abundancias_experimento[is.na(abundancias_experimento)] <- 0

# df_alcohol<-cbind(accesion_alcohol, gene_ID_alcohol, abundancias_experimento)
rownames(abundancias_experimento)<-id
rownames(abundancias_experimento)<-id



prot_control<- sample_metadata %>% 
  filter(alcohol == "abstemio" & stress == "relax") %>% 
  pull("sample.id")


prot_alcohol_relax<- sample_metadata %>% 
  filter(alcohol == "EtOH" & stress == "relax") %>% 
  pull("sample.id")


tabla_control<- abundancias_experimento %>% 
  select(all_of(prot_control))

tabla_alcohol_relax<- abundancias_experimento %>% 
  select(all_of(prot_alcohol_relax))

alcohol<-cbind(tabla_control, tabla_alcohol_relax)



suma<-colSums(alcohol, na.rm = TRUE)


suma_total<-sum(suma)
media_total<-sum(suma)/ncol(alcohol)
factor_correlaccion<-media_total/suma



factor_correlaccion_rep <- lapply(factor_correlaccion, rep, times = nrow(alcohol))
factor_correlaccion_rep <- as.data.frame(do.call(cbind, factor_correlaccion_rep))
ajustado<-factor_correlaccion_rep*alcohol

suma_2<-colSums(ajustado, na.rm = TRUE)



prot_control<- sample_metadata %>% 
  filter(alcohol == "abstemio" & stress == "relax") %>% 
  pull("sample.id")


prot_alcohol_relax<- sample_metadata %>% 
  filter(alcohol == "EtOH" & stress == "relax") %>% 
  pull("sample.id")


tabla_control<- ajustado %>% 
  select(all_of(prot_control))

tabla_alcohol_relax<- ajustado %>% 
  select(all_of(prot_alcohol_relax))

rownames(tabla_control)<-id
rownames(tabla_alcohol_relax)<-id





n<-ncol(tabla_control)
tabla_control_datos<-tabla_control
tabla_control_datos$media<-(rowSums(tabla_control / n, na.rm = TRUE)) 

tabla_alcohol_relax_datos<-tabla_alcohol_relax
tabla_alcohol_relax_datos$media <- rowSums(tabla_alcohol_relax, na.rm = TRUE) / n

nprot<-nrow(tabla_alcohol_relax) 
nratas<-ncol(tabla_alcohol_relax) 


abundancias_alcohol_normalizada<-cbind(tabla_control,tabla_alcohol_relax)
rownames(abundancias_alcohol_normalizada)<-id

p_alcohol <- c()
normalidad_alcohol <- c()
normalidad_control <- c()
cohen_values <- c() 

# Realizar el test y comprobación de normalidad para cada proteína
for (n in 1:nprot) {
  c_alcohol_relax <- as.numeric(tabla_alcohol_relax[n, ])
  c_abstemio_relax <- as.numeric(tabla_control[n, ])
  
  # Calcular la d de Cohen
  resultado_cohen_alcohol <- cohen.d(c_alcohol_relax, c_abstemio_relax)
  cohen <- resultado_cohen_alcohol$estimate
  
  # Almacenar el valor de la d de Cohen
  cohen_values <- c(cohen_values, cohen)
  
  
  # Comprobar normalidad para el grupo "alcohol"
  if (length(unique(c_alcohol_relax)) == 1) {
    p_valor_alcohol <- 1
  } else {
    p_valor_alcohol <- shapiro.test(c_alcohol_relax)$p.value
  }
  normalidad_alcohol <- c(normalidad_alcohol, p_valor_alcohol)
  
  # Comprobar normalidad para el grupo control
  if (length(unique(c_abstemio_relax)) == 1) {
    p_valor_control <- 1
  } else {
    p_valor_control <- shapiro.test(c_abstemio_relax)$p.value
  }
  normalidad_control <- c(normalidad_control, p_valor_control)
  
  # Verificar si los valores de ambos grupos son idénticos
  if (all(c_alcohol_relax == c_abstemio_relax)) {
    p_valor <- 1
  } else {
    # Elegir el test estadístico en función de la normalidad de ambos grupos
    if (p_valor_alcohol > 0.05 && p_valor_control > 0.05) {
      # Ambos grupos son normales, realizar t-test
      resultado_alcohol <- t.test(c_alcohol_relax, c_abstemio_relax, var.equal = FALSE )
      
    } else {
      # Al menos uno de los grupos no es normal, realizar Wilcoxon test
      resultado_alcohol <- wilcox.test(c_alcohol_relax, c_abstemio_relax, var.equal = FALSE)
    }
    
    # Almacenar el p-value
    p_valor <- resultado_alcohol$p.value
    
  }
  
  
  # Almacenar el p-value
  p_alcohol <- c(p_alcohol, p_valor)
}

# Crear un data frame con los resultados
result_alcohol <- data.frame(Protein = rownames(tabla_alcohol_relax), 
                             P_Value = p_alcohol, 
                             Cohen_D = cohen_values, 
                             Normality_P_Value = normalidad_alcohol, 
                             row.names = NULL)



# Filtrar las de proteínas significativas
significativa_alcohol <- result_alcohol %>% filter(P_Value < 0.05)



write.xlsx(significativa_alcohol, ""datos/estadística_alcohol_BR16_machos.xlsx", sheetName = "Desrreguladas")





tabla_alcohol_relax_datos<-tabla_alcohol_relax
tabla_alcohol_relax_datos$media<- rowSums(tabla_alcohol_relax) / nratas

tabla_control_datos<-tabla_control
tabla_control_datos$media<- rowSums(tabla_control) / nratas



cociente_alcohol<-tabla_alcohol_relax_datos$media / tabla_control_datos$media
cociente_alcohol<-data.frame(cociente_alcohol)
log2cociente_alcohol<-apply(cociente_alcohol, 1, log2)
log2cociente_alcohol<-data_frame(log2cociente_alcohol)

p_valor_alcohol<-result_alcohol$P_Value
p_valor_alcohol<-data.frame(p_valor_alcohol)

p_valor_alcohol<-data.frame(p_valor_alcohol)
log10pvalor_alcohol<-apply(p_valor_alcohol, 1, log10)
log10pvalor_alcohol<-data.frame(log10pvalor_alcohol)
log10pvalor_alcohol<-apply(log10pvalor_alcohol, 1, function(x) -x)
log10pvalor_alcohol<-data.frame(log10pvalor_alcohol)

vulcanot<-data.frame(Accession = result_alcohol$Protein,
                     Gene.Symbol = entera_alcohol$`Gene Symbol`,
                     log2cociente_alcohol, 
                     log10pvalor_alcohol, 
                     P_Value = p_alcohol,
                     Cohen_D = cohen_values, 
                     Normality_P_Value = normalidad_alcohol, 
                     tabla_alcohol_relax_datos$media, 
                     tabla_control_datos$media,
                     entera_alcohol$Description)
rata_control<-ncol(tabla_control)
rata_alcohol<-ncol(tabla_alcohol_relax)
p_valor<-1.30103
LFC<-0.58
# LFC<-0.32192809488
a<-  rata_control + rata_alcohol
ref <- 1.35
vulcanot_etiquetas<- vulcanot %>% filter(log2cociente_alcohol < -LFC | log2cociente_alcohol > LFC)
vulcanot_etiquetas<- vulcanot_etiquetas %>% filter(log10pvalor_alcohol>p_valor)



# mas altos sobrexpresados


tops<- vulcanot_etiquetas %>% filter(log2cociente_alcohol>LFC)
tops <- tops %>% arrange(desc(log10pvalor_alcohol)) %>%  # Ordenar el DataFrame por 'log10pvalor_alcohol' en orden descendente
  head(3)                         # Seleccionar las primeras tres filas

top_rows<-tops


# Obtener los valores de la columna 'Accesion' de las tres filas
sobreexpresadas_evidentes_accesion <- top_rows$Accession

# mas altos subrexpresados


high_negative_rows <- vulcanot_etiquetas %>%
  filter(log2cociente_alcohol < -LFC) %>%   # Filtrar filas con 'log2cociente_alcohol' negativo
  arrange(desc(log10pvalor_alcohol)) %>%  # Ordenar el DataFrame resultante por 'log10pvalor_alcohol' en orden descendente
  head(3)                         # Seleccionar lasrownames(abundancias) <- id



# Obtener los valores de la columna 'Accesion' de las tres filas
subrrepresentadas_evidentes_accesion <- high_negative_rows$Accession

vulcanot_etiquetas_evidente<- vulcanot_etiquetas  %>% filter(log10pvalor_alcohol>p_valor)

far_rows <- vulcanot_etiquetas_evidente %>%
  arrange(desc(log2cociente_alcohol)) %>%  # Ordenar el DataFrame por 'log10pvalor_alcohol' en orden descendente
  head(3)                         # Seleccionar las primeras tres filas

# Obtener los valores de la columna 'Accesion' de las tres filas
mas_sobreexpresadas_accesion <- far_rows$Accession



# Obtener los tres valores máximos de 'log10pvalor_alcohol' para filas con 'log2cociente_alcohol' negativo
far_negative_rows <- vulcanot_etiquetas %>%
  filter(log2cociente_alcohol < (-LFC) & log10pvalor_alcohol > p_valor) %>%   # Filtrar filas con 'log2cociente_alcohol' negativo
  arrange(log2cociente_alcohol) %>%  # Ordenar el DataFrame resultante por 'log10pvalor_alcohol' en orden descendente
  head(3)                         # Seleccionar las primeras tres filas

mas_subrrepresentadas_accesion <- far_negative_rows$Accession

accession_values <- c(mas_sobreexpresadas_accesion, mas_subrrepresentadas_accesion, 
                      subrrepresentadas_evidentes_accesion, sobreexpresadas_evidentes_accesion)
accession_values <- data.frame("Accession" = accession_values)
subset_desrreguladas<- subset(vulcanot, Accession %in% accession_values$Accession)

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
sobrexpresadas<-filter(vulcanot, log2cociente_alcohol>LFC & log10pvalor_alcohol>p_valor)
subexpresadas <- filter(vulcanot, log2cociente_alcohol < -LFC & log10pvalor_alcohol > p_valor )

resto<-anti_join(vulcanot, sobrexpresadas) %>% anti_join(subexpresadas)
sobrexpresadas$estado<-"Sobrexpresado"
subexpresadas$estado<-"Subexpresado"


df_desrreguladas<-rbind(subexpresadas, sobrexpresadas)
write.xlsx(df_desrreguladas, ""datos/efecto_alcohol_BR16_machos.xlsx", sheetName = "Desrreguladas")
write.xlsx(vulcanot, ""datos/estadistica_entera_alcohol_BR16_machos.xlsx", sheetName = "todas")




# Calcula el número de puntos en tus datos
numero_de_puntos <- sum(
  nrow(sobrexpresadas),
  nrow(subexpresadas),
  nrow(resto)
)

# Crea el gráfico sin capas geométrimio
tu_grafico <- ggplot() +
  labs(
    title = "Volcano efecto de la proteómica alcohol BR-16 mio",
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
  geom_point(data = sobrexpresadas, aes(x = log2cociente_alcohol, y = log10pvalor_alcohol, color = estado), size = 0.5) +
  geom_point(data = subexpresadas, aes(x = log2cociente_alcohol, y = log10pvalor_alcohol, color = estado), size = 0.5) +
  geom_point(data = resto, aes(x = log2cociente_alcohol, y = log10pvalor_alcohol), color = "grey", size = 0.5) +
  geom_vline(xintercept = -LFC, linetype = "dashed", color = "black") +
  geom_vline(xintercept = LFC, linetype = "dashed", color = "black") +
  geom_hline(yintercept = p_valor, linetype = "dashed", color = "black") +
  geom_label_repel(data = subset_desrreguladas,
                   aes(x = log2cociente_alcohol, y = log10pvalor_alcohol, label = `Gene.Symbol`),
                   color = "black", box.padding = 0.2, point.padding = 0.1,
                   force = 10, nudge_x = 0.000000000000001, nudge_y = 0.0000000000001)

# Agrega el texto con el número de puntos
tu_grafico <- tu_grafico +
  geom_text(aes(x = 2, y = ref - 0.3, label = paste("p-valor(0.05)= 1.30103")), size = 3, hjust = 0)+
  geom_text(aes(x = 2, y = ref - 0.55, label = paste("N (de proteínas=", numero_de_puntos)), size = 3, hjust = 0) +
  geom_text(aes(x = 2, y = ref - 0.80, label = paste("N (ratas)  =", a)), size = 3, hjust = 0) +
  geom_text(aes(x = 2, y = ref - 1.05 , label = paste("subexpresadas =", nrow(subexpresadas))), size = 3, hjust = 0) +
  geom_text(aes(x = 2, y = ref - 1.30 , label = paste("sobrexpresadas =", nrow(sobrexpresadas))), size = 3, hjust = 0) +
  geom_text(aes(x = -5, y = 5, label = paste("Ratas consumo alcohol vs control")), size = 3, hjust = 0)

print(tu_grafico)  
ruta_archivo <- ""datos/efecto_alcohol_BR16_machos.png"
ggsave(filename = ruta_archivo, plot = last_plot(), device = "png", width = 7.5, height = 7.5)



lista_red <- paste(df_desrreguladas$Accession, collapse = ", ")

lista_red <- strsplit(lista_red, ", ")[[1]]

# Agregar comillas a cada elemento
lista_red <- paste0('"', lista_red, '"')

# Cadena formateada
lista_red <- paste("[", paste(lista_red, collapse = ", "), "]", sep = "")

cat(lista_red)
#Aquí obtenemos un fichero RED_ALCOHOL.XLSX utilizando estos códigos para obtener un fichero de interacción de STRING

br_16_alcohol <- read_excel(""datos/RED_ALCOHOL.xlsx")



br_16_interacciones <- br_16_alcohol[, (3:4)]

br_16_interacciones <- replace(br_16_interacciones, br_16_interacciones == "Rps10l1", "Rps10")

br_16_interacciones <- replace(br_16_interacciones, br_16_interacciones == "LOC102555453", "Rpl12")

br_16_interacciones <- replace(br_16_interacciones, br_16_interacciones == "Qk", "Qki")


# # Plot ----
# g <- graph_from_data_frame(br_16_interacciones, directed = FALSE) %>% simplify
# plot.igraph(g, vertex.label.dist = 3.5, layout = layout_nicely, edge.arrow.size = 0.5)

nodos <- unique(br_16_interacciones[, 1])

# Supongo que aquí estás utilizando el dataframe correcto (df_desrreguladas)
proteinas_all <-
  tibble(
    Gene.Symbol = df_desrreguladas$Gene.Symbol,
    P_Value = df_desrreguladas$log10pvalor_alcohol,
    Cociente = df_desrreguladas$log2cociente_alcohol  ,
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
  
  ggtitle("Red interacción de proteínas desreguladas por el consumo de alcohol (mio)")


ruta_archivo <- ""datos/red_alcohol.png"
ggsave(filename = ruta_archivo, plot = last_plot(), device = "png", width = 10, height = 10)









sample_metadata <- read.csv("alcohol/sample-metadata.2.csv")


prot_control<- sample_metadata %>% 
  filter(alcohol == "abstemio") %>% 
  pull("sample.id")

tabla_control_relleno<- abundancias_alcohol_normalizada %>% 
  select(all_of(prot_control))



# Calcular la media de cada fila excluyendo los valores vacíos
tabla_control_relleno$media <- rowMeans(tabla_control_relleno, na.rm = TRUE)

# Reemplazar los valores vacíos en cada fila con la media calculada
for (i in 1:nrow(tabla_control_relleno)) {
  tabla_control_relleno[i, is.na(tabla_control_relleno[i, ])] <- tabla_control_relleno$media[i]
}

# Eliminar la columna de medias si no la necesitas
tabla_control_relleno$media <- NULL


# alcohol tranqui
prot_alcohol_relax<- sample_metadata %>% 
  filter(alcohol == "EtOH") %>% 
  pull("sample.id")

tabla_alcohol_relax_relleno<- abundancias_alcohol_normalizada %>% 
  select(all_of(prot_alcohol_relax))


# Calcular la media de cada fila excluyendo los valores vacíos
tabla_alcohol_relax_relleno$media <- rowMeans(tabla_alcohol_relax_relleno, na.rm = TRUE)

# Reemplazar los valores vacíos en cada fila con la media calculada
for (i in 1:nrow(tabla_alcohol_relax_relleno)) {
  tabla_alcohol_relax_relleno[i, is.na(tabla_alcohol_relax_relleno[i, ])] <- tabla_alcohol_relax_relleno$media[i]
}

# Eliminar la columna de medias si no la necesitas
tabla_alcohol_relax_relleno$media <- NULL




abundancias_alcohol_relleno<-cbind(tabla_alcohol_relax_relleno, tabla_control_relleno)

rownames(abundancias_alcohol_relleno)<-id

suma<-colSums(abundancias_alcohol_relleno, na.rm = TRUE)
print(suma)
max<-which.max(suma)
figura<-barplot(suma, cex.names = 0.6, las = 2, xlab = "Muestras", ylab = "Expresión")
title("Diagrama de barras del total de lecturas por muestras normalizado con HK1" )
# abline(h=median(colSums(abundancias_alcohol_relleno)),col="blue")


n<-(abundancias_alcohol_relleno)

medias_gen<-colSums(abundancias_alcohol_relleno, na.rm = TRUE)

cpm<-cpm(abundancias_alcohol_relleno)
logcpm <- cpm(abundancias_alcohol_relleno, log=TRUE)
plotDensities(logcpm, legend = "topright")




boxplot(logcpm, xlab="", ylab="Recuento en Log2 por millón de lecturas",las=2)
abline(h=median(logcpm),col="blue")
title("Boxplots de los log(CPM)")



sample_metadata$alcohol <- factor(sample_metadata$alcohol, levels = c("abstemio", "EtOH"))
# asociamos un color a cada tipo de muestra
sample.color <- c("blue","red")[sample_metadata$alcohol]
# usamos ese código de colores para pintar el MDSplot
plotMDS(logcpm, col=sample.color)
# añadimos leyenda y título
legend("bottomright", fill=c("blue", "red"), legend=levels(sample_metadata$alcohol))
title("MDSplot con muestras coloreadas según estado")




log_abundancias_alcohol_relleno<-apply(abundancias_alcohol_relleno, 1, log10 )
log_abundancias_alcohol_relleno<-t(log_abundancias_alcohol_relleno)
seleccionadas<-names(sort(varianza, decreasing = TRUE))[1:50]
desrreguladas<-df_desrreguladas$Accession

matriz_elegidas<-log_abundancias_alcohol_relleno[desrreguladas,]
desrreguladas_nombres <- indice_proteina[match(rownames(matriz_elegidas), indice_proteina$Accession),]
rownames(matriz_elegidas) <- desrreguladas_nombres$`Gene Symbol`

rownames(sample_metadata) <- sample_metadata$sample.id
pheatmap(matriz_elegidas,  annotation_col = sample_metadata)
