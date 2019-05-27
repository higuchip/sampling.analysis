#====================================================================================================================================
# sampling.analysis                           
# Função para análise do processo de amostragem em levantamentos fitossociológicos, em função do número de indivíduos e da área basal  
# Autor: Pedro Higuchi, 26/05/2019
# Como citar:Higuchi, P. sampling.analysis: Função em linguagem de programação estatística R para análise do processo amostragem de 
# levantamentos fitossociógicos em função do número de indivíduos e da área basal. 2019. #Disponvel em https://github.com/higuchip/sampling.analysis
#-----------------------------------------------------------------------------------------------------------------------------------
#										                      
# Observações:											                      
# - a) O argumento x (planilha de dados) terá que conter as colunas parc (identificação das parcelas), 
#   spp (id. espécies),dap e estratos (obrigatório apenas no caso de amostragem estratificada);
# - b) arquivo exemplo de entrada, disponível em 
#   https://raw.githubusercontent.com/higuchip/sampling.analysis/master/dados_exemplo_amostragem.csv
# - c) O argumento sys, representa o sistema de amostragem (AS = Amostragem Simples; EST = Estratificada e SIS = Sistemática);
# - d) O argumento plot_size representa o tamanho de cada parcela em m2;
# - e) O argumento forest_area representa o tamanho da área florestal em ha;
# - f) O argumento strata_area deve ser usado apenas quando sys = EST e representa um vetor numérico com as áreas de cada estrato em ha;
# - g) o argumento alfa representa o nível de significância pela tabela t de Student;
# - h) o argumento LE representa o Limite de erro admissível.
#---------------------------------------------------------------------------------------------------------------------------------------
#
#Aplicação
#
#dados_exemplo_amostragem <- read.table("https://raw.githubusercontent.com/higuchip/sampling.analysis/master/dados_exemplo_amostragem.csv", header=T, sep = ";", dec=",")
#source("https://raw.githubusercontent.com/higuchip/sampling.analysis/master/sampling.analysis.R")
#sampling.analysis(dados_exemplo_amostragem,sys = "AS", plot_size = 200,  forest_area = 150, alfa = 0.05, LE = 0.1)
#
# Modificações
# Data: 27/05/2019
# Add: Adicionado a possibilidade de troncos múltiplos para dap
#====================================================================================================================================

sampling.analysis<-function(x, sys, plot_size, forest_area, strata_area=NULL, alfa, LE){
  
  x[is.na(x)] <- 0
  
  matrix.data<-table(x$parc, x$spp)
  n<-dim(matrix.data)[1]  
  sample_size<-(plot_size*n)/10000
  (f <-sample_size/forest_area) 
  
 ndaps<-grep('dap', colnames(x))
   length(ndaps)
  

   if (length(ndaps) > 1){
  
     dbhs <-x[,c(grep('dap', colnames(x)))]
     dbhs
     sas<-(pi*dbhs^2)/40000
     sas
     x$SA<-apply(sas, 1, sum)
  
   } else {
    
  x$SA<-(pi*x$dap^2)/40000 #determinacao de AS para as árvores
   }
  
  
  
  abund<-apply(matrix.data,1,sum) 
  BA<-tapply(x$SA, x$parc, sum) 
  sub_veg <- subset(x, select = c(estratos,parc)) 
  matrix_ind<-table(sub_veg$parc, sub_veg$estratos) 
  matrix_ind<-as.data.frame(as.table(matrix_ind))
  matrix_ind<-matrix_ind[matrix_ind$Freq>0,]
  strata<-matrix_ind$Var2
  structure<-as.data.frame(cbind(strata,abund, BA))
  
  
  
  N<-forest_area/(plot_size/10000)
  ttab<-qt(1-alfa,df=n-1) 
  
  
  if(sys=="AS") {
    
    
    mean.structure<-apply(structure[,2:3],2,mean)
    
    E<-LE* mean.structure #Erro adimitido absoluto
    
    var.structure<-apply(structure[,2:3],2,var) #variancia
    sd.structure<-apply(structure[,2:3],2,sd) #desvio-padrao
    
    
    sde.structure.inf <- sd.structure/sqrt(n)*sqrt((1-f))
    sde.structure.fin <- sd.structure/sqrt(n)
    
    absulute.error.inf<-ttab*sde.structure.inf
    relative.error.inf<-((absulute.error.inf/mean.structure))*100
    
    absulute.error.fin<-(ttab*sde.structure.fin)
    relative.error.fin<-((absulute.error.fin/mean.structure))*100
    
    if(f>0.02) {
      n.min.fin<-(N*(ttab^2)*var.structure)/(N* (E^2)+((ttab^2)*var.structure))
      
      cat("--------------------------------------------")
      cat("\nERRO DE AMOSTRAGEM RELATIVO", fill=TRUE)
      cat("Abundância", round(relative.error.fin[1],1), "%", fill=TRUE )
      cat("Área basal", round(relative.error.fin[2],1), "%", fill=TRUE )
      cat("--------------------------------------------")
      cat("\nNÚMERO NECESSÁRIO DE PARCELAS", fill=TRUE)
      cat("Abundância:",ceiling(n.min.fin[[1]]),fill=TRUE)
      cat("Área Basal:",ceiling(n.min.fin[[2]]),fill=TRUE)
      cat("--------------------------------------------")
      
    } else {
      n.min.inf<-(ttab^2*var.structure)/E^2
      cat("--------------------------------------------")
      cat("\nERRO DE AMOSTRAGEM RELATIVO", fill=TRUE)
      cat("Abundância", round(relative.error.inf[1],1), "%", fill=TRUE )
      cat("Área basal", round(relative.error.inf[2],1), "%", fill=TRUE )
      cat("--------------------------------------------")
      cat("\nNÚMERO NECESSÁRIO DE PARCELAS", fill=TRUE)
      cat("Abundância:",ceiling(n.min.inf[[1]]),fill=TRUE)
      cat("Área Basal:",ceiling(n.min.inf[[2]]),fill=TRUE)
      cat("--------------------------------------------")
      
    }
    
    
  } else if (sys=="EST") {
    
    Wh<-strata_area/forest_area
    wh<-table(structure$strata)/n
    
    
    mean_strata<-aggregate(structure[,2:3], list(structure$strata), mean)[,2:3]
    var_strata<-aggregate(structure[,2:3], list(structure$strata), var) [,2:3]
    
    stratified_mean<-apply((mean_strata[,1:2]*Wh),2,sum) #media estratificada
    stratified_var<-apply((var_strata[,1:2]*Wh),2,sum) #variancia estratificada
    
    
    Wh2_x_var_by_n<-((Wh^2)*(var_strata/table(structure$strata)))
    sum_Wh2_x_var_by_n<-apply(Wh2_x_var_by_n, 2,sum)
    
    
    Wh_x_var_by_N<-((Wh*var_strata)/N)
    sum_Wh_x_var_by_N<-apply(Wh_x_var_by_N, 2,sum)
    
    
    var_mean_stratied<-sum_Wh2_x_var_by_n-sum_Wh_x_var_by_N # Variancia da media estratificada
    
    
    
    erro.padrao.media.estratificada.inf<-sqrt(var_mean_stratied)
    erro.padrao.media.estratificada.fin<-sqrt(var_mean_stratied)*sqrt((1-f))
    
    sum_Wh_x_sh<-apply(Wh*var_strata,2,sum) #somatorio de variancia x Wh
    
    sum_Wh_x_sh_N<-apply((Wh*var_strata)/N,2,sum) #somatorio (variancia x Wh)/N
    
    mean.structure<-apply(structure[,2:3],2,mean)
    
    E<-LE* stratified_mean #Erro adimitido absoluto
    
    
    absulute.error.inf.est<-ttab*erro.padrao.media.estratificada.inf
    relative.error.inf.est<-(absulute.error.inf.est/stratified_mean)*100
    
    absulute.error.fin.est<-ttab*erro.padrao.media.estratificada.fin
    relative.error.fin.est<-(absulute.error.fin.est/stratified_mean)*100
    
    
    if(f>0.02) {
      n.min.fin.est<-((ttab^2)* sum_Wh_x_sh)/((E^2)+((ttab^2)*sum_Wh_x_sh_N))
      cat("--------------------------------------------")
      cat("\nERRO DE AMOSTRAGEM RELATIVO", fill=TRUE)
      cat("Abundância:",round(relative.error.fin.est[[1]],1),"%", fill=TRUE)
      cat("Área Basal:",round(relative.error.fin.est[[2]],1),"%",fill=TRUE)
      cat("--------------------------------------------")
      cat("\nNÚMERO NECESSÁRIO DE PARCELAS", fill=TRUE)
      cat("Abundância:",ceiling(n.min.fin.est[[1]]),fill=TRUE)
      cat("Área Basal:",ceiling(n.min.fin.est[[2]]),fill=TRUE)
      cat("--------------------------------------------")
      
    } else {
      n.min.inf.est<-(ttab^2*sum_Wh_x_sh)/E^2
      cat("--------------------------------------------")
      cat("\nERRO DE AMOSTRAGEM RELATIVO", fill=TRUE)
      cat("Abundância:",round(relative.error.inf.est[[1]],1),"%", fill=TRUE)
      cat("Área Basal:",round(relative.error.inf.est[[2]],1),"%",fill=TRUE)
      cat("--------------------------------------------")
      cat("\nNÚMERO NECESSÁRIO DE PARCELAS", fill=TRUE)
      cat("Abundância:",ceiling(n.min.inf.est[[1]]),fill=TRUE)
      cat("Área Basal:",ceiling(n.min.inf.est[[2]]),fill=TRUE)
      
    }  
    
  }else if (sys=="SIS") {
    
    mean.structure<-apply(structure[,2:3],2,mean)
    
    E.abund<-LE* mean.structure[1] #Erro adimitido absoluto
    E.ba<-LE* mean.structure[2]
    
    # var.structure<-apply(structure[,2:3],2,var) #variancia
    # sd.structure<-apply(structure[,2:3],2,sd) #desvio-padrao
    # 
    
    sde.structure.abund.sis <- sqrt(sum(diff(structure$abund, 1)^2)/(2*n*(n-1))*((N-n)/N)) #erro padrao da media por meio do metodo das diferencas sucessivas
    sde.structure.ba.sis <- sqrt(sum(diff(structure$BA, 1)^2)/(2*n*(n-1))*((N-n)/N)) #erro padrao da media por meio do metodo das diferencas sucessivas
    
    
    absulute.error.abund.sis<-ttab*sde.structure.abund.sis
    relative.error.abund.sis<-((absulute.error.abund.sis/mean.structure[1]))*100
    
    absulute.error.ba.sis<-ttab*sde.structure.ba.sis
    relative.error.ba.sis<-((absulute.error.ba.sis/mean.structure[2]))*100
    
    
    var.structure.abund.sis<-(sde.structure.abund.sis*sqrt(n))^2
    var.structure.abund.sis<-(sde.structure.abund.sis*sqrt(n))^2
    
    var.structure.ba.sis<-(sde.structure.ba.sis*sqrt(n))^2
    var.structure.ba.sis<-(sde.structure.ba.sis*sqrt(n))^2
    
    n.min.abund.sis<-(ttab^2*var.structure.abund.sis)/E.abund^2
    n.min.ba.sis<-(ttab^2*var.structure.ba.sis)/E.ba^2
    
    cat("--------------------------------------------")
    cat("\nERRO DE AMOSTRAGEM RELATIVO", fill=TRUE)
    cat("Abundância:",round(relative.error.abund.sis,1),"%", fill=TRUE)
    cat("Área Basal:",round(relative.error.ba.sis,1),"%",fill=TRUE)
    cat("--------------------------------------------")
    cat("\nNÚMERO NECESSÁRIO DE PARCELAS", fill=TRUE)
    cat("Abundância:",ceiling(n.min.abund.sis),fill=TRUE)
    cat("Área Basal:",ceiling(n.min.ba.sis),fill=TRUE)
    cat("--------------------------------------------")
    
    
  }
}
