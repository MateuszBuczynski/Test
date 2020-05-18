

setwd("C:\\Users\\SilentiumPC\\Desktop")
dane <- read.table("dane_scoring.csv", sep = ';', header = TRUE)

library(pROC)
library(stringr)
library(factorMerger)
head(dane)



# 1) zmienna ma byc podzielona na tyle przedziłow ile daje quantile
# 2) szukaj optmalnego podziału według giniiego
# 3) zwroc zmienna poprwanie zbinowana 




decimalplaces <- function(x){
  if (abs(x-round(x)) > .Machine$double.eps^0.5){
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

dl_osobnika <- function(min_x,max_x){
  ile_po_przecinku <- max(decimalplaces(min_x), decimalplaces(max_x))
  
  ile_wszystkich <- (max_x - min_x)*10^ile_po_przecinku
  
  z = 0
  while(ile_wszystkich >2^z){
    z <- z+1
  }
  return(z)
}

to_real <- function(x, min_x, max_x, dl_osob){
  y <- min_x + strtoi(x, base = 2)*(max_x-min_x)/(2^dl_osob-1)
  return(y)
}




algorytm_genetyczny <- function(zm, flaga, Liczba_osob = 100, Liczba_gen = 300, mutacja = 0.3, krzyzowanie = 0.8){

  ile_arg <- length(unique(quantile(zm, probs = seq(0,1,0.1)))) -2
  
  min_x <- min(zm)
  max_x <- max(zm)
  
  dl_osob <- dl_osobnika(min_x = min_x, max_x = max_x)
  

  
  populacja <- matrix(0,Liczba_osob,dl_osob*ile_arg)
  
  for(i in 1:Liczba_osob){
    populacja[i,] <- sample(c(0,1), dl_osob*ile_arg, replace = TRUE)
    
  }
  
  mut <- mutacja
  Cross <- krzyzowanie
  
  

  ### najlpesi
  
  BESTY <- matrix(0,Liczba_gen,ile_arg*dl_osob)
  
  max_fit <- c()
  
  for(i in 1:Liczba_gen){
    
    #### Sprawdzam wartosc funkcji dla moich argumentow ###
    mic <- c()
    
    
    for(k in 1:Liczba_osob){
      
      rozdzielenie <- matrix(0,ile_arg,dl_osob)
      
      for(j in 1:ile_arg){
        
        rozdzielenie[j,] <- populacja[k,(1+(dl_osob*(j-1))):(dl_osob*j)]
        
      }
      
      roz_char <- apply(rozdzielenie,1,paste0, collapse = '')
      
      xx <- apply(as.matrix(roz_char), 1, to_real, dl_osob = dl_osob, max_x = max_x, min_x = min_x)
      
      xx[which.min(xx)] <- min(zm)
      xx[which.max(xx)] <- max(zm)
      
      zm_x <- cut(zm, breaks = unique( sort(xx)), include.lowest = TRUE) 
      zm_x <-factor(zm_x, ordered = TRUE, levels = levels(zm_x))
      
      mic[k] <- (2*roc(flaga ~ zm_x ,quiet = TRUE)$auc) -1
    }
    
    
    ############## Zapisuje najlepszych gosci w kazdym pokoleniu
    
    
    best <- populacja[which.max(mic), ]
    
    BESTY[i, ] <- best
    
    
    ########## Selekcja  ###########
    zw <- c()
    
    for(l in 1:Liczba_osob){
      x <- sample(1:Liczba_osob, 2)
      if(mic[x[1]] < mic[x[2]]){
        zw[l] <- x[2]
      } else {
        zw[l] <- x[1]
      }
    }
    
    ######## Mutacja #############
    mut_alg <- runif(Liczba_osob)
    
    for(m in 1:Liczba_osob){
      
      if(mut > mut_alg[m]){
        
        ile_mut <- rpois(1,ile_arg)
        
        if(ile_mut > ile_arg*dl_osob){
          
          ile_mut <- ile_arg*dl_osob
          
        } else if(ile_mut == 0){
          
          ile_mut <- 1
          
        }
        
        ktore_miejsce <- sample(1:(ile_arg*dl_osob),ile_mut)
        
        populacja[zw[m],][ktore_miejsce] <- sample(c(0,1),ile_mut, replace = TRUE)
        
      }
    }
    
    Dzieci <- populacja
    
    ###### krzyzowanie #############
    npok <- 0
    pokolenie <- rep(1, ile_arg*dl_osob)
    
    while(npok<Liczba_osob) {
      
      ktore.r <- sample(1:Liczba_osob, 2, replace=FALSE)
      r1 <- Dzieci[ktore.r[1], ]
      r2 <- Dzieci[ktore.r[2], ]
      
      if(runif(1) < Cross){
        
        criss <- sample(1:((ile_arg*dl_osob)-1), 1)
        
        pot1 <- c(r1[1:criss], r2[(criss+1):(ile_arg*dl_osob)])
        pot2 <- c(r2[1:criss], r1[(criss+1):(ile_arg*dl_osob)])
        
        
        pokolenie <- rbind(pokolenie, pot1,pot2)
        
        
      } else {
        pokolenie <- rbind(pokolenie, r1, r2)
      }
      
      npok <- npok+2
    }
    pokolenie <- pokolenie[2:(Liczba_osob+1),]
    pokolenie[1, ] <- best
    populacja <- pokolenie
    
    
    #max_fit[i] <- max(mic)
    
    #if(i %% 100 ==0){
      
     # plot(max_fit)
    #}
    
    
    
  }
  
  return(BESTY)
}

wynik <- algorytm_genetyczny(flaga = flaga, zm = zm)

flaga = dane$SeriousDlqin2yrs[1:10000]
zm = dane$age[1:10000]

### Porównuje 150 zmiennych
gini_org <- matrix(0,15,10)
gini_alg <- matrix(0,15,10)

for(i in 1:15){
  
  if(i ==1){
    dane_alg <- dane[1:10000,]
  } else {
  
    dane_alg <- dane[((i-1)*10000):(10000*i),]
  }

  for(j in 1:10){
    
    wynik <- algorytm_genetyczny(flaga = dane_alg$SeriousDlqin2yrs, zm = dane_alg[,j+1])
    
    ile_arg <- length(unique(quantile(dane_alg[,j+1], probs = seq(0,1,0.1)))) -2
    dl_osob <- dl_osobnika(min_x = min(dane_alg[,j+1]), max_x = max(dane_alg[,j+1]))
    
    
    rozdzielenie <- matrix(0,ile_arg,dl_osob)
    
    for(k in 1:ile_arg){
      
      rozdzielenie[k,] <- wynik[300,(1+(dl_osob*(k-1))):(dl_osob*k)]
      
    }
    
    roz_char <- apply(rozdzielenie,1,paste0, collapse = '')
    
    xx <- apply(as.matrix(roz_char), 1, to_real, dl_osob = dl_osob, max_x = max(dane_alg[,j+1]), min_x = min(dane_alg[,j+1]))
    
    zm_x_org <- cut(dane_alg[,j+1], breaks = unique(c(min(dane_alg[,j+1]), sort(xx), max(dane_alg[,j+1]))), include.lowest = TRUE) 
    zm_x_org <-factor(zm_x_org, ordered = TRUE, levels = levels(zm_x_org))
    
    
    
    
    #Wyliczam Gini zmiennej podzielonej na decyle
    
    zm_x <- cut(dane_alg[,j+1], breaks = unique(quantile(dane_alg[,j+1], probs = seq(0,1,0.1))), include.lowest = TRUE)
    zm_x2<-factor(zm_x, ordered = TRUE, levels = levels(zm_x))
    
    levels(zm_x_org) <- podmiana_litery(levels(zm_x_org))
    levels(zm_x2) <- podmiana_litery(levels(zm_x2))
    
    wynik1 <- mergeFunc(zm_x_org,  flaga = dane_alg$SeriousDlqin2yrs)
    wynik2 <- mergeFunc(zm_x2, flaga = dane_alg$SeriousDlqin2yrs)
    
    
    gini_alg[i,j] <- (2*roc(dane_alg$SeriousDlqin2yrs ~ wynik1 ,quiet = TRUE)$auc) -1
    
    gini_org[i,j] <- (2*roc(dane_alg$SeriousDlqin2yrs ~ wynik2 ,quiet = TRUE)$auc) -1
    
    print(j)
  }
  
  
  print(i)
}





zm_x2 <- cut(zm, breaks = quantile(zm, probs = seq(0,1,0.1)), include.lowest = TRUE)
zm_x3_org <-factor(zm_x2, ordered = TRUE, levels = levels(zm_x2))

(2*roc(flaga ~ zm_x3_org ,quiet = TRUE)$auc) -1 # 0.263346
  

  rozdzielenie <- matrix(0,ile_arg,dl_osob)
  
  for(j in 1:ile_arg){
    
    rozdzielenie[j,] <- wynik[500,(1+(dl_osob*(j-1))):(dl_osob*j)]
    
  }
  
  roz_char <- apply(rozdzielenie,1,paste0, collapse = '')
  
  xx <- apply(as.matrix(roz_char), 1, to_real, dl_osob = dl_osob, max_x = max_x, min_x = min_x)
  
  zm_x_org <- cut(zm, breaks = unique(c(min(zm), sort(xx), max(zm))), include.lowest = TRUE) 
  zm_x_org <-factor(zm_x_org, ordered = TRUE, levels = levels(zm_x_org))
  
 (2*roc(flaga ~ zm_x_org ,quiet = TRUE)$auc) -1 # 0.2744257


  dane_wykres  <- table(zm_x_org)
  dane_wykres1 <- table(zm_x3_org)
  
  ### jak wyglda bad rate w kazdym binie ##
  
  dane_wyk <- table(zm_x_org, flaga) 
  dane_wyk1 <- table(zm_x3_org, flaga)
  

  bad_rate <- dane_wyk[,2]/(dane_wyk[,1]+dane_wyk[,2])
  bad_rate1 <- dane_wyk1[,2]/(dane_wyk1[,1]+dane_wyk1[,2])
  
  
  par(mar = c(4,5,4,3))
  cent <- barplot(dane_wykres, ann = FALSE, yaxt = 'n', las = 2)
  axis(2,las = 2)
  title(ylab = "Bad rate", line = 3.5, cex.axis = 4)
  par(new=TRUE)
  plot(cent,bad_rate, xaxt='n', yaxt='n', ylim = c(0,max(bad_rate)), ann = FALSE)
  lines(cent,bad_rate, xaxt='n', yaxt='n', ylim = c(0,max(bad_rate)))
  axis(4, at = seq(0,max(bad_rate),0.02))  
  
  par(mar = c(4,5,4,3))
  cent <- barplot(dane_wykres1, ann = FALSE, yaxt = 'n', las = 2)
  axis(2,las = 2)
  title(ylab = "Bad rate", line = 3.5, cex.axis = 4)
  par(new=TRUE)
  plot(cent,bad_rate1, xaxt='n', yaxt='n', ylim = c(0,max(bad_rate1)), ann = FALSE)
  lines(cent,bad_rate1, xaxt='n', yaxt='n', ylim = c(0,max(bad_rate1)))
  axis(4, at = seq(0,max(bad_rate1),0.02))   
  
  
#install.packages("factorMerger")
library(factorMerger)


# funcja do mergowania
podmiana_litery <- function(x){
    
  y <- replace(x,1:length(x), letters[1:length(x)])
  
  return(y)
}




zamianaNA <- function(x){
  
  levels(x[is.na(levels(x))]) <- 'n'
  
  return(x)
}




  
mergeFunc <- function(flaga,zm){
  
  x_merg <- mergeFactors(response = flaga, zm, family = 'binomial', method = 'fast-adaptive')
  
  z <- getOptimalPartition(x_merg)
  
  x <- x_merg[[2]]
  
  for(i in 1:length(z)){
    
    if(nchar(z[i] > 3)){
      
      levels(x)[levels(x) %in% substring(z[i], seq(1,nchar(z[i])-1, 3), seq(3,nchar(z[i]),3))] <- z[i]
    }
  }
  return(x)
}

levels(zm_x) <- podmiana_litery(levels(zm_x))
levels(zm_x3) <- podmiana_litery(levels(zm_x3))

wynik1 <- mergeFunc(zm_x,flaga = flaga)
wynik2 <- mergeFunc(zm_x3, flaga = flaga)


(2*roc(flaga ~ wynik1 ,quiet = TRUE)$auc) -1 # 0.2728639

(2*roc(flaga ~ wynik2 ,quiet = TRUE)$auc) -1 # 0.2671494



dane_wykres <-table(wynik1)
dane_wykres1 <-table(wynik2)

### jak wyglda bad rate w kazdym binie ##

dane_wyk <- table(wynik1, flaga) 
dane_wyk1 <- table(wynik2, flaga)


bad_rate <- dane_wyk[,2]/(dane_wyk[,1]+dane_wyk[,2])
bad_rate1 <- dane_wyk1[,2]/(dane_wyk1[,1]+dane_wyk1[,2])


par(mar = c(4,5,4,3))
cent <- barplot(dane_wykres, ann = FALSE, yaxt = 'n', las = 2)
axis(2,las = 2)
title(ylab = "Bad rate", line = 3.5, cex.axis = 4)
par(new=TRUE)
plot(cent,bad_rate, xaxt='n', yaxt='n', ylim = c(0,max(bad_rate)), ann = FALSE)
lines(cent,bad_rate, xaxt='n', yaxt='n', ylim = c(0,max(bad_rate)))
axis(4, at = seq(0,max(bad_rate),0.02))  

par(mar = c(4,5,4,3))
cent <- barplot(dane_wykres1, ann = FALSE, yaxt = 'n', las = 2)
axis(2,las = 2)
title(ylab = "Bad rate", line = 3.5, cex.axis = 4)
par(new=TRUE)
plot(cent,bad_rate1, xaxt='n', yaxt='n', ylim = c(0,max(bad_rate1)), ann = FALSE)
lines(cent,bad_rate1, xaxt='n', yaxt='n', ylim = c(0,max(bad_rate1)))
axis(4, at = seq(0,max(bad_rate1),0.02))  


#quanitative consulting, finalyze fair, isaac, crif 

runif(10) %>% reduce(~sum(.))
runif(10) %>% accumulate(~sum(.))



















