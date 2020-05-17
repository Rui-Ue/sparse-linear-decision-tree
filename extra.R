#*************************************************************************
# gini()
# ==> called by tunePLDA_SLDT(), calculate Gini index
#*************************************************************************



gini = function(p){
  sum(p*(1-p))
}





#*************************************************************************
# clean()
# ==> called by tunePLDA_SLDT(), remove useless variables
#*************************************************************************



clean = function(dat,numcol){
  for(l in numcol){
    jadge = F
    if(var(dat[,l])==0){
      jadge = T
    }
    for(m in unique(dat$Y)){
      if(var(dat[dat$Y==m,l])==0){
        jadge = T
      }
    }
    jadge_mean=T
    meantmp = mean(dat[dat$Y==unique(dat$Y)[1],l])
    for(m in unique(dat$Y)){
      if(meantmp != mean(dat[dat$Y==m,l])){
        jadge_mean=F
      }
    }
    if(jadge_mean==T){
      jadge=T
    }
    if(jadge==T){
      numcol = numcol[-which(numcol==l)]
    }
  }
  return(numcol)
}





#*************************************************************************
# PenalizedLDA_my()
# ==> PenalizedLDA for SLDT
#*************************************************************************



PenalizedLDA_my = function(x,y,K,type,lambda){
  label_table_2 = data.frame(moto=unique(y),kari=1:length(unique(y)))
  ykari = numeric(length(y))
  for(i in 1:length(y)){
    for(j in 1:nrow(label_table_2)){
      if(y[i]==label_table_2[j,1]){
        ykari[i] = label_table_2[j,2]
      }
    }
  }
  res = PenalizedLDA(x,ykari,K=K,type=type,lambda=lambda)
  return(list(res=res,label_table_2=label_table_2))
}





#*************************************************************************
# predict_my()
# ==> prediction for SLDT
#*************************************************************************



predict_my = function(lda,x,label_table_2){
  res = predict(lda,x)
  for(i in 1:ncol(res$ypred)){   # ypredの列を走査．nvec=1の時は1列だけだが，nvec=2の時は2列あるから．
    for(j in 1:nrow(res$ypred)){ # ypredの行を操作．一個一個見て，label_table_2を使って元のクラスラベルに変換．
      for(k in 1:nrow(label_table_2)){
        if(res$ypred[j,i]==label_table_2[k,2]){
          res$ypred[j,i] = label_table_2[k,1]
        }
      }
    }
  }
  return(res)
}

