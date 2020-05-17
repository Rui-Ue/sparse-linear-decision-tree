#*************************************************************************
# selection_SLDT()
# ==> select the optimal subtree
#*************************************************************************

# Arguments
# validation: 部分木選択用データが格納されたデータフレーム．
# subtrees: subtrees_SLDT() が出力する list．
# maxSLDT: max_SLDT() が出力する list．



selection_SLDT <- function(validation, subtrees,maxSLDT) {
  
  col_class = maxSLDT$col_class       # クラス列の名前
  
  ncol_class = which(colnames(validation)==col_class)     # クラス列の番号
  validation = validation[,c(ncol_class,(1:ncol(validation))[-ncol_class])]   # クラス列が先頭になるように入れ替え．
  names(validation) = "Y"                    # クラス列の名を"Y"に変更して，それ以外の列の名前はナシにする．
  validation$Y = as.character(validation$Y)
  
  label_table = maxSLDT$label_table
  
  for(i in 1:nrow(label_table)){            # クラス"名"をクラス"番号"に変換．その際の対応表がlabel_table.
    for(j in 1:dim(validation)[1]){
      if(validation$Y[j] == as.character(label_table$class_name)[i]){
        validation$Y[j] = as.numeric(label_table$value)[i]
      }
    }
  }
  validation$Y = as.numeric(validation$Y)
  
  K = nrow(label_table)
  
  N_tree <- length(subtrees)   # 候補となるsubtreeの個数．
  N <- dim(validation)[1]      # 部分木選択用データのサンプルサイズ．
  pred <- list()               # 各subtreeを用いたときの部分木選択用データに対する予測ラベル，を格納するためのリスト．
  sum_nodes <- list()         # 各subtreeの，各terminal nodesごとの情報のまとめ(summary)，を格納するためのリスト．
  misclass_rate <- rep(NA, N_tree)
  subt_predictions = list()
  
  for (k in 1:N_tree) {         # kでsubtreeを走査．
    
    predictions = Prediction_SLDT(new=validation,tree=as.data.frame(subtrees[[k]]),ldas=maxSLDT$best_ldas,nums=maxSLDT$nums)
    
    nn <- as.numeric(table(as.numeric(as.character(predictions$goal_nodes))))

    nn_class = list()
    pn_class = list()
    
    for(j in 1:K){
      nn_class[[j]] = tapply((as.numeric(as.character(predictions$Y)))==j,as.numeric(as.character(predictions$goal_nodes)),sum)
      pn_class[[j]] = nn_class[[j]]/nn
    }

    predictions$error_pred <-as.numeric(as.character(apply(predictions[, c("hat.Y", "Y")], 1, function(x)ifelse(x[1] != x[2], "1", "0"))))
    # predictionsに新たにerror_pred列を追加．この列には，予測が間違っていたら1を，当たっていたら0が，格納される．
    
    pred[[k]]<-predictions
    # 部分木kにおけるpred_pLDA()の結果，つまり部分木kのテストデータ検証の結果詳細，を，格納．
    
    misclass <-tapply(as.numeric(as.character(predictions$error_pred)), as.numeric(as.character(predictions$goal_nodes)), sum)
    # 各terminal nodesにおいて誤判別された個体の数，を，misclassに格納している．
    
    
    summaryNodes <-as.data.frame(as.numeric(as.character(names(table(as.numeric(as.character(predictions$goal_nodes)))))))  # summaryNodesの1列目は，ターミナルノードのノード番号．列名は後でつける．
    summaryNodes = cbind(summaryNodes,nn)      # 2列目は，1列目に示された各ターミナルノードに，流れ着いたデータの個数．                     
    for(j in 1:K){
      summaryNodes = cbind(summaryNodes,nn_class[[j]]) # 3列目は，2列目のうちクラス1のもの，4列目は，2列目のうちクラス2のもの...
    }
    for(j in 1:K){
      summaryNodes = cbind(summaryNodes,pn_class[[j]]) # 列として，2列目のうちクラス1の比率，クラス2の比率...を列として追加．
    }
    summaryNodes = cbind(summaryNodes,misclass)      # 列として，1列目に示された各ターミナルノードにおける，誤判別個数．
    sumnames = c("index_nodes","N")                      # これ以降で，SummaryNodesに名前づけ．
    for(j in 1:K){
      sumnames = c(sumnames,paste("N[Y=",as.character(j),"]",sep=""))
    }
    for(j in 1:K){
      sumnames = c(sumnames,paste("P[Y=",as.character(j),"]",sep=""))
    }
    sumnames = c(sumnames,"N[hat.Y!=Y]")
    names(summaryNodes) = sumnames
    
    misclass_rate[k] <-sum(summaryNodes[,"N[hat.Y!=Y]"]) / (N)

    sum_nodes[[k]]<-summaryNodes
    # 部分木kにおけるsummaryNoeuds，つまり，部分木kの各terminal nodeごとの情報，を格納．
    
    subt_predictions[[k]] = predictions # k番目の部分木を使ったときの予測結果をリストに格納．
    predictions<-NULL
    summaryNodes<-NULL
  }
  
  misclass_rate<-as.data.frame(misclass_rate)
  #  names(misclass_rate)<-c("Misclass")
  
  selec_tree = subtrees[[which.min(misclass_rate$misclass_rate)]]
  split_nodes = which((selec_tree$leave)!="*")  # 選択された部分木において，分割される(PLDAが行われる)ノードの番号．
  selec_ldas = list()
  for(i in 1:length(split_nodes)){
    selec_ldas[[i]] = maxSLDT$best_ldas[[split_nodes[i]]]
    names(selec_ldas)[i] = paste("node",as.character(split_nodes[i]),sep="_")
  }
  
  
  selected = list(tree=selec_tree,ldas=selec_ldas,label_table=label_table)
  
  res = list(selected = selected,misclass_rate = misclass_rate,pred = pred,summary_nodes = sum_nodes,subt_predictions = subt_predictions) 
  return(res)
}





#*************************************************************************
# predFunction_SLDT
# ==> called by selection_SLDT(), predict from a subtree model
#*************************************************************************



Prediction_SLDT <- function(new, tree, ldas, nums) {
  
  npred = dim(new)[1]               # クラスを予測しようとしているデータ(説明変数値)のサイズ．
  pred = rep(NA, npred)             # 予測値を格納するためのベクトル．
  # scoreは使わないので廃止．
  noeuds = rep(NA, npred)           # 最終的に行き着いたノード番号を格納するためのベクトル．
  
  res = sapply(1:npred,predFunction_SLDT,new=new,ldas=ldas,nums=nums,tree=tree) # sapplyを使うと，ある関数を(第1引数として)複数の値に一気に適用することが出来る．ここでは，predFunctionを，第1引数を1:Pとして一気に実行している．
  
  return(as.data.frame(cbind(
    hat.Y=as.numeric(as.character(res["pred",])),   # 予測クラス
    
    Y=as.numeric(as.character(new$Y)),              # 正解クラス
    
    goal_nodes=as.numeric(as.character(res["goal_node",])))))     # 最終地点ノード番号
}





#*************************************************************************
# predFunction_SLDT
# ==> called by prediction_SLDT()
#*************************************************************************



predFunction_SLDT = function(x, new,ldas, nums,tree){
  
  ind_node = 1  # 今いるノード番号を保持し続けるオブジェクト．
  
  while (as.numeric(as.character(tree$action[ind_node])) == 1) {   # ind_node番目のノードのactionが1,すなわち「ターミナルノード以外」である間，whileを続ける．
    
    
    if(ldas[[ind_node]]$plda$K == 2){  # 今いるnodeにおけるPLDAでクラスを予測．
      temp_pred <-as.numeric((predict_my(ldas[[ind_node]]$plda, new[x, unlist(nums[[ind_node]])],label_table_2=ldas[[ind_node]]$label_table_2)$ypred)[,2])
    } else{
      temp_pred <-as.numeric((predict_my(ldas[[ind_node]]$plda, new[x, unlist(nums[[ind_node]])],label_table_2=ldas[[ind_node]]$label_table_2)$ypred)[,1])
    }
    
    
    ind_node<-as.numeric(as.character(which(as.numeric(as.character(tree$parent))==as.numeric(as.character(ind_node)) & as.numeric(as.character(tree$yhat))==temp_pred)))  # 次の子ノードに進む(indice_nodeをずらす)，というイメージ．ノードiでクラスkと予測された際に次に進むべきノードは，親がiかつ割り当てクラスがkのノード
    
  }
  
  
  pred <- as.numeric(as.character(tree$yhat[ind_node]))  # 最終的な(ターミナルノードまで進んだうえでの)予測クラス．
  
  return(c(pred=as.numeric(as.character(pred)),goal_node=ind_node)) # predに最終的な予測クラス，goal_nodeに最終的に行き着いたノード番号を格納．
}



