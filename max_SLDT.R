#*************************************************************************
# max_SLDT()
# ==> build a maximal tree
#*************************************************************************

# Arguments
# data: 最大木構築用データが格納されたデータフレーム．クラスラベル列以外の全ての列が説明変数として使われる.
# crit: 不純度の指標として，1 の時はクロスバリデーション Gini 係数，2 の時は Gini係数，3 の時はクロスバリデーション誤判別率が使われる．
# case_min: 最大木の成長の停止基準．サンプルサイズが case_min 以上のクラスが複数存在しなくなった時，分割を止める．
# kfold: 不純度の指標の算出時に行われるクロスバリデーションの fold 数．
# lambda_candidates: チューニングパラメータ \lambda の候補値を格納したベクトル.
# col_class: クラスラベルが記録されている列の名前．
# nvec: 用いる判別ベクトルの本数．1 あるいは 2 を指定．



max_SLDT = function(data,crit=1,case_min=3,kfold=3,lambda_candidates=c(1e-4,1e-3,0.01,0.1,1,10),col_class,
                    nvec = 1){
  
  ncol_class = which(colnames(data)==col_class)
  data = data[,c(ncol_class,(1:ncol(data))[-ncol_class])]
  names(data)[1] = "Y"
  
  data$Y = as.character(data$Y)
  tmpY = data$Y
  K = length(unique(data$Y))
  for(j in 1:K){
    data$Y[which(tmpY==unique(tmpY)[j])] = j
  }
  data$Y = as.numeric(data$Y)
  
  label_table = data.frame(class_name = as.character(unique(tmpY)),value = 1:K)
  
  tree = NA
  parent = c()
  depth = c()
  best_ldas = list()
  n = c()
  n_class = list()
  #  yval = c()
  yhat = c()
  prop = list()
  pop = list()
  splits = list()
  action = c()
  best_lambdas = c()
  n_zero1 = c()
  n_zero2 = c()
  best_impurity_downs = c()
  nums = list()
  parent[1] = NA
  depth[1] = 0
  n[1] = dim(data)[1]
  pop[[1]] = rownames(data)
  yhat[1] = NA
  K = length(unique(data$Y))
  
  
  
  i = 1
  
  while(i <= length(n)){
    
    
    
    node = data[pop[[i]],]
    
    nodeYfac = factor(node$Y,levels=1:K)

    prop[[i]] = numeric(K) # ノードiにおけるクラス比率
    for(j in 1:K){
      prop[[i]][j] = round((length(which(node$Y==j))/n[i]),4)
    }
    
    n_class[[i]] = numeric(K) # ノードiにおけるクラス度数分布
    for(j in 1:K){
      n_class[[i]][j] = length(which(node$Y==j))
    }

    cnt = 0                     # 新終了基準。case_minより大きいサンプルサイズを持つクラスが複数あれば、split、なければstop。
    for(j in 1:K){
      if(n_class[[i]][j] >= case_min){
        cnt = cnt + 1
      }
    }
    if(cnt<2){
      stopsplit = T
    } else{
      stopsplit = F
    }
    
    
    if(stopsplit == FALSE){
      
      res_tunePLDA = tunePLDA_SLDT(node = node,nfold = kfold,crit = crit,candidates = lambda_candidates,K=K,nvec=nvec,case_min=case_min)
      
      if(class(res_tunePLDA)!="character"){ # numの長さがゼロ，すなわちPLDAで使える(質の悪くない)変数が1つもなかった場合．分割を止める
        
        best_lambdas[i] = res_tunePLDA$best_lambda
        best_impurity_downs[i] = res_tunePLDA$best_impurity_down # 第i要素に,ノードiを最適に分割した際のimpurity減少量． 
        
        if(best_impurity_downs[i] > 0){    # impurityが減少する分割だったら，実際にその分割を行う．
          
          action[i] = 1                    # ノードiを分割する，というサイン．

          if(res_tunePLDA$best_nvec == 1){  # 
            pred = (predict_my(res_tunePLDA$best_lda$res,as.matrix(node[,res_tunePLDA$num]),label_table_2=res_tunePLDA$best_lda$label_table_2)$ypred)[,1]   # predに,ノードiの分割結果(予測結果)を格納．
          } else{
            pred = (predict_my(res_tunePLDA$best_lda$res,as.matrix(node[,res_tunePLDA$num]),label_table_2=res_tunePLDA$best_lda$label_table_2)$ypred)[,2]
          }
          
          parent = c(parent,rep(i,K))
          depth = c(depth,rep(depth[i]+1,K))
          
          yhat = c(yhat,1:K)
          
          for(j in 1:K){
            n = c(n,length(which(pred==j)))
          }
          
          for(j in 1:K){
            pop = c(pop,list(rownames(node[which(pred==j),])))
          }
          
          best_ldas[[i]] = list(plda=res_tunePLDA$best_lda$res,label_table_2=res_tunePLDA$best_lda$label_table_2)
          
          nums[[i]] = res_tunePLDA$num
          
          if(res_tunePLDA$best_nvec == 2){
            disvec1 = best_ldas[[i]]$plda$discrim[,1]
            disvec2 = best_ldas[[i]]$plda$discrim[,2]
            numvec1 = c()
            numvec2 = c()
            for(j in 1:length(disvec1)){
              if(disvec1[j]!=0){
                numvec1 = c(numvec1,j) 
              }
            }
            for(j in 1:length(disvec2)){
              if(disvec2[j]!=0){
                numvec2 = c(numvec2,j)
              }
            }
            n_zero1[i] = length(disvec1) - length(numvec1)
            n_zero2[i] = length(disvec2) - length(numvec2)
          } else{
            disvec1 = best_ldas[[i]]$plda$discrim[,1]
            numvec1 = c()
            for(j in 1:length(disvec1)){
              if(disvec1[j]!=0){
                numvec1 = c(numvec1,j)
              }
            }
            n_zero1[i] = length(disvec1) - length(numvec1)
            n_zero2[i] = NA
          }
  
        }else{                     #  どうやってもimpurityを減少させられない場合．
          action[i] = -2           #  ↑の理由で分割をstopした，というサイン．
          n_zero1[i] = NA
          n_zero2[i] = NA
          best_lambdas[i] = NA
        }
      }else{                       # numの長さがゼロ，すなわちPLDAで使える(質の悪くない)変数が1つもなかった場合．分割を止める．
        action[i] = -3             # サイン．
        best_impurity_downs[i] = NA
        n_zero1[i] = NA
        n_zero2[i] = NA
        best_lambdas[i] = NA
      }
    }else{                       # case_minより大きいサンプルサイズを持つクラスが複数なかった場合(新停止基準が満たされた場合)．つまり，ノードが完全に純粋になっていた場合
      action[i] = -1             # この理由で分割をstopした，というサイン．
      best_impurity_downs[i] = NA
      n_zero1[i] = NA
      n_zero2[i] = NA
      best_lambdas[i] = NA
    }
    
    
    i = i+1
    
    pred = NA
    
    
  }
  
  
  leave = ifelse(action<0,"*","")
  
  node_index = 1:length(depth)      # 各ノードに割り当てる番号．別にdepth以外でも良いが．length(depth)が総ノード数になってる．
  
  
  n_class_m = matrix(numeric(length(depth)*K),ncol=K)
  for(i in 1:length(depth)){
    n_class_m[i,] = n_class[[i]]
  }
  n_class_m_names = c()
  for(i in 1:K){
    n_class_m_names = c(n_class_m_names,paste("n[",as.character(label_table$class_name[(label_table$value)==i]),"]",sep=""))
  }
  colnames(n_class_m) = n_class_m_names
  
  prop_m = matrix(numeric(length(depth)*K),ncol=K)
  for(i in 1:length(depth)){
    prop_m[i,] = prop[[i]]
  }
  prop_m_names = c()
  for(i in 1:K){
    prop_m_names = c(prop_m_names,paste("p[",as.character(label_table$class_name[(label_table$value)==i]),"]",sep=""))
  }
  colnames(prop_m) = prop_m_names
  
  
  tree = as.data.frame(cbind(action,best_lambdas=round(best_lambdas,3),n_zero1,n_zero2,depth,parent,n,n_class_m,prop_m,yhat,leave,node_index))
  
  
  out = list(tree=tree,best_ldas=best_ldas,best_impurity_downs=best_impurity_downs,nums=nums,label_table=label_table,col_class=col_class)
  
  return(out)
  
}





#*************************************************************************
# tunePLDA_SLDT()
# ==> called by max_SLDT(), perform cross validation and tuning for PLDA.
#*************************************************************************



tunePLDA_SLDT = function(node,nfold,candidates,crit,K,nvec,case_min=case_min){
  
  node_ex = node
  for(i in 1:K){
    if((0<nrow(node_ex[node_ex$Y==i,]))&&(nrow(node_ex[node_ex$Y==i,])<case_min)){
      node_ex = node_ex[-1*which(node_ex$Y==i),]
    }
  }
  
  
  
  if(crit==2){              # ただのGiniを使う場合，すなわちCVしない場合．
    num = 2:ncol(node)
    num = clean(dat = node_ex,numcol = num) # numにpldaで使う変数(質の悪い変数以外の変数)の番号を格納
    if(length(num)>0){                   # pldaで使える変数(質の悪くない変数)が存在する場合。
      if(length(unique(node_ex$Y))==2){     # 2クラスしかなかったら，nvecは必ず1にせざるを得ない．
        nvec = 1
      }
      
      if(nvec==1){  # 第1判別ベクトル1本だけを使う場合
        val_crit = numeric(length(candidates))
        for(i in 1:length(candidates)){
          lda_my = PenalizedLDA_my(as.matrix(node_ex[,num]),as.numeric(as.character(node_ex$Y)),K=1,type="standard",lambda=candidates[i]) # PLDAを実行。# abstに書いた問題点防ぐための自作関数．
          y_fac = as.factor(node$Y)                  # 各個体の所属クラスYをファクターで取得。
          levels(y_fac) = as.character(1:K)          # こうすることで、個体数0のクラスも水準に入れて考慮できる。
          pred_fac = as.factor(predict_my(lda_my$res,as.matrix(node[,num]),label_table_2=lda_my$label_table_2)$ypred[,1]) # 各個体の予測クラスをファクターで取得。abstに書いた問題点防ぐための自作関数．
          levels(pred_fac) = as.character(1:K)
          prop_ko_node = prop.table(table(pred_fac)) # 子ノードの大きさ(サンプルサイズ)の比率
          kondo = table(pred_fac,y_fac)              # PLDAの予測に関する混同行列
          oya_gini = gini(prop.table(table(y_fac)))  # 親ノード(分割前)のGini係数。
          ko_gini = numeric(K)                       # 分割して得られる各子ノードのGini係数。子ノードがK個もない場合があるが、処理的に便宜上K要素作っておく。どうせサンプルサイズ比率をかけるときに0がかけれれて消えるので気にしない。
          for(j in 1:K){
            if(sum(kondo[j,])!=0){                        # 所属個体0の子ノードに対してGiniなんて定義計算できないし，する必要もない．なので適当に0を入れておく．どうせあとで0がかかるから，適当な値入れておけば良い．
              ko_gini[j] = gini(prop.table(kondo[j,])) # 混同行列の第j行：クラスjと予測された子ノードにおける所属クラス度数分布
            } else{                                  # 所属個体が存在する子ノード，
              ko_gini[j] = 0
            }
          }
          val_crit[i] = oya_gini - prop_ko_node %*% ko_gini  # チューニング基準値：Gini減少量(親Gini - 各サイズで重みづけした子ノードGiniの和)
        }
      }
      
      if(nvec==2){  # 判別ベクトルを2本使う場合．
        val_crit = numeric(length(candidates))
        for(i in 1:length(candidates)){
          lda_my = PenalizedLDA_my(as.matrix(node_ex[,num]),as.numeric(as.character(node_ex$Y)),K=2,type="standard",lambda=candidates[i]) # 判別ベクトル2本使うのでK=2．
          y_fac = as.factor(node$Y)                  
          levels(y_fac) = as.character(1:K)          
          pred_fac = as.factor(predict_my(lda_my$res,as.matrix(node[,num]),label_table_2=lda_my$label_table_2)$ypred[,2]) # 判別ベクトル2本使うので$ypred[,2]
          levels(pred_fac) = as.character(1:K)
          prop_ko_node = prop.table(table(pred_fac)) 
          kondo = table(pred_fac,y_fac)              
          oya_gini = gini(prop.table(table(y_fac)))  
          ko_gini = numeric(K)                       
          for(j in 1:K){
            if(sum(kondo[j,])!=0){
              ko_gini[j] = gini(prop.table(kondo[j,])) # 混同行列の第j行：クラスjと予測された子ノードにおける所属クラス度数分布
            } else{                                  # 所属個体が存在する子ノード，すなわち実在する?子ノード．prop.tableでNaNが発生してしまうことに対処するためのエラー処理
              ko_gini[j] = 0
            }
          }
          val_crit[i] = oya_gini - prop_ko_node %*% ko_gini  # チューニング基準値：Gini減少量(親Gini - 各サイズで重みづけした子ノードGiniの和)
        }
      }
      
      if(nvec==0){ # 1本使うべきか2本使うべきかを，critをもとに選択する場合．
        val_crit_mat = matrix(numeric(2*length(candidates)),ncol=2) # tuning基準の値を格納する行列(第1列が1本使用時、第2列が2本使用時)
        for(i in 1:length(candidates)){
          lda1_my = PenalizedLDA_my(as.matrix(node_ex[,num]),as.numeric(as.character(node_ex$Y)),K=1,type="standard",lambda=candidates[i]) # PLDAを実行。
          lda2_my = PenalizedLDA_my(as.matrix(node_ex[,num]),as.numeric(as.character(node_ex$Y)),K=2,type="standard",lambda=candidates[i]) # PLDAを実行。
          y_fac = as.factor(node$Y)                  
          levels(y_fac) = as.character(1:K)          
          pred_fac_1 = as.factor(predict_my(lda1_my$res,as.matrix(node[,num]),label_table_2=lda1_my$label_table_2)$ypred[,1]) # 判別ベクトル1本使った場合
          pred_fac_2 = as.factor(predict_my(lda2_my$res,as.matrix(node[,num]),label_table_2=lda2_my$label_table_2)$ypred[,2]) # 判別ベクトル2本使った場合
          levels(pred_fac_1) = as.character(1:K)
          levels(pred_fac_2) = as.character(1:K)
          prop_ko_node_1 = prop.table(table(pred_fac_1)) 
          prop_ko_node_2 = prop.table(table(pred_fac_2))
          kondo_1 = table(pred_fac_1,y_fac)
          kondo_2 = table(pred_fac_2,y_fac)
          oya_gini = gini(prop.table(table(y_fac)))  
          ko_gini_1 = numeric(K)
          ko_gini_2 = numeric(K)
          for(j in 1:K){
            if(sum(kondo_1[j,])!=0){
              ko_gini_1[j] = gini(prop.table(kondo_1[j,])) # 混同行列の第j行：クラスjと予測された子ノードにおける所属クラス度数分布
            } else{                                  # 所属個体が存在する子ノード，すなわち実在する?子ノード．prop.tableでNaNが発生してしまうことに対処するためのエラー処理
              ko_gini_1[j] = 0
            }
            if(sum(kondo_2[j,])!=0){
              ko_gini_2[j] = gini(prop.table(kondo_2[j,])) # 混同行列の第j行：クラスjと予測された子ノードにおける所属クラス度数分布
            } else{                                  # 所属個体が存在する子ノード，すなわち実在する?子ノード．prop.tableでNaNが発生してしまうことに対処するためのエラー処理
              ko_gini_2[j] = 0
            }
          }
          val_crit_mat[i,1] = oya_gini - prop_ko_node_1 %*% ko_gini_1  # 1本の場合のチューニング基準値：Gini減少量
          val_crit_mat[i,2] = oya_gini - prop_ko_node_2 %*% ko_gini_2  # 2本の場合のチューニング基準値：Gini減少量
        }
        if(max(val_crit_mat[,1])<max(val_crit_mat[,2])){    # nvecに，選択された本数を格納しておく．
          nvec = 2
          val_crit = val_crit_mat[,2]                       # val_critに，選択された本数における基準の値たちを代入．こうすることで，後で統一的に扱えるようにする．
        } else{
          nvec = 1
          val_crit = val_crit_mat[,1]
        }
      } 
    } else{
      return("cannot_execute_plda_for_numzero")
    }
  }
  
  if(crit==1||crit==3){   # CVする場合（CVGiniもしくはCV誤判別率を使う場合）
    
    w = list()            # w[[i]]に，クラスiの個体の行番号が格納される．
    for(i in 1:K){
      w[[i]] = which(node_ex$Y==i)
    }
    
    samples = list()      # samples[[i]]に，w[[i]]の何番目をどのfoldに割り当てるかの表が格納される．
    for(i in 1:K){
      if(length(w[[i]])>0){  # nodeからnode_exにする過程で捨てられたクラスiはlength(w[[i]])が0になっていてcvFolds()でエラー起きる．その対処
        samples[[i]] = cvFolds(length(w[[i]]),K=nfold,R=1,type="random")
      } else{
        samples[[i]] = -1    # クラスiは捨てられているよというサイン．↓のclass()判定で活きる．
      }
    }
    
    index_CVtrains = list() # [[i]]に，第ifoldを抜いた(テストデータとした)時のCV学習データの(node_exにおける)行番号を格納．
    for(i in 1:nfold){
      index_CVtrain = c()
      for(j in 1:K){
        if(class(samples[[j]])=="cvFolds"){
          index_CVtrain = c(index_CVtrain,w[[j]][samples[[j]]$subset[samples[[j]]$which!=i]])
        }
      }
      index_CVtrains[[i]] = index_CVtrain
    }
    
    CVtrains = list()  # [[i]]に，第ifoldを抜いた(テストデータとした)時のCV学習データを格納．
    for(i in 1:nfold){
      CVtrains[[i]] = node_ex[index_CVtrains[[i]],]
    }
    
    CVtests = list()   # [[i]]に，第ifoldを抜いた(テストデータとした)時のCVテストデータ(つまり第ifoldのデータ)を格納．
    for(i in 1:nfold){
      CVtests[[i]] = node_ex[-1*(index_CVtrains[[i]]),]
    }
    
    nums = list() # [[i]]に，第ifoldを抜いた時のCV学習データにおいて，pldaで使う(質の悪くない)変数の列番号を格納．
    for(i in 1:nfold){
      nums[[i]] = clean(dat=CVtrains[[i]],numcol = 2:ncol(CVtrains[[i]]))
    }
    
    numlengths = c() # [i]に，nums[[i]]の長さ(変数の個数)を格納．
    for(i in 1:nfold){
      numlengths[i] = length(nums[[i]])
    }
    
    if(all(numlengths>0)==TRUE){     # CV学習データのうち1つでも「pldaに使える(質の悪くない)変数が0」だったら，分割を止める．そうでなければ分割を行う．
      if(length(unique(node_ex$Y))==2){     # 2クラスしかなかったら，nvecは必ず1にせざるを得ない．
        nvec = 1
      }
      if(nvec==1||nvec==2){  # 必ず1本または必ず2本使う場合．predictとpldaの時に1,2じゃなくてnvecを打てばこうやってまとられた．
        val_crit = numeric(length(candidates)) # [i]にi番目のチューニングパラ値を使った時のチューニング基準の値(gini,誤判別率等)を格納する．
        for(i in 1:length(candidates)){
          kondos = list()        # [[j]]に第j(foldを抜いた)CV学習データを使った時の混同行列を格納．
          for(j in 1:nfold){
            lda_my = PenalizedLDA_my(as.matrix(CVtrains[[j]][,nums[[j]]]),as.numeric(as.character(CVtrains[[j]]$Y)),K=nvec,type="standard",lambda=candidates[i]) # PLDAを実行。# abstに書いた問題点防ぐための自作関数．
            y_fac = as.factor(CVtests[[j]]$Y)                  # 各個体の所属クラスYをファクターで取得。
            levels(y_fac) = as.character(1:K)          # こうすることで、個体数0のクラスも水準に入れて考慮できる。
            pred_fac = as.factor(predict_my(lda_my$res,as.matrix(CVtests[[j]][,nums[[j]]]),label_table_2=lda_my$label_table_2)$ypred[,nvec]) # 各個体の予測クラスをファクターで取得。abstに書いた問題点防ぐための自作関数．
            levels(pred_fac) = as.character(1:K)
            kondos[[j]] = table(pred_fac,y_fac)              # PLDAの予測に関する混同行列
          }
          kondo = kondos[[1]]                         # 各foldでの混同行列を1個にまとめる．卒論でやってたのと同じこと．↑のabst見て．
          for(j in 2:nfold){
            kondo = kondo + kondos[[j]]
          }
          prop_ko_node = prop.table(apply(kondo,1,sum)) # 子ノードの大きさ(サンプルサイズ)の比率
          if(crit==1){ # CVGiniを使う場合．
            oya_gini = gini(prop.table(table(factor(node_ex$Y,levels=as.character(1:K))))) # 親ノードのGini係数
            ko_gini = numeric(K) # [j]に，第j子ノードのGini係数(のCV推定値)                  
            for(j in 1:K){
              if(sum(kondo[j,])!=0){
                ko_gini[j] = gini(prop.table(kondo[j,])) # 混同行列の第j行：クラスjと予測された子ノードにおける所属クラス度数分布
              } else{                                  # 所属個体が存在する子ノード，すなわち実在する?子ノード．prop.tableでNaNが発生してしまうことに対処するためのエラー処理
                ko_gini[j] = 0
              }
            }
            val_crit[i] = oya_gini - prop_ko_node %*% ko_gini  # チューニング基準値：Gini減少量(親Gini - 各サイズで重みづけした子ノードGiniの和)のCVバージョン．
            #            print(oya_gini);print(ko_gini);print(prop_ko_node)
          }
          if(crit==3){ # CV誤判別率を使う場合．
            oya_miscla = 1 - (max(apply(kondo,2,sum))/(sum(kondo))) # !!!!abst参照．仮の誤判別率．
            ko_miscla = numeric(K) # [j]に，第j子ノードの誤判別率(のCV推定値)
            for(j in 1:K){
              if(sum(kondo[j,])!=0){ # このエラー処理はCVGiniの時と同じ．そっちのコメント参照．
                ko_miscla[j] = 1 - kondo[j,j]/sum(kondo[j,])
              } else{
                ko_miscla[j] = 0
              }
            }
            val_crit[i] = oya_miscla - prop_ko_node %*% ko_miscla
          }
        }
      }
    } else{
      return("cannot_execute_plda_for_numzero")             # ↑の条件が満たされなかった，というサインを返却する
    }
  }
  

  if(crit==1||crit==3){ # CVした場合はnum(使える変数)を再設定する必要がある．
    num = 2:ncol(node)
    num = clean(dat = node_ex,numcol = num)
  }
  
  best_lambda = mean(candidates[which.max(val_crit)]) # 選択されたチューニングパラ \lambda の値
  
  if(nvec==1){   # 判別ベクトル1本使う(ことが選択された)場合
    best_lda = PenalizedLDA_my(as.matrix(node_ex[,num]),as.numeric(as.character(node_ex$Y)),type="standard",lambda=best_lambda,K=1)  # 選択された \lambda でのPLDAを再実行して格納
    row.names(best_lda$res$discrim) = names(node[,num])    # 判別ベクトルの要素(係数)に名前を付けておく
    colnames(best_lda$res$discrim) = c("coef_disvec1")
  }
  if(nvec==2){   # 判別ベクトル2本使う(ことが選択された)場合
    best_lda = PenalizedLDA_my(as.matrix(node_ex[,num]),as.numeric(as.character(node_ex$Y)),type="standard",lambda=best_lambda,K=2)  # 選択された \lambda でのPLDAを再実行して格納
    row.names(best_lda$res$discrim) = names(node[,num])  # 判別ベクトルの要素(係数)に名前を付けておく
    colnames(best_lda$res$discrim) = c("coef_disvec1","coef_disvec2")
  }
  
  
  
  if(crit==1||crit==3){ # abst見て．CVした場合は，best_impurity_downを計算し直す必要がある．
    y_fac = as.factor(node$Y)                  
    levels(y_fac) = as.character(1:K)          
    pred_fac = as.factor(predict_my(best_lda$res,as.matrix(node[,num]),label_table_2=best_lda$label_table_2)$ypred[,nvec])
    levels(pred_fac) = as.character(1:K)
    prop_ko_node = prop.table(table(pred_fac)) 
    kondo = table(pred_fac,y_fac)              
    if(crit==1){
      oya_gini = gini(prop.table(table(y_fac)))  
      ko_gini = numeric(K)                       
      for(j in 1:K){
        if(sum(kondo[j,])!=0){
          ko_gini[j] = gini(prop.table(kondo[j,])) # 混同行列の第j行：クラスjと予測された子ノードにおける所属クラス度数分布
        } else{                                  # 所属個体が存在する子ノード，すなわち実在する?子ノード．prop.tableでNaNが発生してしまうことに対処するためのエラー処理
          ko_gini[j] = 0
        }
      }
      best_impurity_down = oya_gini - prop_ko_node %*% ko_gini  # チューニング基準値：Gini減少量(親Gini - 各サイズで重みづけした子ノードGiniの和)
    }
    if(crit==3){
      oya_miscla = 1 - (max(apply(kondo,2,sum))/(sum(kondo))) # !!!!abst参照．仮の誤判別率．
      ko_miscla = numeric(K) # [j]に，第j子ノードの誤判別率(のCV推定値)
      for(j in 1:K){
        if(sum(kondo[j,])!=0){ # このエラー処理はCVGiniの時と同じ．そっちのコメント参照．
          ko_miscla[j] = 1 - kondo[j,j]/sum(kondo[j,])
        } else{
          ko_miscla[j] = 0
        }
      }
      best_impurity_down = oya_miscla - prop_ko_node %*% ko_miscla
    }
  }else{  # CVしなかった場合は，max(val_crit)がbest_impurity_downとなる．abst見て．
    best_impurity_down = max(val_crit)
  }
  
  print(best_impurity_down)  
  return(list(best_lambda=best_lambda,val_crit=val_crit,best_impurity_down=best_impurity_down,best_lda=best_lda,best_nvec=nvec,num=num))
}



















# projection_SLDT_3 = function(plda,label_table,shape_type){ # 背景を白にしたver．
# 
#   for(i in 1:nrow(label_table)){
#     for(j in 1:length(plda$y)){
#       if(plda$y[j]==as.numeric(label_table$value)[i]){
#         plda$y[j] = as.character(label_table$class_name)[i]
#       }
#     }
#   }
#   plda$y = factor(plda$y)
# 
#   dat = data.frame(plda$xproj,class=as.factor(plda$y))
#   print(ggplot(dat,aes(x=dat[,1],y=dat[,2],color=dat$class,shape=dat$class)) +
#           labs(color="class",shape="class") +
#           geom_point() +
#           scale_shape_manual(values = shape_type)+
#           theme_bw()
#   )
# }
# 
