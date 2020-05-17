
# import library


source("max_SLDT.R")
source("subtrees_SLDT.R")
source("selection_SLDT.R")
source("extra.R")

library(penalizedLDA)
library(cvTools)
library(e1071)





# read data


dat = iris[sample(1:nrow(iris), 600, replace = TRUE), ]
head(dat); dim(dat)

train_size = 200
validate_size = 200
predict_size = nrow(dat) - train_size - validate_size

set.seed(1)

idx = sample(1:nrow(dat),nrow(dat),replace=F)   # 1:nrow(dat)をランダムに並び替えてindに格納．

train = dat[idx[1:train_size],]
valid = dat[idx[(train_size+1):(train_size+validate_size)],]
test = dat[idx[(train_size+validate_size+1):nrow(dat)],]





# SLDT の実行 (モデル構築)


## 最大木を成長させる．
max = max_SLDT(data=train,crit=1,case_min=3,col_class="Species",nvec=1,lambda_candidates = seq(0.1,3,by=0.05))  # SLDT
max$tree

## 最大木を刈り込んで部分木群を取得する．
subtrees = subtrees_SLDT(maxtree = max$tree)
subtrees

## 最適な部分木モデルを選択する．
res = selection_SLDT(validation = valid, subtrees = subtrees,maxSLDT = max)
res$misclass_rate
selected_tree = res$selected
selected_index = which.min(res$misclass_rate$misclass_rate)





# 得られたモデルについて解釈と精度評価を行う．


# 得られた SLDT モデルの木構造を確認
selected_tree$tree

# テストデータに対する精度を測り汎化性能を評価
res = selection_SLDT(validation = test, subtrees = subtrees,maxSLDT = max)
res$misclass_rate$misclass_rate[selected_index]
kondo = table(truth=res$subt_predictions[[selected_index]][,2],pred=res$subt_predictions[[selected_index]][,1])
kondo  # confusion matrix (混同行列)
1-sum(diag(kondo))/sum(kondo)  # Accuracy

# モデルを解釈
idx = selected_tree$tree$node_index[selected_tree$tree$leave!="*"]  # 分割が行われたノードリスト．
print(selected_tree$ldas$node_1$plda$discrim)  # 1層目の分割時の判別ベクトル
# 一部の変数の係数が 0 に推定されており解釈が容易．

