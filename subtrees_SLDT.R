#*************************************************************************
# subtree_SLDT()
# ==> get pruned subtrees
#*************************************************************************

# Arguments
# maxtree: max_SLDT() が出力する list の tree 要素.



subtrees_SLDT<-function(maxtree){
  d_max<-max(as.numeric(as.character(unlist(maxtree$depth))))            # maximal treeのdepth．depthの最大値．
  SLDT_seq<-list()
  for (i in 0:d_max) {
    tree <- as.data.frame(maxtree[which(as.numeric(as.character(maxtree$depth)) <= i), ])
    tree$action <-ifelse(as.numeric(as.character(tree$depth)) == i, -1,as.numeric(as.character(tree$action)))
    tree$leave <-ifelse(as.numeric(as.character(tree$action)) < 0, "*", "")
    SLDT_seq[[i+1]] <- tree
  }
  return(SLDT_seq)
}





