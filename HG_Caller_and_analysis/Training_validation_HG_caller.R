library(dplyr)
library(caret)
library(irr)

# Define a function to extract major haplogroup from sub-haplogroup annotation obtained from MITOMAP
HGC=function(x){
  HG=substr(x, 1, 1)
  if (startsWith(x,'L')){HG=substr(x, 1, 2) }
  return(HG)
}

# Load the mitochondrial genotype of patients consisting of 255 sprobes/positions of sufficient quality (189 defining haplogroup marker).
df=read.delim('./patients.tsv',sep = '\t',stringsAsFactors = F,colClasses = "character", header = T, row.names = 1, check.names = F)

# Load genotype found in 61,134 full length sequences obtained from https://www.mitomap.org/cgi-bin/genbank_ids.cgi
dfm=read.delim('./Mitomap.Haplogroup.markers.csv',sep = '\t',stringsAsFactors = F)

# Upload the mitochondrial genotype of the 2535 individuals from the 1000 Genomes Project.
test=read.delim('./1000GP.tsv',sep = '\t',stringsAsFactors = F,colClasses = "character", header = T, row.names = 1, check.names = F)

# Generate a single data-frame with all haplogroup markers and annotate each sequence with its corresponding major haplogroup.
genebank=do.call(rbind,lapply(list.files('./HG'),function(f){read.delim(paste0('./HG/',f),sep = '\t',stringsAsFactors = F)}))
genebank$HG=sapply(genebank$Haplogroup,HGC)

# Obtains unique mitochondrial genomic positions that change with respect to the reference.
positions=unique(unlist(strsplit(split = ', ',fixed = T,x = gsub(pattern = ".","",fixed = T,gsub('[A-Z]','',paste0(genebank$Variants,collapse = ", "))))))

# Define reference matrix for HG markers 
unique.pos=unique(dfm[,c(2,3)])
m=matrix(rep(unique.pos$rCRS,nrow(genebank)),nrow = nrow(genebank),ncol = nrow(unique.pos),byrow = T)
colnames(m)=as.character(unique.pos$Pos)
row.names(m)=genebank$Genbank.ID

# We first limit the genomic positions we are going to work with to those for which we have information in patients. 
m=m[,which(colnames(m) %in% colnames(df))]


# Based on the information contained in the genotypes from the 61,134 complete sequences, we defined the variants of 
# these sequences that affect haplogroup markers, defining matrix of HG markers matrix from 61,134 sequences that 
# will be used to train our random forest model.
for (i in 1:nrow(genebank)){
  p=gsub('.','',fixed = T,gsub('[A-Z]','',fixed = F,unlist(strsplit(genebank$Variants[i],split = ', '))))
  alt=gsub('[0-9]','',fixed = F,unlist(strsplit(genebank$Variants[i],split = ', ')))
  alt[grepl('d.',alt)]='-'
  m[i,p[which(p %in% colnames(m))]]=alt[which(p %in% colnames(m))]
}
# Define training dataframe
m=as.data.frame(m)
m$Label=as.character(genebank$HG)
m <- m %>% mutate_all(as.factor)

# Create train/test partitions to measure model error in 3 fold Cross Validation, defined as Cohen's kappa.
set.seed(1)
trainIndex <- createFolds(m$Label, k = 3, list = TRUE, returnTrain = TRUE)

# Define fine tune random forest parameters
myGrid <- expand.grid(mtry = c(5, 10, 20, 40, 60),
                        splitrule = c("gini", "extratrees"),
                        min.node.size = 1) ## Minimal node size; default 1 for classification

# Define empty vectors where gather model error in each fold and its significance.
k=c()
p=c()

# Calculate model error
for (fold in trainIndex){
  
  set.seed(1)
  model <- train(Label ~ .,
                 data = m[-fold,],
                 method = "ranger",
                 tuneGrid = myGrid,
                 trControl = trainControl(method = "cv",
                                          allowParallel = TRUE,
                                          number = 3,
                                          verboseIter = FALSE))
  
  
  prediction=predict(model,m[fold,])
  results.df=cbind(as.character(as.factor(m$Label[fold])),as.character(as.factor(predict(model,m[fold,]))))
  k <-c(k,kappa2(results.df)$value)
  p <-c(p,kappa2(results.df)$p.value)
}

mean(k)# Mean error Kappa=0.98
mean(p)# Mean p-value=0

# Train our final random forest model using all data
set.seed(1)
myGrid <- expand.grid(mtry = c(5, 10, 20, 40, 60),
                      splitrule = c("gini", "extratrees"),
                      min.node.size = 1) ## Minimal node size; default 1 for classification

model <- train(Label ~ .,
               data = m,
               method = "ranger",
               tuneGrid = myGrid,
               trControl = trainControl(method = "cv",
                                        allowParallel = TRUE,
                                        number = 3,
                                        verboseIter = FALSE))


# Perform external validation

# Identify unique levels in the training set variables
levels_train <- lapply(m[,-which(names(m)=='Label')], function(x) if(is.factor(x)) levels(x) else unique(x))

# Function to identify unclassifiable individuals in the test set
identify_unclassifiable <- function(test, levels_train) {
  unclassifiable <- apply(test, 1, function(individual) {
    any(sapply(names(individual), function(var) {
      !(individual[[var]] %in% levels_train[[var]])
    }))
  })
  return(which(unclassifiable))
}

# Obtain unclassified individuals with our trained model
unclassifiable_indexes <- identify_unclassifiable(test[,which(colnames(test) %in% colnames(m)[-which(colnames(m)=="Label")])], levels_train)

# Computation of model error in external data from the 1000 Genomes Project.
results.df=cbind(as.character(as.factor(test$Label[-unclassifiable_indexes])),as.character(as.factor(predict(model,test[-unclassifiable_indexes,]))))

# Add unclassified individuals in error estimation
results.df.na=cbind(as.character(as.factor(test$Label[unclassifiable_indexes])),rep('-',length(unclassifiable_indexes)))
results.df=rbind(results.df,results.df.na)

kappa2(results.df)# Kappa=0.961, p-value=0


results.patients=data.frame(ID=rownames(df),HG=as.character(as.factor(predict(model,df))))

saveRDS(model,'./HG_caller.rds')
write.table(results.patients,'./HG.patients.tsv',sep = '\t',quote = F, row.names = F)
