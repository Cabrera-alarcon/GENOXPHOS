################################################
### HG risk factors analysis in all patients ###
################################################

library(lme4)

# Load data from patients data
data=read.delim('./Covid_HG_analysis_data.tsv',sep = '\t',stringsAsFactors = F)

# Perform min-max normalization to handle the same scale
data$sex=scale(data$sex,center = min(data$sex),scale = max(data$sex) - min(data$sex))
data$age=scale(data$age,center = min(data$age),scale = max(data$age) - min(data$age))
data$euPC1=scale(data$euPC1,center = min(data$euPC1),scale = max(data$euPC1) - min(data$euPC1))
data$euPC2=scale(data$euPC2,center = min(data$euPC2),scale = max(data$euPC2) - min(data$euPC2))
data$euPC3=scale(data$euPC3,center = min(data$euPC3),scale = max(data$euPC3) - min(data$euPC3))
data$euPC4=scale(data$euPC4,center = min(data$euPC4),scale = max(data$euPC4) - min(data$euPC4))
data$euPC5=scale(data$euPC5,center = min(data$euPC5),scale = max(data$euPC5) - min(data$euPC5))
data$euPC6=scale(data$euPC6,center = min(data$euPC6),scale = max(data$euPC6) - min(data$euPC6))
data$euPC7=scale(data$euPC7,center = min(data$euPC7),scale = max(data$euPC7) - min(data$euPC7))
data$euPC8=scale(data$euPC8,center = min(data$euPC8),scale = max(data$euPC8) - min(data$euPC8))
data$euPC9=scale(data$euPC9,center = min(data$euPC9),scale = max(data$euPC9) - min(data$euPC9))
data$euPC10=scale(data$euPC10,center = min(data$euPC10),scale = max(data$euPC10) - min(data$euPC10))

# Univariate analysis to perform feature selection taking into account population stratification

# Mitochondrial haplogroup assesment
# Fit fit mixed-efect logistic regression model to assess regressor, considerin population stratification (Subgroup)
m.HV=glmer(paste0('Severity~','Branch_HV + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.HV)

# Calculate Odds Ratio and 95% CI
OR.tab <-data.frame(exp(cbind(fixef(m.HV),confint.merMod(m.HV,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.I=glmer(paste0('Severity~','I + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.I)

OR.tab <-data.frame(exp(cbind(fixef(m.I),confint.merMod(m.I,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.X=glmer(paste0('Severity~','X + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.X)

OR.tab <-data.frame(exp(cbind(fixef(m.X),confint.merMod(m.X,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.W=glmer(paste0('Severity~','W + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.W)

OR.tab <-data.frame(exp(cbind(fixef(m.W),confint.merMod(m.W,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.J=glmer(paste0('Severity~','J + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.J)

OR.tab <-data.frame(exp(cbind(fixef(m.J),confint.merMod(m.J,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.T=glmer(paste0('Severity~','T + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.T)

OR.tab <-data.frame(exp(cbind(fixef(m.T),confint.merMod(m.T,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.U=glmer(paste0('Severity~','U + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.U)

OR.tab <-data.frame(exp(cbind(fixef(m.U),confint.merMod(m.U,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.K=glmer(paste0('Severity~','K + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.K)

OR.tab <-data.frame(exp(cbind(fixef(m.K),confint.merMod(m.K,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

# Age assessment
m.age=glmer(paste0('Severity~','age + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.age)

OR.tab <-data.frame(exp(cbind(fixef(m.age),confint.merMod(m.age,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

# Sex assessment
m.sex=glmer(paste0('Severity~','sex + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.sex)

OR.tab <-data.frame(exp(cbind(fixef(m.sex),confint.merMod(m.sex,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

# Nuclear background assessment
m.pc1=glmer(paste0('Severity~','euPC1 + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.pc1)

OR.tab <-data.frame(exp(cbind(fixef(m.pc1),confint.merMod(m.pc1,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.pc2=glmer(paste0('Severity~','euPC2 + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.pc2)

OR.tab <-data.frame(exp(cbind(fixef(m.pc2),confint.merMod(m.pc2,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.pc3=glmer(paste0('Severity~','euPC3 + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.pc3)

OR.tab <-data.frame(exp(cbind(fixef(m.pc3),confint.merMod(m.pc3,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.pc4=glmer(paste0('Severity~','euPC4 + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.pc4)

OR.tab <-data.frame(exp(cbind(fixef(m.pc4),confint.merMod(m.pc4,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.pc5=glmer(paste0('Severity~','euPC5 + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.pc5)

OR.tab <-data.frame(exp(cbind(fixef(m.pc5),confint.merMod(m.pc5,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.pc6=glmer(paste0('Severity~','euPC6 + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.pc6)

OR.tab <-data.frame(exp(cbind(fixef(m.pc6),confint.merMod(m.pc6,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.pc7=glmer(paste0('Severity~','euPC7 + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.pc7)

OR.tab <-data.frame(exp(cbind(fixef(m.pc7),confint.merMod(m.pc7,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.pc8=glmer(paste0('Severity~','euPC8 + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.pc8)

OR.tab <-data.frame(exp(cbind(fixef(m.pc8),confint.merMod(m.pc8,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.pc9=glmer(paste0('Severity~','euPC9 + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.pc9)

OR.tab <-data.frame(exp(cbind(fixef(m.pc9),confint.merMod(m.pc9,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

m.pc10=glmer(paste0('Severity~','euPC10 + (1|Subgroups)'), data = data,family = 'binomial')
summary(m.pc10)

OR.tab <-data.frame(exp(cbind(fixef(m.pc10),confint.merMod(m.pc10,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

# Define a function to go on with feature selection to be considered in multivariate model, using 
# step-wise backward strategy of AIC reduction. It also allows the selection of BIC reduction, but
# we used AIC, to handle a uniform criteria with SCOURGE cohort analysis.
backward_elimination <- function(data, response, fixed_effects, random_effects, criterion = "AIC") {
  # Build the initial formula of the model
  fixed_effects_formula <- paste(fixed_effects, collapse = " + ")
  full_formula <- as.formula(paste(response, "~", fixed_effects_formula, "+", random_effects))
  
  # Adjust the complete model
  model <- glmer(full_formula, data = data, family = binomial)
  
  # Initialize variables
  best_model <- model
  best_formula <- full_formula
  best_criterion <- ifelse(criterion == "AIC", AIC(model), BIC(model))
  improved <- TRUE
  
  while (improved) {
    improved <- FALSE
    current_fixed_effects <- all.vars(update.formula(best_formula, . ~ . - (1|Subgroups)))  # Exclude random effects
    current_fixed_effects <- current_fixed_effects[!current_fixed_effects %in% response]
    
    # Test to eliminate each fixed effect one by one
    for (effect in current_fixed_effects) {
      candidate_fixed_effects <- setdiff(current_fixed_effects, effect)
      candidate_fixed_effects_formula <- paste(candidate_fixed_effects, collapse = " + ")
      candidate_formula <- as.formula(paste(response, "~", candidate_fixed_effects_formula, "+", random_effects))
      candidate_model <- glmer(candidate_formula, data = data, family = binomial)
      candidate_criterion <- ifelse(criterion == "AIC", AIC(candidate_model), BIC(candidate_model))
      
      # If the new model has a lower AIC, we actualize the model
      if (candidate_criterion < best_criterion) {
        best_model <- candidate_model
        best_formula <- candidate_formula
        best_criterion <- candidate_criterion
        improved <- TRUE
      }
    }
  }
  
  # Return the best model
  return(list(model = best_model, formula = best_formula, criterion = best_criterion))
}

# Go on with feature selection strategy
# Define multivariate model to be tested

# 2- For mitochondrial haplogroups & age
response <- "Severity"
fixed_effects <- c("Branch_HV", "I", "J", "U", "K", "age")
random_effects <- "(1|Subgroups)"


result.age <- backward_elimination(data, response, fixed_effects, random_effects, criterion = "AIC")

m0.age=glmer(Severity ~ Branch_HV + I + J + U + K + age + (1 | Subgroups),data = data, family = binomial)
m1.age=glmer(Severity ~ Branch_HV + age + (1 | Subgroups),data = data, family = binomial)
anova(m1.age,m0.age)
summary(m1.age)

# Explore possible interactions between selected features in the step before
m2.age=glmer(Severity ~ Branch_HV * age + (1 | Subgroups),data = data, family = binomial)
anova(m2.age,m1.age)

OR.tab <-data.frame(exp(cbind(fixef(m1.age),confint.merMod(m1.age,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

# 3- For mitochondrial haplogroups & sex

response <- "Severity"
fixed_effects <- c("Branch_HV", "I", "J", "U", "K", "sex")
random_effects <- "(1|Subgroups)"

result.sex <- backward_elimination(data, response, fixed_effects, random_effects, criterion = "AIC")


m0.sex=glmer(Severity ~ Branch_HV + I + J + U + K + sex + (1 | Subgroups),data = data, family = binomial)
m1.sex=glmer(Severity ~ Branch_HV + sex + (1 | Subgroups),data = data, family = binomial)
anova(m1.sex,m0.sex)
summary(m1.sex)

m2.sex=glmer(Severity ~ Branch_HV * sex + (1 | Subgroups),data = data, family = binomial)
anova(m1.sex,m2.sex)

OR.tab <-data.frame(exp(cbind(fixef(m1.sex),confint.merMod(m1.sex,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

# 4- For mitochondrial haplogroups & genetic background represented by the 10 Principal components obtained in PCA combuted over
# GWAS array information.
response <- "Severity"
fixed_effects <- c("Branch_HV", "I", "J", "U", "K", "euPC1", "euPC3", "euPC4", "euPC5", "euPC6", "euPC8", "euPC9", "euPC10")
random_effects <- "(1|Subgroups)"


result.gb <- backward_elimination(data, response, fixed_effects, random_effects, criterion = "AIC")

m0.gb=glmer(Severity ~ Branch_HV + I + J + U + K + euPC1 + euPC3 + euPC4 + euPC5 + euPC6 + euPC8 + euPC9 + euPC10 + (1 | Subgroups),data = data, family = binomial)
m1.gb=glmer(Severity ~ Branch_HV + euPC1 + euPC3 + (1 | Subgroups),data = data, family = binomial)
anova(m1.gb,m0.gb)
summary(m1.gb)

# Check relevance of selected features by comparing the performance of nested models.
m2.gb=glmer(Severity ~ Branch_HV * euPC1 +  Branch_HV * euPC3 + (1 | Subgroups),data = data, family = binomial)

summary(m2.gb)
anova(m1.gb,m2.gb)

m3.gb=glmer(Severity ~ Branch_HV * euPC1 + euPC3 + (1 | Subgroups),data = data, family = binomial)

summary(m3.gb)
anova(m3.gb,m2.gb)
anova(m3.gb,m1.gb)

m4.gb=glmer(Severity ~ Branch_HV : euPC1 + euPC1 + euPC3 + (1 | Subgroups),data = data, family = binomial)

summary(m4.gb)
anova(m4.gb,m3.gb)

m5.gb=glmer(Severity ~ Branch_HV + euPC1 + (1 | Subgroups),data = data, family = binomial)

summary(m5.gb)
anova(m5.gb,m3.gb)

OR.tab <-data.frame(exp(cbind(fixef(m3.gb),confint.merMod(m3.gb,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab

# Assemble general multivariate model, regarding the information obtained in partial multivariate models.


response <- "Severity"
fixed_effects <- c("Branch_HV", "euPC1", "euPC3", "sex", "age","Branch_HV:euPC1")
random_effects <- "(1|Subgroups)"

# Ejecutar la selección de características
result.all <- backward_elimination(data, response, fixed_effects, random_effects, criterion = "AIC")

m0.all=glmer(Severity ~ Branch_HV + euPC1 + Branch_HV:euPC1 + euPC3 + sex + age + (1 | Subgroups),data = data, family = binomial)
summary(m0.all)

m1.all=glmer(Severity ~ Branch_HV:euPC1 + euPC1 + euPC3 + sex + age + (1 | Subgroups),data = data, family = binomial)
anova(m1.all,m0.all)

m1.all=glmer(Severity ~ Branch_HV + euPC1 + euPC3 + sex + age + (1 | Subgroups),data = data, family = binomial)
anova(m1.all,m0.all)

m1.all=glmer(Severity ~ Branch_HV + euPC1 + euPC3 + sex + age + (1 | Subgroups),data = data, family = binomial)
anova(m1.all,m0.all)

m1.all=glmer(Severity ~ Branch_HV:euPC1 + Branch_HV:euPC1 +euPC3 + sex + age + (1 | Subgroups),data = data, family = binomial)
anova(m1.all,m0.all)

m1.all=glmer(Severity ~ Branch_HV:euPC1 + Branch_HV:euPC1 + euPC1 + sex + age + (1 | Subgroups),data = data, family = binomial)
anova(m1.all,m0.all)

m1.all=glmer(Severity ~ Branch_HV:euPC1 + Branch_HV:euPC1 + euPC1 + euPC3 + age + (1 | Subgroups),data = data, family = binomial)
anova(m1.all,m0.all)

m1.all=glmer(Severity ~ Branch_HV:euPC1 + Branch_HV:euPC1 + euPC1 + euPC3 + sex + (1 | Subgroups),data = data, family = binomial)
anova(m1.all,m0.all)


summary(m0.all)
OR.tab <-data.frame(exp(cbind(fixef(m0.all),confint.merMod(m0.all,method = "Wald")[-1,])))
names(OR.tab)=c('OR','Lower','Upper')
OR.tab
