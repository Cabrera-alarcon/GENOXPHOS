####################################################
### HG risk factors analysis in SCOURGE patients ###
####################################################

library(MASS)

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
data$p_EUR=scale(data$p_EUR,center = min(data$p_EUR),scale = max(data$p_EUR) - min(data$p_EUR))
data$p_EAS=scale(data$p_EAS,center = min(data$p_EAS),scale = max(data$p_EAS) - min(data$p_EAS))
data$p_AMR=scale(data$p_AMR,center = min(data$p_AMR),scale = max(data$p_AMR) - min(data$p_AMR))
data$p_AFR=scale(data$p_AFR,center = min(data$p_AFR),scale = max(data$p_AFR) - min(data$p_AFR))

# Select patients from the SCOURGE cohort (cases), who also have all the comorbidities data, 8,788 out of 8,894 patients.
df2=data[data$Subgroups==1 & !is.na(data$AR),]

# Univariate analysis to perform feature selection taking into account population stratification

# Mitochondrial haplogroup assesment

# Fit logistic regression model to assess regressor
m.HV=glm(paste0('Severity~','Branch_HV'), data = df2,family = 'binomial')
summary(m.HV)

# Calculate Odds Ratio and 95% CI
exp(coef(m.HV))
exp(confint(m.HV))

m.U=glm(paste0('Severity~','U'), data = df2,family = 'binomial')
summary(m.U)
exp(coef(m.U))
exp(confint(m.U))

m.K=glm(paste0('Severity~','K'), data = df2,family = 'binomial')
summary(m.K) 

exp(coef(m.K))
exp(confint(m.K))

m.T=glm(paste0('Severity~','T'), data = df2,family = 'binomial')
summary(m.T)

exp(coef(m.T))
exp(confint(m.T))

m.J=glm(paste0('Severity~','J'), data = df2,family = 'binomial')
summary(m.J)
exp(coef(m.J))
exp(confint(m.J))

m.X=glm(paste0('Severity~','X'), data = df2,family = 'binomial')
summary(m.X)
exp(coef(m.X))
exp(confint(m.X))

m.I=glm(paste0('Severity~','I'), data = df2,family = 'binomial')
summary(m.I)
exp(coef(m.I))
exp(confint(m.I))

m.W=glm(paste0('Severity~','W'), data = df2,family = 'binomial')
summary(m.W)
exp(coef(m.W))
exp(confint(m.W))


# Comorbidities assessment
# Vascular history
m.AV=glm(paste0('Severity~','AV'), data = df2,family = 'binomial')
summary(m.AV)
exp(coef(m.AV))
exp(confint(m.AV))

# Cardiac history
m.AC=glm(paste0('Severity~','AC'), data = df2,family = 'binomial')
summary(m.AC)
exp(coef(m.AC))
exp(confint(m.AC))

# Neurological disease history
m.AN=glm(paste0('Severity~','AN'), data = df2,family = 'binomial')
summary(m.AN) 
exp(coef(m.AN))
exp(confint(m.AN))

# Digestive disease history
m.AD=glm(paste0('Severity~','AD'), data = df2,family = 'binomial')
summary(m.AD)
exp(coef(m.AD))
exp(confint(m.AD))

# Respiratory disease history
m.AR=glm(paste0('Severity~','AR'), data = df2,family = 'binomial')
summary(m.AR)
exp(coef(m.AR))
exp(confint(m.AR))

# Onco-hematologic history
m.AOH=glm(paste0('Severity~','AOH'), data = df2,family = 'binomial')
summary(m.AOH)
exp(coef(m.AOH))
exp(confint(m.AOH))


# Nuclear background assessment (All SCOURGE patients)
m.pc1=glm(paste0('Severity~','euPC1'), data = df2,family = 'binomial')
summary(m.pc1)
exp(coef(m.pc1))
exp(confint(m.pc1))

m.pc2=glm(paste0('Severity~','euPC2'), data = df2,family = 'binomial')
summary(m.pc2)
exp(coef(m.pc2))
exp(confint(m.pc2))

m.pc3=glm(paste0('Severity~','euPC3'), data = df2,family = 'binomial')
summary(m.pc3)
exp(coef(m.pc3))
exp(confint(m.pc3))

m.pc4=glm(paste0('Severity~','euPC4'), data = df2,family = 'binomial')
summary(m.pc4)
exp(coef(m.pc4))
exp(confint(m.pc4))

m.pc5=glm(paste0('Severity~','euPC5'), data = df2,family = 'binomial')
summary(m.pc5)
exp(coef(m.pc5))
exp(confint(m.pc5))

m.pc6=glm(paste0('Severity~','euPC6'), data = df2,family = 'binomial')
summary(m.pc6)
exp(coef(m.pc6))
exp(confint(m.pc6))

m.pc7=glm(paste0('Severity~','euPC7'), data = df2,family = 'binomial')
summary(m.pc7)
exp(coef(m.pc7))
exp(confint(m.pc7))

m.pc8=glm(paste0('Severity~','euPC8'), data = df2,family = 'binomial')
summary(m.pc8)
exp(coef(m.pc8))
exp(confint(m.pc8))

m.pc9=glm(paste0('Severity~','euPC9'), data = df2,family = 'binomial')
summary(m.pc9)
exp(coef(m.pc9))
exp(confint(m.pc9))

m.pc10=glm(paste0('Severity~','euPC10'), data = df2,family = 'binomial')
summary(m.pc10)
exp(coef(m.pc10))
exp(confint(m.pc10))

# Age assessment (All SCOURGE patients)
m.age=glm(paste0('Severity~','age'), data = df2,family = 'binomial')
summary(m.age)
exp(coef(m.age))
exp(confint(m.age))

# Sex assessment (All SCOURGE patients)
m.sex=glm(paste0('Severity~','sex'), data = df2,family = 'binomial')
summary(m.sex)
exp(coef(m.sex))
exp(confint(m.sex))

# Go on with feature selection strategy

# 1- For mitochondrial haplogroups & comorbidities

# Define multivariate model to be tested
m.comorbidities=glm(paste0('Severity~','Branch_HV + AV + AOH + AN + AR + AD'), data = df2,family = 'binomial')
summary(m.comorbidities)

# Perform feature selection following a step-wise backward strategy of AIC reduction.
m.comorbidities2=MASS::stepAIC(m.comorbidities)
summary(m.comorbidities2)

# Check result of feature selection by comparing nested models
anova(m.comorbidities2,m.comorbidities)

exp(coef(m.comorbidities2))
exp(confint(m.comorbidities2))

# Delve in feature selection by comparing nested models 

m.comorbidities3=glm(paste0('A1SEV_4~','Branch_HV + AN + AR'), data = df2,family = 'binomial')
anova(m.comorbidities3,m.comorbidities2)
m.comorbidities3=glm(paste0('A1SEV_4~','Branch_HV + AV + AR'), data = df2,family = 'binomial')
anova(m.comorbidities3,m.comorbidities2)
m.comorbidities3=glm(paste0('A1SEV_4~','Branch_HV + AV + AN'), data = df2,family = 'binomial')
anova(m.comorbidities3,m.comorbidities2)
m.comorbidities3=glm(paste0('A1SEV_4~','AR + AV + AN'), data = df2,family = 'binomial')
anova(m.comorbidities3,m.comorbidities2)

# 2- For mitochondrial haplogroups & nuclear background
m.gb=glm(paste0('Severity~','Branch_HV + euPC1 + euPC2'), data = df2,family = 'binomial')
summary(m.gb)

m.gb2=MASS::stepAIC(m.gb)

summary(m.gb2)
anova(m.gb2,m.gb)
exp(coef(m.gb2))
exp(confint(m.gb2))

m.gb3=glm(paste0('Severity~','euPC1'), data = df2,family = 'binomial')
anova(m.gb3,m.gb2)
m.gb3=glm(paste0('Severity~','Branch_HV'), data = df2,family = 'binomial')
anova(m.gb3,m.gb2)
m.gb3=glm(paste0('Severity~','Branch_HV:euPC1'), data = df2,family = 'binomial')
anova(m.gb3,m.gb2)

# 3- For mitochondrial haplogroups & Sex
m.sex1=glm(paste0('Severity~','Branch_HV + sex '), data = df2,family = 'binomial')
summary(m.sex1)
exp(coef(m.sex1))
exp(confint(m.sex1))
# The HV branch was not selected in this multivariate model. Thus we wanted to explore the behavior of the feature Sex in the 
# SCOURGE cohort.

# Study Sex bias
fisher.test(table(df2$sex, df2$Branch_HV)) # Male bias

# Explore HV branch Sex interaction
m.sex2=glm(paste0('Severity~','Branch_HV : sex '), data = df2,family = 'binomial')
summary(m.sex2)
exp(coef(m.sex2))
exp(confint(m.sex2))

anova(m.sex2, m.sex1)

# 4- For mitochondrial haplogroups & Age

m.age=glm(paste0('Severity~','Branch_HV + age '), data = df3,family = 'binomial')
summary(m.age)
exp(coef(m.age))
exp(confint(m.age))

m.age1=glm(paste0('Severity~','Branch_HV + age '), data = df3,family = 'binomial')
m.age2=glm(paste0('Severity~','age'), data = df3,family = 'binomial')
anova(m.age2,m.age1)
m.age2=glm(paste0('Severity~','Branch_HV'), data = df3,family = 'binomial')
anova(m.age2,m.age1)
m.age2=glm(paste0('Severity~','Branch_HV:age'), data = df3,family = 'binomial')
anova(m.age2,m.age1)


# Assemble general multivariate model, regarding the information obtained in partial multivariate models, following the same 
# feature selection strategy.

global.model=glm(paste0('Severity~','Branch_HV + euPC1 + AV + AN + AR+age+ sex +  Branch_HV:euPC1'), data = df2,family = 'binomial')
summary(global.model)
global.model2=MASS::stepAIC(global.model) # Starting model is the best model
summary(global.model)
global.model3=glm(paste0('Severity~','euPC1 + AV + AN + AR+age+ sex + Branch_HV:euPC1'), data = df2,family = 'binomial')
anova(global.model2,global.model3)
global.model3=glm(paste0('Severity~','Branch_HV + euPC1 + AV + AN + AR + age+ sex' ), data = df2,family = 'binomial')
anova(global.model2,global.model3)
global.model3=glm(paste0('Severity~','Branch_HV + euPC1 + AV + AN + AR + age + Branch_HV:euPC1' ), data = df2,family = 'binomial')
anova(global.model2,global.model3)
global.model3=glm(paste0('Severity~','Branch_HV + euPC1 + AV + AN + AR + sex + Branch_HV:euPC1' ), data = df2,family = 'binomial')
anova(global.model2,global.model3)
global.model3=glm(paste0('Severity~','Branch_HV + euPC1 + AV + AN + age + sex + Branch_HV:euPC1' ), data = df2,family = 'binomial')
anova(global.model2,global.model3)
global.model3=glm(paste0('Severity~','Branch_HV + euPC1 + AV + AR + age + sex + Branch_HV:euPC1' ), data = df2,family = 'binomial')
anova(global.model2,global.model3)
global.model3=glm(paste0('Severity~','Branch_HV + euPC1 + AN + AR + age + sex + Branch_HV:euPC1' ), data = df2,family = 'binomial')
anova(global.model2,global.model3)
global.model3=glm(paste0('Severity~','Branch_HV + AR + AV + AN + age + sex + Branch_HV:euPC1' ), data = df2,family = 'binomial')
anova(global.model2,global.model3)


exp(coef(global.model2))
exp(confint(global.model2))


