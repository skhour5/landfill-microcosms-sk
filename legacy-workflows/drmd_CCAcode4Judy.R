Step 12c: Correspondance analysis
  Correspondence analysis is basically the same as doing chi-square distances and then applying PCoA to it. But there's another layer that needs discussion. Chi-square distance is good at retrieving gradients because correspondence analysis as a technique is guaranteed to recover gradient trends if certain conditions are met. Mainly, this condition is when each species has a unimodal response to some environmental gradient, unimodal means each species has a peak somewhere along the gradient of conditions (so temperature...each species will have a peak at a temeperature). If the condition is met then correspondence analysis will figure out what order the samples need to go in, in order to best recover the unimodal response curves of the species. It's putting the samples in order AND putting the species in order (figuring out where the mode is) along that gradient. This results in an arch trend. 
```{r}

### Unconstrained INDIRECT gradient analysis - just ordinating, method doesn't know anything about environ data.

# We will use the vegan package to run correspondence analysis. We can also plot a biplot using vegan by calling plot() on the resulting CA object. [this is the same as what we used for the chi-squared ordination in video 8]
# Use relative abundances for this.
# run CA using vegan command
my.ca <- cca(vst_trans_count_tab.rel) #vst_trans_count_tab.rel is rows= sample names; columns = ASV#s; data= relative counts
plot(my.ca)
#samples are in black, ASVs are in red.

# What fraction of total inertia is explained by each axis? (How much variance is explained by each axis? )
my.ca$CA$eig/my.ca$tot.chi
#this is how people put the % on the axes that say how much of the pattern can be explained by that axis.

### DIRECT gradient analysis: Canonical Correspondence Analysis (CCA) [aka Constrained Correspondence Analysis]. in which we relate species directly to environmental variables. We'll pass in covariates that we think are driving the community structure. We give it environmental variables that we think relate to the species. You can explicitly test variables.
# NOTE: THERE IS ANOTHER CCA = CONONICAL CORRELATION ANALYSIS, WHICH IS SOMETHING TOTALLY DIFFERENT
#1] convert OTUs to chi-square distances
#2] perform CA to get sample scores on each axis
#3] regress these optimal CA coordinates on the environmental variables, which gives us constrained axis scores
#4] use fitted coordinates as CCA coordinates

#There are environmental variables in this example. DAN KNIGHTS RECOMMENDS NOT USING MORE THAN 3-4 VARIABLES AT ONCE. MORE THAN THAT, AND YOU WILL START GETTING ARITFACTS AND JUST OVERFIT THE DATA.
# run CA using vegan command
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear) # CHANGE 'VAR' to your column names in your metadata file. "samdf.pear" is my metadata table

#we're telling it that the OTU table varies as a function of these variables

plot(my.cca)

#HEY. This produces one of those ordination plots with the environmental vectors on it!
#This is a triplot: samples in black, OTUs in red, metadata in blue
#arrows give you the direction and magnitude of the regression coefficient
#An example: if the arrows for chemotaxis and depth go in different directions, it means that chemotaxis genes are higher in shallow samples.

#What fraction of total inertia is explained by each axis in CCA? Compare this to the fraction of total inertia explained by CA (it should be lower than CA).
my.cca$CCA$eig/my.cca$tot.chi
#gives you % that each axis explains

#how to test the significance of this? Well, you can state the fraction of the variance explained by CA that is also explained by CCA. (in other words, what did CA catch that CCA didn't?)
#divide the smaller value by the larger. e.g., if CCA only had axis one as 31%, and CA had axis one as 37%, then 31/37= 84% of the CA analysis was explained by CCA for that axis. 

#run the next 3 lines of code and record the value
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi

# Test what fraction of CA1 variance is explained in CCA1
b[1]/a[1]

#To directly test if CCA variance is significant
#We need to think about the Null hypothesis (which I think here is that 'there is no gradient in the data'). Use a Monte Carlo simulation to decide if CCA gives a significant answer for this dataset. We can simulate random data by shuffling or permuting the metadata values. We will shuffle them together to preserve correlations between metadata variables. If we shuffle them 10,000 times and calculate the variance explained in CCA axis 1 each time, we can compare this to the observed variation explained to get a p-value. What fraction of the time does our randomized dataset expalin the variance better than or equal to our actual ordered dataset?
#The result is a 'CCA with permuted metadata values.'
# store the observed value
obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi

# Perform 999 randomized CCAs
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})

# include the observed value as one of the "null" values to be conservative
mc.vals <- c(mc.vals, obs.val)

# What fraction of the randomized values was greater than the observed value?
# this is the p-value
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])
plot(my.cca)

#interestingly this plot is different than the one made near the begining of this chunk. 


#Playing aorund with making the plots the way I want them
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)

plot(my.cca, display=c("sites")) #this will make a plot with the sample locations listed on it, but it lacks all the other useful features like the ASVs and vectors. But at least you know which point is which!
plot(my.cca, display=c("sites","cn")) #this will make a plot with both sample locations named and the env vectors. No ASVs. In theory, "sp" added to the display term should do this, but when I add it the sample names disappear and I get a plot that's identical to "plot(my.cca)"
plot(my.cca, display=c("sites","bp")) # substrates are vectors here. No ASVs.
plot(my.cca, display=c("sites","bp","sp")) #Better! 3/4 of the substrates are on here, and the DO, and the ASVs! No site names, oddly. Luckily I know where the 4th substrate goes...


#Help at: rfunctions.blogspot.com/2016/11/canonical-correspondence-analysis-cca.html

my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)
plot(my.cca, display=c("sites","bp","sp"))
plot(my.cca, display=c("sites","cn"))

my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])
plot(my.cca, display=c("sites","bp","sp"))
plot(my.cca, display=c("sites","cn"))


###################  BULK TESTING, HERE ARE 10 RUNS YOU CAN EDIT AND JUST GO. It does not plot the second CCA (that explores the randomized data to see how different it is to estimate error) because I couldn't see the value in producing that plot.

## 1
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)  
plot(my.cca, display=c("sites","bp","sp"))
my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])

## 2
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)  
plot(my.cca, display=c("sites","bp","sp"))
my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])

## 3
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)  
plot(my.cca, display=c("sites","bp","sp"))
my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])

## 4
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)  
plot(my.cca, display=c("sites","bp","sp"))
my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])

## 5
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)  
plot(my.cca, display=c("sites","bp","sp"))
my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])

## 6
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)  
plot(my.cca, display=c("sites","bp","sp"))
my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])

## 7
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)  
plot(my.cca, display=c("sites","bp","sp"))
my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])

## 8
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)  
plot(my.cca, display=c("sites","bp","sp"))
my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])

## 9
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)  
plot(my.cca, display=c("sites","bp","sp"))
my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])

## 10
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear)  
plot(my.cca, display=c("sites","bp","sp"))
my.cca$CCA$eig/my.cca$tot.chi
a <- my.ca$CA$eig/my.ca$tot.chi
b <- my.cca$CCA$eig/my.cca$tot.chi
b[1]/a[1]

obs.val <- my.cca$CCA$eig[1]/my.cca$tot.chi
mc.vals <- replicate(999,{my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),]); my.cca$CCA$eig[1]/my.cca$tot.chi})
mc.vals <- c(mc.vals, obs.val)
mean(c(obs.val, mc.vals) >= obs.val)
my.cca <- cca(vst_trans_count_tab.rel ~ VAR1 + VAR2, data=samdf.pear[sample(1:nrow(samdf.pear)),])



```