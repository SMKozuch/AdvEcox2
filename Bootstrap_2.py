
# coding: utf-8

# # Question 1

# In[1]:


#load numpy and set seed
import numpy as np
np.random.seed(1429)

#setup
n = 50
REP = 1000
BOOTREP = 499
m = 0
s = 1.3
mu = np.exp(m + 0.5 * (s ** 2))
beta = np.sin(mu)

#create arrays for storing the variables
xbar = np.zeros([REP, 1])
bhat = np.zeros([REP, 1])
SE = np.zeros([REP, 1])
trat = np.zeros([REP, 1])
LCLasym = np.zeros([REP, 1])
UCLasym = np.zeros([REP, 1])

#Monte Carlo simulation for coverage frequency
for i in range(REP):
    if ((i + 1) % 500) == 0:
        print(i + 1)
    
    #generate X, calculate mean, SE via delta method, t-ratio and bhat
    X = np.exp(np.random.normal(loc=m, scale=s, size=([n, 1])))
    xbar[i] = np.mean(X)
    bhat[i] = np.sin(xbar[i])
    SE[i] = np.sqrt(np.var(X) / n * (np.cos(xbar[i])) ** 2)
    trat[i] = (bhat[i] - beta) / SE[i]
    
    #boundaries of CIs
    LCLasym[i] = bhat[i] - 1.96 * SE[i]
    UCLasym[i] = bhat[i] + 1.96 * SE[i]

#show the results
CoverageFreqasym = np.mean((beta > LCLasym) & (beta < UCLasym)) 
print('Coverage freq.(asym): \n', CoverageFreqasym)


# # Question 2

# In[2]:


#for loop for finding sufficient n
for n in range(50, 10000, 50):
    np.random.seed(1429)
    
    n = n
    REP = 1000
    BOOTREP = 499
    m = 0
    s = 1.3
    mu = np.exp(m + 0.5 * s ** 2)
    beta = np.sin(mu)

    xbar = np.zeros([REP, 1])
    bhat = np.zeros([REP, 1])
    SE = np.zeros([REP, 1])
    trat = np.zeros([REP, 1])
    LCLasym = np.zeros([REP, 1])
    UCLasym = np.zeros([REP, 1])

    for i in range(REP):
        X = np.exp(np.random.normal(loc=m, scale=s, size=([n, 1])))
        xbar[i] = np.mean(X)
        bhat[i] = np.sin(xbar[i])
        SE[i] = np.sqrt(np.var(X) / n * (np.cos(xbar[i])) ** 2)
        trat[i] = (bhat[i] - beta) / SE[i]

        LCLasym[i] = bhat[i] - 1.96 * SE[i]
        UCLasym[i] = bhat[i] + 1.96 * SE[i]

    CoverageFreqasym = np.mean((beta > LCLasym) & (beta < UCLasym)) 
    
    #break the loop when n is sufficiently large
    if ((CoverageFreqasym < 0.9635) & (CoverageFreqasym > 0.9365)):
        print('Minimal sample size: ', n)
        print('Frequency: ', CoverageFreqasym)
        break


# # Question 3

# #### Bias and bias corrected estimates

# In[3]:


#Bootstrap
np.random.seed(1429)
n = 50

#generate sample
X = np.exp(np.random.normal(loc=m, scale=s, size=([n, 1])))

#Bootstrap
bhat_bs = np.zeros([BOOTREP, 1])
for b in range(BOOTREP):
    bhat_bs[b] = np.sin(np.mean(X[np.random.randint(0, n - 1, n)]))
    
#Bootstrap bias
bias_bs = np.mean(bhat_bs) - np.sin(np.mean(X))
bias_bs_correction = 2 * np.sin(np.mean(X)) - np.mean(bhat_bs)

print('Bootstrap bias: ', bias_bs)
print('Bootstrap bias-corrected: ', bias_bs_correction)

#Jackknife
bhat_jk = np.zeros([n, 1])

for i in range(n):
    mask = np.ones(n, dtype=bool)
    mask[i] = 0
    bhat_jk[i] = np.sin(np.mean(X[mask]))
    
bias_jk = (n - 1) * (np.mean(bhat_jk) - np.sin(np.mean(X)))
bias_jk_correction = n * np.sin(np.mean(X)) - (n - 1) * np.mean(bhat_jk)

#Jackknife bias
print('Jackknife bias: ', bias_jk)
print('Jackknife bias-corrected: ', bias_jk_correction)

#evaluation
print('True beta is:', beta)
print('=> JK is closer to the true value.')


# #### Confidence intervals using the respective estimators

# In[4]:


#load numpy and set seed
import numpy as np
np.random.seed(1429)

#setup
n = 50
REP = 1000
BOOTREP = 499
m = 0
s = 1.3
mu = np.exp(m + 0.5 * (s ** 2))
beta = np.sin(mu)

#create arrays for storing the variables
xbar = np.zeros([REP, 1])
bhat_bs = np.zeros([REP, 1])
bhat_jk = np.zeros([REP, 1])
SE = np.zeros([REP, 1])
LCLasym_bs = np.zeros([REP, 1])
UCLasym_bs = np.zeros([REP, 1])
LCLasym_jk = np.zeros([REP, 1])
UCLasym_jk = np.zeros([REP, 1])

for i in range(REP):
    if ((i + 1) % 500) == 0:
        print(i + 1)
    
    #generate X, calculate mean and SE via delta method
    X = np.exp(np.random.normal(loc=m, scale=s, size=([n, 1])))
    xbar[i] = np.mean(X)
    SE[i] = np.sqrt((np.var(X) / n) * (np.cos(xbar[i])) ** 2)
    
    #bootstrap
    bhat_bs_rep = np.zeros([BOOTREP, 1])
    for b in range(BOOTREP):
        bhat_bs_rep[b] = np.sin(np.mean(X[np.random.randint(0, n - 1, n)]))

    # Bias corrected BS bhat
    bhat_bs[i] = 2 * np.sin(np.mean(X)) - np.mean(bhat_bs_rep)
    
    # Bias BS corrected confidence intervals
    LCLasym_bs[i] = bhat_bs[i] - 1.96 * SE[i]
    UCLasym_bs[i] = bhat_bs[i] + 1.96 * SE[i]
  
    #Jack-knife
    bhat_jk_rep = np.zeros([n, 1])
    for b in range(n):
        mask = np.ones(n, dtype=bool)
        mask[b] = 0
        bhat_jk_rep[b] = np.sin(np.mean(X[mask]))
    
    # Bias corrected JK bhat
    bhat_jk[i] = n * np.sin(np.mean(X)) - (n - 1) * np.mean(bhat_jk_rep)
     
    # Bias JK corrected confidence intervals
    LCLasym_jk[i] = bhat_jk[i] - 1.96 * SE[i]
    UCLasym_jk[i] = bhat_jk[i] + 1.96 * SE[i]

#print results    
CoverageFreqasym_bs = np.mean((beta > LCLasym_bs) & (beta < UCLasym_bs)) 
print('Coverage freq.(asym) for bootstrap corrected is: ', CoverageFreqasym_bs)
    
CoverageFreqasym_jk = np.mean((beta > LCLasym_jk) & (beta < UCLasym_jk)) 
print('Coverage freq.(asym) for jackknife corrected is: ', CoverageFreqasym_jk)


# # Question 4

# #### Confidence intervals using the respective SE estimates

# In[5]:


#load numpy and set seed
import numpy as np
np.random.seed(1429)

#setup
n = 50
REP = 1000
BOOTREP = 499
m = 0
s = 1.3
mu = np.exp(m + 0.5 * (s ** 2))
beta = np.sin(mu)

#create arrays for storing the variable
xbar = np.zeros([REP, 1])
bhat_bs = np.zeros([REP, 1])
bhat_jk = np.zeros([REP, 1])
SE_bs = np.zeros([REP, 1])
SE_jk = np.zeros([REP, 1])
LCLasym_bs = np.zeros([REP, 1])
UCLasym_bs = np.zeros([REP, 1])
LCLasym_jk = np.zeros([REP, 1])
UCLasym_jk = np.zeros([REP, 1])

#Monte Carlo
for i in range(REP):
    if ((i + 1) % 500) == 0:
        print(i + 1)
    
    #generate X, obtain mean and bhat
    X = np.exp(np.random.normal(loc=m, scale=s, size=([n, 1])))
    xbar[i] = np.mean(X)
    bhat[i] = np.sin(xbar[i])
    
    #create array for storing the variable
    bhat_jk = np.zeros([n, 1])
    
    #jackknife
    for j in range(n):
        mask = np.ones(n, dtype=bool)
        mask[j] = 0
        bhat_jk[j] = np.sin(np.mean(X[mask]))
    
    #jackknife estimate of SE
    SE_jk[i] = np.sqrt((n - 1) / n * np.sum((bhat_jk - np.mean(bhat_jk)) ** 2))
    
    
    LCLasym_jk[i] = bhat[i] - 1.96 * SE_jk[i]
    UCLasym_jk[i] = bhat[i] + 1.96 * SE_jk[i]
    
    bhat_bs = np.zeros([BOOTREP, 1])

    for b in range(BOOTREP):
        bhat_bs[b] = np.sin(np.mean(X[np.random.randint(0, n - 1, n)]))    
    
    #bootstrap estimate of SE
    SE_bs[i] = np.sqrt((1 / (BOOTREP - 1)) * np.sum((bhat_bs - np.mean(bhat_bs)) ** 2))
    
    
    LCLasym_bs[i] = bhat[i] - 1.96 * SE_bs[i]
    UCLasym_bs[i] = bhat[i] + 1.96 * SE_bs[i]
    
CoverageFreqasym_bs = np.mean((beta > LCLasym_bs) & (beta < UCLasym_bs)) 
print('Coverage freq.(asym) bootstrap corrected SE is: ', CoverageFreqasym_bs)
    
CoverageFreqasym_jk = np.mean((beta > LCLasym_jk) & (beta < UCLasym_jk)) 
print('Coverage freq.(asym) jackknife corrected SE is: ', CoverageFreqasym_jk)


# # Question 5

# #### percentile-$t$ and percentile method

# In[6]:


import numpy as np
np.random.seed(1211)

#setup
n = 50
REP = 1000
BOOTREP = 499
m = 0
s = 1.3
mu = np.exp(m + 0.5 * (s ** 2))
beta = np.sin(mu)

#create arrays for storing the variables
xbar = np.zeros([REP, 1])
bhat = np.zeros([REP, 1])
SE = np.zeros([REP, 1])
trat = np.zeros([BOOTREP, 1])
LCLasym_pert = np.zeros([REP, 1])
UCLasym_pert = np.zeros([REP, 1])
LCLasym_per = np.zeros([REP, 1])
UCLasym_per = np.zeros([REP, 1])
LCL_pert = np.zeros([REP, 1])
UCL_pert = np.zeros([REP, 1])

#Monte Carlo
for i in range(REP):
    if ((i + 1) % 500) == 0:
        print(i + 1)
    
    #generate X, obtain SE, mean and bhat
    X = np.exp(np.random.normal(loc=m, scale=s, size=([n, 1])))
    xbar[i] = np.mean(X)
    bhat[i] = np.sin(xbar[i])
    SE[i] = np.sqrt((np.var(X) / n) * (np.cos(xbar[i])) ** 2)
    
    #create arrays for storing the variables
    bhat_bs = np.zeros([BOOTREP,1])
    SE_bs = np.zeros([BOOTREP,1])
    
    #bootstrap
    for b in range(BOOTREP):
        index = np.random.randint(0, n - 1, n)
        bhat_bs[b] = np.sin(np.mean(X[index]))
        SE_bs[b] =  np.sqrt((np.var(X[index]) / n) * (np.cos(np.mean(X[index]))) ** 2)
        trat[b] = (bhat_bs[b] - bhat[i]) / (SE_bs[b])
    
    #bootstraps confidence intervals for percentile-t statistics
    LCL_pert[i] = np.percentile(trat, q=97.5, interpolation='higher')
    UCL_pert[i] = np.percentile(trat, q=2.5, interpolation='higher')
    
    #asymptotic confidence intervals for percentile-t statistics
    LCLasym_pert[i] = bhat[i] - LCL_pert[i] * SE[i]
    UCLasym_pert[i] = bhat[i] - UCL_pert[i] * SE[i]
    
    #bootstrap confidence intervals for percentile statistics
    LCLasym_per[i] = np.percentile(bhat_bs, q=2.5, interpolation='higher')
    UCLasym_per[i] = np.percentile(bhat_bs, q=97.5, interpolation='higher')

 
CoverageFreqasym_pert = np.mean((beta > LCLasym_pert) & (beta < UCLasym_pert))
CoverageFreqasym_per = np.mean((beta > LCLasym_per) & (beta < UCLasym_per))

#print results
print('Coverage Freq. Percentile-t: ', CoverageFreqasym_pert)
print('Coverage Freq. Percentile: ', CoverageFreqasym_per)


# # Question 6

# #### percentile-$t$ and percentile, moved center

# In[4]:


import numpy as np
np.random.seed(1429)

#setup
n = 50
REP = 1000
BOOTREP = 499
m = 0
s = 1.3
mu = np.exp(m + 0.5 * (s ** 2))
beta = np.sin(mu)

xbar = np.zeros([REP, 1])
bhat = np.zeros([REP, 1])
SE = np.zeros([REP, 1])
trat = np.zeros([BOOTREP, 1])
LCLasym_pert = np.zeros([REP, 1])
UCLasym_pert = np.zeros([REP, 1])
LCLasym_per = np.zeros([REP, 1])
UCLasym_per = np.zeros([REP, 1])
LCL_pert = np.zeros([REP, 1])
UCL_pert = np.zeros([REP, 1])

# Monte Carlo
for i in range(REP):
    if ((i + 1) % 500) == 0:
        print(i + 1)

    #generate X, obtain mean, bhat and SE
    X = np.exp(np.random.normal(loc=m, scale=s, size=([n, 1])))
    xbar[i] = np.mean(X)
    bhat[i] = np.sin(xbar[i])
    SE[i] = np.sqrt((np.var(X) / n) * (np.cos(xbar[i])) ** 2)
       
    bhat_bs = np.zeros([BOOTREP,1])
    SE_bs = np.zeros([BOOTREP,1])
    
    #bootstrap
    for b in range(BOOTREP):
        index = np.random.randint(0, n - 1, n)
        bhat_bs[b] = np.sin(np.mean(X[index]))
        SE_bs[b] =  np.sqrt((np.var(X[index]) / n) * (np.cos(np.mean(X[index]))) ** 2)
        trat[b] = (bhat_bs[b] - bhat[i]) / (SE_bs[b])
        
    #boundaries of confidence intervals for percentile-t statistics
    LCL_pert[i] = np.percentile(trat, q=97.5, interpolation='higher')
    UCL_pert[i] = np.percentile(trat, q=2.5, interpolation='higher')
     
    bias_bs_correction = 2 * bhat[i] - np.mean(bhat_bs)
    
    #bootstrap bias corrected confidence intervals for percentile-t statistics
    LCLasym_pert[i] = bias_bs_correction - LCL_pert[i] * SE[i]
    UCLasym_pert[i] = bias_bs_correction - UCL_pert[i] * SE[i]
    
    #bootstrap bias corrected confidence intervals for percentile statistics
    LCLasym_per[i] = np.percentile(bhat_bs, q=2.5, interpolation='higher')
    UCLasym_per[i] = np.percentile(bhat_bs, q=97.5, interpolation='higher')

 

CoverageFreqasym_pert = np.mean((beta > LCLasym_pert) & (beta < UCLasym_pert))
CoverageFreqasym_per = np.mean((beta > LCLasym_per) & (beta < UCLasym_per))

#print results
print('Coverage Freq. Percentile-t: ', CoverageFreqasym_pert)
print('Coverage Freq. Percentile: ', CoverageFreqasym_per)


# # Question 7

# #### $t$-ratio with bias correction

# In[8]:


import numpy as np
np.random.seed(1429)

#setup
n = 50
REP = 1000
BOOTREP = 499
m = 0
s = 1.3
mu = np.exp(m + 0.5 * (s ** 2))
beta = np.sin(mu)

xbar = np.zeros([REP, 1])
bhat = np.zeros([REP, 1])
SE = np.zeros([REP, 1])
LCLasym_pert = np.zeros([REP, 1])
UCLasym_pert = np.zeros([REP, 1])
LCL_pert = np.zeros([REP, 1])
UCL_pert = np.zeros([REP, 1])
trat = np.zeros([BOOTREP, 1])


for i in range(REP):
    if ((i + 1) % 100) == 0:
        print(i + 1)
    
    #Created the sample + calculated the SE, bhat
    X = np.exp(np.random.normal(loc=m, scale=s, size=([n, 1])))
    xbar[i] = np.mean(X)
    bhat[i] = np.sin(xbar[i])
    SE[i] = np.sqrt((np.var(X) / n) * (np.cos(xbar[i])) ** 2)
    
    bhat_bs = np.zeros([BOOTREP, 1])
    SE_bs = np.zeros([BOOTREP, 1])
    bhat_jk = np.zeros([n, 1])

    for b in range(BOOTREP):
        
        # First, we bootrstrapped the sample on which we calculated the 
        # necessary bhat_bs and the standard error. 
        
        X_bs = X[np.random.randint(0, n - 1, n)]
        bhat_bs[b] = np.sin(np.mean(X_bs))
        SE_bs[b] = np.sqrt((np.var(X_bs) / n) * (np.cos(np.mean(X_bs))) ** 2)
        
        # Using the jackknife bootstrap, we obtain the jackknife bias estimate
        
        for j in range(n):
            
            mask = np.ones(n, dtype=bool)
            mask[j] = 0
            bhat_jk[j] = np.sin(np.mean(X_bs[mask]))

        bias_jk = (n - 1) * (np.mean(bhat_jk) - bhat[i])
        
        # Using the bias estimate and bhat, bhat_bs, we can calculate the t-ratios
        
        trat[b] = (bhat_bs[b] - bias_jk - bhat[i]) / (SE_bs[b])

        
    # Taking the percentile of the t-ratios will give us the lower / upper bounds
    LCL_pert[i] = np.percentile(trat, q=97.5, interpolation='higher')
    UCL_pert[i] = np.percentile(trat, q=2.5, interpolation='higher')
        
    LCLasym_pert[i] = bhat[i] - LCL_pert[i] * SE[i]
    UCLasym_pert[i] = bhat[i] - UCL_pert[i] * SE[i]
 

CoverageFreqasym_pert = np.mean((beta > LCLasym_pert) & (beta < UCLasym_pert))
print('Coverage Freq. Percentile-t: ', CoverageFreqasym_pert)

