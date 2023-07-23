#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import copy
from scipy.stats import truncnorm


# ## generate pools
# 
# for one pool(100 pools) with ns(5000) species, all species are not anomolous species. (consider add anomolous species) 
# number of resources: nr ( nr = 1,2,...,7)
# 
# we need generate the growth rate $g_{\mu i}$. $\sum_i g_{\mu i}^2=1 \mathrm{hr}^{-2}$

# In[2]:


def generate_anomolous_species(ns, nr):
    vector = []
    for _ in range(ns):
        vector = np.random.uniform(0,1,size=nr)
        vector /= np.linalg.norm(vector)
        vectors.append(vector)
    return vectors


# In[3]:


def generate_anomolous_nutrient_rank(ns, nr):
    vectors = []
    for _ in range(ns):
        vector = random.sample(range(nr), nr)
        vectors.append(vector)
    return vectors


# In[4]:


def generate_random_vectors(ns, nr):
    vectors = []
    ranks = []

    for _ in range(ns):
        vector = np.random.uniform(0,1,size=nr)
        vector /= np.linalg.norm(vector)
        vectors.append(vector)

        # Generate ranks
        rank = np.argsort(-vector)
        ranks.append(rank)

    return vectors, ranks

##the rank should be used like maxrate-vectors[rank[0]], secondrate-vectors[rank[1]]...


# In[5]:


# generate_anomolous growth rate & rank
'''
spciesGR = []
rank = []
for i in range(ns):
    spciesGR.append(generate_anomolous_species(ns, nr))
    rank.append(generate_anomolous_nutrient_rank(ns, nr))
'''


# In[6]:


def delta_n(n0,gr,t):
    # n0,gr,t are all numbers
    try:
        delta_n = n0 * math.exp(gr * t) - n0
        #if t<0.000001:
           # print("gr = %s"%(gr))
           # print("t = %s"%(t))
           # print("math.exp(gr * t) = %s"%(math.exp(gr * t)))
            #print("n0 = %s"%(n0))
            #print("delta_n = %s"%(delta_n))
        return delta_n
    except OverflowError:
        print("OverflowError: Input value is too large for exponential calculation.")
        print(gr)
        print(t)
        return None


# In[7]:


def delta_c(n0,gr,boolgr,delta_t,nri ,ns):
    # gr is a matrix for growth rate 
    # nri is the index of resource from 0 to nr-1
    #boolgr is a vector to determine whether the bug is consuming the resource 
    ###### we need to figure how to generate boolgr
    deltaC = 0
    for i in range(ns):
        deltaC = deltaC + boolgr[i]*delta_n(n0[i],gr[i][nri],delta_t)
        #print(gr[i][nri])
    return deltaC


# In[8]:


#####
def find_root(f, low, high, precision):
    while high - low > precision:
        mid = (low + high) / 2.0
        if f(mid) < 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


# ### scipy.optimize.root

# In[9]:


def runouttime(nri, boolgr, n0, s, gr):  
    #calculate the runout time for the (nri+1)th neutrient
    #n0(vec),s(number) remaining bug & neutrient concentration
    #import scipy.optimize as optimize
    if s < 0.000000000001:
        solution = 0 
    elif max(boolgr)==0:
        solution = np.nan
    else:
        def equation(t):
            deltaC = sum(boolgr[i] * delta_n(n0[i], gr[i][nri], t) for i in range(len(n0)))
            return deltaC - s
        #solution = optimize.fsolve(equation, 1.0)
        solution = find_root(equation, 0, 100, 0.000000000001)
        '''
    if solution > 10:
        print("consumption is:")
        print(sum(boolgr[i] * delta_n(n0[i], gr[i][nri], solution) for i in range(len(n0))))
        print("resource left is %s" %(s))
    summ = sum(boolgr[i] * delta_n(n0[i], gr[i][nri], solution) for i in range(len(n0)))
    if summ > 2:
        print("__________________________")
        print("consumption is %s"%(summ))
        print("deltat is %s"%(solution))
        #for i in range(len(n0)):
            #if boolgr[i]==1:
                #print(delta_n(n0[i], gr[i][nri], solution))
    '''
    return solution


# ### get_boolgr for diauxie

# In[10]:


def get_boolgr(rank, nsi, c0,n0,withinNichtime,lag,Nnich):
    boolgr = np.zeros(len(c0))
    # spic_resource[i] = k the i+1_th resource is used by the k+1_th species
    for i in range(len(c0)):
        #print("nsi:%s, i:%s, rank[nsi][i]:%s;"%(nsi,i,rank[nsi][i]))
        if c0[rank[nsi][i]] > 0.000000000001:
            boolgr[rank[nsi][i]]=1
            break

    if n0[nsi] < 0.000000000001:
        boolgr = np.zeros(len(c0))
    print("lag[%s] = %s"%(nsi,lag[nsi]))
    print("withinNichtime = %s"%(withinNichtime))
    print("Nnich = %s"%(Nnich))
    if lag[nsi] > withinNichtime and Nnich > 0:
        print("species %s is frozen!"%(nsi))
        boolgr = np.zeros(len(c0))
    #print(type(boolgr))
    return boolgr


# In[11]:


def sigma(s):
    # Compute the squared differences
    squared_diff = [(si - 1/3)**2 for si in s]
    # Calculate the average of squared differences
    average_squared_diff = np.mean(squared_diff)
    # Calculate the standard deviation
    sigma_RS = np.sqrt(average_squared_diff)
    return sigma_RS


# In[12]:


def PseudoUniformSupply(size = 3, k = 0):
    vector = np.random.uniform(k+0, k+10, size)
    vector /= sum(vector)
    return vector, sigma(vector)


# In[13]:


def normalize_vector(v):
    sum_of_components = np.sum(v)
    if sum_of_components == 0:
        return v
    elif sum_of_components>1000000000:
        max_value = np.max(v)
        max_index = np.argmax(v)
        #print(max_value)
        #print(max_index)
        v = np.zeros(len(v))
        v[max_index] = 1
        return v
    else:
        return v / sum_of_components


# In[14]:


v = np.array([368.3750, 25.6786, 0.4801, 1.4662])
print(normalize_vector(v))


# In[15]:


###test get_boolgr
'''
nr=4
ns=4
rank = np.array([[4, 2, 1, 3],
       [1, 2, 4, 3],
       [3, 2, 4, 1],
       [4, 1, 3, 2]])

rank = rank - 1
#rank = rank.transpose()
print(rank)
c0 = np.array([2.03166546, 1.76250143, 0.0492092 , 0.11662391])
n0 = [0.01,0.01,0.01,0.01]

#res_used_by = {}

#for i in range(0, nr):
#    res_used_by[i] = []  
hhh = []
for nri in range(ns):
    #spic_resource = []
    boolgr = get_boolgr(rank, nri, c0,n0)
    hhh.append(boolgr)
print(np.array(hhh).transpose())
'''


# In[16]:


#print(boolgr.transpose())


# In[17]:


def onecycle(n0new, c0new,spciesGR,round1,lag):
    printround = -1
    ### calculate the dynamic within one cycle
    # n0 vector length ns
    # c0 vector length nr
    Nnich = 0
    n0 = [n0new]
    c0 = [c0new]
    time_node = [0] ## time_node start from 0
    resourceLeft = nr
    while True:
        WithinNichNode = []
        withinNichtime = 0
        lag2 =  list(lag)
        for ii in range(ns+2):
            boolgr = []
            for i in range(ns):
                boolgr_temp = get_boolgr(rank, i, c0[-1],n0[-1],withinNichtime,lag,Nnich)
                boolgr.append(boolgr_temp)
            boolgr = np.array(boolgr).transpose() 
            
            if round1 == printround:
                #print("________________")
                #print("ns = %s"%(ns))
                #print("ii=%s"%(ii))
                #print("rank")
                #print(rank)
                print("c0:")
                print(c0[-1])
                print("n0:")
                print(n0[-1])
                #print("withinNichtime")
                #print(withinNichtime)
                #print("lag2")
                #print(lag2)
                #print("boolgr")
                #print(boolgr)
              
            runoutT = [runouttime(i, boolgr[i], n0[-1],c0[-1][i],spciesGR) for i in range(nr)]
            lagmin = min(lag2)
            min_index = lag2.index(lagmin)
            lag2[min_index] = 100
            try:
                deltaT, runout_index = min((val, idx) for (idx, val) in enumerate(runoutT) if val > 0)
            except ValueError:
                deltaT = lagmin
            
            if round1 == printround:
                print("Tmin:")
                print(deltaT)
                print("lagmin:")
                print(lagmin)
                print("min_index")
                print(min_index)
               
            if deltaT >= lagmin and Nnich != 0:
                deltaT = lagmin
                lag2 = [x - lagmin for x in lag2]
                
            if round1 == printround:
                print("deltaT:")
                print(deltaT)
            
            withinNichtime = withinNichtime + deltaT
            time_node.append(float(deltaT+time_node[Nnich]))
            next_n0 = []
            for i in range(ns):
                summ = sum(delta_n(n0[-1][i],spciesGR[i][j],deltaT)*boolgr[j][i] for j in range(nr))
                '''
                print("deltaT:")
                print(deltaT)
                print("deltan %s %s * boolgr = %s *%s"%(i,0,delta_n(n0[-1][i],spciesGR[i][0],deltaT),boolgr[0][i]))
                print("deltan %s %s * boolgr = %s *%s"%(i,1,delta_n(n0[-1][i],spciesGR[i][1],deltaT),boolgr[1][i]))
                print("deltan %s %s * boolgr = %s *%s"%(i,2,delta_n(n0[-1][i],spciesGR[i][2],deltaT),boolgr[2][i]))
                print("sum")
                print(sum)
                '''
                next_n0.append(n0[Nnich][i] + summ)
            
            n0.append(next_n0) 
            next_c0 = []
            for j in range(nr):
                summ2 = delta_c(n0[-2],spciesGR,boolgr[j],deltaT,j,ns)
                #print("c0[%s][%s] = %s"%(Nnich,j,c0[Nnich][j]))
                #print("deltaC = %s"%(summ2))
                next_c0.append(c0[Nnich][j] - summ2)
            print("_________________")
            #next_c0 = [(c0[Nnich][j]-delta_c(n0[-2],spciesGR,boolgr[j],deltaT,j,ns))for j in range(nr)]
            c0.append(next_c0) 
            #print(next_c0)
            WithinNichNode.append(withinNichtime)
            Nnich += 1
            
            Newcount = sum(1 for i in c0[-1] if i > 0.000000000001)
            if Newcount < resourceLeft:
                resourceLeft = Newcount
                print("Newcount = %s"%(Newcount))
                break
            
        if max(c0[-1]) < 0.000000000001:
            print("bbb")
            break
        
    if round1 == printround:
        plt.plot(time_node, n0)
        print(n0)
        print(time_node)
    #print("______________")
        #print(time_node)

    #print(Nnich)
    #print(time_node[-1])
    return n0, c0, time_node


# ### read pickle file

# In[18]:


import pickle

# Open the .pkl file in binary mode for reading
#with open('3N4R_all_comm.pkl', 'rb') as file:
    # Load the object from the file
    
with open('3N4R.pkl', 'rb') as file:    
    data = pickle.load(file)


# In[19]:


#print(data['communities'][1722])


# In[20]:


def truncated_gaussian(mean, std_dev, lower_bound, upper_bound, size):
    a = (lower_bound - mean) / std_dev
    b = (upper_bound - mean) / std_dev
    samples = truncnorm.rvs(a, b, loc=mean, scale=std_dev, size=size)
    return samples


# In[21]:


mean = 0.5  
std_dev = 1  
lower_bound = 0.1  
upper_bound = 0.8  

#samples = truncated_gaussian(mean, std_dev, lower_bound, upper_bound, size)
#print(samples)


# In[22]:


ns = 3
nr = 3
stat = np.zeros(ns+1)
anomlyIndex = []
onenum = 0
success = 0


# In[23]:


### generate lag
#lag = truncated_gaussian(mean, std_dev, lower_bound, upper_bound, ns).tolist()
lag = [100,100,100]


# In[24]:


print("start")

for iii in range(1):
    data['communities'][iii]['selected']['possible_resources_stablilty'] = []
    print("first check possible_resource_supply")
    ### first check possible_resource_supply
    if data['communities'][iii]['selected']['possible_resources'] == []:
        continue
    else:
        onenum = onenum +  1
        aaa = 0
        #for ww in range(len(data['communities'][iii]['selected']['possible_resources'])):
        for ww in range(1):
            print("setting")
            resourcesele = data['communities'][iii]['selected']['possible_resources'][ww]
            spciesGR = data['communities'][iii]['g'][:, resourcesele]
            
            all_numbers = [0, 1, 2, 3]
            not_in_vec = [num for num in all_numbers if num not in resourcesele]
            dele = not_in_vec[0] +1
            
            matrixx = copy.deepcopy(data['communities'][iii]['pref_list'])

            for i in range(matrixx.shape[0]):
                for j in range(matrixx.shape[1]):
                    if matrixx[i, j] == dele:
                        matrixx[i, j] = 0

            rank = np.array([[value for value in row if value != 0] for row in matrixx])
            rank[rank > dele] -= 1
    
            c0 = []
            #for i in range(nr):
                #c0[0].append(1)
            c0.append(data['communities'][iii]['selected']['possible_resources_supply'][ww])

            print("set up dilution factor")
            ##set up dilution factor 
            DF = 100
            rank = rank - 1

            lendata = 100
            #Ntr = 32

            #start = 100 
            #for Ntr in range(start, 10000):

            #set up initial population
            n0 = [[]]
            mean = 0 
            stddev = 0.0001
            for i in range(ns):
                n0[0].append(0.01)

            time_node = [0]

            #record Population fractions at the end of each growth cycle

            populationEndRelative = []
            populationEndAbs = []

            print("generate the first cycle")

            #generate the first cycle
            ntrial, ctrial, timenodetrial = onecycle(n0[0],c0[0],spciesGR,-1,lag)

            populationEndAbs.append(ntrial[-1])

            log_N = normalize_vector(ntrial[-1])
            #print(log_N)
            populationEndRelative.append(log_N)

            #generate other cycle
            cycleNum = 10

            realcyc = 1
            print("generate other cycle!!!!!!!!!!!!")
            for k in range(cycleNum):
                print("the %s th cycle"%(k))
                n0temp = []
                for i in range(len(ntrial[-1])):
                    if ntrial[-1][i]<0.000000000001:
                        n0temp.append(0)
                    else:
                        n0temp.append(ntrial[-1][i]/ DF)
                n0[0] = n0temp
                #print(n0[0])
                #c0[0] , avg = PseudoUniformSupply(nr)
                # visualize
                #print(c0)
                ntrial, ctrial, timenodetrial = onecycle(n0[0],c0[0],spciesGR,k,lag)

                #print(ntrial[-1])
                if max(ntrial[-1])<0.000000000001:
                    populationEndRelative.append(ntrial[-1])
                    break
                log_N = normalize_vector(ntrial[-1])

                realcyc = realcyc + 1
                populationEndAbs.append(ntrial[-1])
                populationEndRelative.append(log_N)


            populRelative = np.transpose(populationEndRelative)
            populAbs = np.transpose(populationEndAbs)


            ### check survive
            survNum = 0

            if realcyc > cycleNum - 2:
                for i in range(ns):
                    if populRelative[i][-1] > 0.01:
                        survNum = survNum + 1
            else:
                survNum = 1
            stat[int(survNum)] = stat[int(survNum)]+1
            print("surive num of the %s th trial is %f" %(iii,survNum))
            
            if survNum ==3:
                data['communities'][iii]['selected']['possible_resources_stablilty'].append(1) 
                aaa = 1
            else:
                data['communities'][iii]['selected']['possible_resources_stablilty'].append(0) 
            
        if aaa == 1:
            success = success + 1
            
        cycle = np.linspace(0, cycleNum,  cycleNum+1)
        #for i in range(3):
            #plt.plot(cycle, populRelative[i])
            
        plt.show()


# In[25]:


cycle = np.linspace(0, cycleNum,  cycleNum+1)
print(cycle)


# In[26]:


for i in range(3):
    plt.plot(cycle, populRelative[i])
plt.ylim(0, 1)            
plt.show()


# In[ ]:




