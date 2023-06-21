#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import math
import scipy.optimize as optimize
import matplotlib.pyplot as plt


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
        return delta_n
    except OverflowError:
        print("OverflowError: Input value is too large for exponential calculation.")
        print(gr)
        print(t)
        return None


# In[7]:


def delta_c(n0,gr,boolgr,delta_t,nri ,ns = 2):
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
    if s < 0.0001:
        solution = 0 
    elif max(boolgr)==0:
        solution = np.nan
    else:
        def equation(t):
            deltaC = sum(boolgr[i] * delta_n(n0[i], gr[i][nri], t) for i in range(len(n0)))
            return deltaC - s
        #solution = optimize.fsolve(equation, 1.0)
        solution = find_root(equation, 0, 100, 0.000001)
    if solution > 99:
        print(s)
        print(n0)
        print(solution)
        print(boolgr)
    return solution


# ### get_boolgr for diauxie

# In[10]:


def get_boolgr(rank, nri, c0,n0):
    boolgr = []
    # spic_resource[i] = k the i+1_th resource is used by the k+1_th species
    for i in range(ns):
        if n0[i] < 0.0000001:
            boolgr.append(0)
        elif rank[i][nri]<0.0000001 and c0[nri]>0:     #resource_rank1 remains
            boolgr.append(1)
            #spic_resource.append(i)
        elif c0[rank[i][nri]-1] < 0.0000001 and c0[nri]>0:   #resource with higher rank ran out && this resource remains
            boolgr.append(1)
            #spic_resource.append(i)
        else:
            boolgr.append(0)
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


###test get_boolgr
'''
nr=4
ns=3
rank = np.array([[1, 2, 4, 3],
       [3, 1, 2, 4],
       [4, 1, 3, 2]])

rank = rank - 1

print(rank)
c0 = [12,13,1,4,2,3,5]
n0 = [1/3,1/3,1/3]

#res_used_by = {}

#for i in range(0, nr):
#    res_used_by[i] = []  

for nri in range(nr):
    #spic_resource = []
    boolgr = get_boolgr(rank, nri, c0,n0)
    print(boolgr)
#print(res_used_by)
'''


# In[15]:


def onecycle(n0new, c0new,spciesGR):
    ### calculate the dynamic within one cycle
    # n0 vector length ns
    # c0 vector length nr
    Nnich = 0
    n0 = [n0new]
    c0 = [c0new]
    time_node = [0] ## time_node start from 0

    while True:
        boolgr = []
        res_used_by = {}
        for i in range(0, nr):
            res_used_by[i] = [] 
        for i in range(nr):
            boolgr_temp = get_boolgr(rank, i, c0[Nnich],n0[Nnich])
            boolgr.append(boolgr_temp)
            
        runoutT = [runouttime(i, boolgr[i], n0[Nnich],c0[Nnich][i],spciesGR) for i in range(nr)]
        #print(runoutT)
        try:
            deltaT, runout_index = min((val, idx) for (idx, val) in enumerate(runoutT) if val > 0.001)
        except ValueError:
            break
            
        time_node.append(float(deltaT+time_node[Nnich]))

        next_n0 = [n0[Nnich][i] + sum(delta_n(n0[Nnich][i],spciesGR[i][j],deltaT)*boolgr[j][i] for j in range(nr)) for i in range(ns)]
        n0.append(next_n0) 
    
        next_c0 = [(c0[Nnich][j]-delta_c(n0[Nnich],spciesGR,boolgr[j],deltaT,j,ns))for j in range(nr)]
        c0.append(next_c0) 
        #print(next_c0)
        #print(c0[Nnich+1])
        if all(element < 0.001 for element in c0[Nnich+1]):
            break

        Nnich += 1
    #print(Nnich)
    #print(time_node[-1])
    return n0, c0, time_node


# ### read pickle file

# In[16]:


import pickle

# Open the .pkl file in binary mode for reading
with open('3N4R_success_comm.pkl', 'rb') as file:
    # Load the object from the file
    data = pickle.load(file)


# In[18]:


##set up dilution factor 
DF = 100

nr = 4
ns = 3

lendata = 9670
#Ntr = 9393

stat = np.zeros(ns+1)

start = 0 
for Ntr in range(start, start+2000):

    #set up initial population
    n0 = [[]]
    for i in range(ns):
        n0[0].append(1/ns/(DF-1))

    ##set up sources suppliment
    c0 = [[]]
    for i in range(nr):
        c0[0].append(1/nr)
    #s = 1/ns

    extinct_spicies = set()

    time_node = [0]

    #record Population fractions at the end of each growth cycle

    populationEndRelative = []
    populationEndAbs = []

    # generate_unanomolous growth rate & rank
    #spciesGR, rank = generate_random_vectors(ns, nr)
    spciesGR = data[Ntr]["g"]
    rank = data[Ntr]["pref_list"]
    rank = rank - 1

    #generate the first cycle
    ntrial, ctrial, timenodetrial = onecycle(n0[0],c0[0],spciesGR)

    populationEndAbs.append(ntrial[-1])

    log_N = normalize_vector(ntrial[-1])
    #print(log_N)
    populationEndRelative.append(np.round(log_N,12))

    #generate other cycle
    cycleNum = 1000

    realcyc = 1

    for k in range(cycleNum):
        n0[0] = [element / DF for element in ntrial[-1]]
        c0[0] , avg = PseudoUniformSupply(nr)
        # visualize
        ntrial, ctrial, timenodetrial = onecycle(np.round(n0[0],12),c0[0],spciesGR)

        #print(ntrial[-1])
        if max(ntrial[-1])<0.0000001:
            populationEndRelative.append(np.round(ntrial[-1],12))
            break
        log_N = normalize_vector(ntrial[-1])

        realcyc = realcyc + 1
        populationEndAbs.append(ntrial[-1])
        populationEndRelative.append(np.round(log_N,12))


    ### check survive
    populRelative = np.transpose(populationEndRelative)
    populAbs = np.transpose(populationEndAbs)

    survNum = 0

    if realcyc > cycleNum - 2:
        for i in range(ns):
            if populRelative[i][-1] > 0.000001:
                survNum = survNum + 1
    else:
        survNum = 1
    stat[int(survNum)] = stat[int(survNum)]+1
    print("surive num of the %s th trial is %f" %(Ntr,survNum))

    #print(stat)
    if survNum == 0:

        linear_space = np.linspace(1,realcyc ,realcyc)

        for x in range(ns):
                #if x not in extinct_spicies:
                plt.plot(linear_space, populAbs[x])

        #plt.plot(linear_space, popul[1])


        plt.ylim(0, 8)  # set the ylim to be -2 to 2
        plt.show()

        #break


# In[17]:


print(len(data))


# In[26]:


print(stat)


# In[19]:


print(realcyc)


# In[20]:


linear_space = np.linspace(1,realcyc ,realcyc)


for x in range(ns):
        #if x not in extinct_spicies:
        plt.plot(linear_space, populRelative[x])
        
#plt.plot(linear_space, popul[1])


plt.ylim(0, 1.2)  # set the ylim to be -2 to 2
plt.show()


# In[21]:


linear_space = np.linspace(1,realcyc ,realcyc)

for x in range(ns):
        #if x not in extinct_spicies:
        plt.plot(linear_space, populAbs[x])

#plt.plot(linear_space, popul[1])


plt.ylim(0, 200)  # set the ylim to be -2 to 2
plt.show()


# In[22]:


c = [np.nan,np.nan,np.nan]
print(sum(c))


# In[23]:


import numpy as np

x = [[2, 919, 79, 0],
     [3, 914, 83, 0],
     [4, 906, 90, 0],
     [5, 908, 87, 0],
     [4, 898, 98, 0],
     [4, 900, 96, 0],
     [1, 919, 80, 0],
     [7, 916, 77, 0],
     [3, 919, 78, 0],
     [3, 611, 56, 0]]

x_transposed = np.transpose(x)


# In[24]:


print(sum(x_transposed[0]))
print(sum(x_transposed[1]))
print(sum(x_transposed[2]))
print(sum(x_transposed))


# In[27]:


from prettytable import PrettyTable

# Specify the column names while initializing the Table
my_table = PrettyTable(["total trial\survive num", "0", "1", "2","3"])

# Add rows
my_table.add_row(["2000", "2", "1249", "607","142"])


print(my_table)


# In[ ]:




