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


def get_boolgr(rank, nsi, c0,n0):
    boolgr = np.zeros(len(c0))
    # spic_resource[i] = k the i+1_th resource is used by the k+1_th species
    for i in range(nr):
        if c0[rank[nsi][i]] > 0.000000000001:
            boolgr[rank[nsi][i]]=1
            break

    if n0[nsi] < 0.000000000001:
        boolgr = np.zeros(len(c0))

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


def onecycle(n0new, c0new,spciesGR,round1):
    printround = 270
    ### calculate the dynamic within one cycle
    # n0 vector length ns
    # c0 vector length nr
    Nnich = 0
    n0 = [n0new]
    c0 = [c0new]
    time_node = [0] ## time_node start from 0

    while True:
        boolgr = []
        
        for i in range(ns):
            boolgr_temp = get_boolgr(rank, i, c0[Nnich],n0[Nnich])
            boolgr.append(boolgr_temp)
            
        boolgr = np.array(boolgr).transpose() 
        #print(boolgr)
        #print(len(boolgr))
        #print(len(c0[Nnich]))
        runoutT = [runouttime(i, boolgr[i], n0[Nnich],c0[Nnich][i],spciesGR) for i in range(nr)]
        #if round1 == printround:
            #print("t:")
            #print(runoutT)
        #print(runoutT)
        try:
            deltaT, runout_index = min((val, idx) for (idx, val) in enumerate(runoutT) if val > 0)
            #deltaT, runout_index = min((val, idx) for (idx, val) in enumerate(runoutT))
        except ValueError:
            break
            
        #if round1 == printround:
            #print("deltaT:")
            #print(deltaT)
        time_node.append(float(deltaT+time_node[Nnich]))

        next_n0 = [n0[Nnich][i] + sum(delta_n(n0[Nnich][i],spciesGR[i][j],deltaT)*boolgr[j][i] for j in range(nr)) for i in range(ns)]
        '''
        for i in range(ns):
            summ = sum(delta_n(n0[Nnich][i],spciesGR[i][j],deltaT)*boolgr[j][i] for j in range(nr))
            if summ > 2:
                for j in range(nr):
                    if boolgr[j][i] == 1:
                        print("n0[Nnich][i] is %s"%(n0[Nnich][i]))
                        print("deltaT is %s"%(deltaT))
                        print("spciesGR[i][j] is %s"%(spciesGR[i][j]))
                        print("delta_n is %s" %(delta_n(n0[Nnich][i],spciesGR[i][j],deltaT)))
           '''             
                        
        n0.append(next_n0) 
    
        next_c0 = [(c0[Nnich][j]-delta_c(n0[Nnich],spciesGR,boolgr[j],deltaT,j,ns))for j in range(nr)]
        c0.append(next_c0) 
        #print(next_c0)
        #print(c0[Nnich+1])
        '''
        if round1 == printround:
            print("@@@@@@@@@@@@@@@@@@@@@@")
            print("n0:")
            print(n0)
            print(boolgr)
            print("c0:")
            print(next_c0)
            print("@@@@@@@@@@@@@@@@@@@@@@")
        if all(element < 0.000000000001 for element in c0[Nnich+1]):
            break
'''
        Nnich += 1
    #if round1 == printround:
        #plt.plot(time_node, n0)
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
    
with open('1_season_random_species.pkl', 'rb') as file:    
    data = pickle.load(file)


# In[19]:


#print(data['communities'][0])
print(len(data['communities']))


# In[20]:


print(len(data['communities']))
print(data['communities'][1])


# In[21]:


anomlyIndex = [7, 10, 26, 43, 78, 118, 121, 143, 155, 163, 219, 242, 265, 320, 328, 337, 350, 363, 470, 551, 592, 619, 622, 636, 695, 704, 719, 737, 738, 792, 796, 858, 873, 875, 908, 914, 915, 922, 1011, 1074, 1146, 1197, 1220, 1259, 1317, 1346, 1365, 1386, 1425, 1431, 1448, 1481, 1679, 1691, 1723, 1752, 1809, 1819, 1835, 1879, 1881, 1989, 2025, 2050, 2053, 2078, 2125, 2158, 2171, 2173, 2191, 2204, 2257, 2273, 2295, 2317, 2383, 2398, 2453, 2490, 2521, 2605, 2642, 2671, 2683, 2717, 2728, 2738, 2743, 2821, 2822, 2832, 2933, 2983, 3003, 3047, 3065, 3092, 3096, 3153, 3166, 3172, 3233, 3267, 3323, 3338, 3355, 3388, 3396, 3402, 3408, 3411, 3436, 3444, 3473, 3483, 3555, 3568, 3590, 3678, 3725, 3754, 3779, 3784, 3789, 3798, 3812, 3867, 3925, 3965, 4047, 4067, 4103, 4129, 4159, 4200, 4298, 4310, 4334, 4390, 4407, 4434, 4476, 4497, 4548, 4553, 4557, 4609, 4733, 4750, 4781, 4786, 4856, 4937, 4973, 5001, 5056, 5092, 5097, 5134, 5194, 5204, 5340, 5362, 5389, 5398, 5408, 5461, 5488, 5574, 5630, 5672, 5689, 5712, 5731, 5734, 5746, 5798, 5805, 5817, 5840, 5845, 5882, 5897, 5935, 6047, 6072, 6073, 6079, 6158, 6186, 6239, 6285, 6307, 6339, 6346, 6388, 6449, 6450, 6460, 6488, 6491, 6579, 6622, 6655, 6699, 6703, 6712, 6738, 6747, 6788, 6821, 6850, 6878, 6924, 6952, 6955, 6976, 7014, 7021, 7042, 7049, 7059, 7102, 7125, 7187, 7326, 7337, 7373, 7391, 7475, 7487, 7497, 7523, 7609, 7694, 7699, 7749, 7801, 7829, 7834, 7945, 8044, 8069, 8074, 8077, 8084, 8087, 8090, 8114, 8131, 8195, 8257, 8346, 8366, 8380, 8398, 8468, 8471, 8477, 8482, 8493, 8543, 8578, 8583, 8614, 8642, 8650, 8653, 8698, 8894, 8909, 8932, 8940, 8958, 8996, 9015, 9038, 9042, 9067, 9140, 9147, 9188, 9346, 9381, 9384, 9429, 9460, 9489, 9572, 9577, 9594, 9598, 9629, 9674, 9688, 9700, 9742, 9745, 9769, 9797, 9818, 9822, 9824, 9842, 9852, 9868, 9880, 9972, 10029, 10059, 10157, 10193, 10200, 10204, 10221, 10266, 10270, 10271, 10373, 10377, 10467, 10480, 10515, 10519, 10521, 10527, 10581, 10593, 10594, 10630, 10688, 10694, 10695, 10729, 10732, 10747, 10777, 10828, 10867, 10885, 10901, 10906, 10933, 10958, 10975, 11042, 11103, 11127, 11131, 11148, 11153, 11210, 11236, 11239, 11250, 11269, 11321, 11468, 11501, 11506, 11509, 11519, 11557, 11580, 11618, 11707, 11710, 11723, 11737, 11738, 11767, 11819, 11852, 11871, 11890, 11895, 11896, 11898, 11914, 11923, 11969, 12014, 12033, 12039, 12064, 12274, 12289, 12348, 12421, 12444, 12454, 12491, 12496, 12518, 12522, 12573, 12600, 12613, 12642, 12644, 12675, 12692, 12706, 12718, 12752, 12782, 12882, 12903, 12912, 12915, 12925, 12953, 12956, 12968, 12992, 13004, 13012, 13024, 13104, 13115, 13126, 13130, 13251, 13254, 13280, 13295, 13334, 13340, 13360, 13456, 13459, 13462, 13498, 13506, 13617, 13623, 13676, 13724, 13729, 13734, 13755, 13805, 13948, 13971, 13999, 14004, 14064, 14066, 14069, 14115, 14118, 14123, 14135, 14136, 14141, 14173, 14193, 14213, 14218, 14257, 14270, 14327, 14352, 14377, 14401, 14433, 14461, 14470, 14591, 14616, 14624, 14636, 14742, 14764, 14818, 14918, 15037, 15040, 15065, 15116, 15148, 15184, 15194, 15201, 15379, 15399, 15416, 15481, 15501, 15508, 15550, 15561, 15598, 15600, 15628, 15713, 15764, 15809, 15880, 15892, 15938, 15979, 15997, 16005, 16008, 16028, 16067, 16086, 16093, 16102, 16155, 16232, 16336, 16347, 16349, 16351, 16359, 16371, 16374, 16400, 16446, 16523, 16558, 16609, 16621, 16689, 16690, 16709, 16710, 16723, 16769, 16777, 16832, 16878, 16907, 16947, 16972, 17003, 17004, 17136, 17185, 17193, 17220, 17235, 17265, 17369, 17380, 17393, 17399, 17400, 17419, 17461, 17466, 17485, 17498, 17505, 17558, 17561, 17592, 17611, 17650, 17693, 17694, 17717, 17732, 17765, 17858, 17908, 17953, 18014, 18030, 18056, 18097, 18135, 18142, 18174, 18216, 18332, 18450, 18464, 18493, 18519, 18545, 18605, 18644, 18655, 18669, 18814, 18891, 18928, 18941, 19003, 19005, 19012, 19041, 19042, 19044, 19045, 19085, 19090, 19121, 19126, 19131, 19149, 19161, 19246, 19286, 19323, 19370, 19441, 19453, 19472, 19534, 19575, 19637, 19661, 19669, 19690, 19728, 19756, 19769, 19778, 19805, 19821, 19833, 19867, 19876, 19906, 19910, 19921, 19953]


# In[22]:


#anomlyIndex= [1, 9, 22, 32, 44, 49, 50, 53, 81, 140, 144, 157, 159, 163, 164, 174, 190, 225, 226, 245, 248, 255, 261, 271]


# In[23]:


my_dict = {}


# In[24]:


ns = 4
nr=4
stat = np.zeros(ns+1)

onenum = 0


for iii in anomlyIndex:
    ### first check possible_resource_supply
    aaaa = []
    for tt in range(len(data['communities'][iii]['possible_resource_supply'])):
        spciesGR = data['communities'][iii]['g']
        rank = data['communities'][iii]['pref_list']

        c0 = []
        #for i in range(nr):
            #c0[0].append(1/nr)

        c0.append(data['communities'][iii]['possible_resource_supply'][tt])  

        #print(ns)

        ####################################################
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
            n0[0].append(0.0095)

        time_node = [0]

        #record Population fractions at the end of each growth cycle

        populationEndRelative = []
        populationEndAbs = []



        #generate the first cycle
        ntrial, ctrial, timenodetrial = onecycle(n0[0],c0[0],spciesGR,1)

        populationEndAbs.append(ntrial[-1])

        log_N = normalize_vector(ntrial[-1])
        #print(log_N)
        populationEndRelative.append(log_N)

        #generate other cycle
        cycleNum = 1000

        realcyc = 1

        for k in range(cycleNum):
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
            ntrial, ctrial, timenodetrial = onecycle(n0[0],c0[0],spciesGR,k)

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

        linear_space = np.linspace(1,realcyc ,realcyc)
        
        for x in range(ns):
                #if x not in extinct_spicies:
            plt.plot(linear_space, populAbs[x])
        plt.show()
        
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

        #print(stat)
        aaaa.append(survNum)
    print(aaaa)
    my_dict[iii] = aaaa
  


# In[25]:


print(stat)


# In[26]:


print(anomlyIndex)


# In[27]:


kkk=0
for i in anomlyIndex:
    if data['communities'][i]['possible_resource_supply'] == []:
        continue
    else:
        print(len(data['communities'][i]['possible_resource_supply']))
    kkk = kkk+1
print(kkk)


# In[28]:


linear_space = np.linspace(1,realcyc ,realcyc)


for x in range(ns):
        #if x not in extinct_spicies:
    plt.plot(linear_space, populRelative[x])
     
#plt.plot(linear_space, popul[1])

plt.title("4 spcies and 4 resorces")
plt.xlabel("growth cycle")
plt.ylabel("population fraction")

plt.ylim(0, 1.1)  # set the ylim to be -2 to 2
plt.show()


# In[29]:


linear_space = np.linspace(1,realcyc ,realcyc)

for x in range(ns):
        #if x not in extinct_spicies:
    plt.plot(linear_space, populAbs[x])
#plt.plot(linear_space, populAbs[1])
#plt.plot(linear_space, popul[1])


#plt.ylim(0, 1.5)  # set the ylim to be -2 to 2
plt.show()


# In[30]:


print(populAbs)


# In[31]:


print(populAbs[1][273])


# In[32]:


A = np.array([   0,  528, 3572,  900])
B = np.array([   0,  509, 3339,  822])
c = A+B
print(c)


# In[33]:


print(0.01*np.ones(4))


# In[34]:


print(1722/9670)


# In[35]:


print(2311/2331)


# In[37]:


print(my_dict)


# In[39]:


twopossible = np.zeros(3)
threepossible = np.zeros(4)

for ii in anomlyIndex:
    if len(my_dict[ii]) == 2:
        twopossible[my_dict[ii].count(4)] = twopossible[my_dict[ii].count(4)] +1
    elif len(my_dict[ii]) == 3:
        threepossible[my_dict[ii].count(4)] = threepossible[my_dict[ii].count(4)] +1
    else:
        print("index:%s;possible num: %s"%(ii,len(my_dict[ii])))


# In[40]:


print(twopossible)
print(threepossible)


# In[41]:


print(sum(twopossible))
print(sum(threepossible))


# In[54]:


print(4/71)


# In[47]:


fourpossible = np.zeros(5)
Indices = [719, 2743, 4973, 9629, 12675, 14173, 17732, 19575]
for kk in Indices:
    fourpossible[my_dict[ii].count(4)] = fourpossible[my_dict[ii].count(4)] +1
    
print(fourpossible)


# In[ ]:




