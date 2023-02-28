########################################################
#          using Metropolis algorithm                  #
#                     DGlab                            #
########################################################





import numpy as np
import math

#==========================================================
# Defining the system size and system parameters 
#==========================================================
N = 3
j =  1.0
k =  1.0
Temp = 1.0

T = Temp
iteration = [ 40000 ]

#===========================================================
#Generating a ranndom configuration
#===========================================================
def ising(n):   # makes a N*N lattice
    l= []
    for i in range(n):
        l.append([np.random.choice([-1,1]) for i in range(n)])
    return l

#=============================================================
# Calculating energy of a configuration using PBC
#=============================================================
def energyf(l):

    energy = 0
    for i in range(len(l)):
        for k in range(len(l[i])):
            up    = i-1
            down  = i+1
            left  = k-1
            right = k+1


            if(i == 0 ):
                up  = len(l) - 1

            if(i == (len(l) - 1)):
                down = 0

            if (k == 0):
                left = len(l[i]) - 1

            if (k == len(l[i]) - 1):
                right = 0


            energy = energy +  j *l[i][k]* (  l[up][k] + l[down][k]  + l[i][left] + l[i][right])   # 2d energy 
    return energy



#================================================================
#randomly updating spin sites 
#================================================================
def update_config(config, a):
    k =1
    for i in range(k):
        x= np.random.randint(len(config))
        y= np.random.randint(len(config))
        config[x][y]= -1* config[x][y]

        up    = x - 1
        down  = x + 1
        left  = y - 1
        right = y + 1


        if(x == 0 ):
           up  = len(config) - 1

        if(x == (len(config) - 1)):
           down = 0

        if (y == 0):
           left = len(config[x]) - 1

        if (y == len(config[x]) - 1):
           right = 0

    delE= 2 * j * config[x][y] * (  config[up][y] + config[down][y] + config[x][left] + config[x][right] )

    return config, delE



#===========================================================
#Generating the Markov chain using Monte-Carlo sampling 
#===========================================================
def montecarlof(start_config, energy, iter1, T):
    old_energy = energy
    w =0
    # writing the Configs and energies 
    file1 = open("ising_dataset.dat", "w+")      
    file1.write("\n")
    for i in range(len(start_config)):
        for j in range(len(start_config[i])):
            file1.write(str( (start_config[i,j] + 1)/2 )+"  ")
    file1.write(str(old_energy))


    #===========================================================
    #starting monte carlo 
    #===========================================================

    for i in range(iter1):
        new_config, del_E = update_config(start_config.copy(), i)

        if (del_E < 0):
            start_config = new_config.copy()
            old_energy= old_energy + del_E
        else:
            prb = math.exp( (float( -del_E )) / float( T) )
            r = np.random.random()

            if(r < prb):
                start_config = new_config.copy()
                old_energy = old_energy + del_E
                w=w+1

        file1.write("\n")


        for m in range(len(start_config)):
            for j in range(len(start_config[m])):
                # saving configurations and eneegies
                file1.write(str( (start_config[m,j] + 1)/2 )+"  ") 
        file1.write(str(old_energy))
    file1.close()


    print(w/iter1)    # acceptance rate 


    return start_config, old_energy

#========================================================================
# Main function
#========================================================================
if __name__ == "__main__" :
    isng = np.array(ising(N))
    #for i in range(20):
    for iter1 in iteration:
        config , energy =  montecarlof( isng, energyf(isng), iter1, T) 
