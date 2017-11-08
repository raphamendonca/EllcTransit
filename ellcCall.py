#%load_ext Cython
#%pylab inline

import ellc

from datetime import datetime
from numpy import linspace

rho = 0.4   # Stellar density in solar units
period = 1  # Orbital period in days
p_rot  = 23.9 # Rotation period of the star in days
k = 0.1     # R_planet/R_star 
a=(3.75226985055)*((period)**(2./3.))*((rho)**(1./3.))   # semi-major axis in solar radii
print("a = {:.3f} Rsun".format(a))

r_1 = 1/a     # R_star/a
r_2 = k*r_1   # R_planet/a

# Limb darkening parameters for the star
ld_1 = 'quad'
ldc_1 = [0.4,0.3]

t_zero=0.25   # Time of mid-transit

rotfac_1 = period/p_rot

# Spots on star
spots_1 = [[180-360*t_zero*rotfac_1,120-360*t_zero*rotfac_1],
           [0,0],
           [5.739170477266787, 8.11641272572196],
           [0.5,0.5]]

lista_duracoes = []
lista_pontos = []
lista_processos = []

def registraExecucao(curvas, pontos, duracao, processos):
    lista_pontos.append(pontos)
    lista_duracoes.append(duracao)
    lista_processos.append(processos)
    
    #print( lista_pontos, ' | ', lista_processos ,' | ',lista_duracoes)
   
    #df = df.append({'Pontos': pontos, 'Duracao' : duracao, 'Processos' : processos}, ignore_index=True)
    
        
def calcFLUX(curvas,pontos):
    startTime = datetime.now()
    
    time = linspace(0,curvas,pontos)
    
    flux = ellc.lc(time,t_zero=t_zero, period=period, \
        radius_1=r_1, radius_2=r_2,incl=90,sbratio=0, rotfac_1 = rotfac_1, \
        ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',\
        grid_1='sparse',grid_2='sparse',spots_1=spots_1)

    endTime = datetime.now()
    duracao = endTime - startTime
    
    registraExecucao(curvas, pontos, duracao, 1)
    return flux

    
def calcFluxMp(curvas, pontos, processos):
    startTime = datetime.now()
    
    time = linspace(0,curvas,pontos)

    flux = ellc.lcOpenMp(time, t_zero=t_zero, period=period, \
        radius_1=r_1, radius_2=r_2,incl=90,sbratio=0, rotfac_1 = rotfac_1, \
        ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',\
        grid_1='sparse',grid_2='sparse',spots_1=spots_1)

    endTime = datetime.now()
    duracao = endTime - startTime
    
    registraExecucao(curvas, pontos, duracao, processos)  
    return flux