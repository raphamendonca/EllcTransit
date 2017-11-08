# -*- coding: utf-8 -*-
####
import ellcCall as ec
import pandas as pd
import numpy as np 
from numpy import linspace
from matplotlib import pyplot as plt

from datetime import datetime

def saveXls(pts,prcs,durs):
    #print(pts)
    #print(prcs)
    #print(durs)
    
    df = pd.DataFrame({"Pontos": pts, "Processos": prcs, "Duração": durs })
    data = datetime.now()
    filename = "output-"+str(data)+".xlsx"
    #print(filename)
    writer = pd.ExcelWriter(filename,  engine='xlsxwriter')
    #df.to_csv(filename)
    df.to_excel(writer, sheet_name='exexucao')
    writer.save()
    
    print(df)

#pts = np.arange(10000,100001,10000)
def callRange(ini, max, inc):
    for pt in np.arange(ini, max, inc):
        ec.calcFLUX( 25, pt)
        ec.calcFluxMp(25, pt, 2)
        print("--------",pt,"---------")

def callSingle(pt):
    flux = ec.calcFLUX( 25, pt)
    lc_plot(25,pt, flux)
    flux = ec.calcFluxMp(25, pt, 2)
    print("--------",pt,"---------")
    lc_plot(25,pt, flux)

def lc_plot(curvas, pontos, flux):
    time = linspace(0,curvas,pontos)

    fontsize=12
    fig=plt.figure(1,figsize=(12,8))
    plt.subplot(211)
    plt.xlim([0,24])
    plt.plot(time,flux,linewidth=1,color='darkblue')
    plt.xlabel("Time [d]",fontsize=fontsize)
    plt.ylabel("Normalized flux",fontsize=fontsize)
    plt.tick_params(axis='both', labelsize=fontsize)
    plt.subplot(234)
    plt.plot(time-6.25,flux,linewidth=1,color='darkblue')
    plt.locator_params(axis = 'x', nbins = 4)
    plt.xlabel("Time-6.25 [d]",fontsize=fontsize)
    plt.ylabel("Normalized flux",fontsize=fontsize)
    plt.xlim([-0.08,0.08])
    plt.ylim([0.975,0.995])
    plt.tick_params(axis='both', labelsize=fontsize)
    plt.subplot(235)
    plt.plot(time-9.25,flux,linewidth=1,color='darkblue')
    plt.locator_params(axis = 'x', nbins = 4)
    plt.xlabel("Time-9.25 [d]",fontsize=fontsize)
    plt.xlim([-0.08,0.08])
    plt.ylim([0.97,0.99])
    plt.tick_params(axis='both', labelsize=fontsize)
    plt.subplot(236)
    plt.plot(time-14.25,flux,linewidth=1,color='darkblue')
    plt.locator_params(axis = 'x', nbins = 4)
    plt.xlabel("Time-14.25 [d]",fontsize=fontsize)
    plt.xlim([-0.08,0.08])
    plt.ylim([0.98, 1.0])
    plt.tick_params(axis='both', labelsize=fontsize)
    plt.tight_layout()

callSingle(1000000)
#callRange(1000000, 10000001, 1000000)
#for i in np.arange(1, 4, 1):
    #callRange(1000000, 30000001, 1000000)
    #saveXls( ec.lista_pontos,  ec.lista_processos,  ec.lista_duracoes)


#df = pd.DataFrame({"Pontos": ec.lista_pontos, "Processos": ec.lista_processos, "Duração": ec.lista_duracoes })
#saveXls( ec.lista_pontos,  ec.lista_processos,  ec.lista_duracoes)



