# -*- coding: utf-8 -*-
####
import ellcCall as ec
import pandas as pd
import numpy as np

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
    ec.calcFLUX( 25, pt)
    ec.calcFluxMp(25, pt, 2)
    print("--------",pt,"---------")

#callSingle(10000)
#callRange(1000000, 10000001, 1000000)
for i in np.arange(1, 11, 1):
    callRange(1000000, 30000001, 1000000)
    saveXls( ec.lista_pontos,  ec.lista_processos,  ec.lista_duracoes)
    #print(i)


#df = pd.DataFrame({"Pontos": ec.lista_pontos, "Processos": ec.lista_processos, "Duração": ec.lista_duracoes })
#saveXls( ec.lista_pontos,  ec.lista_processos,  ec.lista_duracoes)



