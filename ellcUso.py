import ellcCall as ec
import pandas as pd
import numpy as np
from datetime import datetime

#pts = np.arange(10000,100001,10000)
def callRange():
    for pt in np.arange(10000,20001,10000):
        ec.calcFLUX( 25, pt)
        ec.calcFluxMp(25, pt, 2)
        print("--------",pt,"---------")

def callSingle(pt):
    ec.calcFLUX( 25, pt)
    ec.calcFluxMp(25, pt, 2)
    print("--------",pt,"---------")

callSingle(10000000)

df = pd.DataFrame({"Pontos": ec.lista_pontos, "Processos": ec.lista_processos, "Duração": ec.lista_duracoes })
print(df)


