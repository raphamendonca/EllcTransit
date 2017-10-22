import ellcCall as ec
import pandas as pd
import numpy as np
from datetime import datetime

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

#callSingle()


df = pd.DataFrame({"Pontos": ec.lista_pontos, "Processos": ec.lista_processos, "Duração": ec.lista_duracoes })
print(df)


