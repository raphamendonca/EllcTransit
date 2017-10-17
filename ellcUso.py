import ellcCall as ec
import pandas as pd
import numpy as np
from datetime import datetime

#pts = np.arange(10000,100001,10000)

for pt  in np.arange(10000000,20000001,1000000):
    ec.calcFLUX( 25, pt)
    ec.calcFluxMp(25, pt, 2)
    print("--------",pt,"---------")

df = pd.DataFrame({"Pontos": ec.lista_pontos, "Processos": ec.lista_processos, "Duração": ec.lista_duracoes })
print(df)


