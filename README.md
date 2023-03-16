# Proyecto de investigación aplicada

Elaborado para acreditar:
**2023 1st WBDS LA Camp "Introduction to Bioinformatics and Data Science"**

Por **Yanel Bernardi** 

# TRANSCRIPTÓMICA

###### Podemos definir al transcriptoma de una célula como al conjunto total de transcriptos presentes en ella, para un estadio específico del desarrollo o condición fisiológica. 
###### El objetivo de la transcriptómica es estudiar el ARN total. Incluye determinar la cantidad o concentración de cada molécula de ARN, así como sus identidades moleculares, y por tanto, permite identificar sus niveles cualitativos y cuantitativos así como sus funciones.
###### Actualmente, las tecnologías NGS han revolucionado la transcriptómica al permitir la secuenciación del ADNc, conocida como RNAseq permitiendo un estudio a gran escala del los transcriptomas.

# Sobre el proyecto: 

###### La cresta neural es una población embrionaria que en un momento determinado del desarrollo embrionario mediante cambios en la expresión de genes, entre otros procesos, sufre un cambio en su fenotipo y adquiere la capacidad de migrar para luego diferenciarse en múltiples derivados en diferentes partes del organismo. 
###### Al estado previo a esos cambios se lo conoce como Pre-Migratorio y al posterior como Migratorio. 
###### El objetivo general de este proyecto es encontrar qué genes están expresados diferencialmente entre estas dos poblaciones celulares que corresponden a estos dos estados.  

> ###### *El dataset que vamos a utilizar es una matriz de counts de diferentes genes de pollo. Esta matriz surge del siguiente proceso:*

>> ###### **1-** Control de calidad y filtro de lecturas basado en las medidas de calidad
>>###### **2-** Alineamiento contra el genoma de referencia (gallus gallus)
>>###### **3-** Obtención del número de lecturas por gen/transcrito (Matriz de counts)

#Preparación del entorno

###### Con el siguiente código instalamos las bibliotecas que se van a usar en el proyecto:
```
!pip3 install pandas
!pip3 install seaborn
!pip3 install scipy
!pip3 install matplotlib
!pip3 install biopython
!pip3 install bioinfokit 
```
###### Aclaración: Las bibliotecas `sys`,`subprocess` y `io` forman parte del core de python por lo que no es necesario instalarlas.

###### Con este codigo las cargamos para que todo el código del notebook funcione:
```
import pandas as pd
import numpy  as np
import seaborn as sns
import subprocess
import sys
from scipy import stats
from bioinfokit.analys import norm 
from bioinfokit import analys, visuz
import matplotlib.pyplot as plt
from io import StringIO
```

# Obtención del dataset

###### Fuente: "Reconstruction of the Global Neural Crest Gene Regulatory Network In Vivo". Williams et al., 2019, Developmental Cell 51, 255–276
###### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121526

> ##### El dataset es un archivo gzip, primero hay que descomprimirlo y luego pedimos que sea leído como un archivo CSV indicando que los datos están separados por tabulación. Usamos para esta acción un metodo de pandas a través de este código para obtener los datos.
```
data = pd.read_csv('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE121nnn/GSE121526/suppl/GSE121526_counts_RNA-seq.txt.gz', compression='gzip', sep='\t')
```

> ##### Exploramos la tabla con los codigos siguientes para conocerla y evaluar si está lista para ser usar o hay que hacer un tratamiento previo:
```
data.info()
```
```
data.head(5)
```

##### En dataset se muestra la expresión de diferentes genes (filas) en tres diferentes poblaciones celulares (columnas):
*   non NC: son células que no son cresta neural
*   PM: son células de la cresta neural en un estadio pre migratorio
*   M: son células de la cresta neural en un estadio migratorio

##### Tenemos filas que corresponden a genes de pollo, las dos primeras columnas nos dan información de cada gen y el resto de las columnas corresponden las réplicas biólogicas de cada población celular secuenciada por RNA-seq que contienen las counts registradas para cada gen.

#Pre-procesado del dataset

###### Para trabajar mejor y reducir errores, primero vamos a realizar una copia de la tabla y a partir de la nueva tabla vamos a tratarla para que nos quede lista para el análisis de expresión de genes. 
```
data_copia = data.copy()
```
###### *Aclaración: Las modificaciones tienen que tener efecto en nuestro DataFrame "data_copia" así que vamos a incluir **inplace=True** en los códigos que sean necesarios*.

1.   Eliminar las columnas que corresponden a non NC en un sólo paso porque no las vamos a utilizar:
```
data_copia.drop(
  columns=[
    "./N1_bulk_RNA-seq_non_NC.out.sort.sam",
    "./N2_bulk_RNA-seq_non_NC.out.sort.sam",
    "./N3_bulk_RNA-seq_non_NC.out.sort.sam"],
  inplace=True)
```
2.   Usamos **rename**, para renombrar las columnas.
```
data_copia.rename(columns = {
    "./R1_bulk_RNA-seq_5-6ss.out.sort.sam": "1_PM",
    "./R1_bulk_RNA-seq_8-10ss.out.sort.sam": "1_M",
    "./R2_bulk_RNA-seq_5-6ss.out.sort.sam": "2_PM",
    "./R2_bulk_RNA-seq_8-10ss.out.sort.sam": "2_M",
    "./R3_bulk_RNA-seq_5-6ss.out.sort.sam": "3_PM",
    "./R3_bulk_RNA-seq_8-10ss.out.sort.sam": "3_M",
    "./R4_bulk_RNA-seq_5-6ss.out.sort.sam": "4_PM",
    "./R4_bulk_RNA-seq_8-10ss.out.sort.sam": "4_M"
},inplace=True)
```
3.   Reordenamos las columnas para que las replicas de cada población queden juntas.
```
data_copia = data_copia[['Geneid','Length','1_PM','2_PM','3_PM','4_PM','1_M','2_M','3_M','4_M']] 
```

# Normalización

###### Realizar una transformación de los counts a CPM (Counts Per Million) con el objetivo de lograr una normalización por el tamaño de la librería. Al tener en una misma escala todas las muestras las diferencias observadas entre las poblaciones celulares son debido a la biología y no por factores técnicos.

> CPM   (Counts Per Million) = 10^6 x (counts del gen "a") / número total de counts alineadas

##### Para por hacerlo ejecutamos el siguiente código:
```
# Preparamos los datos
del data_copia["Length"] # Eliminamos la columna Length 
data_copia = data_copia.set_index('Geneid') # Hacemos que la columna Geneid ahora sea el index

# Normalizamos los counts utilizando el método CPM
nm = norm()
nm.cpm(df=data_copia)
cpm_df = nm.cpm_norm # Obtenemos el DataFrame con los datos normalizado en CPM

```

# Clusterización 

###### Análisis de correlación por agrupamiento jerárquico:  Se agrupan las muestras en función de los coeficientes de correlación. 

> ###### Las replicas que son más similares entre sí están juntas y tienen una correlación de 1, es decir, poseen una alta correlación. A menor correlación, el cofieciente se aleja de 1 y lo vemos representado por colores más oscuros.

###### Este agrupamiento es útil para saber si los diferentes tipos de poblaciones celulares se pueden separar. Se espera que las muestras de las diferentes poblaciones sean más diferentes entre sí que las réplicas dentro de la misma población. La correlación entre replicas también es importante para los análisis estadísticos.

###### El clusterizado de muestras lo vamos a realizar utilizando el siguiente código:
```
Agrup_jerarquico_data = sns.clustermap(cpm_df.corr());
```

# Análisis de expresión diferencial

###### Agregamos dos columnas para incluir en nuestra tabla el promedio de las counts por población celular. 


```
cpm_df["PM_mean"] = (cpm_df['1_PM'] + cpm_df['2_PM'] + cpm_df['3_PM'] + cpm_df['4_PM'])/4
cpm_df["M_mean"] = (cpm_df['1_M'] + cpm_df['2_M'] + cpm_df['3_M'] + cpm_df['4_M'])/4
```

#### **Cálculo para obtener el log2Fold-Change**
> ###### Log2FC= log2[(Counts_norm_tratamiento)/(Counts_norm_control)]
Nosotros queremos evaluar cómo es la expresión de los genes en las células migratorias con respecto a las pre  migratorias. Entonces, nuestro "tratamiento" es población migratoria y nuestro "control" es población pre migratoria. 

###### Esta medida se usa para medir el cambio en el nivel de expresión de un gen entre dos muestras, en nuestro caso dos poblaciones celulares diferentes. Los umbrales en el nivel de cambio (“fold change”) de expresión y de valor p (“p-value”) para definir a un transcrito como diferencialmente expresado se ajusta dependiendo del rigor de los análisis.

1.   Eliminamos las filas que tienen counts = 0 en ambas poblaciones
```
cpm_df.drop(cpm_df[(cpm_df['PM_mean'] == 0) & (cpm_df['M_mean'] == 0)].index, inplace=True)
```
2.   La operación para el log2FC es: 
```
log2FC = np.log2((cpm_df["M_mean"])/(cpm_df["PM_mean"]))
```
* Vamos a agregar la columna con los datos del log2FC a nuestra tabla, por practicidad se puede poner en el mismo codigo de agregado de columna a la operación. 
```
cpm_df["log2FC"] = np.log2((cpm_df["M_mean"])/(cpm_df["PM_mean"]))
```

#### **Cálculo de un valor de significancia**

###### El p-value nos sirve para determinar si ese cambio en la expresión de un gen es estadísticamente significativo. 

1.   Creamos dos grupos donde incluimos las replicas para de las poblaciones en cada uno. 
```
PM = [cpm_df['1_PM'], cpm_df['2_PM'], cpm_df['3_PM'], cpm_df['4_PM']]
M = [cpm_df['1_M'], cpm_df['2_M'], cpm_df['3_M'], cpm_df['4_M']]
```
2.   Calculamos la varianza de cada grupo para determinar si asumimos que nuestras dos poblaciones tienen varianzas iguales.  
```
print(np.var(PM), np.var(M))
# El resultado, para este caso, fue:
var(PM): 194210.81928503313 
var(M):222081.1841629431
```
###### Podemos asumir que las poblaciones tienen varianzas iguales si: 
   varianza de la muestra más grande/la varianza de la muestra < 4. 
   **Si esto pasa realizamos realice una prueba t estándar independiente de 2 muestras que asume varianzas poblacionales iguales.**
```
print(np.var(M)/np.var(pM))
# 1.1435057273354379
```
3.   Con este codigo vamos a generar un DataFrame que contenga el valor estadistico y el p-value. Utilizamos *equal_var = True* porque asumimos varianzas poblaciones iguales.
```
ttest_data = pd.DataFrame(stats.ttest_ind(a = PM, b = M, equal_var = True)).transpose()
ttest_data.rename(columns = {
    0: "statistic",
    1: "pvalue",
},inplace=True)
```

###### Finalmente para incorporar estos valores a nuestra tabla, primero restablecemos el index para poder agregar la columna del DataFrame *ttest_data* a nuestro DataFrame *cpm_df*. 

```
cpm_df = cpm_df.reset_index()
cpm_df["p-val"] = ttest_data["pvalue"]
```

# Visualización

###### Hay muchos genes que tienen un p-value mayor a 0.05 por lo tanto no se estan expresando diferencialmente. Por lo tanto, resulta útil hacer un recorte de los datos. Entonces ejecutaremos:
```
viz_df = cpm_df[(cpm_df["p-val"] < 0.075)]
```
###### A continuación vamos a ordenar los datos para que los genes con el p-value más bajo nos queden en las primeras filas. Y a partir de DataFrame vamos a recuperar el Geneid de los 10 primeros para una posterior búsqueda. 
```
viz_sort = viz_df.sort_values("p-val")
list(viz_sort["Geneid"].iloc[0:10])
```
###### Recuento de genes expresados diferencialmente en las células migratorias con respecto a las pre migratorias

```
Genes_up = len(cpm_df[(cpm_df['p-val'] < 0.05) & (cpm_df['log2FC'] > 1)])
Genes_up
#217
Genes_down = len(cpm_df[(cpm_df['p-val'] < 0.05) & (cpm_df['log2FC'] < -1)])
Genes_down
#278
Genes_no = len(cpm_df[(cpm_df['p-val'] < 0.05) & (cpm_df['log2FC'] > -1) & (cpm_df['log2FC'] < 1)])
Genes_no
#577
```
```
Numeros_DE = pd.DataFrame()
Numeros_DE["Genes"] = ["Genes up", "Genes down", "Genes no dif"]
Numeros_DE["Recuento"] = [217, 278, 577]
Numeros_DE.plot.bar(x = "Genes")
```
###### Por último, realizamos un volcano plot para ver la distribución de nuestros datos relacionando el log2FC con el p-value de los genes. 
```
visuz.GeneExpression.volcano(df=viz_sort, lfc='log2FC', pv='p-val',sign_line=True, show=True, geneid="Geneid", genenames= ('ENSGALG00000008732','ENSGALG00000020488','ENSGALG00000008999','ENSGALG00000000620'))
```
# Resultados

###### Despues de identificar los Ensembl Id de los genes que mostraron expresión diferencial, buscamos en la página de Ensembl y de Uniprot más información sobre ellos y si estaban vinculados a algun proceso que ocurre en estas células cuando pasan de un estado  otro.

###### Las paginas usadas fueron estas:

*   https://may2015.archive.ensembl.org/Gallus_gallus/Info/Index
*   https://www.uniprot.org/ 
*   https://www.ncbi.nlm.nih.gov/

###### Genes Up regulados: 
###### ENSGALG00000008732 = OTOR (Otoraplin/Fdp/CDRAP/MIA)
> https://pubmed.ncbi.nlm.nih.gov/10998416/
https://pubmed.ncbi.nlm.nih.gov/10873378/
###### ENSGALG00000000620 = FABP3 
> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3742220/

###### Genes Down regulados: 
###### ENSGALG00000020488 = C11orf24 (DM4E3)
> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6755966/
###### ENSGALG00000008999 = AFAP1L2 
> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9835296/
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2822450/

###### Las células de la cresta neural pasan de un estado pre migratorio a uno migratorio para continuar con el desarrollo normal embrionario a traves de un proceso llamado transición epitelial a mesenquimatosa (EMT). 
###### La EMT implica en las células una pérdida de la polaridad apico-basal, ocurren regulaciones negativas sobre la expresión de proteínas de adhesión y unareorganización del citoesqueleto para lograr los movimientos celulares que darán lugar a las capas germinales embrionarias y posteriormente a los tejidos.
###### Este proceso de transición no solo ocurre en las células embrionarias, como la cresta neural, sino que se ha descripto que las células cancerosas imitan o se valen de estos procesos para favorecer la metastasisi y el establecimientos de tumores secundarios [Nieto et al., 2016; Lamouille et al., 2014; Pastushenko et al., 2019]

###### Los genes OTOR (Fdp, CD-RAP/MIA) y FABP3 tienen su expresión aumentada en las células de la cresta neural migratoria con respecto a las pre migratorias, mientras que la de los genes AFAP1L1 y C11orf24 (DM4E3) se encuentra disminuida. 
###### Esto podría estar en concordancia como algunos descubrimientos que se reportaron de estos genes. A modo resumen, se puede mencionar que:
* MIA/CD-RAP, en una línea celular de melanoma humano causaba una disminución de la capacidad adhesiva celular. Por lo tanto, tenia la capacidad de promover la invasión tumoral.
* La represión de la expresión de FDP producía inhibición de la la condrogénesis en mesénquima durante el desarrollo de la cápsula ótica.
* FABP3 actuaba como gen tumor-supresor. Se vio en muestras de cancer de mama su expresión estaba disminuida pero si era agregado exógenamente por transfección mostraba una actividad anti-proliferativa. acts as a tumor suppressor gene: its expression was downregulated in breast cancer samples;
* Por otra parte, experimentos knockdown de FABP3 en zebra fish demostraron que causaba un deterioro significativo del desarrollo cardíaco.
* AFAP1L1 promueve la proliferación, migración, invasión de células in vitro y crecimiento tumoral, metástasis in vivo al inducir la transición epitelio mesénquima (EMT).
* A pesar de que no hay mucha evidencia científica se puede asociar que la caída de C11orf24 y TMED3, dos proteínas asociadas al retículo endoplasmático y a la misma red secretora, llevaba la represión de la señalización WNT canónica y esto podría ser suficiente para impulsar la metástasis.

#### En conclusión, un aumento o una represión de expresión en los genes correctos conduce al éxito de las regulaciones moleculares y celulares que deben ocurrir en las células para garantizar el desarrollo adecuando. Los genes mencionados podrían tener un rol fi¡uncional en la inducción de la EMT, asegurar la migración y mantener el estado mesenquimal de las células de la cresta neural hasta que llegue el momento de la diferenciación.
