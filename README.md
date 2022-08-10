# Método para el agrupamiento dinámico de celdas - DyTRA 

<!---
Información de celdas resúmenes (contiene información de trayectorias GPS, y los rangos de velocidadd promedio)
-->

## Requisitos
+ Base de datos PostgreSQL con tablas que incluyan datos de trayectorias vehiculares
+ Los datos deben contener los campos: latitud, longitud, velocidad, tiempo, id_trayectoria.
+ RStudio 
+ R 4.1+

## Descripción del método

+ Representación del flujo vehicular

  En base al área geográfica donde se analizarán los datos, se realiza una división del área en zonas mas pequeñas denominadas *celdas* con un tamaño de 200x200 metros, la metodología realiza el análisis de lo que va ocurriendo en cada una de las celdas para considerar la representación de cada área en vez de cada trayectoria vehicular por separado. Este proceso se realiza de manera secuencial, estas celdas consolidan la información de todos los puntos GPS que recaen sobre el área de cada celda y extrae la información resumen de los puntos que se ubican dentro de los límites de las celdas analizadas.

  ![Celdas_Batch](https://github.com/gary-reyes-zambrano/Dynamic-grouping-of-vehicle-trajectories--DyTRA-/blob/master/images/Celdas_Batch.png?raw=true)

+ Agrupamiento basado en distancia

  Las celdas de cada ciclo se agrupan en una etapa basada en distancia que considera la menor diferencia de velocidad existente entre el grupo y la celda analizada para que la celda forme parte de ese grupo y en cuyo caso la información de este grupo se actualiza; si no existiese algún grupo cercano se formará un nuevo microcluster con la celda analizada. 

+ Agrupamiento basado en densidad

  Cuando finaliza la etapa basada en distancia las agrupaciones resultantes son evaluadas en una etapa basada en densidad y según su respectivo valor de densidad se clasifican en grupos densos y en grupos poco densos.

## Pasos para la ejecución
### Preparación del entorno
+ Se debe establecer las parametrizaciones con las que se desea ejecutar el algoritmo en el script principal (`DyTRA.R`):
  + Credenciales del repositorio (PostgreSQL)
  + Elección del Conjunto de datos a recuperar (previamente debe estar creadas las sentencias para la realización de la consulta) (`location`)
  + Mínimos y máximos de la dimensión "Velocidad" a analizar en el agrupamiento (`min_v`,`max_v`)
  + Mínimos y máximos de la dimensión "Velocidad" para la representación de los mapas (`vv_min`,`vv_max`)
  + `relativeSize`: valor porcentual a partir de los maximos y minimos usados en el agrupamiento para definir los rangos de velocidades 
  + `tiempo_a_ejecutar_algoritmo`: Tiempo de duración de cada ciclo (en minutos)
  + `rango`: Tamaño de la celda en largo y ancho (en grados decimales, se debe realizar la conversión si fuera necesario)
  + `PER`: Numero del periodo en el que se realizará el análisis (Cada valor considera 5 ciclos, si se establece `tiempo_a_ejecutar_algoritmo` en 3 minutos la duración de cada periodo durará 15 minutos, cada valor representa un periodo distinto; para PER=1 desde 0:00:00 hasta 0:15:00; para PER=2 desde 0:15:00 hasta 0:30:00...)

Nota: El script contiene valores por defecto para esats variables

### Ejecución del agrupamiento para cada periodo de tiempo
+ Habiendo establecido las parametrizaciones y teniendo seleccionado un numero de periodo válido (no sea inferior ni superior al periodo de tiempo que abarca los datos) se puede realizar la ejecución de todas las sentencias del script principal ( `DyTRA.R`), las cuales incluirán:
  + Cargado en memoria de las funciones a utilizar
  + Cargado de las parametrizaciones
  + Consulta a la base de datos
  + Creación de la ruta de almacenamiento para los resultados (Por defecto `C:/Pruebas-R/Algorimto-Fecha/Dataset/Hora`)
  + Calculo de la información de velocidad para los conjuntos de datos que no cuentan con esta información
  + División del área de procesamiento en celdas
  + Consolidación de las celdas
  + Selección de los datos del periodo seleccionado
  + Agrupamiento basado en distancia
  + Agrupamiento basado en densidad
  + Exportación de los resultados de las agrupaciones
  + Exportación de métricas de rendimiento

### Tratamiento de agrupaciones poco densas para la representación en los mapas

+ Teniendo una sesión de R cargada en memoria, abrir el script secundario `Generacion de mapas.R`
+ Para la visualización de los resultados se procesan las agrupaciones que resultaron poco densas, se utiliza todos los resúmenes de las agrupaciones resultantes y se aplican algunos criterios de unión, que por defecto son:
1. Se los ordena según su velocidad
2. Se excluyen las agrupaciones atípicas con un único punto GPS y que a su vez se encuentren sin mas grupos atípicos contiguos. 
3. Aquellas agrupaciones atípicas que se encuentren con mas agrupaciones atípicas contiguas, se seleccionan para considerarlas como un único grupo, esta unión se debe realizar calculando de manera ponderada según la cantidad de celdas de cada grupo para determinar la velocidad promedio del nuevo grupo y la información de cantidad de puntos, celdas y vehículos se totaliza.  

+ Se pueden considerar otros métodos de tratamientos para agrupaciones poco densas
+ La información de la velocidad así como el número de los grupos poco densos que se agruparán se deben agregar al script secundario de manera manual y este proceso se realiza para cada nueva agrupación (poco densa) que se desee agregar en el mapa

+ Si se omite este paso, la representación de los resultados en el mapa omitirán las agrupaciones poco densas 

### Representación y visualización en los mapas

+ Independientemente si se han tratado las agrupaciones poco densas o no, este proceso tomará la información de las agrupaciones densas (y poco densas según se haya agregado en el proceso anterior) y se representará en el mapa con un color de acuerdo a la escala de su velocidad
+ Se generará un fichero `.html` por periodo (por cada ejecución) donde se puede visualizar el resultado de los agrupamientos con la información acumulada desde el primer ciclo hasta el último ciclo

![img3-info-mapaA](https://github.com/gary-reyes-zambrano/Dynamic-grouping-of-vehicle-trajectories--DyTRA-/blob/master/images/img3-info-mapaA.png?raw=true)

![img3-info-mapaB](https://github.com/gary-reyes-zambrano/Dynamic-grouping-of-vehicle-trajectories--DyTRA-/blob/master/images/img3-info-mapaB.png?raw=true)

