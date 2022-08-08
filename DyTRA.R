#V20.2-RB16.1 - Dyclee por Celdas solo Vel
  #Clear Global Environment.####
  #' Clear Global Environment.
  #'
  #' @param ls() A list.
  #' @function "options" Permita que el usuario establezca y examine una variedad
  #'  de opciones globales que afectan la forma en que R calcula y muestra sus
  #'   resultados.
  #' @param digits  A number. permite modficar la cantidad de decimales que R 
  #' calculara.
  rm(list = ls())
  options( digits = 15 ) ## 
  
  #==================================================================================================
  #Class -  CF####
  #'Vector de rasgo caracteristico CF_k
  #'
  #' @param n A number: numero de elementos.
  #' @param LS A vertor: vector que contiene la suma lineal de cada
  #' caracterÃ­stica sobre los n objetos.
  #' @param SS A vector: la suma cuadrada de caracterÃ­sticas sobre los n objetos.
  #' @param tl A number: hora en la que se asignÃ³ el Ãºltimo objeto al ?? racimo.
  #' @param ts A timeStamp: tiempo cuando el ?? se creÃ³ el clÃºster.
  #' 
  #' @return la creacion de CF en el correspodiente IDk de los Uclusters.
  CF <- function(n, LS, SS, tl, ts){
    n = n
    LS = LS
    SS = SS
    tl = tl
    ts = ts
    self.CF <- list(n, LS, SS, tl, ts)
    names(self.CF) <- c("n", "LS", "SS", "tl", "ts")
    return(self.CF)
  }
  
  #Class -  Timestamp class####
  Timestamp <- function(){
    self.Timestamp <- list(timestamp=0)
    return(self.Timestamp)
  }
  
  #Class -  Micro Cluster####
  MicroCluster <- function(hyperboxSizePerFeature, currTimestamp, point, pos_Tray, data){
    self.MicroCluster <- list(currTimestamp, # timestamp object,
                              hyperboxSizePerFeature,
                              CF=initializeCF(currTimestamp,point),
                              label = -1,
                              previousCentroid = list(),
                              id_celds = pos_Tray,
                              data = data,
                              tramo = point,
                              idMC=idMC,
                              prev_label=0,
                              prev_elem=0,
                              log=data.frame()
    )
    names(self.MicroCluster) <- c("currTimestamp",
                                  "hyperboxSizePerFeature",
                                  "CF",
                                  "label",
                                  "previousCentroid",
                                  "id_celds",
                                  "data",
                                  "tramo",
                                  "idMC",
                                  "prev_label",
                                  "prev_elem",
                                  "log"
    )
    return(self.MicroCluster)
  }
  
  #Encapsulameinto de las dimensiones
  zip <-  function(point1, point2){
    count <- min(length(point1),length(point2))
    unif <-list()
    for (i in 1:count) {
      unif[[i]] <- cbind(point1[,i],point2[,i])
    }
    return(unif)
  }
  
  #Inicializar el vector de rasgos caractesisticos
  initializeCF <- function(currTimestamp, point){
    # we assume point is a list of features
    #create Variable LS = point central  of segment
    LS = data.frame(longitude= mean(point[,c("longitude")]),
                    latitude= mean(point[,c("latitude")]))
    res_zip <- zip(point[,c("longitude","latitude")], 
                   point[,c("longitude","latitude")])
    multi <- data.frame(matrix(data = 0,nrow = 1, ncol = length(res_zip)))
    for (ab in 1:length(res_zip)) {
      multi[,ab] <-  (mean(res_zip[[ab]][,1]) * mean(res_zip[[ab]][,2]))
    }
    # this vector will only have point elements squared
    SS = multi
    now = currTimestamp
    # CF creation
    cf = CF(n=1, LS=LS, SS=SS, tl=now, ts=now)
    return (cf)
  }
  
  #Obtener el centroide del grupo
  getCentroid <- function(microCluster){
    if(length(microCluster)!=0){
      centroid = data.frame(matrix(data = 0,
                                   nrow = 1,
                                   ncol = length(microCluster$CF$LS[,c("longitude","latitude")])))
      colnames(centroid) <- c("longitude","latitude")
      # para cada uC
      for (i in 1:length(centroid)){
        centroid[1,i] <- as.numeric(microCluster$CF$LS[,c("longitude","latitude")][,i]/ microCluster$CF$n)
      }
    }
    return (centroid )
  }
  
  #Adiciona un nuevo elemento y actualiza el grupo
  addElement <- function( microCluster, point, lambd, pos_Tray, data){
    res_decayComponent = decayComponent(microCluster,lambd)
    microCluster <- updateN(microCluster,res_decayComponent, point)
    microCluster <- updateLS(microCluster,res_decayComponent, point)
    microCluster <- updateSS(microCluster,res_decayComponent, point)
    microCluster <- updateTl(microCluster)
    microCluster$id_celds <- c(microCluster$id_celds, pos_Tray)
    microCluster$data <- rbind(microCluster$data, data)
    microCluster$tramo <- rbind(microCluster$tramo, point)
    return(microCluster)
  }
  
  # Componente de olvido
  decayComponent <- function(microCluster, lambd){
    r_dt = as.numeric(microCluster$CF$tl) - as.numeric(microCluster$currTimestamp$timestamp)
    var = 2.71828182845904 ^ (- lambd * r_dt)
    return(var)
  }
  
  # Actualizacion del vector de rasgos caracteristicos (N = Elementos)
  updateN <- function(microCluster, decayComponent, point=NULL){
    N = microCluster$CF$n * decayComponent
    if (!is.null(point)){
      N = N + 1
    }
    microCluster$CF$n = N
    return(microCluster)
  }

  # Actualizacion del vector de rasgos caracteristicos (LS = Suma lineal)
  updateLS <- function(microCluster, res_decayComponent, point=NULL){
    # forget
    for(i in 1:length(microCluster$CF$LS[,c("longitude","latitude")])){
      microCluster$CF$LS[,c("longitude","latitude")][,i] = microCluster$CF$LS[,c("longitude","latitude")][,i] *
        res_decayComponent
    }# add element
    if (!is.null(point)){
      for (i in 1:length(point[,c("longitude","latitude")])){
        microCluster$CF$LS[,c("longitude","latitude")][,i] = get_sum_vectors(microCluster$CF$LS[,c("longitude","latitude")][,i],
                                                                             point[,c("longitude","latitude")][,i])
      }
    }
    return(microCluster)
  }
  
  # Actualizacion del vector de rasgos caracteristicos (SS = Suma cuadratica)
  updateSS <- function(microCluster, res_decayComponent, point=NULL){
    # forget
    for (i in 1:length(microCluster$CF$SS)){
      microCluster$CF$SS[[i]] = mean(microCluster$CF$SS[[i]] * res_decayComponent)
    }# add element
    if (!is.null(point)){
      for (i in 1:length(point[,c("longitude","latitude")])){
        microCluster$CF$SS[,i] = get_sum_vectors(microCluster$CF$SS[,i],(point[,c("longitude","latitude")][,i] **2))
      }
    } 
    return(microCluster)
  }
  
  # Actualizacion del vector de rasgos caracteristicos (Tl = tiempo del elemento que ingresó)
  updateTl <-  function(microCluster){
    microCluster$CF$tl = as.numeric(Sys.time())
    return(microCluster)
  }
  
  # Determina si el elemento se enceuntra dentro del rango del grupo analizado
  isReachableFrom <-  function(microCluster, point){
    myCentroid = getCentroid(microCluster)
    maxDiff = -Inf
    featureIndex = 0
    boolean_point <- c()
    # for each feature
    for (i in 1:length(point[,c("longitude","latitude")])){
      # difference between the element feature and the cluster centroid for that feature
        boolean_point <- c(boolean_point, abs(mean(point[,c("longitude","latitude")][,i]) - myCentroid[i]) < 
        (microCluster$hyperboxSizePerFeature[i] / 2))
    }
    if (sum(boolean_point)==2){
      return (T)
    }# the element fits the u cluster
    return (F)
  }
  

  #Class -  Dyclee####
  
  Dyclee <- function( 
    dataContext, 
    relativeSize = 0.6,
    speed = Inf, 
    uncommonDimensions = 0,
    lambd = 0, 
    periodicRemovalAt = Inf,
    periodicUpdateAt = Inf, 
    timeWindow = 5, 
    findNotDirectlyConnButCloseMicroClusters = T,
    closenessThreshold = 1.5){
    
    # hyper parameters
    self.Dyclee <- list(
      relativeSize = relativeSize,
      uncommonDimensions = uncommonDimensions,
      processingSpeed = speed,
      lambd = lambd,
      periodicUpdateAt = periodicUpdateAt,
      timeWindow = timeWindow,
      periodicRemovalAt = periodicRemovalAt,
      findNotDirectlyConnButCloseMicroClusters = findNotDirectlyConnButCloseMicroClusters,
      closenessThreshold = closenessThreshold,
      dataContext = dataContext, # must be a bounding box instance
      # define hyperboxSizePerFeature
      hyperboxSizePerFeature = getHyperboxSizePerFeature(),
      # internal vis
      aList = list(),
      oList = list(),
      processedElements = 0,
      timestamp = 0,
      currTimestamp = Timestamp(), # initialized at 0
      densityMean = 0,
      densityMedian = 0)
    return(self.Dyclee)
  }
  
  #Stage1####
  
  # Obtener los rangos maximos para el centroide
  getHyperboxSizePerFeature <- function(dyclee){
    hyperboxSizePerFeature = c()
    for( context in 1:length(dataContext)){
      aux = dataContext[[context]]['maximun'] -  dataContext[[context]]['minimun']
      hyperboxSizePerFeature <- c(hyperboxSizePerFeature , (relativeSize * abs(aux)))
    }
    return (hyperboxSizePerFeature)
  }
  
  timeToIncTimestamp <- function(dyclee){
    return((dyclee$processedElements %% dyclee$processingSpeed) == 0)
  }
  
  # modifies reareachableMicroClusters iterating over a given list of u clusters
  # Obtener la lista de microclusters a los que puede alcanzar el elemento analizado
  getReachableMicroClustersFrom <- function(dyclee, microClustersList, point){
    #microClustersList= dyclee$oList
    res <- list()
    if (length(microClustersList)!=0) {
      count <- 1
      for (p_microCluster in 1:length(microClustersList)){
        if (isReachableFrom(microClustersList[[p_microCluster]], point)){
          res[[count]] <- microClustersList[[p_microCluster]]
          count <-  count+1
        }
      }
    }
    return(res)
  }

  # Obtener la lista de microclusters a los que puede alcanzar el elemento analizado
  # Grupos densos obligatorio; grupos poco densos opcional
  findReachableMicroClusters <- function(dyclee, point){
    reachableMicroClusters = getReachableMicroClustersFrom(dyclee, dyclee$aList, point)
    if (length(reachableMicroClusters)==0) {
      reachableMicroClusters = getReachableMicroClustersFrom(dyclee, dyclee$oList, point)
    }
    return(reachableMicroClusters)
  }
  
  # Evaluar el listado de grupos cercanos y determinar el que tnega menor distancia
  findClosestReachableMicroCluster <- function(dyclee, point, reachableMicroClusters){
    closestMicroCluster = NULL
    minDistance = Inf
    for (pos_microCluster in 1:length(reachableMicroClusters)){
      distance = dist(x = rbind(point[,c("longitude","latitude")],getCentroid(reachableMicroClusters[[pos_microCluster]])), method = "manhattan")
      if (distance < minDistance){
        minDistance = distance
        closestMicroCluster = reachableMicroClusters[[pos_microCluster]]
      }
    }
    return(closestMicroCluster)
  }
  
  # Procesar un elemento para determinar a que grupo va a pertenecer o si debe crear un nuevo microcluster
  processPoint <- function(dyclee, point, pos_Tray, info_c){
    # ASSUMPTION: point is a list of floats
    # find reachable u clusters for the new element
    reachableMicroClusters = findReachableMicroClusters(dyclee, point)
    if (length(reachableMicroClusters)==0){
      # empty list -> create u cluster from element
      # the microCluster will have the parametrized relative size, and the Timestamp object to being able to access the
      # current timestamp any atime
      count_uC <- length(dyclee$oList)+1
      idMC <<- idMC + 1
      microCluster = MicroCluster(dyclee$hyperboxSizePerFeature, dyclee$currTimestamp, point,pos_Tray, data = info_c)
      dyclee$oList[[count_uC]] <- microCluster
    }else{
      # find closest reachable u cluster
      closestMicroCluster = findClosestReachableMicroCluster(dyclee, point, reachableMicroClusters)
      if (length(dyclee$oList)>0) {
        for (comp in 1:length(dyclee$oList)) {
          if (isTRUE(compare(model = dyclee$oList[[comp]],comparison = closestMicroCluster))) {
            dyclee$oList[[comp]]  <-  addElement( closestMicroCluster, point=point, lambd=dyclee$lambd, pos_Tray, data = info_c)          
          }
        }
      }
      if (length(dyclee$aList)>0) {
        for (comp in 1:length(dyclee$aList)) {
          if (isTRUE(compare(model = dyclee$aList[[comp]],comparison = closestMicroCluster))) {
            dyclee$aList[[comp]]  <-  addElement( closestMicroCluster, point=point, lambd=dyclee$lambd, pos_Tray, data = info_c)
          }
        }
      }
    }
    # at this point, self self.aList and self.oList are updated
    return(dyclee)
  }
  
  timeToCheckMicroClustersTl <- function(self){
    return (self$processedElements %% (self$periodicUpdateAt * self$processingSpeed) == 0)
  }
  
  timeToPerformPeriodicClusterRemoval <- function(self){
    return (self$processedElements %% (self$periodicRemovalAt * self$processingSpeed) == 0)
  }
    
  trainOnElement <- function(dyclee, info_c){#point, pos_Tray){
    # get an object matching the desired format (to guarantee consistency)
    dyclee$processedElements = dyclee$processedElements + 1
    dyclee$timestamp = as.numeric(Sys.time())
    dyclee$currTimestamp$timestamp = dyclee$timestamp
    point <- data.frame(longitude = info_c$vel_prom,
                        latitude = info_c$vel_prom,
                        unixtime = info_c$unixtime,
                        file_name = info_c$cant_veh,
                        velocidad = info_c$vel_prom)  
    # now, check what to do with the new point
    dyclee <- processPoint(dyclee, point, info_c$proc_id, info_c)
    
    return(dyclee)
  }
  
  #Stage2####
  
  # Etapa 2: basada en densidad
  getClusteringResult <- function(dyclee){
    # update density mean and median values with current ones
    dyclee <- calculateDensityMeanAndMedian(dyclee)
    # rearrange lists according to microClusters density, considering density mean and median limits
    dyclee <- rearrangeLists(dyclee)
    # form final clusters
    dyclee$aList <- formClusters(dyclee)    
    dyclee <- updateMicroClustersPrevCentroid(dyclee)
    return(dyclee)
  }
  
  # Actualiza el centoride previo
  updateMicroClustersPrevCentroid <- function(dyclee){
    if(length(dyclee$aList)>0){
      for (microCluster in 1:length(dyclee$aList)){
        dyclee$aList[[microCluster]]$previousCentroid = getCentroid(dyclee$aList[[microCluster]])
      } 
    }
    return(dyclee)
  }

  # Calcular la densidad media y mediana de todos los grupos densos y poco densos
  calculateDensityMeanAndMedian <- function(dyclee){
    concatenatedLists = append(dyclee$aList, dyclee$oList)
    dyclee$densityMean <- calculateMeanFor(concatenatedLists)
    dyclee$densityMedian <- calculateMedianFor(concatenatedLists)
    return(dyclee)
  }
  
  #Calcular la media de densidades
  calculateMeanFor <- function(concatenatedLists){
    acc_mean <- c()
    for (microCluster in concatenatedLists) {
      acc_mean <- c(acc_mean, getD(microCluster))
    }
    return(mean(acc_mean))
  }
  
  # Calcular la mediana de densidades
  calculateMedianFor <- function(concatenatedLists){
    acc_median <- c()
    for (microCluster in concatenatedLists) {
      acc_median <- c(acc_median, getD(microCluster))
    }
    return(median(acc_median))
  }
  
  #Calcular la densidad de una grupo
  getD <- function(microCluster){
    V <- prod(microCluster$hyperboxSizePerFeature)
    return(microCluster$CF$n / V)
  }

  # Reasigna las etiquetas d elos grupos
  rearrangeLists <- function(dyclee){
    newAList = list()
    newOList = list()
    concatenatedLists = append(dyclee$aList, dyclee$oList)
    concatenatedLists <- concatenatedLists[!duplicated(concatenatedLists, fromLast = TRUE)]
    countAList <- 1
    countOList <- 1
    if(length(concatenatedLists)>0){
      for (microCluster in 1:length(concatenatedLists)){
        if (isTRUE(isOutlier(dyclee, concatenatedLists[[microCluster]]))){
          newOList[[countOList]] <- concatenatedLists[[microCluster]]
          countOList <- countOList + 1
        }else{
          newAList[[countAList]] <- concatenatedLists[[microCluster]]
          countAList <- countAList + 1
        }
      }
    }
    dyclee$aList = newAList
    dyclee$oList = newOList
    return(dyclee)
  }
  
  isOutlier <- function(dyclee, microCluster){
    return(getD(microCluster) < dyclee$densityMean && getD(microCluster) < dyclee$densityMedian)
  }
  
  isDense <- function(dyclee, microCluster){
    # returns true if a given u cluster is considered dense
    return (getD(microCluster) >= dyclee$densityMean && getD(microCluster) >= dyclee$densityMedian)
  }
  
  isSemiDense <- function(dyclee, microCluster){
    # returns true if a given u cluster is considered semi dense
    # xor
    return (getD(microCluster) >= dyclee$densityMean) != (getD(microCluster) >= dyclee$densityMedian)
  }
  
  # Formar grupos mas densos en caso de solapamiento de grupos  
  formClusters <- function(dyclee) {
    #reset microClusters labels as -1
    dyclee$aList <- resetLabelsAsUnclass(dyclee$aList)
    dyclee$oList <- resetLabelsAsUnclass(dyclee$oList)
    
    # start clustering
    #contendra la nueva lista de aList de los grupos finales
    alreadySeen = list()
    # init currentClusterId
    currentClusterId <- 0
    #uCk_Id <- 1
    tam <-  0
    if(length(dyclee$aList)>0){
      for (uCk_Id in 1:length(dyclee$aList)) {
        current_uCk <- dyclee$aList[[uCk_Id]]
        if (!has_element(.x = alreadySeen, .y = current_uCk)) {
          tam <- tam + 1
          alreadySeen[[tam]] <- current_uCk
          if (current_uCk$label == -1) {
            currentClusterId <- currentClusterId + 1
            alreadySeen[[tam]]$label <- currentClusterId
            dyclee$aList[[uCk_Id]]$label <- currentClusterId# <- NULL
          }
          connected_uCk <- findSimilarMicroClustersFor(dyclee, dyclee$aList[[uCk_Id]], dyclee$aList)
          i <-  1
          while (i < length(connected_uCk)) {
            neighbor <- connected_uCk[[i]]
            if (!has_element(.x = alreadySeen, .y = neighbor)) {
              if (isDense(dyclee, neighbor) || isSemiDense(dyclee, neighbor)) {
                #|| isSemiDense(dyclee, neighbor)
                neighbor$label <- currentClusterId
                tam <- tam + 1 # length(alreadySeen) + 1
                alreadySeen[[tam]] <- neighbor
                for (id in 1:length(dyclee$aList)) {
                  if (compare.list(a = list(dyclee$aList[[id]]),
                                   b = list(connected_uCk[[i]]))) {
                    dyclee$aList[[id]]$label <- currentClusterId
                    #print(c(dyclee$aList[[id]]$vehicles, neighbor$vehicles))
                    #dyclee$aList[[id]]$vehicles <- unique(c(dyclee$aList[[id]]$vehicles, neighbor$vehicles))
                  }
                }
                newConnected_uCk <- findSimilarMicroClustersFor(dyclee, neighbor, dyclee$aList)
                if(length(newConnected_uCk)>0){
                  for (newNeighbor in 1:length(newConnected_uCk)) {
                    if (!has_element(.x = alreadySeen, .y = newConnected_uCk[[newNeighbor]])) {
                      if (isDense(dyclee, newConnected_uCk[[newNeighbor]]) || isSemiDense(dyclee, newConnected_uCk[[newNeighbor]])) {
                        connected_uCk[[length(connected_uCk) + 1]] <- newConnected_uCk[[newNeighbor]]
                      }
                    }
                  }}
                connected_uCk <- connected_uCk[!duplicated(connected_uCk, fromLast = TRUE)]
              }
            }
            i <- i + 1
          }
        }
      }
    }
    
    return(alreadySeen)
  }
  
  findSimilarMicroClustersFor <- function(dyclee, microCluster, microClusters){
    directlyConn <-  findDirectlyConnectedMicroClustersFor(dyclee, microCluster, microClusters)
    if (!dyclee$findNotDirectlyConnButCloseMicroClusters){
      return (directlyConn)
    }else{
      notDirectlyConnButClose = findCloseMicroClustersFor(dyclee, microCluster, microClusters)
      return (append(directlyConn, notDirectlyConnButClose))
    }
  }
  
  # Buscar grupos muy cercanos que se esten solapando
  findCloseMicroClustersFor <- function(dyclee, microCluster, microClusters){
    stddevProportion = dyclee$closenessThreshold
    # for encompassing more micro clusters
    encompassing = getAvgDistToMicroClustersFor(dyclee, microCluster, microClusters)
    avgDistToAllMicroClusters <- encompassing[[1]]
    distances <- encompassing[[2]] 
    stdev = stddev(distances, avgDistToAllMicroClusters)
    limit = avgDistToAllMicroClusters - (stdev * stddevProportion)
    # the set of close micro clusters which will be used to expand a macro one
    res = list()
    ac <- 1
    for (mc in microClusters) {
      if (length(mc)!=0 && length(microCluster)!=0) {
        #Boolean
        mcIsClose = distanceTo(microCluster, mc) < limit
        if (!(isDirectlyConnectedWith(microCluster, mc, dyclee$uncommonDimensions)) && mcIsClose){
          res[[ac]] <-  mc
          ac <-  ac + 1
        }
      }
    }
    return(res)
  }
  
  # Calcular desviacion de los elementos con respecto a los centroides de los grupos
  stddev <- function(data, meanData, ddof=0){
    #' """Calculates the population standard deviation
    #'   by default; specify ddof=1 to compute the sample
    #'   standard deviation."""
    n = length(data)
    if( n < 2){
      print('variance requires at least two data points')
    }
    ss = get_ss(data, meanData)
    pvar = ss/(n-ddof)
    return (pvar^0.5)
  }
  
  get_ss <- function(data, meanData){
    #'
    #' """Return sum of square deviations of sequence data."""
    #' 
    ss <- 0
    for (x in data) {
      ss <- ss + ((x - meanData)^2)
    }
    return (ss)
  }
  
  # Distancia promedio entre dos grupos
  getAvgDistToMicroClustersFor <- function(dyclee, microCluster, microClusters){
    sum = 0
    dists = list()
    ac = 1
    
    for (mc in microClusters){
      if (length(microCluster)!=0 && length(mc)!=0 ) {
        dist = distanceTo(microCluster, mc)
        sum = sum + dist
        dists[[ac]] <- dist
        ac = ac + 1
      }
    }
    return (list(sum/length(microClusters), dists))
  }
  
  distanceTo <- function(microCluster, mc){
    distance = dist(x = rbind(getCentroid(microCluster),getCentroid(mc)),
                    method = "manhattan")
    return (distance)
  }
  
  #Buscar grupso directamente conectados
  findDirectlyConnectedMicroClustersFor <- function(dyclee, microCluster, microClusters){
    res = list()
    ac <- 1
    for (mc in microClusters){
      if (!isempty(mc) && !isempty(microCluster)) {
        if(!compare.list(a = list(mc), b = list(microCluster))){          
          if (isDirectlyConnectedWith(microCluster, mc, dyclee$uncommonDimensions)){
            res[[ac]] <- mc
            ac <- ac + 1
          }
        }
      }
    }
    return (res)
  }
  
  isDirectlyConnectedWith <- function(microCluster, mc, uncommonDimensions){
    #uncommonDimensions <- dyclee$uncommonDimensions
    # retunrs true if the microCluster is directly connected to another microCluster
    featuresCount = length(microCluster$CF$LS[,c('longitude','latitude')])
    currentUncommonDimensions = 0
    myCentroid = getCentroid(microCluster)
    microClusterCentroid = getCentroid(mc)
    currentUncommonDimensions <-  dist(x = rbind(myCentroid,microClusterCentroid),
                                       method = "manhattan")
    #currentUncommonDimensions <-get_distance_Hausdorff(tray_obj1 = myCentroid, tray_obj2 = microClusterCentroid)

    if(is.nan(currentUncommonDimensions)){
      currentUncommonDimensions <-  0
    }
    return (currentUncommonDimensions <= uncommonDimensions)
    
  }
  
  # Reestablecer las etiquetas de los grupos
  resetLabelsAsUnclass<-function(microClusters){
    if (length(microClusters) >= 1) {
      for (microCluster in 1:length(microClusters)) {
        microClusters[[microCluster]]$label <- -1
      }
    }
    return(microClusters)
  }
  
  
  
  #Funciones - Manipulacion de datos#####
  
  #Crear directorio
  # C:/Pruebas-R/Dyclee-[Fecha:YYYY-MM-DD]/[Dataset]/[Hora: hh-mm-ss]/
  get_workSpace<- function(carpeta, algoritmo, dataset, angulo_tolerancia, relativeSize, tGlobal){
    algoritmo <- paste(algoritmo,Sys.Date(),sep = "-")
    simbol<-"/";
    root<- "C:";  	##directorio raiz
    setwd(paste("C:",simbol,sep=""));
    
    dir.create(carpeta);
    setwd(paste(getwd(),simbol,carpeta,sep=""));
    
    dir.create(algoritmo);
    setwd(paste(getwd(),simbol,algoritmo,sep=""));
    
    dir.create(dataset);						
    setwd(paste(getwd(),simbol,dataset,sep=""));

    hora <- format(Sys.time(), "%H-%M-%S")
    dir.create(hora);
    setwd(paste(getwd(),simbol,hora,sep=""))

    return(getwd())
  }
  
  #Cargar Librerias y dependencias a la sesion
  ipak <- function(pkg) { #cargar e instalar paquetes de manera dinamica
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
      install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
  }
  
  get_norm_lon_lat <- function(dataset){
    data_norm <-  scale(x = dataset[,1:2])
    colnames(data_norm) <- c('longitude','latitude')
    return(data.frame(data_norm))
  }
  
  
  #Funciones DataSet de DB 
  # Importar los datos desde una conexion PostgreSQL
  get_DatasetsFromDB <- function(location){
    #Crea DATA SOURCE
    print("CARGANDO LOS DATOS DE LA BD...")
    drv <- dbDriver(dbdriver)
    
    con <- dbConnect(drv, host = host, port = port, dbname = dbname, user = user, pass = pass)
    #Seleccion de los datos segun localidad
    if (location == 1) {
      time_in_timestamp = T
      big_int <- T
      dataset <- "Ecuador-Guayaquil"
      tz <- 'America/Guayaquil'
      data <-dbGetQuery(con, "SELECT longitud as longitude, latitud as latitude, fecha as unixtime, id as file_name FROM guayaquil WHERE  longitud BETWEEN  '-79.899' and '-79.876' and latitud BETWEEN  '-2.2' and '-2.13' AND to_char(to_timestamp(fecha/1000),'DD/MM/YY HH24:MI:SS')> '28/10/17 16:00:00' AND to_char(to_timestamp(fecha/1000),'DD/MM/YY HH24:MI:SS')<='28/10/17 18:30:00'")
    }else if (location == 5) { 
      time_in_timestamp = F
      big_int <- F
      dataset <- "Roma"
      tz <- 'Europe/Rome'
      #data <-dbGetQuery(con, "select lat as latitude, lon as longitude, time as unixtime, n_speed_kmh as velocidades, id as file_name from public.roma where time > '14:00:00' and time <= '15:00:00' and problem is null")
      data <-dbGetQuery(con, "select lat as latitude, lon as longitude, time as unixtime, n_speed_kmh as velocidades, id as file_name from public.roma where time > '12:00:00' and time <= '14:00:00' and lat >= 41.88 and lat <= 41.92 and lon >= 12.48 and lon <= 12.52 and problem is null")
      data$unixtime <- paste("2014-02-12",data$unixtime, sep=" ")
      data$unixtime <- as.POSIXct(data$unixtime)
    }else if (location == 6) { 
      time_in_timestamp = F
      big_int <- F
      dataset <- "China-Beijing"
      tz <- 'Asia/Shanghai'
      data <- dbGetQuery(con, "select latitud latitude, longitud longitude, fecha unixtime, substring(id, '[0-9]+') as file_name from public.beijing g WHERE longitud BETWEEN  '116.28' and '116.485' and latitud BETWEEN  '39.825' and '39.99' and fecha between '02/02/08 13:30:00' and '02/02/08 15:30:00' order by fecha asc") 
      
    }else if (location == 7) { 
      print("Read data Simulated")
      time_in_timestamp = F
      big_int <- F
      dataset <- "Data Simulated"
      tz <- 'America/Guayaquil'
      data<- read.delim(file = "E:/Proyectos Trayectorias/algoritmos de agrupacion/dyclee/Dyclee-python MOD Onofri 29-12-2020/DyClee-master/blobs.csv",
                        header = F,
                        sep = ",")
      
      data <- data.frame(longitude = data$V1, latitude = data$V2, unixtime=0, file_name=0)
      data <- get_norm_lon_lat(dataset = data)
      data <- data.frame(longitude = data$longitude, latitude = data$latitude, unixtime=0, file_name=0)
    }
    
    
    #Desconexion de DB
    dbDisconnect(con)
    if (data$unixtime != 0) {
      #Conversion de Col Tiempo si es TimeStamp  y BingInt
      if (isTRUE(time_in_timestamp) && isTRUE(big_int)) {
        data$unixtime <- as.POSIXct(data$unixtime/1000, origin="1970-01-01")
      }else if (isTRUE(time_in_timestamp) && !isTRUE(big_int)) {
        #Conversion de Col Tiempo si es TimeStamp
        data$unixtime <- as.POSIXct(data$unixtime, origin="1970-01-01") 
      }
      data$unixtime <-with_tz(data$unixtime, tzone = tz)
    }
    pto_medio <- matrix(0, nrow=1, ncol=2) #punto central en un plano carteciano X, Y
    min_latitud <- min(data$latitude)
    max_latitud <- max(data$latitude)
    min_longitud <- min(data$longitude)
    max_longitud <- max(data$longitude)
    pto_medio[2] <- (max_latitud + min_latitud)/2   ## valor de x para el centro de la ventana
    pto_medio[1] <- (max_longitud + min_longitud)/2 ## valor de y para el centro de la ventana
    
    result <- list(dataset,nrow(data),pto_medio,summary(data),data)
    names(result) <- c('dataset','size','central_point','summary','data')
    print(T)
    return(result)
  }
  
  # Generar escala de colores
  colors <- function(){
    # rainbow(n = 10)
    # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'seq',]
    # unlist(mapply(brewer.pal,                     # Create vector with all colors
    #               qual_col_pals$maxcolors,
    #               rownames(qual_col_pals)))
    col_vector =rainbow(n = 10)
    #  unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    return(col_vector)
  }
  
  #Funciones de tranformacion de datos
  
  get_list_per_trajectory <- function(dataset){
    #Separarlo en listas de tray
    ac <- 1
    lista_trayectorias <- list()
    #for(count in 1:nrow(dataset$data)){
    for(count in names_tray){
      tray_corresp <-subset(dataset$data, file_name == count)
      if (nrow(tray_corresp)>5) {
        lista_trayectorias[[ac]] <-  tray_corresp #dataset$data[count,]
        ac <- ac + 1
      }
    }
    return(lista_trayectorias)
  }
  
  get_list_trajectory <- function(lista){
    #Separarlo en listas de tray
    ac <- 1
    lista_trayectorias <- list()
    for(count in 1:nrow(lista)){
    #for(count in names_tray){
      #tray_corresp <-subset(dataset$data, file_name == count)
      #if (nrow(tray_corresp)>5) {
        lista_trayectorias[[ac]] <-  lista[count,] #tray_corresp
        ac <- ac + 1
      #}
    }
    return(lista_trayectorias)
  }
  
  get_pre_process_info <- function(trajectories){
    cant_trays <- length(trajectories)
    trayectorias_informacion<-data.frame( No_Trayectoria = c(1:cant_trays),#numero de la trayecotria o sub-tray 
                                          #vehiculo=0,
                                          microCluster = 0,#micro cluster al que a sido asignado - Etapa 1
                                          clusterFinal =0, #grupo final segun la densidad de los microclusters
                                          silhouette = 0,#criterio de calidad 
                                          vecino_cercano=0,# segmento del grupo mas cercano 
                                          dist_centroide=0,#distancia de la tray o seg de tray  con su centroide asignado
                                          dist_promedio_al_grupo = 0,#distancia promedio de la tray o seg tray con respecto a su grupo asignado
                                          velocidad_Prom = 0, #velocidad promedio de tray o seg tray
                                          tiempo_recorrido=0,
                                          cant_puntos = 0, 
                                          estado = 0,
                                          Men_Latitud =0,#menor latitud
                                          Long_Corresp1 = 0,#longitud correspondiente
                                          May_Latitud  =0 ,
                                          Long_Corresp2=0 ,
                                          Men_Longitud=0, 
                                          Lat_Corresp1 = 0,
                                          May_Longitud=0 ,
                                          Lat_Corresp2 = 0 , 
                                          Dist_ancho=0, 
                                          Dist_alto = 0
                                          )
    for (i in 1:cant_trays) {
      tray_actual <- trajectories[[i]]
      cant_puntos <- nrow(tray_actual)
      #trayectorias_informacion$vehiculo[i] <- tray_actual$file_name
      #ordeno de menor a mayor la longitud  xÂ´
      tray_order <- tray_actual[order(tray_actual$longitude,decreasing=FALSE),]
      trayectorias_informacion$Men_Longitud[i] <- tray_order$longitude[1]
      trayectorias_informacion$Lat_Corresp1[i] <- tray_order$latitude[1]
      ##la mayor es la ultima, pues esta ordenado de menor a mayor  
      trayectorias_informacion$May_Longitud[i] <- tray_order$longitude[cant_puntos]
      trayectorias_informacion$Lat_Corresp2[i] <- tray_order$latitude[cant_puntos]
      ##distancia ancho
      trayectorias_informacion$Dist_ancho[i] <- dist(rbind(tray_order[1,c('longitude','latitude')], 
                                                           tray_order[cant_puntos,c('longitude','latitude')]),
                                                     method = "euclidean")
      
      ##ordeno de menor a mayor la latitud  yÂ´
      tray_order <- tray_actual[order(tray_actual$latitude, decreasing=FALSE),]
      trayectorias_informacion$Men_Latitud[i] <- tray_order$latitude[1]
      trayectorias_informacion$Long_Corresp1[i] <- tray_order$longitude[1]
      ## la mayor es la ultima, pues esta ordenado de menor a mayor  
      trayectorias_informacion$May_Latitud[i] <- tray_order$latitude[cant_puntos]
      trayectorias_informacion$Long_Corresp2[i] <- tray_order$longitude[cant_puntos]
      ##distancia alto
      trayectorias_informacion$Dist_alto[i] <- dist(rbind(tray_order[1,c('longitude','latitude')],
                                                          tray_order[cant_puntos,c('longitude','latitude')]),
                                                    method = "euclidean")
      
      ##almaceno la cant de puntos de la tray
      trayectorias_informacion$cant_puntos[i] <- cant_puntos
      #Velocidad promedio de la tray o segmento
      trayectorias_informacion$velocidad_Prom[i] <- tray_actual$velocidades#mean(tray_actual$unixtime)# no tomar en cuenta porque toma el tiempo 
      #tiempo recorrido entre los puntos iniciales y finales
      #trayectorias_informacion$tiempo_recorrido[i] <- get_time_difference(self = tray_actual[cant_puntos,],
      #                                                                 previous = tray_actual[1,])
      
    }
    return(trayectorias_informacion)
  }
  
  get_calculate_angle <-  function (coord1 , coord2){
    ## angulo en grados
    anguloGrados <- (atan2(coord2[1,2] - coord1[1,2] , coord2[1,1] -coord1[1,1]) * 180) / pi
    return (anguloGrados)
  }
  
  get_segmented_trajectories <- function(segmentation_type, trajectories){
    numeracion <- 1
    cant_trayectorias<-length(trajectories)
    lista_trayectorias_segmentadas <- list()
    #evalua las trayectorias y la segementa segun un valor de tolernacia angular
    #se calcula el angulo entre P1 y P2 y se guarda como angulo de referfencia
    #se calcula el angulo entre P2 y P3, y se compara contra el angulo de referencia, si la diferencia absoluta de los angulos comparados supera la tolerancia (umbral) > se crea una subtrayectoria
    #en caso contrario, se continua evaluando los sigt puntos contra P1
    #si se crea una subtrayectoria con los puntos P1 y P2, entonces se procede a calcular de nuevo el angulo de referencia que sera el angulo entre P2 y P3
    if(segmentation_type == 1){
      for(i in 1:cant_trayectorias){
        tray_actual <- trajectories[[i]]
        cant_puntos <- nrow(tray_actual)
        pos_punto_inicial <- 1
        distancia_aux <- 0
        part <- 1
        if(cant_puntos > 2){
          pto_inicial <- tray_actual[1,1:2] 
          pto_sigt <- tray_actual[2,1:2]
          angulo_referencia <- get_calculate_angle(pto_inicial, pto_sigt)  
          for(j in 2: (cant_puntos-1)){
            pto_sigt <- tray_actual[j + 1,1:2]
            angulo_aux <- get_calculate_angle(pto_inicial, pto_sigt)  
            ##si el valor absoluto, de la resta es mayor que la tolerancia entonces se segmenta la trayectoria
            if(abs(angulo_referencia - angulo_aux) > angulo_tolerancia){
              pto_inicial <- tray_actual[j,1:2]
              angulo_referencia <- get_calculate_angle(pto_inicial, pto_sigt)##crear una lista de segmentos 
              if((j + 1) == cant_puntos ){
                tray_actual[pos_punto_inicial: cant_puntos,'file_name'] <- paste(tray_actual[j,'file_name'], "_part_", as.character(part), sep="_")
                lista_trayectorias_segmentadas[[numeracion]] <- tray_actual[ pos_punto_inicial: cant_puntos ,]
                pos_punto_inicial <- cant_puntos
              }else{
                tray_actual[pos_punto_inicial: j,'file_name'] <- paste(tray_actual[j,'file_name'], "_part_", as.character(part), sep="_")
                lista_trayectorias_segmentadas[[numeracion]] <- tray_actual[ pos_punto_inicial: j ,]
                pos_punto_inicial <- j
              }
              part <- part + 1
              numeracion <- numeracion +1 ##pos en la lista,
            }
          }
        }
        ##aniado el ultimo tramo de la trayectoria por si fue segementado
        if(pos_punto_inicial < cant_puntos){
          tray_actual[pos_punto_inicial: cant_puntos,'file_name'] <- paste(tray_actual[pos_punto_inicial,'file_name'], "_part_", as.character(part), sep="_")
          lista_trayectorias_segmentadas[[numeracion]] <- tray_actual[pos_punto_inicial: cant_puntos,]
          numeracion <- numeracion + 1 ##pos en la lista,
        }    
      }
    }
    return(lista_trayectorias_segmentadas) ## lista de trayectoria lo remplazo con la lista segmentada
  }
  
  get_amplitude <- function(geoCoordinates,centralpoint){
    amplitude <- data.frame(longitude=0, latitude=0)
    min_latitude  <-min(geoCoordinates$latitude)
    max_latitude  <-max(geoCoordinates$latitude)
    min_longitude <-min(geoCoordinates$longitude)
    max_longitude <-max(geoCoordinates$longitude)
    amplitude$longitude <-  abs(max_latitude - centralpoint[2]) 				##amplitud de la ventana desde el centro a sus aristas
    amplitude$latitude <-  abs(max_longitude - centralpoint[1])	 			##amplitud de la ventana desde el centro a sus aristas
    return(amplitude)
  }
  
  
  
  get_sum_vectors <- function(vector1, vector2){
    resp <- c(0)
    if (length(vector1)== length(vector2)) {
      resp <-  vector1 + vector2
    }else if(length(vector1) < length(vector2)){
      for (k in 1:length(vector1)) {
        resp[k] <-  vector1[k]+vector2[k]
        if (k ==length(vector1)) {
          resp[k] <- vector1[k]+mean(vector2[k:length(vector2)])
        }
      }
    }
    return(resp)
  }
  
  
  ##area donde enfocar el estudio (se selecciona solos los trayectorias q esten en dicha area)
  create_window <- function(pto_medio_aux, amplitud_x_aux, amplitud_y_aux){
    pto_izq_arr <- matrix(0, nrow=1, ncol=2)
    pot_izq_ab <- matrix(0, nrow=1, ncol=2)
    pto_der_arr <- matrix(0, nrow=1, ncol=2)
    pto_der_ab <- matrix(0, nrow=1, ncol=2)
    
    pto_izq_arr [1] <-  pto_medio_aux[1] - amplitud_x_aux
    pto_izq_arr [2] <-  pto_medio_aux[2] + amplitud_y_aux
    
    pot_izq_ab [1] <-   pto_medio_aux[1] - amplitud_x_aux
    pot_izq_ab [2] <-   pto_medio_aux[2] - amplitud_y_aux
    
    pto_der_arr [1] <-  pto_medio_aux[1] + amplitud_x_aux
    pto_der_arr [2] <-  pto_medio_aux[2] + amplitud_y_aux
    
    pto_der_ab [1] <-   pto_medio_aux[1] + amplitud_x_aux
    pto_der_ab [2] <-   pto_medio_aux[2] - amplitud_y_aux
    
    ventana <- matrix(0,nrow=4, ncol = 2)
    ventana[1,]<-pto_izq_arr
    ventana[2,]<-pot_izq_ab
    ventana[3,]<-pto_der_arr
    ventana[4,]<-pto_der_ab
    
    ##pintar ventana en grafica
    par(lwd=3)
    points(rbind(ventana[1,1:2],ventana[2,1:2]) ,type="b", col = "red")
    points(rbind(ventana[2,1:2],ventana[4,1:2]) ,type="b", col = "red")
    points(rbind(ventana[4,1:2],ventana[3,1:2]) ,type="b", col = "red")
    points(rbind(ventana[3,1:2],ventana[1,1:2]) ,type="b", col = "red")
    return (ventana)
  }
  
  #Funciones Algoritmo
  get_time_difference<- function(self, previous){
    #""" Calcultes the time difference against another point
    
    #     Args:
    #          previous (:obj:`Point`): Point before
    #      Returns:
    #          Time difference in seconds
    #      """
    return(difftime(self[,c('unixtime')],previous[,c('unixtime')],
                    units = "secs"))
  }
  
  get_distance_Hausdorff<- function(tray_obj1, tray_obj2){                       #depura los objetos para la Metrica hausdorff
    if(is.na(tray_obj1) || is.na(tray_obj2) || is.null(tray_obj1) || is.null(tray_obj2)){
      return(0);
    }else{
      #aux1<-as.numeric(as.matrix(tray_obj1[,1:2],nrow = nrow(tray_obj1),ncol = ncol(tray_obj1)))
      #aux2<-as.numeric(as.matrix(tray_obj2[,1:2],nrow = nrow(tray_obj2),ncol = ncol(tray_obj2)))
      #dist<- hausdorff_dist(aux1,aux2)   
      #dist<- hausdorff_distance(aux1,aux2)
      dist <- hausdorff_distance(lista=tray_obj1, point=tray_obj2)
    }
    return(dist)  
  }
  
  hausdorff_distance <- function(lista, point){
    # Formule: D = acos(sin(N1)*sin(N2) + cos(N1)*cos(N2)*cos(E1-E2) )
    #where distance and coordinates are expressed in radians. N1 and N2 is the latitude of origin and destination, and E1 and E2 is longitude.
    radianes = 180/pi
    distances = matrix(0, nrow = nrow(lista), ncol = 1)
    #Expresado en radianes
    # N1 = Latitud del origen   = o_lat
    # N2 = Latitud del destino  = d_lat
    # E1 = Longitud del origen  = o_lon
    # E2 = Longitud del destino = d_lon
    d_lat = as.numeric(point[1, 'latitude'])
    d_lon = as.numeric(point[1, 'longitude'])
    for (i in 1:nrow(lista)) {
      o_lat = as.numeric(lista[i, 'latitude'])
      o_lon = as.numeric(lista[i, 'longitude'])
      
      factor1 = sin(round(o_lat, 12) / radianes)
      factor2 = sin(round(d_lat, 12) / radianes)
      factor3 = cos(round(o_lat, 12) / radianes)
      factor4 = cos(round(d_lat, 12) / radianes)
      factor5 = round(o_lon, 12) / radianes
      factor6 = round(d_lon, 12) / radianes
      
      precalc1 = factor1 * factor2
      precalc2 = factor3 * factor4
      precalc3 = cos(factor5 - factor6)
      precalc4 = precalc2 * precalc3
      precalc5 = precalc1 + precalc4
      calculo = acos(precalc5)
      if(is.nan(calculo)){
        calculo = 0
      }
      #D = 60 * radianes * 
      #  acos( 
      #    sin(round(o_lat, 12) / radianes) * sin(round(d_lat, 12) / radianes) 
      #    + cos(round(o_lat, 12) / radianes) * cos(round(d_lat, 12) / radianes)
      #    * cos((round(o_lon, 12) / radianes) - (round(d_lon, 12) / radianes))
      #  ) * 1.852
      D = 60 * radianes * calculo * 1.852

      if(is.nan(D)){
        D = 0
      }
      distances[i,] <- D
    }
    return(max(distances))
  } 
  
  
  update_uClustersFinals_to_binnacle <- function(dyclee, trayectorias_informacion){
    if (length(dyclee$aList)>0) {
      for (uCk in 1:length(dyclee$aList)) {
        trayectorias_informacion[dyclee$aList[[uCk]]$id_trays,'clusterFinal'] <- dyclee$aList[[uCk]]$label 
      }
    }
    if (length(dyclee$oList)>0) {
      for (uCk in 1:length(dyclee$oList)) {
        trayectorias_informacion[dyclee$oList[[uCk]]$id_trays,'clusterFinal'] <- -1#dyclee$oList[[uCk]]$label
      }
    }
    return(trayectorias_informacion)
  }
  
  update_uClusters_to_binnacle <- function(dyclee, trayectorias_informacion){
    acc <- 0
    
    if (length(dyclee$aList)>0) {
      for (uCk in 1:length(dyclee$aList)) {
        atemp <- dyclee$aList[[uCk]]$id_trays
        trayectorias_informacion[atemp,'microCluster'] <- uCk
      }
      acc <- uCk
    }

    if (length(dyclee$oList)>0) {
      continued <- acc + 1#+ length(dyclee$oList)
      
      for (uCk in 1:length(dyclee$oList)) {
        otemp <- dyclee$oList[[uCk]]$id_trays
        trayectorias_informacion[otemp,'microCluster'] <- continued#-1#dyclee$oList[[uCk]]$label
        continued <- continued + 1
      }
    }

    return(trayectorias_informacion$microCluster)
  }
  
  #Creo copias de los Clusters Finales segun el "tGlobal" para guardar los resultados de la evolucion 
  #'de datos
  #' @param CopyDyclee: list()
  #' @param CopyBinnacle: list()
  #'
  copy_according_to_time_series <-  function(time_series, dyclee, data_proc){
    CopyDyclee[[time_series]] <<-  dyclee
    CopyBinnacle[[time_series]] <<- data_proc
  }
  
  #Graficar los resultados en un plano simple
  simple_plot <- function(dyclee){
    colors <- rainbow(n = length(dyclee$aList))
    par(lwd=2)
    plot(NULL ,ylim=c(min(all_celds$vel_prom), max(all_celds$vel_prom)),
         xlim= c(min(all_celds$vel_prom), max(all_celds$vel_prom)),
         xlab = "Velocidad Promedio", ylab = "Velocidad Promedio", cex.lab=1.5, cex.sub=1.5, main = dataset$dataset)
    if(length(dyclee$oList) > 0){
      for(k in 1:length(dyclee$oList)){
        ol <- dyclee$oList[[k]]$data
        points(x = ol[,"vel_prom"], y = ol[,"vel_prom"], pch=21, col="black", cex=0.5)
      }
    }
    if(length(dyclee$aList) > 0){
      for(k in 1:length(dyclee$aList)){
        al <- dyclee$aList[[k]]$data
        points(x = al[,"vel_prom"], y = al[,"vel_prom"], pch=21, col=colors[k], cex=0.5)
      }
    }
  }

  # Generar bitacora para etapa 1
  generate_binnacle_E1 <- function(dyclee, general_data_celds){
    temp_bitacora <- general_data_celds
    if(nrow(temp_bitacora)>0){
      temp_bitacora$etapa1 <- 0
      if(length(dyclee$aList) > 0){
        for(ak in 1:length(dyclee$aList)){
          label <- dyclee$aList[[ak]]$label
          ids <- dyclee$aList[[ak]]$data$proc_id
          temp_bitacora[temp_bitacora$proc_id %in% ids,"etapa1"] <- label
        }
      }
      if(length(dyclee$oList) > 0){
        for(ok in 1:length(dyclee$oList)){
          ids <- dyclee$oList[[ok]]$data$proc_id
          temp_bitacora[temp_bitacora$proc_id %in% ids,"etapa1"] <- 0
        }
      }
    }
    return(temp_bitacora)
  }
  
  # Generar bitacora para etapa 2
  generate_binnacle_E2 <- function(dyclee, general_data_celds){
    temp_bitacora <- general_data_celds
    if(nrow(temp_bitacora)>0){
      temp_bitacora$etapa2 <- 0
      if(length(dyclee$aList) > 0){
        for(ak in 1:length(dyclee$aList)){
          label <- dyclee$aList[[ak]]$label
          ids <- dyclee$aList[[ak]]$data$proc_id
          temp_bitacora[temp_bitacora$proc_id %in% ids,"etapa2"] <- label
        }
      }
      if(length(dyclee$oList) > 0){
        for(ok in 1:length(dyclee$oList)){
          ids <- dyclee$oList[[ok]]$data$proc_id
          temp_bitacora[temp_bitacora$proc_id %in% ids,"etapa2"] <- 0
        }
      }
    }
    return(temp_bitacora)
  }

  #==================================================================================================
  
  # Actualizar Log de cada grupo
  update_log <- function(dyclee, evolucion, etapa=0){
    if(length(dyclee$aList)>0){
      for(i in 1:length(dyclee$aList)){
        mc <- dyclee$aList[[i]]
        #LS.longitude       LS.latitude
        #SS.X1            SS.X2
        df_log <- data.frame(idMC=mc$idMC,
                             evolucion = evolucion,
                             etapa = etapa,
                             grupo_anterior = mc$prev_label,
                             grupo_actual = mc$label,
                             elementos_anteriores = mc$prev_elem,
                             elementos_actuales = nrow(mc$tramo),
                             elementos_proc = abs(nrow(mc$tramo) - mc$prev_elem),
                             velocidad=mean(mc$tramo[,"velocidad"]),
                             LS = mc$CF$LS,
                             #SS = mc$CF$SS,
                             n = mc$CF$n,
                             V = prod(mc$hyperboxSizePerFeature),
                             Densidad = getD(mc),
                             Actual_dens_mean = dyclee$densityMean,
                             Actual_dens_median = dyclee$densityMedian
        )
        mc$log <- rbind(mc$log, df_log)
        mc$prev_label <- mc$label
        mc$prev_elem <- nrow(mc$tramo)
        dyclee$aList[[i]] <- mc
      }
    }
    if(length(dyclee$oList)>0){
      for(i in 1:length(dyclee$oList)){
        mc <- dyclee$oList[[i]]
        #LS.longitude       LS.latitude
        #SS.X1            SS.X2
        df_log <- data.frame(idMC=mc$idMC,
                             evolucion = evolucion,
                             etapa = etapa,
                             grupo_anterior = mc$prev_label,
                             grupo_actual = -1*i,#mc$label,
                             elementos_anteriores = mc$prev_elem,
                             elementos_actuales = nrow(mc$tramo),
                             elementos_proc = abs(nrow(mc$tramo) - mc$prev_elem),
                             velocidad=mean(mc$tramo[,"velocidad"]),
                             LS = mc$CF$LS,
                             #SS = mc$CF$SS,
                             n = mc$CF$n,
                             V = prod(mc$hyperboxSizePerFeature),
                             Densidad = getD(mc),
                             Actual_dens_mean = dyclee$densityMean,
                             Actual_dens_median = dyclee$densityMedian
        )
        mc$log <- rbind(mc$log, df_log)
        mc$prev_label <- -1*i#mc$label
        mc$prev_elem <- nrow(mc$tramo)
        dyclee$oList[[i]] <- mc
      }
    }
    
    return(dyclee)
  }
  

#########################################################################################
############   INICIO DEL ALGORITMO   ###################################################
#########################################################################################


  #Librerias y Data-Source####
  
  #Librerias
  packages <- c('RPostgreSQL','lubridate', 'schoolmath', 'compare', 'pracma', 'purrr', 'RColorBrewer', 'useful', 'rjson', 'animation','dplyr','purrr','magick',"chron","leaflet")
  ipak(packages) 
  
  #Data Source y Creddenciales del PostgreSQL
  dbdriver <- "PostgreSQL";
  host <- 'localhost';
  port='5432';
  dbname='CPPP2'
  user='postgres';
  pass='postg369';

  #Datos / Exploracion####
  #*Guayaquil = 1  :=  30557
  #* Roma     = 5       36650   entre las 14h00 y 15h00  Zona horaria de Roma, Italia (GMT+2) 
  #* Beijing  = 6
  
  # Establecer el numero del dataset a importar
  location <- 1 
  dataset <- get_DatasetsFromDB(location)
  

  # Eliminacion de datos repetidos y ordenacion por fecha/hora
  print(dataset$size)
  dataset$data <- dataset$data[!duplicated(dataset$data),]
  dataset$data <- dataset$data[order(dataset$data$unixtime, decreasing = FALSE),]
  nrow(dataset$data)
  length(unique(dataset$data$file_name))
  

  # Gereacion del Workspace
  folder <- "Pruebas-R"
  algorithm <- "Dyclee"
  
  #lista de los identificadores de trayectorias
  names_tray <- unique(dataset$data$file_name) 
  boolean_segmentacion <- T #TRUE o FLASE para aplicar segmentacion
  
  ##
  ## PARAMETRIZACIONES
  ##
  
  #' @param relativeSize: Integer = Este parÃ¡metro especifica el tamaÃ±o de las hipercajas que dan forma a los grupos ??C,
  #' tamaÃ±o-relativo es un nÃºmero real en el intervalo [0, 1].
  #' @param speed:  Integer 
  #' @param lambd: 
  #'  No usadas (legacy):
  #' @param periodicUpdateAt:
  #' @param timeWindow:
  #' @param periodicRemovalAt:
  #' @param closenessThreshold:
  #' @param uncommonDimensions:
  #' @param tGlobal: Integer Serie de tiempo o registro para la evolucion de  datos
  #' @param ac: Integer Acumulador  que se evalua con tGlobal  para proceder a calcular lod grupos densos
  #' @param ac_Evolution: Acumulador que guarda el numero de cuantas evoluciones seran
  #' @param CopyDyclee: Lista que contiene las caracteristicas del algortimo por cada evolucion
  #' @param CopyBinnacle: Lista que contiene  la bitacora de los datos conforme a las evoluciones
  
  #===== COMUNMENTE USADAS ===================================================
  relativeSize = 0.1#0.5
  lambd = 0 # 0.75 # if it has a value over 0, when a micro cluster is updated, tl will be checked and the diff with current time will matter
  #Minimo y maximo del rango admitido en el AGRUPAMIENTO:
  min_v <- 0     
  max_v <- 200
  #Minimo y maximo del rango admitido en los MAPAS:
  vv_min <- 0
  vv_max <- 100
  #====== NO USADAS ==========================================================
  speed = 1 #not used
  periodicUpdateAt = 3# float("inf")
  timeWindow = 3# 5# float("inf")
  periodicRemovalAt = 10#3 # float("inf")
  tGlobal <-500
  closenessThreshold = 1.5#1.5
  uncommonDimensions = 0 #.001
  #===========================================================================

  ac <- 0
  ac_Evolution <- 0
  CopyDyclee <- list()
  CopyBinnacle <- list()
  idMC <- 0
  


  #Creacion del directorio para guardar los resultados
  get_workSpace(carpeta = folder, algoritmo = algorithm, dataset = dataset$dataset, angulo_tolerancia = angulo_tolerancia,relativeSize=relativeSize, tGlobal)#Crea el directorio donde se alojaran los resultadso del algortimo
  folder_result <- getwd() 
  save.image(file="RDataSet.RData")

  
  # INICIO ALGORITMO
  print("INICIO ALGORITMO")
  
  # Area Original de procesamiento 
  odataContext <- dataContext <-  list(c(maximun=max(dataset$data$longitude),
                         minimun=min(dataset$data$longitude)),
                       c(maximun=max(dataset$data$latitude),
                         minimun=min(dataset$data$latitude))
  )
  
  dyclee <- Dyclee(dataContext=dataContext, 
                   relativeSize=relativeSize,
                   closenessThreshold=closenessThreshold,
                   uncommonDimensions = uncommonDimensions,
                   speed = speed,
                   lambd=lambd)
  
  

  #===========================================================================
  # Adaptación a lotes por tiempos (MICRO-BATCH)
  # Colocar tiempo en minutos:   10 := 10 minutos; 0.3 := 0.3 minutos 
  # Duracion del periodo de tiempo por cada lote
  #===========================================================================
  # Representa Flujos de cada X minutos:
  tiempo_a_ejecutar_algoritmo <- 3 
  #===========================================================================
  mins_x_dia<- 1440
  
  #Calculo de la brecha de tiempo de la Dataset en minutos
  datos <- dataset[["data"]]
  time_init <- datos[order(datos$unixtime, decreasing=FALSE),][1,]$unixtime 
  time_end  <- datos[order(datos$unixtime, decreasing=TRUE),][1,]$unixtime 
  time_diff_dataset_mins <- as.numeric(time_end - time_init, units = "mins")
  
  
  #Partametros para ejecutar el algortimos por filtro de tiempo
  time_increase1 <- times(tiempo_a_ejecutar_algoritmo/mins_x_dia)
  
  #para este algoritmo em particular se le asigna 1 mas por que el priemro siempre sera el historico mas grande  y los demas dependeras del restante de tiempo
  cant_execut_dataset <- ceiling(time_diff_dataset_mins / tiempo_a_ejecutar_algoritmo) 
  #cant_execut_dataset <- 5
  lista_trayectorias <- get_list_per_trajectory(dataset = dataset)
  # Calcular velocidades en el caso de que no existieran
  if(is.null(dataset[["data"]]$velocidad)){#V no hay velocidades
    for(li in 1:length(lista_trayectorias)){
      trayect <- lista_trayectorias[[li]]
      #print(nrow(trayect))
      if(nrow(trayect)>1){
        velocidades_celda <- c()
        # V * T = D      V = D / T
        delta_time <- 0
        delta_dist <- 0
        velocidades <- c()
        vel <- 0
        for(i_seg in 2:nrow(trayect)){
          delta_time <- difftime(trayect[i_seg,"unixtime"],trayect[i_seg-1,"unixtime"],units = "secs")
          delta_dist <- abs(dist(rbind(trayect[i_seg-1,c('longitude','latitude')], trayect[i_seg,c('longitude','latitude')]), method = "manhattan"))
          if(delta_time > 0){
            seg_to_hour <- as.numeric(delta_time)/(60*60)
            vel = delta_dist / as.numeric(seg_to_hour) * 100
          }else{
            vel = 0
          }
          velocidades <- c(velocidades, round(vel,2))
        }
        velocidades <- c(velocidades[1],velocidades)
      }else if(nrow(trayect)==1){
        velocidades <- c(0)
      }
      lista_trayectorias[[li]]$velocidad <- velocidades
      lista_trayectorias[[li]] <- lista_trayectorias[[li]][lista_trayectorias[[li]]$velocidad < 200,]
    }
  }
  lista_trayectorias <- do.call(rbind.data.frame, lista_trayectorias)
  print(as.numeric(object.size(lista_trayectorias)))
  lista_trayectorias <- lista_trayectorias[order(lista_trayectorias$unixtime, decreasing = FALSE),]
  lista_trayectorias <- get_list_trajectory(lista = lista_trayectorias)
  
  # Guardado de los parametros del experimento
  par <- t(data.frame(relativeSize=relativeSize,speed,lambd,periodicUpdateAt,timeWindow,periodicRemovalAt,closenessThreshold,uncommonDimensions,tiempo_a_ejecutar_algoritmo,cant_execut_dataset))
  write.table(x = par, file = "Parametros.txt", sep = ";", dec = ",", col.names = T, row.names = T)
  
  save.image(file="RData_Informacion.RData")
  
  #===========================================================================
  #Reticulado - Generar Reticulado####
  #===========================================================================
  # Definir N areas del reticulado en base a tama?s de las celdas
  rango <- 0.002  #Tama? de la celda:   rango x rango
  #===========================================================================

  areaA <- dataContext[[1]][["maximun"]]-dataContext[[1]][["minimun"]]
  ngruposlon <- ceiling(areaA / rango)
  
  areaB <- dataContext[[2]][["maximun"]]-dataContext[[2]][["minimun"]]
  ngruposlat <- ceiling(areaB / rango)
  
  ngruposarea <- ngruposlon * ngruposlat
  
  datos_general <- list() #Almacena los datos de todas las evoluciones
  datos_general_filtrado <- list() #Almacena los datos filtrados(sin valore nulos)
  datos_general_puntos <- list()#Almacena los puntos que pertenecen a cada celda
  #timeinicio <- datos[order(datos$unixtime, decreasing=FALSE),][1,]$unixtime 
  timeinicio <- time_init
  timefin <- time_end
  
  rangos_intercuartiles <- data.frame()
  
  dir.create(paste("Reticulado",sep=''))
  dir_ret <- getwd()
  setwd(paste(dir_ret,'//',"Reticulado",sep=''))
  
  #cant_execut_dataset <- 10
  time_ccel <- list()
  time_cel <- list()
  size_buffer_puntos <- list()
  size_buffer_celdas <- list()
  size_puntos <- as.numeric(object.size(do.call(rbind.data.frame,lista_trayectorias)))
  size_celdas <- 0
  
  ini_cel<- as.POSIXct(Sys.time())
  for(evolucion in 1:cant_execut_dataset){
    #iteraciones de las evoluciones
    print(paste("Evolucion ", evolucion,sep=""))
    id_celda <- 0
    cant=0
    tr_x_flujo <- list()
    datos_evolucion <- data.frame()
    info_puntos_celda <- list()
    print("+-Puntos correspondientes")
    ini_ccel<- as.POSIXct(Sys.time())
    for(point in 1:length(lista_trayectorias)){
      if(lista_trayectorias[[point]]$unixtime>= timeinicio 
         &  lista_trayectorias[[point]]$unixtime<timeinicio + duration(tiempo_a_ejecutar_algoritmo, "minutes")){
        cant = cant + 1
        tr_x_flujo[[cant]] <- lista_trayectorias[[point]]
      }
    }
    end_ccel <- as.POSIXct(Sys.time())
    print("+-Informacion de celdas")
    #Recorrer cada celda, orden: menor a mayor longitud y menor a mayor latitud
    ini_cel<- as.POSIXct(Sys.time())
    for(grupo_lat in 1:ngruposlat){
      for(grupo_long in 1:ngruposlon){
        #print(paste(grupo_long," - ",grupo_lat,sep=""))
        id_celda <- id_celda + 1
        #Datos de la celda
        datos_celda <- list(id=id_celda, grupo_long=0, grupo_lat=0,
                            lon_minima=0, lon_maxima=0,
                            lat_minima=0, lat_maxima=0,
                            cant_segmentos=0, cant_veh=0,
                            densidad=0,
                            vel_min=0,vel_max=0,
                            vel_prom=0,angulo=0,
                            unixtime=0)#Almacena los datos de la evolucion actual
        info_celda <- data.frame()
        #Asignacion de datosde la celda
        datos_celda$grupo_long <- grupo_long
        datos_celda$grupo_lat  <- grupo_lat
        datos_celda$lon_minima <- dataContext[[1]][["minimun"]] + (rango * (grupo_long - 1))
        datos_celda$lon_maxima <- dataContext[[1]][["minimun"]] + (rango * (grupo_long))
        datos_celda$lat_minima <- dataContext[[2]][["minimun"]] + (rango * (grupo_lat - 1))
        datos_celda$lat_maxima <- dataContext[[2]][["minimun"]] + (rango * (grupo_lat))
        
        puntos_celda <- 0
        velocidades_celda <- c()
        punto_inicio <- 0
        punto_fin <- 0
        flag <- TRUE
        #identificacion de segmentos que correspondan a la celda
        if(length(tr_x_flujo) > 0){
          for(tam_lista in 1:length(tr_x_flujo)){
            if(tr_x_flujo[[tam_lista]]$longitude  >= datos_celda$lon_minima
               & tr_x_flujo[[tam_lista]]$longitude < datos_celda$lon_maxima
               & tr_x_flujo[[tam_lista]]$latitude >= datos_celda$lat_minima
               & tr_x_flujo[[tam_lista]]$latitude  < datos_celda$lat_maxima
               & tr_x_flujo[[tam_lista]]$velocidad > 0){
              puntos_celda <- puntos_celda + 1
              info_celda <- rbind(info_celda,tr_x_flujo[[tam_lista]])
              if(tr_x_flujo[[tam_lista]]$velocidad <= max_v){
                velocidades_celda <- c(velocidades_celda, tr_x_flujo[[tam_lista]]$velocidad)
              }
              if(flag){
                punto_inicio <- tr_x_flujo[[tam_lista]]
                punto_fin <- tr_x_flujo[[tam_lista]]
                flag <- FALSE
              }else{
                punto_fin <- tr_x_flujo[[tam_lista]]
              }
            }
          }
        }
        
        #Calculo Velocidad y cantidad de elementos
        if(puntos_celda <= 0){
          datos_celda$vel_prom <- 0
          datos_celda$cant_segmentos <- 0
          datos_celda$vel_min <- 0
          datos_celda$vel_max <- 0
          datos_celda$angulo <- 0
        }else{
          datos_celda$cant_veh <- length(unique(info_celda[,"file_name"]))
          datos_celda$cant_segmentos <- puntos_celda
          datos_celda$vel_min <- min(velocidades_celda)
          datos_celda$vel_max <- max(velocidades_celda)
          datos_celda$angulo <- round(get_calculate_angle(punto_inicio[,c("longitude","latitude")],punto_fin[,c("longitude","latitude")]),2) #point[,c("longitude","latitude")]
          a <- as.character(as.POSIXct(mean(info_celda[,"unixtime"]), origin="1970-01-01"))
          #print(a)
          datos_celda$unixtime <- a
          vel_prom <- round(mean(velocidades_celda),2)
          if(is.infinite(vel_prom) || is.na(vel_prom)){
            vel_prom <- 0
          }
          datos_celda$vel_prom <- vel_prom
        }
        
        #Calculo Densidad de la celda
        #V <- prod(abs(datos_celda$lon_maxima - datos_celda$lon_minima),
        #          abs(datos_celda$lat_maxima - datos_celda$lat_minima))
        V <- prod(dyclee$hyperboxSizePerFeature)
        dens_celda <- datos_celda$cant_segmentos / V
        datos_celda$densidad = round(dens_celda,3)
          
        #TODO: Quitar
        #cat <- cut(vel_prom, breaks = c(0,10,20,30,40,50,60,70,80,90,100,Inf),
        #           labels = c(1,2,3,4,5,6,7,8,9,10,-1))
        
        #Unir la celda a datos_evolucion
        datos_evolucion <- rbind(datos_evolucion,as.data.frame(datos_celda))
        #Unir la informacion de los puntos de la celda a puntos_celda
        info_puntos_celda[[id_celda]] <- info_celda
      }#for grupo long
    }#for grupo lat
    end_cel <- as.POSIXct(Sys.time())
    
    time_ccel[[evolucion]] <- c(ini_ccel,end_ccel)
    time_cel[[evolucion]] <- c(ini_cel,end_cel)
    
    #Guardado de la info de las celdas en 
    datos_general_puntos[[evolucion]] <- info_puntos_celda
    
    #Unir datos_evolucion a datos_general
    datos_general[[evolucion]] <- datos_evolucion
    datos_general_filtrado[[evolucion]] <- subset(datos_evolucion, subset = datos_evolucion$vel_prom > 0 & datos_evolucion$cant_segmentos > 0) # 
    
    #Guardar dataframe
    #write.table(x = datos_evolucion, file = paste("Datos Reticula Ev",evolucion,".csv"), sep = ";", dec = ",", col.names = T, row.names = F)
    #write.table(x = datos_general_filtrado[[evolucion]], file = paste("Datos Reticula toProc Ev",evolucion,".csv"), sep = ";", dec = ",", col.names = T, row.names = F)
    
    if(FALSE){
      #Reticulado
      #dev.list()
      print("+-Guardado reticula Velocidad")
      #Graficar la evolucion con el conjunto de celdas
      nombre <- paste("Evolucion ", evolucion," Velocidad.png", sep="")
      jpeg(nombre, width = 1080, height = 1080, units = 'px', pointsize = 15 ,res = 100)
      titulo <- paste(format(timeinicio, "%H:%M:%S")," a ",
                      format(timeinicio + duration(tiempo_a_ejecutar_algoritmo, "minutes"),"%H:%M:%S"), 
                      " -> ", ngruposarea, " Celdas",
                      " - Velocidades",sep="")
      # reticula <- ggplot(datos_general_filtrado[[evolucion]],aes(lon_minima, lat_minima, fill=vel_prom)) +
      #  geom_tile(aes(fill =vel_prom)) +
      #  geom_text(aes(label = cant_segmentos),color="white",size=2.5) +
      #  #scale_fill_gradient(low="white",high = "red") +
      #  scale_fill_continuous(type = "viridis",name = "Velocidades") +
      #  #scale_fill_viridis_b()+
      #  labs(title = titulo, x = "Longitude", y = "Latitude" ) +
      #  theme(plot.title = element_text(size = rel(1.5)))
      reticula <- ggplot(datos_general_filtrado[[evolucion]],aes(x = lon_minima, y = lat_minima)) +
        lims(x = c(dataContext[[1]][["minimun"]]-(rango/2), dataContext[[1]][["maximun"]]),
             y = c(dataContext[[2]][["minimun"]]-(rango/2), dataContext[[2]][["maximun"]])) +
        #xlim(dataContext[[1]][["minimun"]]-(rango/2),dataContext[[1]][["maximun"]]+(rango/2)) +
        #ylim(dataContext[[2]][["minimun"]]-(rango/2),dataContext[[2]][["maximun"]]+(rango/2)) +
        #xlab("Longitude") + ylab("Latitude") + ggtitle(titulo) +
        labs(title = titulo, x = "Longitude", y = "Latitude") +
        geom_tile(aes(fill = vel_prom)) +
        scale_fill_continuous(type = "viridis",name = "Velocidad",limits=c(0,100)) +
        #scale_fill_viridis_b() +
        geom_text(aes(label = cant_segmentos),color="white",size=2.5) +
        theme(plot.title = element_text(size = rel(1.5)))
      print(reticula)
      dev.off()
      
      print("+-Guardado reticula Cantidad de Puntos")
      nombre <- paste("Evolucion ", evolucion," Puntos.png", sep="")
      jpeg(nombre, width = 1080, height = 1080, units = 'px', pointsize = 15 ,res = 100)
      
      titulo2 <- paste(format(timeinicio, "%H:%M:%S")," a ",
                       format(timeinicio + duration(tiempo_a_ejecutar_algoritmo, "minutes"),"%H:%M:%S"), 
                       " -> ", ngruposarea, " Celdas",
                       " - Cantidad de Puntos",sep="")
      
      #reticula2 <- ggplot(datos_general_filtrado[[evolucion]],aes(lon_minima, lat_minima)) +
      #  geom_tile(aes(fill =cant_segmentos)) +
      #  scale_fill_continuous(type = "viridis",name = "Puntos") +
      #  geom_text(aes(label = cant_segmentos),color="white",size=2.5) +
      #  labs(title = titulo2, x = "Longitude", y = "Latitude" ) +
      #  theme(plot.title = element_text(size = rel(1.5)))
      reticula2 <- ggplot(datos_general_filtrado[[evolucion]],aes(x = lon_minima, y = lat_minima)) +
        lims(x = c(dataContext[[1]][["minimun"]]-(rango/2), dataContext[[1]][["maximun"]]),
             y = c(dataContext[[2]][["minimun"]]-(rango/2), dataContext[[2]][["maximun"]])) +
        #xlim(dataContext[[1]][["minimun"]]-(rango/2),dataContext[[1]][["maximun"]]+(rango/2)) +
        #ylim(dataContext[[2]][["minimun"]]-(rango/2),dataContext[[2]][["maximun"]]+(rango/2)) +
        #xlab("Longitude") + ylab("Latitude") + ggtitle(titulo) +
        labs(title = titulo2, x = "Longitude", y = "Latitude") +
        geom_tile(aes(fill = cant_segmentos)) +
        scale_fill_continuous(type = "viridis",name = "Puntos") +
        #scale_fill_viridis_b() +
        geom_text(aes(label = cant_segmentos),color="white",size=2.5) +
        theme(plot.title = element_text(size = rel(1.5)))
      print(reticula2)
      dev.off()
      
    }
    
    print("+-Rangos Intercuartiles")
    #Rangos intercuartiles
    r_interc <- quantile(datos_general_filtrado[[evolucion]]$vel_prom)
    rangos_intercuartiles <- rbind(rangos_intercuartiles,cbind(evolucion,t(r_interc)))
    
    #Actualizr tiempos
    timeinicio <- timeinicio + duration(tiempo_a_ejecutar_algoritmo, "minutes")
  }
  end_cel <- as.POSIXct(Sys.time())
  
  write.table(x = rangos_intercuartiles, file = paste("Rangos_Intercuartiles.csv"), sep = ";", dec = ",", col.names = T, row.names = F)
  
  size_celdas <- as.numeric(object.size(do.call(rbind.data.frame,datos_general_filtrado)))
  
  #DyClee####
  save.image(file="Reticula.RData")
  
  print(" DYCLEE ")
  setwd(folder_result)
  dir.create(paste("Procesamiento",sep=''))
  dir_ret <- getwd()
  setwd(paste(dir_ret,'//',"Procesamiento",sep=''))
  all_celds <- do.call(rbind.data.frame,datos_general_filtrado)
  all_celds <- cbind(proc_id=c(1:nrow(all_celds)),all_celds)
  general_data_celds <- data.frame()
  
  dataContext <-  list(c(maximun=max(max_v),
                         minimun=min(min_v)),
                       c(maximun=max(max_v),
                         minimun=min(min_v)) )
  
  dyclee$dataContext <- dataContext
  dyclee$hyperboxSizePerFeature <- getHyperboxSizePerFeature()
  save.image(file="Bck-PreProc.RData")
  
  #===============================================================================
  
  time_e1 <- list()
  time_e2 <- list()
  time_proc <- list()
  DCE1 <- list()
  DCE2 <- list()
  save.image(file="DC.EX-00.RData")
  
  #===============================================================================
  
  
  
  
  
  
  
  
  
  
#===============================================================================
# Agrupamiento
#===============================================================================
  ipak(c("xlsx"))
  folder_result <- getwd() 
  save.image(file="EX-00.RData")
  
  #Establecer el peridoo a ejecutar
  PER <- 8
  #1-5/6-10/11-15/16-20/21-25/26-30/31-35/36-40
  ev_inicio <- (5 * (PER - 1)) + 1
  ev_fin <- 5 * (PER)
  
  #===============================================================================
  
  
  for(pos_Tray in ev_inicio:ev_fin){
    print(paste("== Flujo # ", pos_Tray," ==", sep=""))
    ini_iter<- as.POSIXct(Sys.time()) 

    print("Etapa: 1 - Inicio")
    #1er Etapa - Creacion de los MicroCLusters (UClusters) - 
    # Actualizo la bitacora de forma Global  por medio de pos_Tray
    #Primera etapa: Algoritmo basaado en distancia
    celdas_con_datos <- datos_general_filtrado[[pos_Tray]]
    print(nrow(celdas_con_datos))
    ini_e1<- as.POSIXct(Sys.time())
    if(nrow(celdas_con_datos) > 0){
      proc_inicio <- nrow(general_data_celds)+1
      proc_fin <- nrow(general_data_celds)+nrow(celdas_con_datos)
      celdas_con_datos <- cbind(proc_id=c(proc_inicio:proc_fin),celdas_con_datos)
      general_data_celds <- rbind(general_data_celds, celdas_con_datos)
      
      for(li in 1:nrow(celdas_con_datos)){
        uniq_celd <- celdas_con_datos[li,]
        dyclee <- trainOnElement(dyclee, info_c = uniq_celd) #lista_trayectorias[[punto$idor]],punto$idor,punto$celda)
      }
      
    }
    end_e1 <- as.POSIXct(Sys.time())
    
    
    #contador de la evolucion de los datos
    ac_Evolution <- ac_Evolution + 1
    print("Etapa: 1 - Fin")

    #Actualizacion y guardado de metricas de rendimiento
    dyclee <- update_log(dyclee, evolucion = ac_Evolution, etapa = 1)
    DCE1[[ac_Evolution]] <- dyclee
    B_E1 <- generate_binnacle_E1(dyclee, general_data_celds)
    #write.table(x = B_E1, file = paste("MicroClusters ",pos_Tray,".csv"), sep = ";", dec = ",", col.names = T, row.names = F)
    
    
    print("Etapa: 2 - Inicio")
    #' 2da Etapa - Creacion de las densidades (Aclusters) a partir de los Uclusters
    #   al finalizar esta estapa de densidad aClust contrendra los grupos(uC) finales y oClust los uC atipicos
    ini_e2<- as.POSIXct(Sys.time())
    if((length(dyclee$aList)+length(dyclee$oList))>0){
      dyclee <- getClusteringResult(dyclee)
    }
    end_e2 <- as.POSIXct(Sys.time())
    time_e1[[ac_Evolution]] <- c(ini_e1,end_e1)
    time_e2[[ac_Evolution]] <- c(ini_e2,end_e2)
    print(length(dyclee$aList))
    #actualizar la bitacora con los grupos densos formados
    #trayectorias_informacion <- update_uClustersFinals_to_binnacle(dyclee = dyclee, trayectorias_informacion = trayectorias_informacion)
    print("Etapa: 2 - Fin")
    #Actualizacion y guardado de metricas de rendimiento
    dyclee <- update_log(dyclee, evolucion = ac_Evolution, etapa = 2)
    DCE2[[ac_Evolution]] <- dyclee
    B_E2 <- generate_binnacle_E2(dyclee, B_E1)
    #write.table(x = B_E2, file = paste("MicroClustersFinals ",pos_Tray,".csv"), sep = ";", dec = ",", col.names = T, row.names = F)
    
    #guardo cada evolucion  en una super lista  tanto uC  como la bitacora
    copy_according_to_time_series(time_series = pos_Tray,dyclee = dyclee,data_proc = B_E2)
    
    end_iter <- as.POSIXct(Sys.time())
    time_proc[[ac_Evolution]] <- c(ini_iter,end_iter)
    
    n_ev <- paste("RData_Ev-",pos_Tray,".RData",sep="")
    #save.image(file=n_ev)
  }
  
  save.image(file="RDataResult.RData")
  setwd(folder_result)
  
  #save_information( CopyBinnacle = CopyBinnacle)

##
## METRICAS E INDICADORES EN GENERAL
##

  dir.create(paste("ResumenIndicadores",sep=''))
  dir_ret <- getwd()
  setwd(paste(dir_ret,'//',"ResumenIndicadores",sep=''))
  
  res_evs <- data.frame()
  rango <- dyclee$hyperboxSizePerFeature[["maximun"]]
  rango_2 <- rango / 2
  for(k in ev_inicio:ev_fin){
    #Informacion de Indicadores por evolucion
    ult_evol <- CopyDyclee[[k]]
    df_indicadores <- data.frame()
    if(length(ult_evol$oList)>0){
      atipico <- data.frame()
      for (at in 1:length(ult_evol$oList)){
        atipico <- rbind(atipico, ult_evol$oList[[at]]$data)
      }
      a_df <- data.frame(id_grupo=0,
                         total_celdas=length(atipico$id),
                         total_puntos=sum(atipico$cant_segmentos),
                         min_puntos=min(atipico$cant_segmentos),
                         max_puntos=max(atipico$cant_segmentos),
                         prom_puntos=round(mean(atipico$cant_segmentos),2),
                         desv_puntos=sd(atipico$cant_segmentos),
                         cant_veh=sum(atipico$cant_veh),
                         vel_min=min(atipico$vel_prom),
                         vel_max=max(atipico$vel_prom),
                         vel_prom=round(mean(atipico$vel_prom),2),
                         desv_vel=round(sd(atipico$vel_prom),2),
                         vel_rango_min = round(mean(atipico$vel_prom),2) - rango_2,
                         vel_rango = rango,
                         vel_rango_max = round(mean(atipico$vel_prom),2) + rango_2
                         )
      df_indicadores <- rbind(df_indicadores, a_df)
    }
    
    if(length(ult_evol$aList)>0){
      for(i in 1:length(ult_evol$aList)){
        g <- ult_evol$aList[[i]]
        g_df <- data.frame(id_grupo=i,
                           total_celdas=length(g$id_celds),
                           total_puntos=sum(g$data$cant_segmentos),
                           min_puntos=min(g$data$cant_segmentos),
                           max_puntos=max(g$data$cant_segmentos),
                           prom_puntos=round(mean(g$data$cant_segmentos),2),
                           desv_puntos=sd(g$data$cant_segmentos),
                           cant_veh=sum(g$data$cant_veh),
                           vel_min=min(g$data$vel_prom),
                           vel_max=max(g$data$vel_prom),
                           vel_prom=round(mean(g$data$vel_prom),2),
                           desv_vel=round(sd(g$data$vel_prom),2),
                           vel_rango_min = round(mean(g$data$vel_prom),2) - rango_2,
                           vel_rango = rango,
                           vel_rango_max = round(mean(g$data$vel_prom),2) + rango_2)
        df_indicadores <- rbind(df_indicadores, g_df)
      }
    }
    
    if(nrow(df_indicadores)>0){
      r_ev <- data.frame(evolucion=k,
                         cant_celdas=sum(df_indicadores$total_celdas),
                         cant_puntos=sum(df_indicadores$total_puntos),
                         total_puntos=sum(df_indicadores$total_puntos),
                         min_puntos=min(df_indicadores$total_puntos),
                         max_puntos=max(df_indicadores$total_puntos),
                         prom_puntos=round(mean(df_indicadores$total_puntos),2),
                         desv_puntos=sd(df_indicadores$total_puntos),
                         cant_veh=sum(df_indicadores$cant_veh),
                         vel_min=min(df_indicadores$vel_prom),
                         vel_max=max(df_indicadores$vel_prom),
                         vel_prom=round(mean(df_indicadores$vel_prom),2),
                         desv_vel=round(sd(df_indicadores$vel_prom),2),
                         vel_rango_min = round(mean(df_indicadores$vel_prom),2) - rango_2,
                         vel_rango = rango,
                         vel_rango_max = round(mean(df_indicadores$vel_prom),2) + rango_2
      )
      res_evs <- rbind(res_evs, r_ev)
      write.table(x = df_indicadores, file = paste("Resumen_Ev-",k,".csv"), sep = ";", dec = ",", col.names = T, row.names = F)
      
    }
  }
  write.table(x = res_evs, file = paste("Resumen_General.csv"), sep = ";", dec = ",", col.names = T, row.names = F)
  save.image(file="Indicadores.RData")
  
  setwd(folder_result)
  
  print("FIN ALGORITMO")
  
  setwd(folder_result)
  save.image(file="RData.RData")
  #save_information( CopyBinnacle = CopyBinnacle)
  
  print("FIN ALGORITMO")
  
  #Consolidar y eportar los resultados (celdas) de las agrupaciones
  
  ult_evol <- DCE2[[length(DCE2)]]
  df_indicadores <- data.frame()
  
  if(length(ult_evol$oList)>0){
    for(i in 1:length(ult_evol$oList)){
      mc <- dyclee$oList[[i]]$tramo
      mc_df <- data.frame(id_grupo=(-1*i),
                          celdas=length(unique(dyclee$oList[[i]]$data[,"id"])),
                          total_puntos=nrow(mc),
                          cant_veh=length(unique(mc$file_name)),
                          vel_min=min(mc$velocidad),
                          vel_max=max(mc$velocidad),
                          vel_prom=round(mean(mc$velocidad),2),
                          desv_vel=round(sd(mc$velocidad),2)
      )
      df_indicadores <- rbind(df_indicadores, mc_df)
      rm(mc)
      rm(mc_df)
      gc()
    }
  }
  
  if(length(ult_evol$aList)>0){
    for(i in 1:length(ult_evol$aList)){
      mc <- dyclee$aList[[i]]$tramo
      mc_df <- data.frame(id_grupo=i,
                          celdas=length(unique(dyclee$aList[[i]]$data[,"id"])),
                          total_puntos=nrow(mc),
                          cant_veh=length(unique(mc$file_name)),
                          vel_min=min(mc$velocidad),
                          vel_max=max(mc$velocidad),
                          vel_prom=round(mean(mc$velocidad),2),
                          desv_vel=round(sd(mc$velocidad),2)
      )
      df_indicadores <- rbind(df_indicadores, mc_df)
      rm(mc)
      rm(mc_df)
      gc()
    }
  }
  
  write.table(x = df_indicadores, file = paste("Resumen_Ev-",ev_fin,".csv"), sep = ";", dec = ",", col.names = T, row.names = F)
  save.image(file="Indicadores.RData")
  

  
  
  
  #Exportar los resultados d elas agrupaciones (Resumen)
  Ejecucion <- PER
  
  
  ult_evol <- DCE2[[length(DCE2)]]
  df_indicadores <- data.frame()
  
  if(length(ult_evol$oList)>0){
    for(i in 1:length(ult_evol$oList)){
      mc <- dyclee$oList[[i]]$tramo
      mc_df <- data.frame(id_grupo=(-1*i),
                          celdas=length(unique(dyclee$oList[[i]]$data[,"id"])),
                          total_puntos=sum(dyclee$oList[[i]]$data[,"cant_segmentos"]),#nrow(mc),
                          cant_veh=length(unique(mc$file_name)),
                          vel_min=min(mc$velocidad),
                          vel_max=max(mc$velocidad),
                          vel_prom=round(mean(mc$velocidad),2),
                          desv_vel=round(sd(mc$velocidad),2)
      )
      df_indicadores <- rbind(df_indicadores, mc_df)
      rm(mc)
      rm(mc_df)
      gc()
    }
  }
  
  if(length(ult_evol$aList)>0){
    for(i in 1:length(ult_evol$aList)){
      mc <- dyclee$aList[[i]]$tramo
      mc_df <- data.frame(id_grupo=i,
                          celdas=length(unique(dyclee$aList[[i]]$data[,"id"])),
                          total_puntos=sum(dyclee$aList[[i]]$data[,"cant_segmentos"]),#nrow(mc),
                          cant_veh=length(unique(mc$file_name)),
                          vel_min=min(mc$velocidad),
                          vel_max=max(mc$velocidad),
                          vel_prom=round(mean(mc$velocidad),2),
                          desv_vel=round(sd(mc$velocidad),2)
      )
      df_indicadores <- rbind(df_indicadores, mc_df)
      rm(mc)
      rm(mc_df)
      gc()
    }
  }
  
  print(sum(df_indicadores[,"total_puntos"]))
  
  df_indicadores[is.na(df_indicadores[,"desv_vel"]),"desv_vel"] <- 0
  write.table(x = df_indicadores, file = paste("Resumen_Ev-",Ejecucion, "--", length(DCE2),".csv"), sep = ";", dec = ",", col.names = T, row.names = F)
  save.image(file="Indicadores.RData")
  
  if(FALSE){
    df_indicadores <- df_indicadores[order(df_indicadores[,"vel_prom"]),]
    
    bffr <- data.frame()
    flag <- FALSE
    for(k in 1:nrow(df_indicadores)){
      if(df_indicadores[k, "id_grupo"] < 0){
        bffr <- rbind(bffr, df_indicadores[k, ])
        
      }else if(df_indicadores[k, "id_grupo"] > 0){
        if(nrow(bffr) > 0){
          print(bffr)
        }
        bffr <- data.frame()
      }
    }
  }
  
  setwd(folder_result)
  
  print("FIN ALGORITMO")
  
  
  #===============================================================================
  # Exportacion de tiempos de ejecucion
  
  setwd(folder_result)
  ipak(packages)
  
  tt_dc <- 0 #time_proc
  tt_celdas <- 0 #time_cel
  tt_cc <- 0
  tt_e1 <- 0 #time_e1
  tt_e2 <- 0 #time_e2
  tt_map <- 0
  
  
  
  A <- list()
  for(i in 1:length(time_proc)){
    if(!is.null(time_proc[[i]])){
      tt_dc <- tt_dc + (time_proc[[i]][2]-time_proc[[i]][1])
      A <- A + time_proc[[i]][2]-time_proc[[i]][1]
    }
  }
  B <- list()
  for(i in 1:length(time_cel)){
    if(!is.null(time_cel[[i]])){
      #print(time_cel[[i]][2]-time_cel[[i]][1])
      tt_celdas <- tt_celdas + (time_cel[[i]][2]-time_cel[[i]][1])
      B[[i]] <- time_cel[[i]][2]-time_cel[[i]][1]
    }
  }
  C <- list()
  for(i in 1:length(time_ccel)){
    if(!is.null(time_ccel[[i]])){
      tt_cc <- tt_cc + (time_ccel[[i]][2]-time_ccel[[i]][1])
      C[[i]] <- time_ccel[[i]][2]-time_ccel[[i]][1]
    }
  }
  BC <- list()
  for(i in 1:length(B)){
    BC[[i]] <- as.numeric(B[[i]])+as.numeric(C[[i]])
  }
  D <- list()
  for(i in 1:length(time_e1)){
    if(!is.null(time_e1[[i]])){
      tt_e1 <- tt_e1 + (time_e1[[i]][2]-time_e1[[i]][1])
      D[[i]] <- time_e1[[i]][2]-time_e1[[i]][1]
    }
  }
  E <- list()
  for(i in 1:length(time_e2)){
    if(!is.null(time_e2[[i]])){
      tt_e2 <- tt_e2 + (time_e2[[i]][2]-time_e2[[i]][1])
      E[[i]] <- time_e2[[i]][2]-time_e2[[i]][1]
    }
  }
  F <- list()
  for(i in 1:length(time_map)){
    if(!is.null(time_map[[i]])){
      tt_map <- tt_map + (time_map[[i]][2]-time_map[[i]][1])
      F[[i]] <- time_map[[i]][2]-time_map[[i]][1]
    }
  }
  
  print(paste("Ejecucion:",Ejecucion))
  print(tt_dc)
  print(tt_celdas)
  print(tt_cc)
  print(tt_celdas+tt_cc)
  print(tt_e1)
  print(tt_e2)
  print(tt_map)
  
  tiempos <- list(DC=tt_dc, 
                  Est_celd=end_cel-ini_cel, 
                  Cr_Celdas=tt_celdas, 
                  Celdas_T=tt_celdas+tt_cc,
                  E1=tt_e1, 
                  E2=tt_e2,
                  Mapas=tt_map)
  ciclo <- ev_fin
  write.table(t(as.data.frame(tiempos)), file = "Tiempos.txt", sep = ";", row.names = TRUE, col.names = FALSE)
  write.xlsx2(t(as.data.frame(tiempos)), file = paste("Tiempos 2.xlsx",sep=""), sheetName = "Times", append = TRUE, row.names = TRUE, col.names = TRUE)
  write.xlsx2(as.data.frame(as.numeric(tiempos)), file = paste("Tiempos 3.xlsx",sep=""), sheetName = "Times", append = TRUE, row.names = TRUE, col.names = TRUE)
  
  write.xlsx(as.data.frame(t(rbind(size_puntos))), file = paste("Size_Points Ex_",Ejecucion,"-Ciclo_",ciclo, ".xlsx",sep=""), sheetName = "Size_Point", append = TRUE, row.names = TRUE, col.names = TRUE)
  write.xlsx(as.data.frame(t(rbind(size_celdas))), file = paste("Size_Celds Ex_",Ejecucion,"-Ciclo_",ciclo, ".xlsx",sep=""), sheetName = "Size_Celds", append = TRUE, row.names = TRUE, col.names = TRUE)
  
  
  #write.xlsx2(as.data.frame(as.numeric(rbind(A))), file = paste("T_DC Ex_",Ejecucion,"-Ciclo_",ciclo, ".xlsx",sep=""), sheetName = "t_agrup", append = TRUE, row.names = TRUE, col.names = TRUE)
  write.xlsx2(as.data.frame(as.numeric(rbind(B))), file = paste("T_Celdas Ex_",Ejecucion,"-Ciclo_",ciclo, ".xlsx",sep=""), sheetName = "t_celd", append = TRUE, row.names = TRUE, col.names = TRUE)
  write.xlsx2(as.data.frame(as.numeric(rbind(D))), file = paste("T_E1 Ex_",Ejecucion,"-Ciclo_",ciclo, ".xlsx",sep=""), sheetName = "t_e1", append = TRUE, row.names = TRUE, col.names = TRUE)
  write.xlsx2(as.data.frame(as.numeric(rbind(E))), file = paste("T_E2 Ex_",Ejecucion,"-Ciclo_",ciclo, ".xlsx",sep=""), sheetName = "t_e2", append = TRUE, row.names = TRUE, col.names = TRUE)
  #write.xlsx2(as.data.frame(as.numeric(rbind(F))), file = paste("T_Mapas Ex_",Ejecucion,"-Ciclo_",ciclo, ".xlsx",sep=""), sheetName = "t_map", append = TRUE, row.names = TRUE, col.names = TRUE)
  
  
##
## METRICAS DE VALIDACION (SILHOUETTE)
##

  
  ipak(c("cluster","clusterSim","xlsx","leaflet"))
  options( digits = 15 )
  setwd(folder_result)
  
  #SILHOUETE
  Silhouette <- function(DCE, etapa="0", meth = "manhattan"){
    library(cluster)
    library(xlsx)
    setwd(folder_result)
    dir.create(paste("Silhouette",sep=''))
    folder_Result <- getwd()
    setwd(paste(folder_Result,'/',"Silhouette",sep=''))
    
    calidad_SC_evolucion <- data.frame(matrix(0, nrow=length(DCE), ncol=5))
    colnames(calidad_SC_evolucion) <- c("Evolucion", "Grupos", "Puntos", "SC_Silhouette_Puntos", "SC_Silhouette_Grupos")
    
    for(m in 1:length(DCE)){
      dc <- DCE[[m]]
      cant_grupos <- length(dc$aList)
      calidad_SC_evolucion[m,"Evolucion"] <- m
      calidad_SC_evolucion[m,"Grupos"] <- cant_grupos
      
      if(cant_grupos > 1){
        silhouette <- data.frame()
        for(i in 1:cant_grupos){
          silhouette <- rbind(silhouette, cbind(dc$aList[[i]]$tramo,cluster=i))
        }
        
        distancias <- dist(silhouette[,c("latitude","longitude")],method = meth)
        coef.silueta <- silhouette(silhouette$cluster, distancias)
        summ.coef.silueta <- summary(coef.silueta)
        
        silhouette[,"cluster_cercano"] <- coef.silueta[,"neighbor"]
        silhouette[,"sc_silhouette"] <- coef.silueta[,"sil_width"]
        resumen_grupos <- data.frame(Cluster=c(1:cant_grupos), Elementos=t(summ.coef.silueta$clus.sizes)[1,], SC_Silhouette=summ.coef.silueta$clus.avg.widths)
        
        write.xlsx(silhouette, file = paste("SC_Puntos Ev ", m, ".xlsx"), sheetName = paste("Etapa",etapa), append = TRUE, row.names = FALSE, col.names = TRUE)
        write.xlsx(resumen_grupos, file = paste("SC_Grupos Ev ", m, ".xlsx"), sheetName = paste("Etapa",etapa), append = TRUE, row.names = FALSE, col.names = TRUE)
        #write.table(do.call(rbind,summ.coef.silueta), file = paste("SummarySC Ev_", m, ".txt"), sep = ";", row.names = TRUE, col.names = TRUE)
        
        calidad_SC_evolucion[m,"Puntos"] <- nrow(silhouette)
        calidad_SC_evolucion[m,"SC_Silhouette_Puntos"] <- summ.coef.silueta$avg.width
        calidad_SC_evolucion[m,"SC_Silhouette_Grupos"] <- mean(summ.coef.silueta$clus.avg.widths)
        
        
        rm(distancias)
        rm(coef.silueta)
        rm(summ.coef.silueta)
        #gc()
      }
    }
    s_res <- data.frame(Evolucion="General",
                        Grupos="-",
                        Puntos="-",
                        SC_Silhouette_Puntos=mean(calidad_SC_evolucion[,"SC_Silhouette_Puntos"]),
                        SC_Silhouette_Grupos=mean(calidad_SC_evolucion[,"SC_Silhouette_Grupos"]))
    s_sd <- data.frame(Evolucion="Desv.",
                       Grupos="-",
                       Puntos="-",
                       SC_Silhouette_Puntos=sd(calidad_SC_evolucion[,"SC_Silhouette_Puntos"]),
                       SC_Silhouette_Grupos=sd(calidad_SC_evolucion[,"SC_Silhouette_Grupos"]))
    print(paste("Silhouette Score: ", mean(calidad_SC_evolucion[,"SC_Silhouette_Puntos"])))
    print(paste("Silhouette Score: ", mean(calidad_SC_evolucion[,"SC_Silhouette_Grupos"])))
    calidad_SC_evolucion <- rbind(calidad_SC_evolucion,s_res)
    calidad_SC_evolucion <- rbind(calidad_SC_evolucion,s_sd)
    write.xlsx(calidad_SC_evolucion, file = paste("SC_General.xlsx"), sheetName = paste("Etapa",etapa), append = TRUE, row.names = FALSE, col.names = TRUE)
    
    rm(calidad_SC_evolucion)
    rm(dc)
    gc()
  }
  setwd(folder_result)
  Silhouette(DCE1, "1", meth = "manhattan")
  Silhouette(DCE2, "2", meth = "manhattan")
  
  