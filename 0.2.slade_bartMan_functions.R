####################
## Description:
##  - Functions taken from bartMan. Used for plotting trees from bartMachine
####################

########
## Extract Tree Data

extractTreeData <- function(model, data){
  trees <- extractTrees(model, data)
  
  hideHelper <- function(df){
    class(df) <- c("hideHelper", class(df))
    df
  }
  trees$structure <- hideHelper(trees$structure)
  return(trees)
}

# Main function:
extractTrees <- function(model, data) {
  UseMethod("extractTrees")
}


extractTrees.bartMachine <- function(model, data){
  # library needed
  require(tidyverse)
  
  # Get variable names
  varNames <- model$training_data_features
  
  # Get No of iterations after burn in
  iter <- model$num_iterations_after_burn_in
  
  # function to extract tree data
  nodeData <- vector("list", iter)
  nodeData <- lapply(1:iter,  function(i){
    extract_raw_node_dataSP(model, g = i, iter = iter)
  })
  
  # Melting the tree data into useable format
  df <- rrapply::rrapply(nodeData, how = 'melt')
  nCol <- ncol(df)
  
  suppressMessages(
    df <- df %>%
      pivot_longer(cols = 3:(nCol-1), values_drop_na = TRUE, names_repair = "unique") %>%
      filter(grepl('depth|isLeaf|is_stump|string_location|y_pred|splitValue|splitAttributeM', value...5)) %>%
      select(-name) %>%
      mutate(rn = data.table::rowid(L1, L2, value...5)) %>%
      pivot_wider(names_from = value...5, values_from = value...3) %>%
      select(-rn) %>% as_tibble()
  )
  
  # convert to correct types
  df$depth    <- as.numeric(df$depth)
  df$isLeaf   <- as.logical(df$isLeaf)
  #df$n_eta    <- as.numeric(df$n_eta)
  df$is_stump <- as.logical(df$is_stump)
  df$string_location <- as.character(df$string_location)
  df$splitAttributeM <- as.numeric(df$splitAttributeM)
  df$splitValue <- as.numeric(df$splitValue)
  df$y_pred <- as.numeric(df$y_pred)
  
  
  # match var number to varName
  names(varNames) <- c(0:(length(varNames)-1))
  df$var <- varNames[as.character(df$splitAttributeM)]
  
  # Add tree number (ignoring iteration)
  df$treeNumID <-cumsum(df$string_location=="P")
  
  # define node number sequentially
  df <- df %>%
    group_by(treeNumID) %>%
    mutate(node = row_number())
  
  # get the parent node
  df <- df %>%
    group_by(treeNumID) %>%
    mutate(parentNode = substr(string_location, 0, nchar(string_location)-1))
  
  # # define parent node of P as NA
  df$parentNode[df$string_location == "P"] <- NA
  
  # # define parent nodes nodeID with still empty parent node name as "P"
  df$parentNode[df$parentNode==""] <- "P"
  
  # Match parent node names to node numbers
  df <- df %>%
    group_by(treeNumID) %>%
    mutate(parentNodeNo = match(parentNode, string_location))
  
  # round values
  df$splitValue <- round(as.numeric(df$splitValue),4)
  df$y_pred <- round(as.numeric(df$y_pred),  4)
  
  # add value column
  df <-  df %>%
    mutate(value = coalesce(splitValue, y_pred))
  
  # add label column
  df <- transform(df, label = ifelse(is.na(splitValue), value, paste(var, value, sep = " ≤ ")))
  
  # add new column defining the 'to', for the nodes 'from-to'
  df <- df %>%
    mutate(to = node)
  df <- transform(df, to = ifelse(is.na(parentNode), NA, to))
  
  # turn into tibble
  df <- as_tibble(df)
  
  # rename columns and reorder/remove cols
  names(df) <- c("iteration", "treeNum", "depth", "isLeaf", "isStump", "direction",
                 "splitAtt", "splitValue", "leafValue", "var", "treeNumID",
                 "node", "parentNode", "from", "value", "label", "to")
  df <- df %>%
    select(
      "var",
      "splitValue",
      "node",
      "isLeaf",
      "leafValue",
      "iteration",
      "treeNum",
      "label",
      "value",
      "depth",
      "isStump",
      "direction",
      "parentNode",
      "treeNumID",
      'from',
      "to")
  
  df$iteration <- as.numeric(df$iteration)
  df$treeNum <- as.numeric(df$treeNum)
  
  # add depth column
  df <- rename(df, c('depthAll'= 'depth'))
  
  df <- df %>%
    ungroup() %>%
    group_by(iteration, treeNum) %>%
    mutate(depthMax = max(depthAll))
  
  # get which observations
  
  #dat <- as.data.frame(model$X)
  dat <- model$model_matrix_training_data
  dat <- as.data.frame(dat)
  dat <- dat[,-(length(dat))]
  
  dfObs <-  df %>%
    group_by(iteration, treeNum) %>%
    mutate(obsList = evalNode(dat, var, splitValue))
  
  obsIndex <- lapply(dfObs$obsList, function(x) {
    lapply(x, row.names)
  })
  
  whichObs <- lapply(obsIndex, rapply, f = c)
  whichObs <- lapply(whichObs, as.numeric)
  
  df$obsNode <- whichObs
  
  # get number of observation
  noObser <- NULL
  for(i in 1:nrow(dfObs)){
    noObser[[i]] <- lapply(dfObs$obsList[[i]], dim)
  }
  
  df$noObs <- sapply(noObser, function(y) sum(do.call(rbind, y)[, 1]))
  
  
  
  trees <- list()
  trees$structure <- df
  trees$nMCMC <- iter
  trees$nTree <- model$num_trees
  trees$nVar <- model$p
  trees$varName <- model$training_data_features
  
  
  class(trees) <- c("bartMach", "list")
  
  return(trees)
}

# Function to find obs in each node ---------------------------------------

evalNode <- function(df, x, v) {
  
  out <- vector("list", length(x))
  stk <- vector("list", sum(is.na(x)))
  pos <- 1L
  stk[[pos]] <- df
  
  for (i in seq_along(x)) {
    if (!is.na(x[[i]])) {
      
      subs <- pos + c(0L, 1L)
      stk[subs] <- split(stk[[pos]], stk[[pos]][[x[[i]]]] <= v[[i]])
      
      names(stk)[subs] <- trimws(paste0(
        names(stk[pos]), ",", x[[i]], c(">", "<="), v[[i]]
      ), "left", ",")
      
      out[[i]] <- rev(stk[subs])
      pos <- pos + 1L
    } else {
      
      out[[i]] <- stk[pos]
      stk[[pos]] <- NULL
      pos <- pos - 1L
    }
  }
  return(out)
}


# Function to improve bartMachine speeds ----------------------------------

extract_raw_node_dataSP <- function (bart_machine, g = 1, iter)
{
  
  raw_data_java = .jcall(bart_machine$java_bart_machine, "[LbartMachine/bartMachineTreeNode;",
                         "extractRawNodeInformation", as.integer(g - 1), simplify = TRUE)
  
  raw_data <- vector('list', iter)
  raw_data <- lapply(raw_data_java, bMachineNode)
  raw_data
}


# recursivly go through java object
bMachineNode <- function (node_java)
{
  
  BAD_FLAG_INT = -2147483647
  BAD_FLAG_DOUBLE = -1.7976931348623157e+308
  
  
  node_data = list()
  #node_data = vector("list", 19)
  node_data$java_obj = node_java
  node_data$depth = node_java$depth
  node_data$isLeaf = node_java$isLeaf
  node_data$n_eta = node_java$n_eta
  node_data$is_stump = node_java$isStump()
  node_data$string_location = node_java$stringLocation()
  
  
  if (node_java$splitAttributeM == BAD_FLAG_INT) {
    node_data$splitAttributeM = NA
  }
  else {
    node_data$splitAttributeM = node_java$splitAttributeM
  }
  
  
  if (node_java$splitValue == BAD_FLAG_DOUBLE) {
    node_data$splitValue = NA
  }
  else {
    node_data$splitValue = node_java$splitValue
  }
  
  
  if (node_java$y_pred == BAD_FLAG_DOUBLE) {
    node_data$y_pred = NA
  }
  else {
    node_data$y_pred = node_java$y_pred
  }
  
  
  if (!is.jnull(node_java$left)) {
    node_data$left = bMachineNode(node_java$left)
  }
  else {
    node_data$left = NA
  }
  if (!is.jnull(node_java$right)) {
    node_data$right = bMachineNode(node_java$right)
  }
  else {
    node_data$right = NA
  }
  node_data
}


print.hideHelper <- function(x, ...) {
  x$childLeft <- NULL
  x$childRight <- NULL
  x$parent <- NULL
  x$isStump <- NULL
  x$depthAll <- NULL
  x$direction <- NULL
  x$parentNode <- NULL
  x$treeNumID <- NULL
  x$from <- NULL
  x$to <- NULL
  NextMethod(x, ...)
}


########
## Acceptence of Trees

acceptRate <- function(treeData) {
  df <- treeData$structure
  
  maxIter <- max(df$iteration)
  
  acceptance <- df %>%
    filter(!is.na(var)) %>%
    group_by(iteration, treeNum) %>%
    summarise(values = paste0(sort(unique(label)), collapse = ",")) %>%
    group_by(treeNum) %>%
    mutate(changed = values != lag(values)) %>%
    replace_na(list(changed = TRUE)) %>%
    group_by(iteration) %>%
    summarise(percent_change = mean(changed))
  
  p <- ggplot(acceptance, aes(x = iteration, y = percent_change)) +
    geom_point(alpha = 0.5, colour = 'blue') +
    geom_smooth(formula = y ~ x, data = acceptance,
                method = "lm", colour = "black", se = F) +
    theme_bw() +
    xlab("Iteration") +
    ylab("% Acceptence Rate of Trees")
  
  
  return(p)
}

########
## Tree Depth

treeDepth <- function(treeData) {
  maxIter <- max(treeData$structure$iteration)
  
  newTrees <- treeData$structure %>%
    select(iteration, treeNum, depthMax) %>%
    as_tibble() %>%
    group_by(iteration, treeNum) %>%
    summarize(maxDepth = max(depthMax)) %>%
    ungroup() %>%
    group_by(iteration) %>%
    summarize(avgDepth = mean(maxDepth))
  
  ylimMax <- max(newTrees$avgDepth)
  
  p <-  ggplot(newTrees, aes(iteration, avgDepth)) +
    geom_point(alpha = 0.5, colour = 'blue') +
    #geom_line(alpha = 0.5, colour = 'blue') +
    geom_smooth(formula = y ~ x, method = "loess", colour = "black", se = F) +
    theme_bw() +
    xlab("Iteration") +
    ylab("Average Tree Depth")
  

  return(p)
}

########
## Tree Nodes

treeNodes <- function(treeData) {
  df <- treeData$structure
  maxIter <- max(df$iteration)
  
  df <- df %>%
    group_by(iteration, treeNum) %>%
    summarize(count = n()) %>%
    group_by(iteration) %>%
    summarize(new = mean(count))
  
  p <- ggplot(df, aes(iteration, new)) +
    geom_point(alpha = 0.5, colour = 'blue') +
    # geom_line(alpha = 0.5, colour = 'blue') +
    geom_smooth(formula = y ~ x, method = "loess", colour = "black", se = F) +
    theme_bw() +
    xlab("Iteration") +
    ylab("Average Tree Nodes")
  
  return(p)
}


########
## Split in Densities

splitDensity <- function(treeData, data, colBy = NULL, display = "histogram") {
  
  if (!(display %in% c("histogram", "ridge", "density", 'both'))) {
    stop("display must be \"histogram\", \"ridge\", \"density\", or \"both\"")
  }
  
  # get just the variable and split value
  tt <- treeData$structure %>%
    ungroup() %>%
    select(var, splitValue) %>%
    na.omit()
  
  # create plotting order to match order of data
  nam <- treeData$varName
  tt$var <- factor(tt$var, levels = nam)
  
  varNames <- unique(tt$var)
  
  # create plot
  
  if (display == "density") {
    dPlot <- tt %>%
      ggplot(aes(x = splitValue)) +
      geom_density(aes(colour = var, fill = var)) +
      facet_wrap(~var) +
      theme_bw() +
      ylab("Density") +
      xlab("Split value") +
      theme(legend.position = "none")
  } else if (display == "ridge") {
    dPlot <- tt %>%
      ggplot(aes(x = splitValue, y = var, fill = stat(x))) +
      geom_density_ridges(aes(fill = var, alpha = 0.1)) +
      ylab("Variable") +
      xlab("Split value") +
      theme_bw() +
      theme(legend.position = "none")
  } else if(display == "histogram") {
    dPlot <- tt %>%
      ggplot(aes(x = splitValue)) +
      geom_histogram(aes(colour = var, fill = var), bins = 30) +
      facet_wrap(~var) +
      theme_bw() +
      ylab("Density") +
      xlab("Split value") +
      theme(legend.position = "none")
  }else if(display == 'both'){
    
    dataIdx <- which((names(data) %in% varNames))
    dat <- data[, dataIdx]
    
    meltDat <- melt(dat)
    names(tt) <- c('variable', 'value')
    
    dataList <- list(meltDat, tt)
    names(dataList)  <- c('dat', 'sv')
    dfList <- plyr::ldply(dataList)
    
    dPlot <- ggplot(dfList) +
      geom_density(aes(x = value, fill = .id), alpha = 0.5) +
      facet_wrap(~variable) +
      scale_fill_discrete(name = "", labels = c("Data", "Split Value")) +
      ylab('Density') +
      xlab("Split value") +
      theme_bw()
  }
  
  suppressMessages(print(dPlot))
  #return(dPlot)
}


########
## Proximity Matrix

proximityMatrix <- function(treeData, data, nRows, normalize = TRUE, reorder = TRUE, iter = NULL) {
  
  
  if (!is.null(iter)) {
    treeData$structure <- treeData$structure %>%
      filter(iteration == iter)
    treeTotal <- max(treeData$structure$treeNum)
  } else {
    treeNumber <- treeData$nTree
    iterNumber <- treeData$nMCMC
    treeTotal <- treeNumber*iterNumber
  }
  
  dfObs <- treeData$structure %>%
    select(var, splitValue, iteration, treeNum, value, noObs)
  whichObs <- treeData$structure$obsNode
  
  # get all indicies
  idxLength <- length(whichObs)
  
  # rename the list elements to work with stack function
  myLetters <- function(n) {
    unlist(Reduce(paste0,
                  replicate(n %/% length(letters), letters, simplify = FALSE),
                  init = letters,
                  accumulate = TRUE
    ))[1:n]
  }
  
  myNames <- myLetters(idxLength)
  allIndex <- setNames(whichObs, myNames)
  
  #add back to dataframe to find terminal nodes
  dfObs$whichNode <- allIndex
  
  # filter terminal nodes
  dfTerm <- dfObs %>%
    filter(is.na(var))
  
  # get observations
  allIndex <- dfTerm$whichNode
  
  # turn into matrix
  resMat <- tcrossprod(table(stack(allIndex)))
  diag(resMat) <- 0
  
  # normalize
  if (normalize) {
    resMat <- resMat / treeTotal
  }
  
  if (reorder) {
    resMat <- proxReorder(resMat)
  }
  
  
  return(resMat)
}

# reorder proximity matrix ------------------------------------------------

proxReorder <- function(d) {
  prox <- as.dist(d)
  score <- apply(as.matrix(prox), 1, max)
  o <- DendSer::dser(-prox, -score, cost = DendSer::costLS)
  res <- d[o, o]
  class(res) <- class(d)
  res
}
########
## Plot Proximity

plotProximity <- function(matrix,
                          pal = rev(colorspace::sequential_hcl(palette = "Blues 2", n = 100)),
                          limit = NULL){
  
  # Set the limits
  if (is.null(limit)) {
    limit <- range(as.dist(matrix))
    limit <- range(labeling::rpretty(limit[1], limit[2]))
  }
  
  # melt the matrox
  suppressWarnings(
    df <- reshape::melt(matrix)
  )
  
  # order axis names
  varNames <- unique(df$values)
  df$values   <- factor(df$values,   levels = varNames)
  df$values.1 <- factor(df$values.1, levels = varNames)
  
  p <- ggplot(df, aes(values, values.1)) +
    geom_point(aes(size = value), show.legend = F) +
    geom_tile(aes(fill = value)) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev(levels(df$values.1))) +
    scale_fill_gradientn(
      colors = pal, limits = limit, name = "Proximity",
      guide = guide_colorbar(
        order = 1,
        frame.colour = "black",
        ticks.colour = "black"
      ), oob = scales::squish
    ) +
    xlab("") +
    ylab("") +
    theme_bw() +
    theme(aspect.ratio = 1)
  
  return(p)
  
}


########
## Plot Tree

plotTree <- function(treeData,
                     iter = 1,
                     treeNo = 1,
                     plotType = c("dendrogram", "icicle")) {
  
  p <- plotAnyTree(treeData, iter = iter, treeNo = treeNo)
  
  if (plotType == "dendrogram") {
    gp <- ggraph::ggraph(p, "dendrogram") +
      ggraph::geom_edge_elbow() +
      ggraph::geom_node_label(aes(label = label, color = label)) +
      ggraph::theme_graph() +
      theme(legend.position = "none")
  } else if (plotType == "icicle") {
    gp <- ggraph::ggraph(p, "partition") +
      ggraph::geom_node_tile(aes(fill = var), size = 0.25) +
      ggraph::geom_node_label(aes(label = label, color = var)) +
      scale_y_reverse() +
      theme_void() +
      theme(legend.position = "none")
  }
  
  return(gp)
}



# -------------------------------------------------------------------------


# Main plot function:
plotAnyTree <- function(treeData, iter = 1, treeNo = 1) {
  UseMethod("plotAnyTree")
}

# BARTMACHINE -------------------------------------------------------------

plotAnyTree.bartMach <- function(treeData,
                                 treeNo = 1,
                                 iter = 1) {
  df <- treeData$structure
  
  # select tree num from iteration
  df <- df %>%
    filter(iteration == iter, treeNum == treeNo)
  
  # create dataframe of edges
  dfOfEdges <- data.frame(
    from = df$from,
    to = df$to
  )
  
  dfOfEdges <- na.omit(dfOfEdges)
  
  # remove unnecessary columns
  df <- df %>%
    select(-isStump, -to, -from, -node, -parentNode, -treeNumID)
  
  
  # Turn into a table graph object
  singleTree <- tidygraph::tbl_graph(nodes = df, edges = dfOfEdges)
  
  return(singleTree)
}

########
## Plot All trees

plotAllTrees <- function(treeData,
                         iter = NULL,
                         treeNo = NULL,
                         sampleSize = NULL,
                         cluster = NULL,
                         sizeNode = TRUE,
                         pal = RColorBrewer::brewer.pal(9, "Purples"),
                         fillBy = NULL,
                         selectedVars = NULL,
                         removeStump = FALSE
) {
  
  if(length(selectedVars) > length(treeData$varName)){
    message("SelectedVars is longer than number of available variables. Selecting all variables")
    selectedVars <- c(1:length(treeData$varName))
  }
  
  allTrees <- plotAll(treeData, iter = iter, treeNo = treeNo, cluster = cluster)
  
  suppressWarnings(
    p <- plotAllTreesPlotFn(allTrees,
                            sampleSize = sampleSize,
                            sizeNode = sizeNode,
                            pal = pal,
                            fillBy = fillBy,
                            name = treeData$varName,
                            selectedVars = selectedVars,
                            removeStump = removeStump
    )
    
  )
  return(p)
}


plotAll <- function(treeData, iter = NULL, treeNo = NULL, cluster = NULL) {
  UseMethod("plotAll")
}


# bartMachine -------------------------------------------------------------

plotAll.bartMach <- function(treeData, iter = NULL, treeNo = NULL, cluster = NULL) {
  
  df <- treeData$structure
  maxIter <- treeData$nMCMC
  noObservations <- treeData$structure$noObs[1]
  
  if (is.null(iter) & is.null(treeNo)) {
    df <- df %>%
      filter(iteration == maxIter)
  } else if (is.null(iter) & !is.null(treeNo)) {
    df <- df %>%
      filter(treeNum == treeNo)
  } else if (!is.null(iter) & is.null(treeNo)) {
    df <- df %>%
      filter(iteration == iter)
  } else {
    df <- df %>%
      filter(iteration == iter, treeNum == treeNo)
  }
  
  df <- df %>%
    group_by(iteration, treeNum)
  
  # -------------------------------------------------------------------------
  
  # add mean response per node:
  respNode <- df$obsNode #apply(df, 1, function(x) {y[x$obsNode]})
  respNode <- lapply(respNode, mean)
  df$respNode <- unlist(respNode)
  
  
  
  # cluster trees
  suppressWarnings(
    if (!is.null(cluster)) {
      if (cluster == "depth") {
        df <- df[with(df, order(-depthMax)), ]
        
        # split the dataframe into a list of dfs, one for each tree
        list_edges <- df %>%
          ungroup() %>%
          group_split(cumsum(noObs == noObservations))
      }
    }
  )
  
  suppressWarnings(
    if(cluster == 'var' || is.null(cluster)){
      # split the dataframe into a list of dfs, one for each tree
      list_edges <- df %>%
        group_split(cumsum(noObs == noObservations))
    }
  )
  
  
  # remove unnecessary columns
  treesSplit <- lapply(list_edges, function(x) {
    x["isStump"] <- x["to"] <- x["from"] <- x["node"] <- x["parentNode"] <- x["treeNumID"] <- NULL
    x
  })
  
  # create dataframe of edges
  dfOfEdges <- lapply(list_edges, function(df_tree) {
    res <- data.frame(
      from = df_tree$from,
      to = df_tree$to
    )
    # delete NAs from result
    res <- na.omit(res)
    return(res)
  })
  
  
  # Turn into a table graph object
  eachTree <- list()
  for (i in 1:length(treesSplit)) {
    eachTree[[i]] <- tidygraph::tbl_graph(treesSplit[[i]], dfOfEdges[[i]])
  }
  
  if (!is.null(cluster)) {
    if (cluster == "var") {
      eachTree <- clusterTrees(eachTree)
    }
  }
  
  return(eachTree)
}


# -------------------------------------------------------------------------
# cluster function by variable --------------------------------------------
# -------------------------------------------------------------------------

clusterTrees <- function(treeList) {
  
  #df <- data
  indIDS <- map(treeList, function(x) {
    x %>%
      pull(var) %>%
      replace_na("a") %>%
      paste0(collapse = "")
  }) %>%
    unlist(use.names = F) %>%
    as_tibble() %>%
    mutate(ids = 1:n()) %>%
    group_by(value) %>%
    mutate(count = n():1) %>%
    arrange(value)
  
  ind <- indIDS %>%
    group_by(value) %>%
    mutate(valrank = max(count)) %>%
    ungroup() %>%
    arrange(-valrank, value, -count) %>%
    pull(ids)
  
  
  treeList <- treeList[ind]
  
  return(treeList)
}


plotAllTreesPlotFn <- function(treeList,
                               sampleSize = NULL,
                               sizeNode = TRUE,
                               pal =  rev(colorspace::sequential_hcl(palette = "Purples 3", n = 100)),
                               fillBy = NULL,
                               selectedVars = NULL,
                               removeStump = FALSE,
                               name) {
  
  # remove stumps
  if(removeStump){
    treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)
    stumpIdx <- NULL
  }else{
    # get the stump index
    whichStump = NULL
    for(i in 1:length(treeList)){
      whichStump[[i]] <-  which(igraph::gsize(treeList[[i]]) == 0)
    }
    stumpIdx <- which(whichStump == 1)
    
    if(length(stumpIdx >=1)){
      # create new tree list
      newTrees <- treeList[stumpIdx]
      
      # create df of tree stumps
      newTreesDF <- NULL
      for(i in 1:length(newTrees)){
        newTreesDF[[i]] <- newTrees[[i]] %>%
          tidygraph::activate(nodes) %>%
          data.frame()
        newTreesDF[[i]]$var <- "Stump"
        
      }
      
      
      # create edge data for stumps
      newDF_Nodes <- newDF_Edges <- NULL
      for(i in 1:length(newTreesDF)){
        newDF_Nodes[[i]] <- rbind(newTreesDF[[i]], newTreesDF[[i]][rep(1), ])
        newDF_Edges[[i]] <- data.frame(from = c(1,1), to = c(1,1))
      }
      
      # turn into tidygraph trees
      newTree <- NULL
      for(i in 1:length(newDF_Nodes)){
        newTree[[i]] <- tbl_graph(nodes = newDF_Nodes[[i]], edges = newDF_Edges[[i]])
      }
      # replace stumps with new stumps
      treeList[stumpIdx] <- newTree
    }
  }
  
  
  # get limits
  if(!is.null(fillBy)){
    if(fillBy == 'response'){
      lims <- range(unlist(lapply(treeList, . %>% tidygraph::activate(nodes) %>% pull(respNode))))
      lims <- labeling::rpretty(lims[1], lims[2])
      lims <- c(min(lims), max(lims))
      nam <- 'Mean \nResponse'
    } else if(fillBy == "mu"){
      lims <- range(unlist(lapply(treeList, . %>% tidygraph::activate(nodes) %>% filter(is.na(var) | var == "Stump") %>% pull(value))))
      lims <- c(-max(abs(lims)), max(abs(lims)))
      nam <- 'Mu'
    }
  }else{
    nam <- 'Variable'
  }
  
  
  # set node colours
  if(!is.null(selectedVars)){
    nodeNamesImp <- name[selectedVars]
    nodeNamesOthers <- name[-selectedVars]
    nodeColorsImp <- setNames(scales::hue_pal(c(0, 360) + 15, 100, 64, 0, 1)(length(nodeNamesImp)), nodeNamesImp)
    nodeColorsOth <- setNames(rep('#e6e6e6', length(nodeNamesOthers)), nodeNamesOthers)
    nodecolors <- c(nodeColorsImp, nodeColorsOth)
    
    namedOthers <- setNames('#e6e6e6',  "Others")
    legColours <- c(nodeColorsImp, namedOthers)
    dfLegend <- reshape::melt(legColours) %>%
      tibble::rownames_to_column(var = 'varName')
    dfLegend$val <- rep(1, times = length(dfLegend$varName))
    dfLegend$varName <- factor(dfLegend$varName, levels = names(legColours))
    pLeg <- ggplot(dfLegend, aes(x = varName, y = val, fill = varName)) +
      geom_bar(stat = 'identity') +
      scale_fill_manual(values = dfLegend$value, name = 'Variable')
    if(length(stumpIdx) >=  1){
      nodecolors[["Stump"]] <- '#e6e6e6'
    }
  }else{
    nodeNames <- unique(na.omit(unlist(lapply(treeList, . %>% tidygraph::activate(nodes) %>% pull(var)))))
    nodeNames <- sort(nodeNames)
    nodecolors <- setNames(scales::hue_pal(c(0, 360) + 15, 100, 64, 0, 1)(length(nodeNames)), nodeNames)
    
    # colour stumps grey
    if(length(stumpIdx) >=  1){
      if(is.null(fillBy)){
        nodecolors[["Stump"]] <- '#808080'
      }else if(fillBy == 'response'){
        nodecolors[["Stump"]] <- min(pal)
      }else if(fillBy == 'mu'){
        # stumpVal <- lapply(treeList, . %>% tidygraph::activate(nodes) %>% filter(is.na(var) | var == "Stump") %>% pull(value))
        # sv  <- as.data.frame(stumpVal[stumpIdx])
        # cdfLims <- ecdf(lims)
        # cdfVal <- cdfLims(sv[1,1])
        # nodecolors[["Stump"]] <-  pal[abs(length(pal) * cdfVal)]
        nodecolors[["Stump"]] <- pal[ceiling(length(pal)/2)]
      }
    }
  }
  
  suppressMessages(
    allPlots <- lapply(treeList,
                       plotFun,
                       n = length(treeList),
                       color = nodecolors,
                       sizeNode = sizeNode,
                       pal = pal,
                       range = lims,
                       name = nam,
                       fill = fillBy)
  )
  
  
  # get legend
  if(!is.null(selectedVars) & !is.null(fillBy)){
    themeMargin <- theme(legend.box.margin = margin(100, 15, 80, 20))
    legend1 <- cowplot::get_legend(pLeg + themeMargin)
    legend <- cowplot::get_legend(allPlots[[1]] + themeMargin)
    legend$grobs[[2]] <- legend1
  }else if(!is.null(selectedVars) & is.null(fillBy)){
    themeMargin <- theme(legend.box.margin = margin(100, 15, 80, 20))
    legend <- cowplot::get_legend(pLeg + themeMargin)
  }else{
    legend <- cowplot::get_legend(allPlots[[1]])
  }
  
  
  # remove legends from individual plots
  allPlots <- lapply(allPlots, function(x) x + theme(legend.position = "none"))
  if(removeStump == FALSE){
    for(i in stumpIdx){
      allPlots[[i]]$data <- allPlots[[i]]$data[-2, ]
    }
  }
  n <- length(allPlots)
  nRow <- floor(sqrt(n))
  allTreesPlot <- gridExtra::arrangeGrob(grobs=allPlots, nrow=nRow)
  
  cowplot::plot_grid(allTreesPlot, legend, rel_widths = c(0.9, 0.13), ncol = 2)
  
  
}

# Plot function -----------------------------------------------------------

plotFun <- function(List,
                    color = NULL,
                    n,
                    sizeNode = TRUE,
                    pal = rev(colorspace::sequential_hcl(palette = "Purples 3", n = 100)),
                    range,
                    name,
                    fill) {
  
  if(is.null(pal)){
    pal = "grey"
  }
  
  if(!is.null(fill)){
    if (fill == "response") {
      if (sizeNode) {
        plot <- ggraph::ggraph(List, "partition", weight = noObs) +
          ggraph::geom_node_tile(aes(fill = var), size = 0.25) +
          ggraph::geom_node_text(aes(label = ""), size = 4) +
          scale_y_reverse() +
          theme_void()
        if (!is.null(colors)) {
          plot <- plot + scale_fill_manual(values = color, name = "Variable") +
            ggnewscale::new_scale_fill() +
            ggnewscale::new_scale_color() +
            ggraph::geom_node_tile(
              size = 0.15,
              data = . %>% filter(is.na(var)),
              aes(fill = respNode)
            ) +
            scale_fill_gradientn(
              colours = pal,
              limits = range,
              name = name,
              guide = guide_colorbar(
                frame.colour = "black",
                ticks.colour = "black",
                order = 2
              )
            )
        }
      } else {
        plot <- ggraph::ggraph(List, "partition") +
          ggraph::geom_node_tile(aes(fill = var), size = 0.25) +
          ggraph::geom_node_text(aes(label = ""), size = 4) +
          scale_y_reverse() +
          theme_void()
        if (!is.null(colors)) {
          plot <- plot + scale_fill_manual(values = color, name = "Variable") +
            ggnewscale::new_scale_fill() +
            ggnewscale::new_scale_color() +
            ggraph::geom_node_tile(
              size = 0.15,
              data = . %>% filter(is.na(var)),
              aes(fill = respNode)
            ) +
            scale_fill_gradientn(
              colours = pal,
              limits = range,
              name = name,
              guide = guide_colorbar(
                frame.colour = "black",
                ticks.colour = "black",
                order = 2
              )
            )
        }
      }
    }
  }
  
  if(!is.null(fill)){
    if (fill == "mu") {
      if (sizeNode) {
        plot <- ggraph::ggraph(List, "partition", weight = noObs) +
          ggraph::geom_node_tile(aes(fill = var), size = 0.25) +
          ggraph::geom_node_text(aes(label = ""), size = 4) +
          scale_y_reverse() +
          theme_void()
        if (!is.null(colors)) {
          plot <- plot + scale_fill_manual(values = color, name = "Variable") +
            ggnewscale::new_scale_fill() +
            ggnewscale::new_scale_color() +
            ggraph::geom_node_tile(
              size = 0.15,
              data = . %>% filter(is.na(var)),
              aes(fill = value)
            ) +
            scale_fill_gradientn(
              colours = pal,
              limits = range,
              name = name,
              guide = guide_colorbar(
                frame.colour = "black",
                ticks.colour = "black",
                order = 2
              )
            )
        }
      } else {
        plot <- ggraph::ggraph(List, "partition") +
          ggraph::geom_node_tile(aes(fill = var), size = 0.25) +
          ggraph::geom_node_text(aes(label = ""), size = 4) +
          scale_y_reverse() +
          theme_void()
        if (!is.null(colors)) {
          plot <- plot + scale_fill_manual(values = color, name = "Variable") +
            ggnewscale::new_scale_fill() +
            ggnewscale::new_scale_color() +
            ggraph::geom_node_tile(
              size = 0.15,
              data = . %>% filter(is.na(var)),
              aes(fill = value)
            ) +
            scale_fill_gradientn(
              colours = pal,
              limits = range,
              name = name,
              guide = guide_colorbar(
                frame.colour = "black",
                ticks.colour = "black",
                order = 2
              )
            )
        }
      }
    }
  }
  
  if(is.null(fill)){
    if (sizeNode) {
      plot <- ggraph::ggraph(List, "partition", weight = noObs) +
        ggraph::geom_node_tile(aes(fill = var), size = 0.25) +
        ggraph::geom_node_text(aes(label = ""), size = 4) +
        scale_y_reverse() +
        theme_void()
      if (!is.null(colors)) {
        plot <- plot + scale_fill_manual(values = color, name = "Variable")
      }
    } else {
      plot <- ggraph::ggraph(List, "partition") +
        ggraph::geom_node_tile(aes(fill = var), size = 0.25) +
        ggraph::geom_node_text(aes(label = ""), size = 4) +
        scale_y_reverse() +
        theme_void()
      if (!is.null(colors)) {
        plot <- plot + scale_fill_manual(values = color, name = "Variable")
      }
    }
  }
  
  return(plot)
}



########
## Tree Bar Plot

treeBarPlot <- function(treeData,
                        iter = NULL,
                        treeNo = NULL,
                        topTrees = NULL,
                        removeStump = FALSE){
  
  library(patchwork)
  
  treeList <- plotAll(treeData, iter = iter, treeNo = treeNo, cluster = NULL)
  
  # remove stumps
  treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)
  
  if(removeStump){
    treeList <- Filter(function(x) igraph::gsize(x) > 0, treeList)
    stumpIdx <- NULL
  }else{
    # get the stump index
    whichStump = NULL
    for(i in 1:length(treeList)){
      whichStump[[i]] <-  which(igraph::gsize(treeList[[i]]) == 0)
    }
    stumpIdx <- which(whichStump == 1)
    
    if(length(stumpIdx >=1)){
      # create new tree list
      newTrees <- treeList[stumpIdx]
      
      # create df of tree stumps
      newTreesDF <- NULL
      for(i in 1:length(newTrees)){
        newTreesDF[[i]] <- newTrees[[i]] %>%
          tidygraph::activate(nodes) %>%
          data.frame()
        newTreesDF[[i]]$var <- "Stump"
      }
      # create edge data for stumps
      newDF_Nodes <- newDF_Edges <- NULL
      for(i in 1:length(newTreesDF)){
        newDF_Nodes[[i]] <- rbind(newTreesDF[[i]], newTreesDF[[i]][rep(1), ])
        newDF_Edges[[i]] <- data.frame(from = c(1,1), to = c(1,1))
      }
      # turn into tidygraph trees
      newTree <- NULL
      for(i in 1:length(newDF_Nodes)){
        newTree[[i]] <- tbl_graph(nodes = newDF_Nodes[[i]], edges = newDF_Edges[[i]])
      }
      # replace stumps with new stumps
      treeList[stumpIdx] <- newTree
    }
  }
  
  
  
  # get the frequency of similar trees:
  freqs <- map(treeList, function(x){
    x %>%
      pull(var) %>%
      tidyr::replace_na("..") %>%
      paste0(collapse = "")
  }) %>%
    unlist(use.names = F) %>%
    as_tibble() %>%
    mutate(ids = 1:n()) %>%
    group_by(value) %>%
    mutate(val = n():1)
  
  
  freqDf <-  freqs %>%
    slice(1) %>%
    arrange(-val) %>%
    rename(frequency = val)  # frequency tibble
  freqDf$treeNum <- seq(1:nrow(freqDf)) # add tree number
  
  
  if(!is.null(topTrees)){
    freqDf <- freqDf[1:topTrees,]
  }
  
  lengthFreq <- length(freqDf$value)
  ids <- freqDf$ids
  #ids <- freqs %>% slice(1) %>% pull(ids) # remove duplicates
  freqs <- freqs[ids,] %>% pull(val) # get frequencies
  
  
  treeList <- treeList[ids]
  treeListNew <- purrr::imap(treeList, ~.x %>%
                               mutate(frequency = freqs[.y]) %>%
                               select(var, frequency))
  
  # # return new list of trees
  # treeList <- treeListNew[sort(ids)]
  treeList <- treeListNew
  
  # add plot name as number
  # for(i in 1:(length(treeList))){
  #   treeList[[i]] <- treeList[[i]] %>%
  #     tidygraph::activate(nodes) %>%
  #     mutate(name = c(i, rep("", length.out = igraph::gsize(treeList[[i]]))))
  # }
  
  
  # Create barplot of frequency ---------------------------------------------
  names <- factor(freqDf$value, levels = freqDf$value)
  
  bp <- freqDf %>%
    ggplot() +
    geom_bar(aes(x = value, y = frequency), fill = 'steelblue', stat = "identity") +
    scale_x_discrete(limits = rev(levels(names))) +
    ggtitle("") +
    ylab("Count") +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none") +
    coord_flip()
  
  
  # ggraph plotting funtion -------------------------------------------------
  
  
  # set node colours
  nodenames <- unique(na.omit(unlist(lapply(treeList, .%>% tidygraph::activate(nodes) %>% pull(var)))))
  nodenames <- sort(nodenames)
  nodecolors <- setNames(scales::hue_pal(c(0,360)+15, 100, 64, 0, 1)(length(nodenames)), nodenames)
  
  if(length(stumpIdx) >=  1){
    nodecolors[["Stump"]] <- '#808080'
  }
  
  
  
  # plotting function
  plotFun <- function(List, colors = NULL, n) {
    
    plot <- ggraph::ggraph(List, "partition") +
      ggraph::geom_node_tile(aes(fill = var), size = 0.25) +
      ggraph::geom_node_text(aes(label = ''), size = 4) +
      theme(legend.position = "bottom") +
      scale_y_reverse() +
      theme_void()
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = colors, name = "Variable") +
        scale_color_manual(values = colors, na.value = "grey")  +
        theme(legend.position = "bottom")
    }
  }
  
  
  allPlots <- lapply(treeList, plotFun, n = length(treeList), color = nodecolors)
  
  # get legend
  legend <- cowplot::get_legend(allPlots[[1]])
  
  # remove legends from individual plots
  allPlots <- lapply(allPlots, function(x) x + theme(legend.position = "none"))
  
  if(removeStump == FALSE){
    
    whichStumpData = NULL
    for(i in 1:length(allPlots)){
      whichStumpData[[i]] <-  which(allPlots[[i]]$data$leaf[1] == TRUE)
    }
    
    stumpIdxNew <- which(whichStumpData == 1)
    for(i in stumpIdxNew){
      allPlots[[i]]$data <- allPlots[[i]]$data[-2, ]
    }
  }
  
  # filter top X% of plots
  if(!is.null(topTrees)){
    allPlots <- allPlots[1:topTrees]
  }
  
  # Create final barplot ----------------------------------------------------
  width = 1
  p_axis <- ggplot(freqDf) +
    geom_blank(aes(y = value)) +
    purrr::map2(allPlots,
                rev(seq_along(allPlots)),
                ~ annotation_custom(ggplotGrob(.x),
                                    ymin = .y - width / 2,
                                    ymax = .y + width / 2,
                )) +
    theme_void()
  
  bp1 <- bp + theme(axis.text.y = element_blank())
  ppp <- p_axis + theme(aspect.ratio = 6)
  px <- ppp|bp1
  
  bpFinal <- cowplot::plot_grid(px, legend, rel_heights = c(1, .1), ncol = 1)
  
  return(bpFinal)
  
}
