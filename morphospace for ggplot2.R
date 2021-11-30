#For PCA based on Elliptic Fourier Analysis(not rfourier, dfourier and Procrustes etc).
#Visualize "morphospace" which calculated from PC plane on ggplot2 package.  

#preparation 
#define some essential functions(from source codes of Momocs. "https://rdrr.io/cran/Momocs/src/R/gr-morphospaces.R")
#PCA2shp_efourier(), to extract coo object for morphospace.
#domestic,work under the function "morphospacePCA_custom()". Non-independent.

#================================================================================
PCA2shp_efourier <- function(pos, rot, mshape, amp.shp = 1, pts.shp = 60) {
  if (ncol(pos) != ncol(rot))
    stop("'rot' and 'pos' must have the same ncol")
  if (length(mshape) != nrow(rot))
    stop("'mshape' and ncol(rot) lengths differ")
  nb.h <- length(mshape)/4
  n <- nrow(pos)
  # we prepare the array
  res <- list()
  for (i in 1:n) {
    ax.contrib <- .mprod(rot, pos[i, ]) * amp.shp
    coe <- mshape + apply(ax.contrib, 1, sum)
    xf <- coeff_split(coe)
    coo <- efourier_i(xf, nb.h = nb.h, nb.pts = pts.shp)
    # reconstructed shapes are translated on their centroid if
    # (trans) {
    dx1 <- pos[i, 1] - coo_centpos(coo)[1]
    dy1 <- pos[i, 2] - coo_centpos(coo)[2]
    coo <- coo_trans(coo, dx1, dy1)
    # }
    res[[i]] <- coo
  }
  return(res)
}

#================================================================================
#domestic, too. For rendering list type data to data.frame.
Unzip.list<-function(...) rbind(data.frame(), ...)

#================================================================================
.mprod <- function(m, s) {
  res <- m
  for (i in 1:ncol(m)) {
    res[, i] <- m[, i] * s[i]
  }
  return(res)
}
#================================================================================
.wdw <- function() {
  wdw <- par("usr")
  x <- wdw[2] - wdw[1]
  y <- wdw[4] - wdw[3]
  return(c(x, y))
}
#================================================================================

#visualize "morphospace" 
#based on "morphospacePCA()"
#set PCA object (PCA<-PCA(efourier(data)))
#nr.shp * nc.shp must be = nb.shp(ex:nb.shp=18,nr.shp=6,nc.shp=3)
#ex: xax="PC1",yax="PC2" (PCA$x[,PC1]) 
#pos.shp="range", if you try other arguments, see "?morphospace_positions" 

morphospacePCA_custom<-function(PCA, xax, yax, pos.shp="range", nb.shp,
                                nr.shp, nc.shp, amp.shp = 1,
                                rotate.shp = 0,
                                flipx.shp = FALSE,
                                flipy.shp = FALSE,
                                size.shp = 0.5, wdw=max(.wdw()),
                                pts.shp = 60,
                                col.shp=col_alpha("#000000", 0.95), border.shp = col_alpha("#000000", 0.5),
                                lwd.shp = 1, plot = TRUE) {
  # we check here, though it should have been before
  if (length(PCA$method)>4 | is.null(PCA$method)) {
    stop("morphospacePCA needs a $method of length <= 5")}
  if (nb.shp!=nr.shp*nc.shp) {
    stop("nb.shp should be equal to the product of nr and nc.shp")}
  
  # we retrive the values corresponding to the two plotted axes and the meanshape
  xy <- PCA$x[, c(xax, yax)]
  rot <- PCA$rotation[, c(xax, yax)]
  mshape <- PCA$mshape
  # we define the position of shapes
  pos <- morphospace_positions(xy, pos.shp = pos.shp, nb.shp = nb.shp,
                               nr.shp = nr.shp, nc.shp = nc.shp)
  # according to the type of morphometrics applied, we switch the method
  # and the way we plot reconstruct shapes (polygon, lines, points for Out, Opn, Ldk)
  # when the object combines different morphometric approaches (up to 4)
  # their size is divided by 2 and the shapes and set (of d) around the (x; y) coordinates of pos.shp
  method <- PCA$method
  lm <- length(method)
  if (length(rotate.shp) != lm){
    rotate.shp <- rep(rotate.shp, lm)
  }
  if (length(flipx.shp) != lm){
    flipx.shp <- rep(flipx.shp, lm)
  }
  if (length(flipy.shp) != lm){
    flipy.shp <- rep(flipy.shp, lm)
  }
  
  if (length(size.shp)!=lm) size.shp <- rep(size.shp[1], lm)
  size.shp.final <- (size.shp*wdw/14) / ifelse(lm<2, 1, 2)
  d <- mean(size.shp.final) / 2
  # here we define the translation x and y for every sub-morphoshape
  # and the coe to retrieve
  dx <- 0
  dy <- 0
  # indices of successive coe to select
  col.start <- 1
  col.end   <- length(mshape)
  # not very satisfactory...
  # hack in case of multi
  if (!plot) SHP <- list()
  for (i in seq(along=method)){
    shp <- NULL
    plot.method <- NULL
    ids <- col.start[i]:col.end[i]}
    # outlines
    # efourier
    if (method[i] == "efourier") {
      shp <- PCA2shp_efourier(pos = pos, rot = rot[ids, ], mshape = mshape[ids],
                              amp.shp = amp.shp, pts.shp = pts.shp)
      
      shp <- lapply(shp, coo_close)
      plot.method <- "poly"}
    
    ###Here we will ignore other methods(for speed sake), if you need, see "morphospace for ggplot2.R". 
    
    ### Then...
    # we template shapes
    shp <- lapply(shp, coo_template, size = size.shp.final[i])
    # since coo_template does not center shapes but the bounding box
    shp <- lapply(shp, coo_center)
    # we rotate shapes
    shp <- lapply(shp, coo_rotate, rotate.shp[i])
    # we flip (if required)
    if (flipx.shp[i]) shp <- lapply(shp, coo_flipx)
    if (flipy.shp[i]) shp <- lapply(shp, coo_flipy)
    
    # we translate shapes
    if (plot) { # to translate only for morphospace PCA, not PCcontrib, etc.
      for (s in 1:length(shp)) {
        shp[[s]] <- coo_trans(shp[[s]], pos[s, 1] + dx[i], pos[s, 2] + dy[i])}
    } 
    shp<<-shp
    
    # otherwise, we plot the morphospace
    #===========================================================================================
    #for visualization in ggplot2
        Unzip.list<-function(...) rbind(data.frame(), ...)
        submorphospace <- rep(c(1:nb.shp), times = c(rep(pts.shp+1,nb.shp)))
        shp.frame <- do.call(Unzip.list, shp)
        shp.frame <- data.frame(shp.frame, submorphospace)
        
        
        #call ggplot2 
       
        gg<-ggplot()+
          geom_polygon(aes(x=shp.frame[,1],y=shp.frame[,2],subgroup=submorphospace),
                       fill=col.shp,
                       colour=border.shp,
                       data=shp.frame)+
          theme(legend.text = element_blank())
        
        return(list(gg=gg,shape=shp))
      
}
      #==============================================================================================
      
    
#================================================================================



