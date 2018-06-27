
cppFunction( "
  NumericVector density_points(const NumericMatrix& x, double sigma=0.2 ){

    int n=x.nrow();
    NumericVector h(n);
    for(int i = 0; i<n; i++){
      double s = 0;
      R_CheckUserInterrupt();
      for (int j = 0; j<n; j++) {
        double dx = x(i,0) - x(j,0);
        double dy = x(i,1) - x(j,1);
        double d2 = ( dx*dx + dy*dy ) / (2 * sigma*sigma );
        if( d2 < 30 )
           s += exp( -d2 );
      }
      h[i] = 1/(sqrt(2*sigma*sigma)) * s;
    }
    return h;
  }"
)

pept_norm_pool=peptides %>%
  select(starts_with("Reporter intensity corrected")) %>%
  select(contains("TMT")) %>%
  rename_all(str_replace,"Reporter intensity corrected ", "")
pept_norm_pool[pept_norm_pool==0]=NA
colsA <- (2:ncol(pept_norm_pool))[(-1+(2:ncol(pept_norm_pool)))%%10!=0]
colsB <- seq(1, ncol(pept_norm_pool), by= 10)
pept_norm_pool[,colsA] <- pept_norm_pool[,colsA] / pept_norm_pool[,rep(colsB,each=9)]
pept_norm_pool <- pept_norm_pool %>%
  select(-contains("0")) %>%
  as.matrix()

mf <- assignment_mice[assignment_mice$sample!="Pool",]
row.names(mf) <- mf$pept_tbl_col
pept_msn <- MSnSet(pept_norm_pool, as.data.frame(peptides[, c("index",
                                                              "Sequence","Proteins","Leading razor protein","Gene names","Protein names")]), 
                   as.data.frame(mf))
pept_msn <- filterNA(pept_msn, pNA = 0)
pept_msn <- normalise(pept_msn, "quantiles.robust")

pept_msn_norm <- log2(exprs(pept_msn)) %>%
  as.tibble()
colnames(pept_msn_norm) <- sprintf( "%02d %02d %05s",  pData(pept_msn)$weeksWD, pData(pept_msn)$age, pData(pept_msn)$pept_tbl_col )

samples <- c("3 TMT2", "5 TMT7", "2 TMT5", "4 TMT2", "2 TMT9", "1 TMT5")
idx <- match(samples, mf$pept_tbl_col)

mf$processing <- rep("01/24",length(mf$tag_idx))
mf$processing[idx] <- "01/10"
mf$processing[mf$age==34] <- "02/12"
mf$processing[mf$age==32] <- "02/08"
mf$processing[mf$age==30] <- "02/08"
mf$processing[mf$age==4] <- "01/22"
mf$processing[mf$age==6] <- "01/22"
mf$processing[mf$age==8] <- "01/22"
mf$processing[mf$age==10] <- "01/22"
mf$processing[mf$age==12] <- "01/22"
mf$processing[mf$age==14 & mf$weeksWD==0] <- "01/22"
mf$processing[mf$age==16 & mf$weeksWD==0] <- "01/22"
mf$processing[mf$age==18 & mf$weeksWD==0] <- "01/22"
mf$processing[mf$age==20 & mf$weeksWD==0] <- "01/22"
mf$processing[mf$age==22 & mf$weeksWD==0] <- "01/22"
mf$processing[mf$age==12 & mf$weeksWD==8] <- "01/22"
mf$processing[mf$age==16 & mf$weeksWD==12] <- "01/22"
mf$processing[mf$age==20 & mf$weeksWD==16] <- "01/22"

fit <- lmFit( pept_msn_norm, model.matrix( ~ processing + age + factor(weeksWD), mf ))

coeffs<-fit$coefficients[,7:11]
colnames(coeffs) <- c("WD 08", "WD 12", "WD 16", "WD 20", "WD 26")

data <- 
  map_dfr( 1:(ncol(coeffs)-1), function(i)
    map_dfr( (i+1):ncol(coeffs), function(j)
      tibble( xCol = colnames(coeffs)[i], yCol = colnames(coeffs)[j], x=coeffs[,i], y=coeffs[,j] ) ) )
data <-  mutate( data, xCol=str_replace( data$xCol,"factor\\(weeksWD\\)","WD ") ) %>%
  mutate( yCol=str_replace( data$yCol,"factor\\(weeksWD\\)","WD ") ) %>%
  mutate( xCol=str_replace( data$xCol,"8","08") ) %>%
  mutate( yCol=str_replace( data$yCol,"8","08") ) %>%
  group_by( xCol, yCol ) %>%
  mutate( dens = density_points( cbind( x, y ), 0.1) )

ggplot(data) +
  geom_point( aes( x=x, y=y,col=dens ), size=0.5 ) +
  scale_color_gradientn( colours=cubeHelix(11, start = 0.5, r = -1.5 )[9:1]) +
  geom_abline( alpha=.3 ) + geom_hline(yintercept=0,alpha=.3) + geom_vline(xintercept = 0,alpha=.3) +
  facet_grid( xCol ~ yCol )
