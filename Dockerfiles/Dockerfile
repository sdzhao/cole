FROM rocker/verse
RUN R -e "install.packages('Rcpp', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('RcppArmadillo', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('glmnet', repos = 'http://cran.us.r-project.org')"
ADD . .
RUN R -e "install.packages('cole_1.1.0.tar.gz')"

