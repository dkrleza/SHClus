#include "SigmaIndex_R.cpp"
#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;

RCPP_EXPOSED_AS(SigmaIndex_R)
RCPP_MODULE(SigmaIndexModule) {
  Rcpp::class_<SigmaIndex_R>("SigmaIndex_R")
  .constructor()
  .constructor<double,double,bool>()
  .method("addPopulation",&SigmaIndex_R::addPopulation,"Add a new population to the sigma-index")
  .method("addDataPoint",&SigmaIndex_R::addDataPoint,"Add a data point to an existing population in the sigma-index")
  .method("queryDataPoint",&SigmaIndex_R::queryDataPoint,"Query a data point in the sigma-index")
  .method("getTotalPopulationNumber",&SigmaIndex_R::getTotalPopulationNumber,"Get total population number")
  .method("print",&SigmaIndex_R::print,"Printout the sigma-index")
  .method("getPopulations",&SigmaIndex_R::getPopulations,"Get all populations in the sigma-index")
  .method("removePopulation",&SigmaIndex_R::removePopulation,"Remove population")
  .method("resetStatistics",&SigmaIndex_R::resetStatistics,"Reset statistics")
  .method("getHistogram",&SigmaIndex_R::getHistogram,"Calculate and return the sigma-index histogram")
  .method("getStatistics",&SigmaIndex_R::getStatistics,"Get sigma-index collected statistics");
}