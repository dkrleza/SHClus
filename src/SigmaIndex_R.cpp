#include <SigmaIndex>
#include <SigmaIndexProxy>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <sstream>
using namespace Rcpp;
using namespace std;
using namespace Eigen;

class SigmaIndex_R : public SigmaIndex<ManualPopulation*> {
public:
  SigmaIndex_R(double theta, double neighborhood_theta, bool precision_switch=true) : 
    SigmaIndex<ManualPopulation*>(theta, neighborhood_theta, precision_switch) {
  }
  SigmaIndex_R() : SigmaIndex<ManualPopulation*>(3, 3, false) {
  }
  
  void addPopulation(string id, NumericVector mean, NumericMatrix covariance, long elements) {
    VectorXd *mean_eigen=new VectorXd(as<VectorXd>(mean));
    MatrixXd *covariance_eigen=new MatrixXd(as<MatrixXd>(covariance));
    ManualPopulation *pop=new ManualPopulation(id, mean_eigen, covariance_eigen, elements,
                                               [&](string id){completeUpdate(nodes[id]->getPopulation());});
    create(pop, true);
    //completeUpdate(pop); // no need, this is done in create
    delete mean_eigen;delete covariance_eigen;
  }
  
  void removePopulation(string id) {
    if(nodes.find(id)!=nodes.end()) {
      ManualPopulation *pop=nodes[id]->getPopulation();
      remove(pop);
    } else {
      stringstream ss;
      ss << "Cannot find the population with id=" << id;
      throw AbstractPopulation_Exception(ss.str().c_str());
    }
  }
  
  void addDataPoint(string id, NumericVector data_point, List neighborhood) {
    VectorXd *data_point_eigen=new VectorXd(as<VectorXd>(data_point));
    if(nodes.find(id)!=nodes.end()) {
      ManualPopulation *pop=nodes[id]->getPopulation();
      if(neighborhood.size()==0) pop->addNewElement(data_point_eigen);
      else {
        unordered_map<ManualPopulation*,double> *n_map=new unordered_map<ManualPopulation*,double>();
        for(int i=0;i<neighborhood.size();++i) {
          List l=neighborhood[i];
          String nid=l["id"];
          if(nodes.find(nid)!=nodes.end()) {
            ManualPopulation *pop_n=nodes[nid]->getPopulation();
            (*n_map)[pop_n]=(double)l["distance"];
          }
        }
        pop->addNewElement(data_point_eigen,false);
        update(pop,n_map);
      }
    } else {
      delete data_point_eigen;
      stringstream ss;
      ss << "Cannot find the population with id=" << id;
      throw AbstractPopulation_Exception(ss.str().c_str());
    }
    delete data_point_eigen;
  }
  
  List queryDataPoint(NumericVector data_point) {
    VectorXd *data_point_eigen=new VectorXd(as<VectorXd>(data_point));
    SigmaIndexQueryResults<ManualPopulation*> *q_res=query(data_point_eigen);
    List res=List::create();
    List classified=List::create();
    for(pair<ManualPopulation*,double> it:*q_res->classified) 
      classified[it.first->getId()]=it.second;
    res["outcome"]=classified;
    List neigh=List::create();
    for(pair<ManualPopulation*,double> it:*q_res->neighborhood) 
      neigh[it.first->getId()]=it.second;
    res["neighborhood"]=neigh;
    delete data_point_eigen;delete q_res;
    return res;
  }
  
  int getTotalPopulations() {
    return nodes.size();
  }
  
  void print() {
    set<string> visited;
    SigmaIndex<ManualPopulation*>::print(ROOT, Rcout, &visited);
  }
  
  List getPopulations() {
    List res=List::create();
    for(pair<string,SigmaIndexNode<ManualPopulation*>*> it:nodes) {
      if(it.first!=ROOT) {
        List node=List::create();
        VectorXd *mean=it.second->getPopulation()->getMean();
        if(mean) node["mean"]=wrap(*mean);
        MatrixXd *icov=it.second->getPopulation()->getICovariance();
        if(icov) node["icovariance"]=wrap(*icov);
        node["elements"]=it.second->getPopulation()->getElements();
        res[it.first]=node;
      }
    }
    return res;
  }
  
  void resetStatistics() {
    SigmaIndex<ManualPopulation*>::resetStatistics();
  }
  
  NumericMatrix getHistogram() {
    unordered_map<int,long> *hist=calculateHistogram();
    IntegerVector xv,yv;
    for(int it=0;it<=100;it++) {
      xv.push_back(it);
      if(hist->find(it)!=hist->end()) yv.push_back((*hist)[it]);
      else yv.push_back(0);
    }
    NumericMatrix res(1,101,yv.begin());
    colnames(res)=xv;
    delete hist;
    return res;
  }
  
  List getStatistics() {
    SigmaIndexStatistics *stats=SigmaIndex<ManualPopulation*>::getStatistics();
    List res=List::create();
    res["totalCount"]=stats->totalCount;
    res["classifiedNodes"]=stats->tsCount;
    res["missedNodes"]=stats->tpCount;
    res["sequentialNodes"]=stats->totalSequential;
    res["computationCostReduction"]=stats->R;
    delete stats;
    return res;
  }
};

