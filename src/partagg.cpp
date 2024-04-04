#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List partagg(List grpl1, DataFrame corrbin, int B) {
  CharacterVector v1 = corrbin["Var1"];
  CharacterVector v2 = corrbin["Var2"];
  CharacterVector ftmp = {""}, ftmp2 = {"", ""};
  List grpl = clone(grpl1);
  int nsize = v1.size();
  int msize = grpl.size();  
  int ind1, ind2, lsize, n1, n2;

  Rcout << "nrow(corrfile): " << nsize << std::endl;
  Rcout << "n starting modules: " << msize << std::endl;

  for(int i = 0; i < nsize; i++){
      ind1 = -1;
      ind2 = -1;
      for(int j = 0; j < msize; j++){
        Rcpp::CharacterVector ftmp = grpl[j];
        lsize = ftmp.size();
        for(int k = 0; k < lsize; k++){
          if(v1[i] == ftmp[k]){
            ind1 = j;
          }
          if(v2[i] == ftmp[k]){
            ind2 = j;
          }
          if(ind1 > -1 && ind2 > -1){
            break;
          }
        } // for k
        if(ind1 > -1 && ind2 > -1){
          break;
        }
      } // for j
            
      if(ind1 < 0 && ind2 < 0){
        ftmp = {v1[i], v2[i]};
        grpl.push_back(ftmp);
        msize = grpl.size();
      }
      if(ind1 > -1 && ind2 < 0){
        ftmp = grpl[ind1];
        ftmp.push_back(v2[i]);
        if(ftmp.size() <= B){
          grpl[ind1] = ftmp;
        }
      }
      if(ind1 < 0 && ind2 > -1){
        ftmp = grpl[ind2];
        ftmp.push_back(v1[i]);
        if(ftmp.size() <= B){
      	  grpl[ind2] = ftmp;
        }
      }
      if(ind1 > -1 && ind2 > -1){
        if(ind1 != ind2){
          Rcpp::CharacterVector ftmp1 = grpl[ind1];
          Rcpp::CharacterVector ftmp2 = grpl[ind2];
          n1 = ftmp1.size();
          n2 = ftmp2.size();
          if((n1 + n2) <= B){
        	  Rcpp::CharacterVector ftmp3(n1 + n2);
            for(int j = 0; j < n1; j++){
              ftmp3[j] = ftmp1[j];
            } // for j
            for(int j = 0; j < n2; j++){
              ftmp3[n1 + j] = ftmp2[j];
            } // for j
            grpl[ind1] = ftmp3;
            grpl.erase(ind2);
            msize = grpl.size();
          }
        }
      }
    } // for i
  
  Rcout << "n ending modules: " << msize << std::endl;        
  return(grpl);
}
