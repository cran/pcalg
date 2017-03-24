/*   This program is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*/

#include "pcalg/cpdag.hpp"

//Constructor
cpdag::cpdag(arma::sp_imat adjmat_input, arma::mat covmat_input): 
           adjmat(adjmat_input), covmat(covmat_input) {
   nodes = adjmat.n_rows;
   skeleton = (adjmat + adjmat.t());
   unique = (adjmat - adjmat.t());
   for(arma::uword i=0; i<nodes;i++){
      for(arma::uword j=0; j<nodes;j++){
         if(skeleton(i,j)>1) skeleton(i,j)=1;
         if(unique(i,j)<0) unique(i,j)=0;
      }
   }
   ambiguous = adjmat - unique ;
   // Sparse matrices are sorted by columns, so we 
   // have to take the counterintuitive convention
   // of working with the transpose matrices.
   unique_t = unique.t();
}



// idaFast as it is in the R function
arma::mat  cpdag::idaFast(uint32_t x ){
   //std::cout << "idaFast: computation of causal effects.\n" ; 
   // Las sparse matrices estan ordenada por columnas
   // asi que hemos tomado la traspuesta arriba. 
   arma::mat beta;
      arma::irowvec uniq_pa, amb_pa; 
      uniq_pa.set_size(nodes);
      amb_pa.set_size(nodes);
      arma::uword p = 0;
      for (uint32_t col=0; col<nodes; col++){
         if (ambiguous(x,col) == 1){
            amb_pa(p) = col;
            p++;
         }
      }
      amb_pa.resize(p);
      p = 0;
      for (uint32_t col=0; col<nodes; col++){
         if (unique(x,col) == 1){
            uniq_pa(p) = col;
            p++;
         }
      }
      uniq_pa.resize(p);
      beta.set_size(0,nodes);
      arma::uword cbeta=0;
      // if no ambiguous parents, construct the effects with the unique parents. 
      // No collider checks needed.
      if ( amb_pa.size() == 0){
         beta.insert_rows(cbeta++, getCausal(uniq_pa,amb_pa,x));
      } 
      else{
         if(amb_pa.size() > 64){
            printf ("The number of ambiguous parents is too high");
            exit (EXIT_FAILURE);
         }
        arma::irowvec amb_pa_1, amb_pa_2;
        arma::uword n = amb_pa.size();
        uint64_t  L  = (1 << (n-1)); 
        for (uint64_t bits = 0; bits < L; bits++){
            directArrows(amb_pa_1,amb_pa_2,amb_pa,bits,n);
            if( !hasNewCollider(uniq_pa,amb_pa_1,amb_pa_2,x) ){
                beta.insert_rows(cbeta++,getCausal(uniq_pa,amb_pa_1,x));
            }
            if( !hasNewCollider(uniq_pa,amb_pa_2,amb_pa_1,x) ){
                beta.insert_rows(cbeta++,getCausal(uniq_pa,amb_pa_2,x));
            }
        }
      }
      return (beta.t());
}


bool cpdag::hasNewCollider(arma::irowvec uniq_pa, arma::irowvec amb_pa_t, 
                         arma::irowvec amb_pa_f, arma::uword x){
   arma::irowvec::iterator uniqst, uniqen, amb_tst,amb_ten, amb_fst, amb_fen;
   uniqst = uniq_pa.begin();   uniqen = uniq_pa.end();
   amb_tst = amb_pa_t.begin(); amb_ten = amb_pa_t.end(); 
   amb_fst = amb_pa_f.begin(); amb_fen = amb_pa_f.end();
   if ( amb_pa_t.size() > 0 ){
      // Check whether all amb_pa_t and uniq_pa are connected. If not --> collider.
      if ( uniq_pa.size() > 0 ){
         for (arma::irowvec::iterator it1 = uniqst; it1 < uniqen; it1++){
           for (arma::irowvec::iterator it2=amb_tst; it2 < amb_ten; it2++){
              if ( skeleton((*it1),(*it2)) == 0 ) return true;
           }
         }
      }
      // Check whether all amb_pa_t are connected among each other. If not --> collider.
      if ( amb_pa_t.size() > 1 ){
        for ( arma::irowvec::iterator it1=amb_tst; it1 < amb_ten; it1++ ){ 
          for ( arma::irowvec::iterator it2=it1+1; it2 < amb_ten; it2++ ){ 
              if ( skeleton((*it1),(*it2)) == 0 ) return true;
          }
        }
      }
   }
   // if any of the uniq parents of amb_pa_f is NOT directly connected to x --> collider.
   if ( amb_pa_f.size() > 0 ){
      arma::irowvec papa;
      papa.set_size(nodes);
      arma::uword p=0;
      for (arma::uword i = 0; i < nodes; i++ ){
         if( i != x ){
            for (arma::irowvec::iterator it=amb_fst; it < amb_fen; it++){
               if( unique((*it),i) ){
                  papa[p++] = i; break;
               }
            }
         }
      }
      papa.resize(p);
      if  (p > 0){
         for (arma::irowvec::iterator it=papa.begin(); it < papa.end(); it++){
            if( skeleton(x,(*it)) == 0  ) return true;
         } 
      }
   }
   return false;
}

arma::rowvec cpdag::getCausal(arma::irowvec uniq_pa, arma::irowvec amb_pa_t, arma::uword x){
    arma::uword n_uniq, n_amb, n_all; 
    n_uniq = uniq_pa.size(); 
    n_amb = amb_pa_t.size(); 
    n_all = n_uniq + n_amb + 1; 
    arma::irowvec all_pa; 
    all_pa.set_size(n_all);
    all_pa(0) = x;
    if ( n_uniq > 0)  memcpy(all_pa.begin()+ 1, uniq_pa.begin(),n_uniq*sizeof(long));
    if ( n_amb > 0)  memcpy(all_pa.begin()+ n_uniq + 1, amb_pa_t.begin(),n_amb*sizeof(long)); 
    arma::mat A,Y,M;
    A.set_size(n_all,n_all); 
    Y.set_size(n_all,nodes);
    arma::irowvec::iterator start = all_pa.begin();
    arma::irowvec::iterator end = all_pa.end();
    for (arma::irowvec::iterator it1 = start ; it1 < end; it1++){
       arma::uword row = it1-start;
       Y.row(row) = covmat.row((*it1));
       for (arma::irowvec::iterator it2 = start ; it2 < end; it2++){
          arma::uword col = it2 - start;
          A(row,col) = covmat((*it1),(*it2));
       }
    }
    M = solve(A,Y);
    for (arma::irowvec::iterator it1 = start+1 ; it1 < end; it1++){
       M(0,(*it1)) = 0.0;
    }
    return M.row(0);
}


// Divide the vector amb_pa into two vectors, following the partition k
void cpdag::directArrows(arma::irowvec &pa_1, arma::irowvec &pa_2, 
                         arma::irowvec amb_pa, uint64_t bits, arma::uword n){
      uint8_t  i, idx1, idx2;
      uint8_t  k; 
      idx1 = 0; idx2=0;
      k = __builtin_popcount(bits);
      pa_1.resize(k); 
      pa_2.resize(n-k); 
      for (i=0; i < n ; i++ ){
         if ( bits  &  (1 << i) ){
            pa_1(idx1++) = amb_pa(i);
         }  else {
            pa_2(idx2++) = amb_pa(i);
         }  
      }
}

