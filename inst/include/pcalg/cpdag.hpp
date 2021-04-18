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

#ifndef DIFFEXP_H
#define DIFFEXP_H

#include <armadillo>
#include <iostream>
#include <stdint.h>


class cpdag{
   private: 
      arma::uword nodes;
      arma::mat covmat;           // Covariance matrix 
      arma::sp_imat adjmat;       // Adjacency matrix CPDAG
      arma::sp_imat skeleton;     // Undirected adjacency matrix
      arma::sp_imat ambiguous;    // Ambiguous parents
      arma::sp_imat unique;       // Unique parents
      arma::sp_imat unique_t;     // Unique parents transposed (used for computations)
   public: 
      cpdag(arma::sp_imat adjmat_input, arma::mat covmat_input); 
      ~cpdag() {}
      arma::mat idaFast(uint32_t x);
      bool hasNewCollider(arma::irowvec uniq_pa, arma::irowvec amb_pa_t, 
                          arma::irowvec amb_pa_f, arma::uword x);
      arma::rowvec getCausal(arma::irowvec uniq_pa, arma::irowvec amb_pa_t, arma::uword x);
      void directArrows(arma::irowvec &pa_1, arma::irowvec &pa_2, 
                        arma::irowvec amb_pa, uint64_t k, arma::uword n);

};
#endif //cpdag.h 

