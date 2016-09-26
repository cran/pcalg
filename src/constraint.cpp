/**
 * $Id: $
 */

#ifndef CONSTRAINT_HPP_
#define CONSTRAINT_HPP_

#include "pcalg/constraint.hpp"
#include "pcalg/gies_debug.hpp"

#include <algorithm>
#include <utility>
#include <iterator>
#include <limits>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/dynamic_bitset.hpp>

double IndepTestRFunction::test(uint u, uint v, std::vector<uint> S) const
{
	// Adapt indices to R convention
	std::vector<uint> shiftS;
	shiftS.reserve(S.size());
	std::vector<uint>::iterator vi;
	for (vi = S.begin(); vi != S.end(); ++vi)
		shiftS.push_back(*vi + 1);

	// Call R function to perform independence test
	return Rcpp::as<double>(_testFunction(u + 1, v + 1, shiftS, *_suffStat));
}

double IndepTestGauss::test(uint u, uint v, std::vector<uint> S) const
{
	// Return NaN if any of the correlation coefficients needed for calculation is NaN
	arma::mat C_sub;
	arma::uvec ind(S.size() + 2);
	ind(0) = u;
	ind(1) = v;
	uint i, j;
	for (i = 0; i < S.size(); ++i) ind(i + 2) = S[i];
	C_sub = _correlation.submat(ind, ind);
	for (i = 0; i < C_sub.n_rows; ++i)
		for (j = 0; j < C_sub.n_cols; ++j)
			if ((boost::math::isnan)(C_sub(i, j)))
				return std::numeric_limits<double>::quiet_NaN();

	// Calculate (absolute value of) z statistic
	#define CUT_THR 0.9999999
	double r, absz;
	//dout.level(3) << " Performing independence test for conditioning set of size " << S.size() << std::endl;
	if (S.empty())
		r = _correlation(u, v);
	else if (S.size() == 1)
		r = (C_sub(0, 1) - C_sub(0, 2) * C_sub(1, 2))/sqrt((1 - C_sub(1, 2)*C_sub(1, 2)) * (1 - C_sub(0, 2)*C_sub(0, 2)));
	else {
		arma::mat PM;
		pinv(PM, C_sub);
		// TODO include error handling
		r = - PM(0, 1)/sqrt(PM(0, 0) * PM(1, 1));
	}
	// Absolute value of r, respect cut threshold
	r = std::min(CUT_THR, std::abs(r));

	// Absolute value of z statistic
	// Note: log1p for more numerical stability, see "Aaux.R"; log1p is also available in
	// header <cmath>, but probably only on quite up to date headers (C++11)?
	absz = sqrt(_sampleSize - S.size() - 3.0) * 0.5 * boost::math::log1p(2*r/(1 - r));

	// Calculate p-value to z statistic (based on standard normal distribution)
	boost::math::normal distN;
	return (2*boost::math::cdf(boost::math::complement(distN, absz)));
}

void Skeleton::addFixedEdge(const uint a, const uint b)
{
	boost::add_edge(a, b, _fixedEdges);
	addEdge(a, b);
}

bool Skeleton::isFixed(const uint a, const uint b) const
{
	bool result;
	boost::tie(boost::tuples::ignore, result) = boost::edge(a, b, _fixedEdges);
	return result;
}

void Skeleton::removeEdge(const uint a, const uint b)
{
	if (!isFixed(a, b))
		boost::remove_edge(a, b, _graph);
}

bool Skeleton::hasEdge(const uint a, const uint b) const
{
	bool result;
	boost::tie(boost::tuples::ignore, result) = boost::edge(a, b, _graph);
	return result;
}


std::set<uint> Skeleton::getNeighbors(const uint vertex) const
{
	std::set<uint> result;
	UndirOutEdgeIter outIter, outLast;

	for (boost::tie(outIter, outLast) = boost::out_edges(vertex, _graph); outIter != outLast; outIter++)
			result.insert(boost::target(*outIter, _graph));

	return result;
}

Rcpp::LogicalMatrix Skeleton::getAdjacencyMatrix()
{
	Rcpp::LogicalMatrix result(getVertexCount(), getVertexCount());
	UndirEdgeIter ei, eiLast;
	for (boost::tie(ei, eiLast) = boost::edges(_graph); ei != eiLast; ei++) {
		dout.level(3) << "  Edge {" << boost::source(*ei, _graph) <<
				", " << boost::target(*ei, _graph) << "}\n";
		result(boost::source(*ei, _graph), boost::target(*ei, _graph)) = true;
		result(boost::target(*ei, _graph), boost::source(*ei, _graph)) = true;
	}

	return result;
}

void Skeleton::fitCondInd(
		const double alpha,
		Rcpp::NumericMatrix& pMax,
		SepSets& sepSet,
		std::vector<int>& edgeTests,
		int maxCondSize,
		const bool NAdelete) {
	if (maxCondSize < 0)
		maxCondSize = getVertexCount();

	dout.level(2) << "Significance level " << alpha << std::endl;
	dout.level(2) << "Maximum order: " << maxCondSize << std::endl;
	bool found = true;

	UndirEdgeIter ei, eiLast;

	// edgeTests lists the number of edge tests that have already been done; its size
	// corresponds to the size of conditioning sets that have already been checked
	// TODO: improve handling of check_interrupt, see e.g.
	// https://github.com/jjallaire/Rcpp/blob/master/inst/examples/OpenMP/piWithInterrupts.cpp
	for (uint condSize = edgeTests.size();
			!check_interrupt() && found && (int)condSize <= maxCondSize;
			++condSize) {
		dout.level(1) << "Order = " << condSize << "; remaining edges: " << getEdgeCount() << std::endl;

		// Make a list of edges in the graph; this is needed for OpenMP
		std::vector<uint> u, v;
		u.reserve(getEdgeCount());
		v.reserve(getEdgeCount());
		for (boost::tie(ei, eiLast) = boost::edges(_graph); ei != eiLast && !check_interrupt(); ei++) {
			uint node1 = boost::source(*ei, _graph);
			uint node2 = boost::target(*ei, _graph);
			if (node1 > node2)
				std::swap(node1, node2);
			if (std::max(getDegree(node1), getDegree(node2)) > condSize && !isFixed(node1, node2)) {
				u.push_back(node1);
				v.push_back(node2);
			}
		}
		boost::dynamic_bitset<> deleteEdges(u.size());
		arma::ivec localEdgeTests(u.size(), arma::fill::zeros);

		// There is a conditioning set of size "condSize" if u is not empty
		found = u.size() > 0;

		// Iterate over all edges in the graph
		#pragma omp parallel for
		for (std::size_t l = 0; l < u.size(); l++) {
			bool edgeDone = false;

			int k;
			UndirOutEdgeIter outIter, outLast;
			std::vector<uint> condSet(condSize);
			std::vector<std::vector<uint>::iterator> si(condSize);

			// Check neighborhood of u
			if (getDegree(u[l]) > condSize) {
				// Get neighbors of u (except v)
				std::vector<uint> neighbors(0);
				neighbors.reserve(getDegree(u[l]) - 1);
				for (boost::tie(outIter, outLast) = boost::out_edges(u[l], _graph); outIter != outLast; outIter++)
					if (boost::target(*outIter, _graph) != v[l])
						neighbors.push_back(boost::target(*outIter, _graph));

				// Initialize first conditioning set
				for (std::size_t i = 0; i < condSize; ++i)
					si[i] = neighbors.begin() + i;

				// Iterate over conditioning sets
				do {
					for (std::size_t i = 0; i < condSize; ++i)
						condSet[i] = *(si[i]);

					// Test of u and v are conditionally independent given condSet
					double pval = _indepTest->test(u[l], v[l], condSet);
					localEdgeTests(l)++;
					dout.level(1) << "  x = " << u[l] << ", y = " << v[l] << ", S = " <<
							condSet << " : pval = " << pval << std::endl;
					if ((boost::math::isnan)(pval))
						pval = (NAdelete ? 1. : 0.);
					if (pval > pMax(u[l], v[l]))
						pMax(u[l], v[l]) = pval;
					if (pval >= alpha) {
						deleteEdges.set(l);
						// arma::ivec condSetR(condSet.size());
						sepSet[v[l]][u[l]].set_size(condSet.size());
						for (std::size_t j = 0; j < condSet.size(); ++j)
							sepSet[v[l]][u[l]][j] = condSet[j] + 1;
						edgeDone = true;
						break; // Leave do-while-loop
					}

					// Proceed to next conditioning set
					for (k = condSize - 1;
							k >= 0 && si[k] == neighbors.begin() + (neighbors.size() - condSize + k);
							--k);
					if (k >= 0) {
						si[k]++;
						for (k++; k < (int)condSize; ++k)
							si[k] = si[k - 1] + 1;
					}
				} while(k >= 0);
			} // IF getDegree(u[l])

			// Check neighborhood of v
			if (!edgeDone && getDegree(v[l]) > condSize) {
				// Get neighbors of u (except v); common neighbors of u and v are listed in the end
				std::vector<uint> neighbors(0);
				std::vector<uint> commNeighbors(0);
				neighbors.reserve(getDegree(v[l]) - 1);
				commNeighbors.reserve(getDegree(v[l]) - 1);
				uint a;
				for (boost::tie(outIter, outLast) = boost::out_edges(v[l], _graph);
						outIter != outLast; outIter++) {
					a = boost::target(*outIter, _graph);
					if (a != u[l]) {
						if (hasEdge(u[l], a))
							commNeighbors.push_back(a);
						else
							neighbors.push_back(a);
					}
				}

				// m: number of neighbors of v that are not neighbors of u
				uint m = neighbors.size();
				neighbors.insert(neighbors.end(), commNeighbors.begin(), commNeighbors.end());
				dout.level(2) << "  v: " << v << "; neighbors: " << neighbors << " (m = " << m << ")\n";

				// If all neighbors of v are also adjacent to u: already checked all conditioning sets
				if (m > 0) {
					// Initialize first conditioning set
					for (std::size_t i = 0; i < condSize; ++i)
						si[i] = neighbors.begin() + i;

					// Iterate over conditioning sets
					do {
						for (std::size_t i = 0; i < condSize; ++i)
							condSet[i] = *(si[i]);

						// Test of u and v are conditionally independent given condSet
						double pval = _indepTest->test(v[l], u[l], condSet);
						localEdgeTests(l)++;
						dout.level(1) << "  x = " << v[l] << ", y = " << u[l] << ", S = " <<
								condSet << " : pval = " << pval << std::endl;
						if ((boost::math::isnan)(pval))
							pval = (NAdelete ? 1. : 0.);
						if (pval > pMax(u[l], v[l]))
							pMax(u[l], v[l]) = pval;
						if (pval >= alpha) {
							deleteEdges.set(l);
							// arma::ivec condSetR(condSet.size());
							sepSet[v[l]][u[l]].set_size(condSet.size());
							for (std::size_t j = 0; j < condSet.size(); ++j)
								sepSet[v[l]][u[l]][j] = condSet[j] + 1;
							edgeDone = true;
							break; // Leave do-while-loop
						}

						// Proceed to next conditioning set
						for (k = condSize - 1;
								k >= 0 && si[k] == neighbors.begin() + (neighbors.size() - condSize + k);
								--k);
						// Make sure first element does not belong to neighborhood of u: otherwise
						// we would redo a test already performed
						if (k == 0 && si[0] == neighbors.begin() + (m - 1))
							k = -1;
						if (k >= 0) {
							si[k]++;
							for (k++; k < (int)condSize; ++k)
								si[k] = si[k - 1] + 1;
						}
					} while(k >= 0);
				} // IF m
			} // IF getDegree(v[l])
		} // FOR l

		// Delete edges marked for deletion
		for (std::size_t l = deleteEdges.find_first(); l < deleteEdges.size(); l = deleteEdges.find_next(l))
			removeEdge(u[l], v[l]);

		// Calculate total number of edge tests
		if (found)
			edgeTests.push_back(arma::accu(localEdgeTests));
	} // FOR condSize
}

#endif /* CONSTRAINT_HPP_ */
