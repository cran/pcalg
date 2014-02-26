/**
 * $Id: $
 */

#ifndef CONSTRAINT_HPP_
#define CONSTRAINT_HPP_

#include "constraint.hpp"
#include "gies_debug.hpp"

#include <algorithm>
#include <utility>
#include <iterator>
#include <limits>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

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
	absz = sqrt(_sampleSize - S.size() - 3) * 0.5 * boost::math::log1p(2*r/(1 - r));

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
	uint condSize, u, v, a, m, u_min, u_max;
	int i, j;
	double pval;
	bool found = true;
	bool edgeDone;
	UndirOutEdgeIter outIter, outLast;
	std::set<uint>::iterator vi;
	std::vector<std::vector<uint>::iterator> si;
	std::vector<uint> condSet, neighbors, commNeighbors;
	std::vector<uint>::iterator vnext;
	arma::ivec condSetR;

	UndirEdgeIter ei, eiLast;
	std::set< std::pair<uint, uint> > deleteEdges;
	std::set< std::pair<uint, uint> >::iterator di;

	// edgeTests lists the number of edge tests that have already been done; its size
	// corresponds to the size of conditioning sets that have already been checked
	for (condSize = edgeTests.size(); found && condSize <= maxCondSize; ++condSize) {
		dout.level(1) << "Order = " << condSize << "; remaining edges: " << getEdgeCount() << std::endl;
		// oldSkel._graph = _graph;
		condSet.resize(condSize);
		si.resize(condSize);
		deleteEdges.clear();
		found = false;
		edgeTests.push_back(0);

		// Iterate over all edges in the graph
		for (boost::tie(ei, eiLast) = boost::edges(_graph); ei != eiLast && !check_interrupt(); ei++) {
			// Get endpoints u, v of edge; make sure that deg(u) >= deg(v)
			u = boost::source(*ei, _graph);
			v = boost::target(*ei, _graph);
			u_min = std::min(u, v);
			u_max = std::max(u, v);
			if (getDegree(u) < getDegree(v))
				std::swap(u, v);

			// There is a conditioning set of size "condSize" if deg(u) > condSize
			if (getDegree(v) > condSize)
				found = true;
			edgeDone = false;

			// Check neighborhood of u, if edge is not fixed
			if (!isFixed(u, v) && getDegree(u) > condSize) {
				// Get neighbors of u (except v)
				neighbors.clear();
				neighbors.reserve(getDegree(u) - 1);
				for (boost::tie(outIter, outLast) = boost::out_edges(u, _graph); outIter != outLast; outIter++)
					if (boost::target(*outIter, _graph) != v)
						neighbors.push_back(boost::target(*outIter, _graph));

				// Initialize first conditioning set
				for (i = 0; i < condSize; ++i)
					si[i] = neighbors.begin() + i;

				// Iterate over conditioning sets
				do {
					for (i = 0; i < condSize; ++i)
						condSet[i] = *(si[i]);

					// Test of u and v are conditionally independent given condSet
					pval = _indepTest->test(u, v, condSet);
					edgeTests.back()++;
					dout.level(1) << "  x = " << u << ", y = " << v << ", S = " <<
							condSet << " : pval = " << pval << std::endl;
					if ((boost::math::isnan)(pval))
						pval = (NAdelete ? 1. : 0.);
					if (pval > pMax(u_min, u_max))
						pMax(u_min, u_max) = pval;
					if (pval >= alpha) {
						deleteEdges.insert(std::make_pair(u, v));
						condSetR.set_size(condSet.size());
						for (j = 0; j < condSet.size(); ++j)
							condSetR[j] = condSet[j] + 1;
						sepSet[u_max][u_min] = condSetR;
						edgeDone = true;
						break; // Leave do-while-loop
					}

					// Proceed to next conditioning set
					for (i = condSize - 1;
							i >= 0 && si[i] == neighbors.begin() + (neighbors.size() - condSize + i);
							--i);
					if (i >= 0) {
						si[i]++;
						for (i++; i < condSize; ++i)
							si[i] = si[i - 1] + 1;
					}
				} while(i >= 0);
			}

			// Check neighborhood of v, if edge is not fixed
			if (!edgeDone && !isFixed(u, v) && getDegree(v) > condSize) {
				// Get neighbors of u (except v); common neighbors of u and v are listed in the end
				neighbors.clear();
				commNeighbors.clear();
				neighbors.reserve(getDegree(v) - 1);
				commNeighbors.reserve(getDegree(v) - 1);
				for (boost::tie(outIter, outLast) = boost::out_edges(v, _graph); outIter != outLast; outIter++) {
					a = boost::target(*outIter, _graph);
					if (a != u) {
						if (hasEdge(u, a))
							commNeighbors.push_back(a);
						else
							neighbors.push_back(a);
					}
				}

				// m: number of neighbors of v that are not neighbors of u
				m = neighbors.size();
				neighbors.insert(neighbors.end(), commNeighbors.begin(), commNeighbors.end());
				dout.level(2) << "  v: " << v << "; neighbors: " << neighbors << " (m = " << m << ")\n";

				// If all neighbors of v are also adjacent to u: already checked all conditioning sets
				if (m > 0) {
					// Initialize first conditioning set
					for (i = 0; i < condSize; ++i)
						si[i] = neighbors.begin() + i;

					// Iterate over conditioning sets
					do {
						for (i = 0; i < condSize; ++i)
							condSet[i] = *(si[i]);

						// Test of u and v are conditionally independent given condSet
						pval = _indepTest->test(v, u, condSet);
						edgeTests.back()++;
						dout.level(1) << "  x = " << v << ", y = " << u << ", S = " <<
								condSet << " : pval = " << pval << std::endl;
						if ((boost::math::isnan)(pval))
							pval = (NAdelete ? 1. : 0.);
						if (pval > pMax(u_min, u_max))
							pMax(u_min, u_max) = pval;
						if (pval >= alpha) {
							deleteEdges.insert(std::make_pair(u, v));
							condSetR.set_size(condSet.size());
							for (j = 0; j < condSet.size(); ++j)
								condSetR[j] = condSet[j] + 1;
							sepSet[u_max][u_min] = condSetR;
							edgeDone = true;
							break; // Leave do-while-loop
						}

						// Proceed to next conditioning set
						for (i = condSize - 1;
								i >= 0 && si[i] == neighbors.begin() + (neighbors.size() - condSize + i);
								--i);
						// Make sure first element does not belong to neighborhood of u: otherwise
						// we would redo a test already performed
						if (i == 0 && si[0] == neighbors.begin() + (m - 1))
							i = -1;
						if (i >= 0) {
							si[i]++;
							for (i++; i < condSize; ++i)
								si[i] = si[i - 1] + 1;
						}
					} while(i >= 0);
				}
			}
		}

		// Delete edges marked for deletion
		for (di = deleteEdges.begin(); di != deleteEdges.end(); ++di)
			removeEdge(di->first, di->second);
	} // FOR condSize

	// Adjust vector of edge test numbers
	if (edgeTests.back() == 0)
		edgeTests.pop_back();
}

#endif /* CONSTRAINT_HPP_ */
