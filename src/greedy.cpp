/**
 * greedy.cpp
 *
 * @author Alain Hauser
 * $Id: greedy.cpp 500 2019-11-20 13:46:14Z alhauser $
 */

#include "pcalg/greedy.hpp"
#include "pcalg/gies_debug.hpp"

#include <algorithm>
#include <map>
#include <stack>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>

EssentialGraph::EssentialGraph(const uint vertexCount) :
	_graph(vertexCount),
	_fixedGaps(vertexCount),
	_gapsInverted(false),
	_maxVertexDegree(vertexCount, vertexCount),
	_childrenOnly(vertexCount),
	_loggers()
{
	disableCaching();
}

double EssentialGraph::_minScoreDiff = sqrt(std::numeric_limits<double>::epsilon());

void EssentialGraph::clear()
{
	// Clear graph
	boost::graph_traits<InternalEssentialGraph>::edge_iterator ei, ei_end, next;
	boost::tie(ei, ei_end) = boost::edges(_graph);
	for (next = ei; ei != ei_end; ei = next) {
		++next;
		boost::remove_edge(*ei, _graph);
	}
}

void EssentialGraph::addEdge(const uint a, const uint b, bool undirected)
{
	if (!hasEdge(a, b)) {
		// Add edge and log it
		boost::add_edge(a, b, _graph);
		for (std::set<GraphOperationLogger*>::iterator logger = _loggers.begin();
				logger != _loggers.end(); ++logger) {
			(*logger)->log(GO_ADD_EDGE, a, b);
		}
	}

	if (undirected && !hasEdge(b, a)) {
		// Add edge and log it
		boost::add_edge(b, a, _graph);
		for (std::set<GraphOperationLogger*>::iterator logger = _loggers.begin();
				logger != _loggers.end(); ++logger) {
			(*logger)->log(GO_ADD_EDGE, b, a);
		}
	}
}

void EssentialGraph::removeEdge(const uint a, const uint b, bool bothDirections)
{
	if (hasEdge(a, b)) {
		// Remove edge and log it
		boost::remove_edge(a, b, _graph);
		for (std::set<GraphOperationLogger*>::iterator logger = _loggers.begin();
				logger != _loggers.end(); ++logger) {
			(*logger)->log(GO_REMOVE_EDGE, a, b);
		}
	}

	if (bothDirections && hasEdge(b, a)) {
		// Remove edge and log it
		boost::remove_edge(b, a, _graph);
		for (std::set<GraphOperationLogger*>::iterator logger = _loggers.begin();
				logger != _loggers.end(); ++logger) {
			(*logger)->log(GO_REMOVE_EDGE, b, a);
		}
	}
}

bool EssentialGraph::gapFixed(const uint a, const uint b) const
{
	bool result;
	boost::tie(boost::tuples::ignore, result) = boost::edge(a, b, _fixedGaps);
	return result ^ _gapsInverted;
}

void EssentialGraph::setFixedGap(const uint a, const uint b, const bool fixed)
{
	if (fixed ^ _gapsInverted) {
		boost::add_edge(a, b, _fixedGaps);
	} else {
		boost::remove_edge(a, b, _fixedGaps);
	}
}

bool EssentialGraph::existsPath(const uint a, const uint b, const std::set<uint>& C, const bool undirected)
{
	// Mark "forbidden" vertices as visited
	boost::dynamic_bitset<> visited(getVertexCount());
	std::set<uint>::iterator si;
	for (si = C.begin(); si != C.end(); ++si)
		visited.set(*si);

	// Trivial cases: if a or b are in C, return false
	if (visited.test(a) || visited.test(b))
		return false;

	// If (a, b) is an edge in the graph, remove it -- and add it again in the end
	bool restore = hasEdge(a, b);
	if (restore)
		removeEdge(a, b);

	// Check with depth-first search whether b is reachable from a without
	// using vertices in C
	std::stack<uint> nbhd;
	nbhd.push(a);
	visited.set(a);
	uint v;
	AdjacencyIter vi, vi_last;
	while (!nbhd.empty()) {
		v = nbhd.top();
		nbhd.pop();
		for (boost::tie(vi, vi_last) = boost::adjacent_vertices(v, _graph); vi != vi_last; ++vi)
			if (!undirected || hasEdge(*vi, v)) {
				if (*vi == b) {
					if (restore)
						addEdge(a, b);
					return true;
				}
				if (!visited.test(*vi)) {
					nbhd.push(*vi);
					visited.set(*vi);
				}
			}
	}

	if (restore)
		addEdge(a, b);
	return false;
}

bool EssentialGraph::existsPath(const uint a, const std::set<uint>& B, const std::set<uint>& C, const bool undirected) const
{
	// Exclude trivial cases
	if (B.empty() || std::includes(C.begin(), C.end(), B.begin(), B.end()) || C.find(a) != C.end())
		return false;

	// Mark "forbidden" vertices as visited
	boost::dynamic_bitset<> visited(getVertexCount());
	std::set<uint>::iterator si;
	for (si = C.begin(); si != C.end(); ++si)
		visited.set(*si);

	// Check with depth-first search whether any vertex in B is reachable from a without
	// using vertices in C
	std::stack<uint> nbhd;
	nbhd.push(a);
	visited.set(a);
	uint v;
	AdjacencyIter vi, vi_last;
	while (!nbhd.empty()) {
		v = nbhd.top();
		nbhd.pop();
		for (boost::tie(vi, vi_last) = boost::adjacent_vertices(v, _graph); vi != vi_last; ++vi)
			if (!undirected || hasEdge(*vi, v)) {
				if (B.find(*vi) != B.end())
					return true;
				if (!visited.test(*vi)) {
					nbhd.push(*vi);
					visited.set(*vi);
				}
			}
	}

	return false;
}

bool EssentialGraph::existsPath(const std::set<uint>& C, const uint a, const std::set<uint>& B)
{
	// Copy set of allowed vertices to bitset for faster lookup
	boost::dynamic_bitset<> allowed(getVertexCount());
	std::set<uint>::iterator si;
	for (si = C.begin(); si != C.end(); ++si)
		allowed.set(*si);

	// Exclude trivial cases
	std::set<uint> T = set_intersection(B, C);
	if (T.empty() || !allowed.test(a))
		return false;

	// Check with depth-first search whether any vertex in B is reachable from a
	// by only using vertices in C
	std::stack<uint> nbhd;
	boost::dynamic_bitset<> visited(getVertexCount());
	nbhd.push(a);
	visited.set(a);
	uint v;
	AdjacencyIter vi, vi_last;
	while (!nbhd.empty()) {
		v = nbhd.top();
		nbhd.pop();
		for (boost::tie(vi, vi_last) = boost::adjacent_vertices(v, _graph); vi != vi_last; ++vi)
			if (allowed.test(*vi)) {
				if (T.find(*vi) != T.end())
					return true;
				if (!visited.test(*vi)) {
					nbhd.push(*vi);
					visited.set(*vi);
				}
			}
	}

	return false;
}

std::vector<uint> EssentialGraph::greedyColoring(std::vector<uint> vertices)
{
	// Initialize coloring: vector of same length as list of vertices,
	// color of first vertex is 0
	std::vector<uint> coloring(vertices.size());
	boost::dynamic_bitset<> available;
	std::set<uint> adjacent;

	for (std::size_t i = 1; i < vertices.size(); ++i){
		// Assign vertex i the smallest color that has not yet been
		// used among its neighbors with smaller index
		adjacent = getAdjacent(vertices[i]);
		available.resize(adjacent.size());
		available.set();
		for (std::size_t j = 0; j < i; ++j)
			if (isAdjacent(vertices[j], vertices[i]) && coloring[j] < adjacent.size())
				available.reset(coloring[j]);
		coloring[i] = (available.any() ? available.find_first() : adjacent.size());
	}

	return coloring;
}

bool EssentialGraph::hasEdge(const uint a, const uint b) const
{
	bool result;
	boost::tie(boost::tuples::ignore, result) = boost::edge(a, b, _graph);
	return result;
}

std::set<uint> EssentialGraph::getParents(const uint vertex) const
{
	std::set<uint> result;
	InEdgeIter ei, ei_end;
	Edge edge;

	for (boost::tie(ei, ei_end) = boost::in_edges(vertex, _graph); ei != ei_end; ei++) {
		edge = Edge(*ei, _graph);
		if (!hasEdge(edge.target, edge.source))
			result.insert(edge.source);
	}

	return result;
}

std::set<uint> EssentialGraph::getChildren(const uint vertex) const
{
	std::set<uint> result;
	OutEdgeIter ei, ei_end;
	Edge edge;

	for (boost::tie(ei, ei_end) = boost::out_edges(vertex, _graph); ei != ei_end; ei++) {
		edge = Edge(*ei, _graph);
		if (!hasEdge(edge.target, edge.source))
			result.insert(edge.target);
	}

	return result;
}

std::set<uint> EssentialGraph::getNeighbors(const uint vertex) const
{
	std::set<uint> result;
	InEdgeIter ei, ei_end;
	Edge edge;

	for (boost::tie(ei, ei_end) = boost::in_edges(vertex, _graph); ei != ei_end; ei++) {
		edge = Edge(*ei, _graph);
		if (hasEdge(edge.target, edge.source))
			result.insert(edge.source);
	}

	return result;
}

std::set<uint> EssentialGraph::getAdjacent(const uint vertex) const
{
	std::set<uint> result;
	InEdgeIter inIter, ei_end;
	OutEdgeIter outIter, outLast;

	for (boost::tie(inIter, ei_end) = boost::in_edges(vertex, _graph); inIter != ei_end; inIter++)
		result.insert(boost::source(*inIter, _graph));
	for (boost::tie(outIter, outLast) = boost::out_edges(vertex, _graph); outIter != outLast; outIter++)
		result.insert(boost::target(*outIter, _graph));

	return result;
}

std::set<uint> EssentialGraph::getInEdges(const uint vertex) const
{
	std::set<uint> result;
	InEdgeIter inIter, ei_end;

	for (boost::tie(inIter, ei_end) = boost::in_edges(vertex, _graph); inIter != ei_end; inIter++)
		result.insert(boost::source(*inIter, _graph));

	return result;
}

uint EssentialGraph::getDegree(const uint vertex) const
{
	return getAdjacent(vertex).size();
}

void EssentialGraph::setFixedGaps(const EssentialGraph& fixedGaps, const bool inverted)
{
	_fixedGaps = fixedGaps._graph;
	_gapsInverted = inverted;
}


void EssentialGraph::limitVertexDegree(const std::vector<uint>& maxVertexDegree)
{
	if (maxVertexDegree.size() != getVertexCount())
		throw std::runtime_error("Number of vertex degrees must coincide with number of vertices");
	std::copy(maxVertexDegree.begin(), maxVertexDegree.end(), _maxVertexDegree.begin());
}

void EssentialGraph::limitVertexDegree(const uint maxVertexDegree)
{
	std::fill(_maxVertexDegree.begin(), _maxVertexDegree.end(), maxVertexDegree);
}

void EssentialGraph::limitVertexDegree(const double maxRelativeDegree)
{
	// TODO: note that getDataCount() is not implemented at the moment in the class
	// ScoreRFunction. Either implement it there, or change the behaviour here in
	// this function: limiting the vertex degree based on the number of data points
	// does not necessarily make sense since the number of parameters used to
	// describe the distribution of the child does not necessarily scale linearly
	// with the number of parents, as it is the case in the Gaussian distribution.
	for (uint i = 0; i < getVertexCount(); ++i)
		_maxVertexDegree[i] = static_cast<uint>(maxRelativeDegree * (_score->getDataCount(i)));
}

boost::dynamic_bitset<> EssentialGraph::getAnteriorSet(const std::set<uint>& A)
{
	boost::dynamic_bitset<> result(getVertexCount());
    InEdgeIter ei, ei_end;
    Edge edge;
    std::set<uint>::iterator vi;
    std::stack<uint> nbhd;
    uint a;

    for (vi = A.begin(); vi != A.end(); ++vi) {
    	// Start a DFS at vertex *vi (one of the vertices in the start set)
    	nbhd.push(*vi);
    	result.set(*vi);
    	while (!nbhd.empty()) {
    		a = nbhd.top();
    		nbhd.pop();
    		// Move through the graph "backward", i.e. from edge targets to edge sources
    		for (boost::tie(ei, ei_end) = boost::in_edges(a, _graph); ei != ei_end; ++ei) {
    			edge = Edge(*ei, _graph);
    			// If newly detected source was not yet visited, add it to the neighborhood
    			// and to the result set (= anterior set)
    			if (!result.test(edge.source)) {
    				nbhd.push(edge.source);
    				result.set(edge.source);
    			}
    		}

    	}
    }

	return result;
}

boost::dynamic_bitset<> EssentialGraph::getPosteriorSet(const std::set<uint>& A)
{
	boost::dynamic_bitset<> result(getVertexCount());
    OutEdgeIter ei, ei_end;
    Edge edge;
    std::set<uint>::iterator vi;
    std::stack<uint> nbhd;
    uint a;

    for (vi = A.begin(); vi != A.end(); ++vi) {
    	// Start a DFS at vertex *vi (one of the vertices in the start set)
    	nbhd.push(*vi);
    	result.set(*vi);
    	while (!nbhd.empty()) {
    		a = nbhd.top();
    		nbhd.pop();
    		// Move through the graph "forward", i.e. from edge soruces to edge targets
    		for (boost::tie(ei, ei_end) = boost::out_edges(a, _graph); ei != ei_end; ++ei) {
    			edge = Edge(*ei, _graph);
    			// If newly detected target was not yet visited, add it to the neighborhood
    			// and to the result set (= posterior set)
    			if (!result.test(edge.target)) {
    				nbhd.push(edge.target);
    				result.set(edge.target);
    			}
    		}

    	}
    }

	return result;

}

//arma::umat EssentialGraph::getAdjacencyMatrix() const
//{
//	arma::umat adjacency(getVertexCount(), getVertexCount());
//	adjacency.zeros();
//	boost::graph_traits<InternalEssentialGraph>::edge_iterator ei, ei_end;
//	for (boost::tie(ei, ei_end) = boost::edges(_graph); ei != ei_end; ++ei)
//		adjacency(boost::source(*ei, _graph), boost::target(*ei, _graph)) = 1;
//
//	return adjacency;
//}

EssentialGraph EssentialGraph::getRepresentative() const
{
	EssentialGraph representative;
	representative._graph = _graph;

	// Orient all edges in the chain components according to a LexBFS-ordering
	boost::dynamic_bitset<> notVisited(getVertexCount());
	notVisited.set();
	std::set<uint> chainComp;
	std::set<uint>::iterator vi;
	uint v;
	while ((v = notVisited.find_first()) < getVertexCount()) {
		chainComp = representative.getChainComponent(v);
		representative.lexBFS(chainComp.begin(), chainComp.end(), true);
		for (vi = chainComp.begin(); vi != chainComp.end(); ++vi)
			notVisited.reset(*vi);
	}

	return representative;
}

void EssentialGraph::enableCaching()
{
	if (!_doCaching) {
		_doCaching = true;
		_actualPhase = SD_NONE;
		_scoreCache = std::vector<ArrowChange>(getVertexCount());
	}
}

void EssentialGraph::disableCaching()
{
	_doCaching = false;
	_actualPhase = SD_NONE;
	_scoreCache.clear();
}

std::set<uint> EssentialGraph::getChainComponent(const uint v) const
{
	std::vector<uint> nbhd(1, v);
	std::set<uint> chainComp;
	uint a;
	AdjacencyIter vi, vi_end;
	while (!nbhd.empty()) {
		a = nbhd.back();
		nbhd.pop_back();
		chainComp.insert(a);
		for (boost::tie(vi, vi_end) = boost::adjacent_vertices(a, _graph); vi != vi_end; vi++)
			if (hasEdge(*vi, a) && std::find(nbhd.begin(), nbhd.end(), *vi) == nbhd.end() && chainComp.find(*vi) == chainComp.end())
				nbhd.push_back(*vi);
	}
	return chainComp;
}

bool EssentialGraph::addLogger(GraphOperationLogger* logger)
{
	bool result;
	boost::tie(boost::tuples::ignore, result) = _loggers.insert(logger);
	return result;
}

bool EssentialGraph::removeLogger(GraphOperationLogger* logger)
{
	return _loggers.erase(logger) == 0;
}

ArrowChange EssentialGraph::getOptimalArrowInsertion(const uint v)
{
	// For DEBUGGING purposes: print vertex being processed
	dout.level(2) << "Calculating optimal arrow insertion for vertex " << v << "\n";

	// Initialize optimal score gain
	ArrowChange result;
	result.score = 0.;

	// Respect maximum vertex degree and vertices that can only have children
	if (_childrenOnly.test(v) || (_maxVertexDegree[v] < getVertexCount() && getDegree(v) >= _maxVertexDegree[v]))
		return result;
	else
		result.score = _minScoreDiff;

	double diffScore;
	boost::unordered_map<std::set<uint>, double > localScore;
	boost::unordered_map<std::set<uint>, double >::iterator hmi;

	// Find maximal cliques in the neighborhood of v
	std::set<uint> neighbors = getNeighbors(v);
	std::vector<std::set<uint> > maxCliques = getMaxCliques(neighbors.begin(), neighbors.end());

	// Get parents of v (used for calculation of partial scores later on)
	std::set<uint> parents = getParents(v);

	// Exclude forbidden sources:
	// - vertices reachable from children of v
	boost::dynamic_bitset<> forbidden = getPosteriorSet(getChildren(v));
	// - vertices adjacent to v
	std::set<uint> tempSet = getAdjacent(v);
	for (std::set<uint>::iterator si = tempSet.begin(); si != tempSet.end(); ++si)
		forbidden.set(*si);
	// - v itself :-)
	forbidden.set(v);
	// - vertices which have reached the maximum degree, or which have a fixed
	//   gap to v
	for (uint u = 0; u < getVertexCount(); ++u)
		if (getDegree(u) >= _maxVertexDegree[u] || gapFixed(u, v))
			forbidden.set(u);

	// Calculate vertices not reachable from v: for those, the "path condition"
	// for the clique C does not have to be checked later
	tempSet = std::set<uint>();
	tempSet.insert(v);
	boost::dynamic_bitset<> posterior = getPosteriorSet(tempSet);

	for (uint u = 0; u < getVertexCount(); ++u)
		if (!forbidden[u]) {
			// Calculate ne(v) \cap ad(u)
			std::set<uint> N = set_intersection(neighbors, getAdjacent(u));

			// Add N as a set to check, and at the same time as a stop set.
			// Actually, N will be checked _only_ if it is a clique, i.e. subset
			// of a maximal clique
			CliqueStack cliqueStack;
			cliqueStack.push_back(N);
			cliqueStack.stop_sets.insert(N);

			for (std::size_t i = 0; i < maxCliques.size(); ++i) {
				// Only consider maximal cliques that contain N
				if (std::includes(maxCliques[i].begin(), maxCliques[i].end(), N.begin(), N.end())) {
					// Check all subsets of the actual maximal clique
					cliqueStack.append(maxCliques[i]);
					while(!cliqueStack.empty()) {
						std::set<uint> C = cliqueStack.back();
						cliqueStack.pop_back();

						// Check whether there is a v-u-path that does not go through C
						if (!posterior.test(u) || !existsPath(v, u, C)) {
							// Calculate BIC score difference for current clique C
							// Note: if calculation is not possible (too low rank of
							// submatrices), local should return NaN; then the
							// test below fails
							// Use "localScore" as (additional) cache
							std::set<uint> C_par = set_union(C, parents);
							hmi = localScore.find(C_par);
							if (hmi == localScore.end()) {
								dout.level(3) << "calculating partial score for vertex " << v << ", parents " << C_par << "...\n";
								diffScore = - _score->local(v, C_par);
								localScore[C_par] = diffScore;
							}
							else
								diffScore = hmi->second;
							dout.level(3) << "partial score for vertex " << v << ", parents " << C_par << ": " << -diffScore << "\n";
							C_par.insert(u);
							diffScore += _score->local(v, C_par);
							dout.level(3) << "partial score for vertex " << v << ", parents " << C_par << ": " << _score->local(v, C_par) << "\n";

							// If new score is better than previous optimum, and there is no
							// v-u-path that does not go through C, store (u, v, C) as new optimum
							if (diffScore > result.score) {
								result.source = u;
								result.clique = C;
								result.score = diffScore;
							}
						}

						// Add all subsets of C that differ from C in only one vertex to the stack
						for (std::set<uint>::iterator si = C.begin(); si != C.end(); ++si) {
							std::set<uint> C_sub = C;
							C_sub.erase(*si);
							cliqueStack.append(C_sub);
						}
					}

					// Mark the actual maximal set as checked, i.e. add it to the stop sets
					cliqueStack.stop_sets.insert(maxCliques[i]);
				}
			}
		}

	if (result.score == _minScoreDiff)
		result.score = 0.;
	return result;
}

ArrowChange EssentialGraph::getOptimalArrowDeletion(const uint v)
{
	std::set<uint> neighbors, parents, candidates;
	std::vector<std::set<uint> > maxCliques;
	std::set<uint> C, C_par, C_sub, N;
	std::set<uint>::iterator iter, ui;
	double diffScore;
	CliqueStack cliqueStack;
	boost::unordered_map<std::set<uint>, double > localScore;
	boost::unordered_map<std::set<uint>, double >::iterator hmi;

	// For DEBUGGING purposes: print vertex being processed
	dout.level(2) << "Calculating optimal arrow deletion for vertex " << v << "\n";

	// Initialize optimal score gain
	ArrowChange result;
	result.score = _minScoreDiff;

	// Get neighbors and parents of v (used for calculation of BIC score later on)
	neighbors = getNeighbors(v);
	parents = getParents(v);
	candidates = set_union(neighbors, parents);
	for (ui = candidates.begin(); ui != candidates.end(); ++ui) {
		N = set_intersection(neighbors, getAdjacent(*ui));

		// Find all maximal cliques on N
		maxCliques = getMaxCliques(N.begin(), N.end());
		cliqueStack.clear_all();

		// Calculate the score difference for all cliques in N
		for (std::size_t i = 0; i < maxCliques.size(); ++i) {
			// Check all subsets of the actual maximal clique
			cliqueStack.append(maxCliques[i]);
			while(!cliqueStack.empty()) {
				C = cliqueStack.back();
				cliqueStack.pop_back();

				// Calculate BIC score difference for current clique C
				// Use "localScore" as (additional) cache
				C_par = set_union(C, parents);
				C_par.insert(*ui);
				hmi = localScore.find(C_par);
				if (hmi == localScore.end()) {
					dout.level(3) << "calculating partial score for vertex " << v << ", parents " << C_par << "...\n";
					diffScore = - _score->local(v, C_par);
					localScore[C_par] = diffScore;
				}
				else
					diffScore = hmi->second;
				C_par.erase(*ui);
				diffScore += _score->local(v, C_par);

				// If new score is better than previous optimum, store (u, v, C) as new optimum
				if (diffScore > result.score) {
					result.source = *ui;
					result.clique = C;
					result.score = diffScore;
				}

				// Add all subsets of C that differ from C in only one vertex to the stack
				for (iter = C.begin(); iter != C.end(); ++iter) {
					C_sub = C;
					C_sub.erase(*iter);
					cliqueStack.append(C_sub);
				}
			} // while !cliqueStack.empty()

			// Mark the actual maximal set as checked, i.e. add it to the stop sets
			cliqueStack.stop_sets.insert(maxCliques[i]);
		} // for i
	} // for ui

	if (result.score == _minScoreDiff)
		result.score = 0.;
	return result;
}

ArrowChange EssentialGraph::getOptimalArrowTurning(const uint v)
{
	std::set<uint> children, neighbors, neighbors2, parents, C, C_par, C_opt, C_sub, CminN, N;
	std::vector<std::set<uint> > maxCliques;
	std::set<uint>::iterator iter, ui;
	double diffScore;
	CliqueStack cliqueStack;

	// For DEBUGGING purposes: print vertex being processed
	dout.level(2) << "Calculating optimal arrow turning for vertex " << v << "\n";

	// Initialize optimal score gain
	ArrowChange result;
	result.score = _minScoreDiff;

	// Respect vertices that are not allowed to have parents, but only children
	if (!_childrenOnly.test(v))  {
		dout.level(3) << "  checking edges incident to v = " << v << "...\n";

		// Store parents and neighbors of v
		neighbors = getNeighbors(v);
		parents = getParents(v);

		// Search over all neighbors of v; turning non-essential arrows
		for (ui = neighbors.begin(); ui != neighbors.end(); ++ui) {
			N = set_intersection(neighbors, getNeighbors(*ui));

			// Find all maximal cliques in the neighborhood of v (without u)
			neighbors2 = neighbors;
			neighbors2.erase(*ui);
			maxCliques = getMaxCliques(neighbors2.begin(), neighbors2.end());
			cliqueStack.clear_all();

			// Calculate the score difference for all (admissible) cliques in the neighborhood of v
			for (std::size_t i = 0; i < maxCliques.size(); ++i) {
				// Check all subsets of the actual maximal clique
				cliqueStack.append(maxCliques[i]);
				while(!cliqueStack.empty()) {
					C = cliqueStack.back();
					cliqueStack.pop_back();

					// Check if C \ N is not empty, and if C \cap N separates C \ N and N \ C in
					// the neighborhood of v
					CminN = set_difference(C, N);
					if (!(CminN.empty()) && !existsPath(set_difference(neighbors, set_intersection(C, N)), *(CminN.begin()), set_difference(N, C))) {
						// Calculate BIC score difference for current clique C
						C_par = set_union(C, parents);
						diffScore = - _score->local(v, C_par);
						C_par.insert(*ui);
						diffScore += _score->local(v, C_par);
						C_par = set_union(set_intersection(C, N), getParents(*ui));
						diffScore += _score->local(*ui, C_par);
						C_par.insert(v);
						diffScore -= _score->local(*ui, C_par);

						// If new score is better than previous optimum, store (u, v, C) as new optimum
						if (diffScore > result.score) {
							result.source = *ui;
							result.clique = C;
							result.score = diffScore;
						}
					}

					// Add all subsets of C that differ from C in only one vertex to the stack
					for (iter = C.begin(); iter != C.end(); ++iter) {
						C_sub = C;
						C_sub.erase(*iter);
						cliqueStack.append(C_sub);
					}
				} // while !cliqueStack.empty()

				// Mark the actual maximal set as checked, i.e. add it to the stop sets
				cliqueStack.stop_sets.insert(maxCliques[i]);
			} // for i
		} // for ui

		// Find all maximal cliques in the neighborhood of v
		maxCliques = getMaxCliques(neighbors.begin(), neighbors.end());

		// Search over all children of v; turning essential arrows
		children = getChildren(v);
		for (ui = children.begin(); ui != children.end(); ++ui) {
			dout.level(3) << "  trying to turn arrow (" << v << ", " << *ui << ")...\n";
			N = set_intersection(neighbors, getParents(*ui));

			// Add N as a set to check, and at the same time as a stop set.
			// Actually, N will be checked _only_ if it is a clique, i.e. subset
			// of a maximal clique
			cliqueStack.clear_all();
			cliqueStack.push_back(N);
			cliqueStack.stop_sets.insert(N);

			for (std::size_t i = 0; i < maxCliques.size(); ++i) {
				// Only consider maximal cliques that contain N
				if (std::includes(maxCliques[i].begin(), maxCliques[i].end(), N.begin(), N.end())) {
					// Check all subsets of the actual maximal clique
					cliqueStack.append(maxCliques[i]);
					while(!cliqueStack.empty()) {
						C = cliqueStack.back();
						cliqueStack.pop_back();

						// Check if there is no v-u-path (except (v, u)) that does not visit C \cup ne(u)
						if (!existsPath(v, *ui, set_union(C, getNeighbors(*ui)))) {
							// Calculate BIC score difference for current clique C
							C_par = set_union(C, parents);
							diffScore = - _score->local(v, C_par);
							C_par.insert(*ui);
							diffScore += _score->local(v, C_par);
							C_par = getParents(*ui);
							diffScore -= _score->local(*ui, C_par);
							C_par.erase(v);
							diffScore += _score->local(*ui, C_par);
							dout.level(3) << "  score difference for (u, v, C) = (" << *ui << ", " << v << ", " << C << "): " << diffScore << "\n";


							// If new score is better than previous optimum, and there is no
							// v-u-path that does not go through C, store (u, v, C) as new optimum
							if (diffScore > result.score) {
								result.source = *ui;
								result.clique = C;
								result.score = diffScore;
							}
						}

						// Add all subsets of C that differ from C in only one vertex to the stack
						for (iter = C.begin(); iter != C.end(); ++iter) {
							C_sub = C;
							C_sub.erase(*iter);
							cliqueStack.append(C_sub);
						}
					}

					// Mark the actual maximal set as checked, i.e. add it to the stop sets
					cliqueStack.stop_sets.insert(maxCliques[i]);
				}
			}
		} // for ui
	} // IF !_childrenOnly.test(v)

	if (result.score == _minScoreDiff)
		result.score = 0.;
	return result;
}

std::set<uint> EssentialGraph::_bitsToParents(const uint vertex, const uint32_t bits)
{
	std::set<uint> parents;
	uint32_t pattern = 1;
	for (uint i = 0; i < getVertexCount(); i++) {
		if (i != vertex) {
			if (bits & pattern)
				parents.insert(i);
			pattern *= 2;
		}
	}

	return parents;
}

std::set<uint> EssentialGraph::_getOptimalUnrestrTarget()
{
	std::set<uint> target;

	// Get a greedy coloring of each chain component; return all vertices
	// with a color in the lower half
	boost::dynamic_bitset<> notVisited(getVertexCount());
	notVisited.set();
	std::set<uint> chainComp;
	std::vector<uint> ordering, coloring;
	std::vector<uint>::iterator vi;
	uint v, i, chromaticNumber;
	while ((v = notVisited.find_first()) < getVertexCount()) {
		// Get LexBFS-ordering of chain component of v
		chainComp = getChainComponent(v);
		ordering = lexBFS(chainComp.begin(), chainComp.end());

		// Get greedy coloring w.r.t. LexBFS-ordering
		coloring = greedyColoring(ordering);

		// Add vertices of lower half of coloring to the target
		chromaticNumber = *(std::max_element(coloring.begin(), coloring.end())) + 1;
		for (i = 0; i < ordering.size(); ++i) {
			if (coloring[i] < chromaticNumber/2)
				target.insert(ordering[i]);
			notVisited.reset(ordering[i]);
		}
	}

	return target;
}

uint EssentialGraph::_getOptimalSingleVertexTarget()
{
	uint u, v_opt, v_ind, i;
	uint eta, eta_opt, eta_min;
	std::set<uint> chainComp, chainCompVert, C, C_sub, neighbors;
	std::set<uint>::iterator si, sj;
	std::vector<uint> startOrder;
	std::vector<std::set<uint> > maxCliques;
	CliqueStack cliqueStack;
	EssentialGraph chainCompGraph, chainCompGraphCopy;
	TargetFamily targets;

	// Initialize optimal number of orientable edges
	eta_opt = 0;
	v_opt = getVertexCount();  // NO intervention helps anymore

	// Vertices that have not been checked yet
	boost::dynamic_bitset<> notChecked(getVertexCount());
	notChecked.set();

	// Search over all chain components
	while ((u = notChecked.find_first()) < getVertexCount()) {
		dout.level(3) << "  checking chain component of vertex " << u << "...\n";
		// Find chain component of v, and check whether it is non-trivial
		chainComp = getChainComponent(u);
		notChecked.reset(u);

		if (chainComp.size() > 1) {
			// Extract subgraph corresponding to chain component of v
			chainCompGraph = inducedSubgraph(chainComp.begin(), chainComp.end());

			v_ind = 0;
			for (si = chainComp.begin(); si != chainComp.end(); ++si, ++v_ind) {
				dout.level(3) << "  checking vertex " << *si << "...\n";
				notChecked.reset(*si);
				chainCompVert.clear();
				for (i = 0; i < chainComp.size(); ++i)
					chainCompVert.insert(i);

				// From now on, work in the induced subgraph

				// Set target family
				targets.clear();
				C.clear();
				targets.push_back(C);
				C.insert(v_ind);
				targets.push_back(C);
				chainCompGraph.setTargets(&targets);

				// Find maximal cliques in the neighborhood of v
				neighbors = chainCompGraph.getNeighbors(v_ind);
				maxCliques = chainCompGraph.getMaxCliques(neighbors.begin(), neighbors.end());

				// Initialize minimal number of orientable edges for given vertex
				eta_min = getVertexCount()*getVertexCount();

				// Check all cliques in the neighborhood of v...
				cliqueStack.clear_all();
				for (i = 0; i < maxCliques.size(); ++i) {
					cliqueStack.append(maxCliques[i]);
					while (!cliqueStack.empty()) {
						// ... and orient the edges of the chain component accordingly
						chainCompGraphCopy = chainCompGraph;
						C = cliqueStack.back();
						cliqueStack.pop_back();
						startOrder = std::vector<uint>(C.begin(), C.end());
						startOrder.push_back(v_ind);
						std::set_difference(chainCompVert.begin(), chainCompVert.end(), C.begin(), C.end(), std::inserter(startOrder, startOrder.end()));
						chainCompGraphCopy.lexBFS(startOrder.begin(), startOrder.end(), true);

						// Replace unprotected arrows (w.r.t. target v) in chain component
						chainCompGraphCopy.replaceUnprotected();
						eta = chainCompGraph.getEdgeCount() - chainCompGraphCopy.getEdgeCount();
						if (eta < eta_min)
							eta_min = eta;

						// Add all subsets of C that differ from C in only one vertex to the stack
						for (sj = C.begin(); sj != C.end(); ++sj) {
							C_sub = C;
							C_sub.erase(*sj);
							cliqueStack.append(C_sub);
						}
					} // WHILE !cliqueStack.empty()
				} // FOR i

				// Find maximal "eta_min"
				if (eta_min > eta_opt) {
					eta_opt = eta_min;
					v_opt = *si;
				}
			} // FOR si
		} // IF chainComp.size() > 1
	}

	return v_opt;
}

std::set<Edge, EdgeCmp> EssentialGraph::replaceUnprotected()
{
	// Map of arrow labels
	std::map<Edge, edge_flag, EdgeCmp> arrowFlags;
	// Set of undecidable arrows
	std::set<Edge, EdgeCmp> undecidableArrows;
	// Set of "affected" vertices, that is, vertices that are sources or targets
	// of an unprotected arrow
	std::set<Edge, EdgeCmp> result;

	dout.level(2) << "  replacing unprotected arrows...\n";

	Edge edge;

	// Find all arrows in the graph. Mark them as "protected", if they are
	// protected by an intervention target, and "undecidable" otherwise
	boost::graph_traits<InternalEssentialGraph>::edge_iterator edgeIter, edgeLast;
	for (boost::tie(edgeIter, edgeLast) = boost::edges(_graph); edgeIter != edgeLast; edgeIter++) {
		edge = Edge(*edgeIter, _graph);
		if (!hasEdge(edge.target, edge.source)) {
			if (_targets->protects(edge.source, edge.target)) arrowFlags[edge] = PROTECTED;
			else {
				undecidableArrows.insert(edge);
				arrowFlags[edge] = UNDECIDABLE;
			}
		}
	}

	// Check whether the arrows are part of a v-structure; if yes, mark them as "protected".
	std::map<Edge, edge_flag, EdgeCmp>::iterator arrIter1, arrIter2;
	uint v;
	for (v = 0; v < getVertexCount(); v++) {
		for (arrIter1 = arrowFlags.lower_bound(Edge(0, v)); arrIter1 != arrowFlags.upper_bound(Edge(getVertexCount(), v)); arrIter1++) {
			arrIter2 = arrIter1;
			for (arrIter2++; arrIter2 != arrowFlags.upper_bound(Edge(getVertexCount(), v)); arrIter2++)
				if (!isAdjacent((arrIter1->first).source, (arrIter2->first).source)) {
					arrIter1->second = PROTECTED;
					arrIter2->second = PROTECTED;
					undecidableArrows.erase(arrIter1->first);
					undecidableArrows.erase(arrIter2->first);
				}
		}
	}

	// Successively check all undecidable arrows, until no one remains
	std::set<Edge, EdgeCmp>::iterator undIter;
	edge_flag flag;
	// If the graph is in a valid state in the beginning, the following loop
	// finally flags all undecidable arrows as protected or unprotected. In
	// case the graph is in an invalid state in the beginning, it might happen
	// that the loop does not terminate; to avoid this, we also check that the
	// loop indeed labels undecidable arrows (as PROTECTED or NOT_PROTECTED) in
	// every run, and throw an error otherwise.
	int labeledArrows = 1;
	while (!undecidableArrows.empty() && labeledArrows > 0) {
		// Find unprotected and protected arrows
		for (undIter = undecidableArrows.begin(); undIter != undecidableArrows.end(); undIter++) {
			edge = *undIter;
			flag = NOT_PROTECTED;

			// Check whether the arrow is in configuration (a)
			for (arrIter1 = arrowFlags.lower_bound(Edge(0, edge.source));
					flag != PROTECTED && arrIter1 != arrowFlags.upper_bound(Edge(getVertexCount(), edge.source));
					arrIter1++)
				if (!isAdjacent((arrIter1->first).source, edge.target))
					flag = (arrIter1->second == PROTECTED ? PROTECTED : UNDECIDABLE);

			// Check whether the arrow is in configuration (c)
			for (arrIter1 = arrowFlags.lower_bound(Edge(0, edge.target));
					flag != PROTECTED && arrIter1 != arrowFlags.upper_bound(Edge(getVertexCount(), edge.target));
					arrIter1++)
				if (isParent(edge.source, (arrIter1->first).source))
					flag = (arrIter1->second == PROTECTED && arrowFlags[Edge(edge.source, (arrIter1->first).source)] == PROTECTED ? PROTECTED : UNDECIDABLE);

			// Check whether the arrow is in configuration (d)
			for (arrIter1 = arrowFlags.lower_bound(Edge(0, edge.target));
					flag != PROTECTED && arrIter1 != arrowFlags.upper_bound(Edge(getVertexCount(), edge.target));
					arrIter1++) {
				arrIter2 = arrIter1;
				for (arrIter2++;
						flag != PROTECTED && arrIter2 != arrowFlags.upper_bound(Edge(getVertexCount(), edge.target));
						arrIter2++)
					if (isNeighbor(edge.source, (arrIter1->first).source) && isNeighbor(edge.source, (arrIter2->first).source) && !isAdjacent((arrIter1->first).source, (arrIter2->first).source))
						flag = (arrIter1->second == PROTECTED && arrIter2->second == PROTECTED ? PROTECTED : UNDECIDABLE);
			}

			// Store flag
			arrowFlags[edge] = flag;
		}

		// Replace unprotected arrows by lines; store affected edges in result set
		labeledArrows = undecidableArrows.size();
		for (arrIter1 = arrowFlags.begin(); arrIter1 != arrowFlags.end(); ) {
			arrIter2 = arrIter1;
			arrIter1++;
			if (arrIter2->second != UNDECIDABLE) {
				undecidableArrows.erase(arrIter2->first);
			}
			if (arrIter2->second == NOT_PROTECTED) {
				addEdge((arrIter2->first).target, (arrIter2->first).source);
				result.insert(arrIter2->first);
				arrowFlags.erase(arrIter2);
			}
		}
		labeledArrows = labeledArrows - undecidableArrows.size();
		dout.level(3) << "  Labeled " << labeledArrows << " undecidable arrows\n";
	}

	if (labeledArrows == 0 && !undecidableArrows.empty()) {
		throw std::runtime_error("Invalid graph passed to replaceUnprotected().");
	}
	dout.level(2) << "  Done.\n";

	return result;
}

void EssentialGraph::insert(const uint u, const uint v, const std::set<uint> C)
{
	// Get a LexBFS-ordering on the chain component of v, in which all edges of C
	// point toward v, and all other edges point away from v, and orient the edges
	// of the chain component accordingly
	std::set<uint> chainComp = getChainComponent(v);
	std::vector<uint> startOrder(C.begin(), C.end());
	startOrder.push_back(v);
	chainComp.erase(v);
	std::set_difference(chainComp.begin(), chainComp.end(), C.begin(), C.end(), std::inserter(startOrder, startOrder.end()));
	lexBFS(startOrder.begin(), startOrder.end(), true);

	// Add new arrow
	addEdge(u, v);

	// Successively replace unprotected arrows by lines
	replaceUnprotected();

	/* MOVED TO greedyForward()

	// If caching is enabled, recalculate the optimal arrow insertions where
	// necessary
	if (_doCaching) {
		// Genereate set of vertices whose anterior set is the set of vertices
		// whose cache has to be refreshed:
		// u, if there was no path from u to v before
		// TODO check conditions!!!
		// if (!oldGraph.existsPath(u, v))
			recalcAnt.insert(u);
		recalc.insert(u);
		// v, if the arrow was undirected and there was no path from v to u before
		// TODO check conditions!!
		// if (hasEdge(v, u) && !oldGraph.existsPath(v, u))
		if (hasEdge(v, u))
			recalcAnt.insert(v);
		recalc.insert(v);
		// the target of any newly directed edge
		diffSet = set_difference(directed, undirected);
		for (ei = diffSet.begin(); ei != diffSet.end(); ++ei) {
			recalcAnt.insert(ei->target);
			recalc.insert(ei->source);
		}
		// the source of any newly undirected edge
		diffSet = set_difference(undirected, directed);
		for (ei = diffSet.begin(); ei != diffSet.end(); ++ei) {
			recalcAnt.insert(ei->source);
			recalc.insert(ei->target);
		}

		// Calculate anterior set of that candidate set, and add vertices that
		// have to be recalculated without the complete anterior set
		refreshCache = getAnteriorSet(recalcAnt);
		for (si = recalc.begin(); si != recalc.end(); ++si)
			refreshCache.set(*si);

		// If v or u have reached the maximum degree, recalculate the optimal
		// arrow insertion for all vertices for which an insertion with new
		// parent u or v is proposed by the cache
		if (getDegree(u) >= _maxVertexDegree[u])
			for (a = 0; a < getVertexCount(); ++a)
				if (_scoreCache[a].source == u)
					refreshCache.set(a);
		if (getDegree(v) >= _maxVertexDegree[v])
			for (a = 0; a < getVertexCount(); ++a)
				if (_scoreCache[a].source == v)
					refreshCache.set(a);

		// Refresh cache: recalculate arrow insertions
		for (a = refreshCache.find_first(); a < getVertexCount(); a = refreshCache.find_next(a))
			_scoreCache[a] = getOptimalArrowInsertion(a);
	}
	*/
}

void EssentialGraph::remove(const uint u, const uint v, const std::set<uint> C)
{
	// Get a LexBFS-ordering on the chain component of v, in which all edges of C
	// (or C \cup \{u\}, if u lies in the chain component)
	// point toward v, and all other edges point away from v, and orient the edges
	// of the chain component accordingly
	std::set<uint> chainComp = getChainComponent(v);
	std::vector<uint> startOrder(C.begin(), C.end());
	if (chainComp.find(u) != chainComp.end())
		startOrder.push_back(u);
	startOrder.push_back(v);
	chainComp.erase(v);
	std::set_difference(chainComp.begin(), chainComp.end(), C.begin(), C.end(), std::inserter(startOrder, startOrder.end()));
	lexBFS(startOrder.begin(), startOrder.end(), true);

	// Remove the edge between u and v
	removeEdge(u, v, true);

	// Successively replace unprotected arrows by lines
	replaceUnprotected();
}

void EssentialGraph::turn(const uint u, const uint v, const std::set<uint> C)
{
	std::set<uint> chainComp;
	std::vector<uint> startOrder;

	// If u and v lie in different chain components, first order the edges in the
	// chain component of u
	if (!hasEdge(u, v)) {
		chainComp = getChainComponent(u);
		chainComp.erase(u);
		startOrder.push_back(u);
		startOrder.insert(startOrder.end(), chainComp.begin(), chainComp.end());
		lexBFS(startOrder.begin(), startOrder.end(), true);
		startOrder.clear();
	}

	// Get a LexBFS-ordering on the chain component of v in which all edges of C
	// point towards v, and all other edges point away from v
	chainComp = getChainComponent(v);
	startOrder.insert(startOrder.end(), C.begin(), C.end());
	startOrder.push_back(v);
	chainComp.erase(v);
	if (hasEdge(u, v)) {
		startOrder.push_back(u);
		chainComp.erase(u);
	}
	std::set_difference(chainComp.begin(), chainComp.end(), C.begin(), C.end(), std::inserter(startOrder, startOrder.end()));
	lexBFS(startOrder.begin(), startOrder.end(), true);

	// Turn the edge (v, u)
	removeEdge(v, u);
	addEdge(u, v);

	// Successively replace unprotected arrows by lines
	replaceUnprotected();
}

bool EssentialGraph::greedyForward(const ForwardAdaptiveFlag adaptive)
{
	uint v_opt = 0;
	ArrowChangeCmp comp;
	ArrowChange insertion, optInsertion;

	// For DEBUGGING purposes: print phase
	dout.level(1) << "== starting forward phase ("
			<< (adaptive ? "" : "not ") << "adaptive)...\n";

	// Initialize optimal score gain
	optInsertion.score = _minScoreDiff;

	// If caching is disabled calculate the score differences for all possible edges
	if (!_doCaching) {
		for (uint v = 0; v < getVertexCount(); v++) {
			// Calculate optimal arrow insertion for given target vertex v
			insertion = getOptimalArrowInsertion(v);

			// Look for optimal score
			if (insertion.score > optInsertion.score) {
				optInsertion = insertion;
				v_opt = v;
			}
		}
	}

	// If caching is enabled, search the cache for the best score
	dout.level(3) << "vertex count: " << getVertexCount() << "\n";
	if (_doCaching) {
		// If score has to be initialized (current phase is not forward), do it
		if (_actualPhase != SD_FORWARD)
			for (uint v = 0; v < getVertexCount(); v++)
				_scoreCache[v] = getOptimalArrowInsertion(v);

		// Find optimal arrow insertion from cache
		std::vector<ArrowChange>::iterator si = std::max_element(_scoreCache.begin(), _scoreCache.end(), comp);
		v_opt = std::distance(_scoreCache.begin(), si);
		optInsertion = *si;
		_actualPhase = SD_FORWARD;
	}

	// If the score can be augmented, do it
	if (!check_interrupt() && optInsertion.score > _minScoreDiff) {
		// For DEBUGGING purposes: print inserted arrow
		dout.level(1) << "  inserting edge (" << optInsertion.source << ", " << v_opt << ") with C = "
				<< optInsertion.clique << ", S = " << optInsertion.score << "\n";

		uint u_opt = optInsertion.source;
		EdgeOperationLogger edgeLogger;
		if (_doCaching) {
			addLogger(&edgeLogger);
		}
		insert(u_opt, v_opt, optInsertion.clique);

		// Adapt fixed gaps if requested (cf. "ARGES")
		if (adaptive == VSTRUCTURES && !hasEdge(v_opt, u_opt)) {
			std::set<uint> sources = set_difference(getParents(v_opt), getAdjacent(u_opt));
			sources.erase(u_opt);
			for (std::set<uint>::iterator si = sources.begin(); si != sources.end(); ++si) {
				setFixedGap(*si, u_opt, false);
				setFixedGap(u_opt, *si, false);
			}
		} // IF VSTRUCTURES
		else if (adaptive == TRIPLES) {
			// Adjacent sets of u_opt and v_opt
			std::vector< std::set<uint> > adjacentSets(2);
			adjacentSets[0] = getAdjacent(u_opt);
			adjacentSets[1] = getAdjacent(v_opt);
			std::vector<uint> edgeVertices(2);
			edgeVertices[0] = u_opt;
			edgeVertices[1] = v_opt;

			// Vertices adjacent to u, but not to v (without v itself) (and vice versa)
			// build new unshielded triples
			std::set<uint> triples;
			for (uint j = 0; j <= 1; j++) {
				triples = set_difference(adjacentSets[j % 2], adjacentSets[(j + 1) % 2]);
				triples.erase(edgeVertices[(j + 1) % 2]);
				for (std::set<uint>::iterator si = triples.begin(); si != triples.end(); ++si) {
					setFixedGap(*si, edgeVertices[(j + 1) % 2], false);
					setFixedGap(edgeVertices[(j + 1) % 2], *si, false);
				} // FOR si
			} // FOR j
		} // IF TRIPLES

		// If caching is enabled, recalculate the optimal arrow insertions where
		// necessary
		if (_doCaching) {
			std::set<uint> recalc, recalcAnt;

			// Genereate set of vertices whose anterior set is the set of vertices
			// whose cache has to be refreshed:
			// u, if there was no path from u to v before
			// TODO check conditions!!!
			// if (!oldGraph.existsPath(u, v))
			recalcAnt.insert(u_opt);
			recalc.insert(u_opt);
			// v, if the arrow was undirected and there was no path from v to u before
			// TODO check conditions!!
			// if (hasEdge(v, u) && !oldGraph.existsPath(v, u))
			if (hasEdge(v_opt, u_opt))
				recalcAnt.insert(v_opt);
			recalc.insert(v_opt);
			// the target of any newly directed edge
			for (std::set<Edge, EdgeCmp>::iterator ei = edgeLogger.removedEdges().begin();
					ei != edgeLogger.removedEdges().end(); ++ei) {
				dout.level(3) << "New directed edge: (" << ei-> source << ", " << ei->target << ")\n";
				recalcAnt.insert(ei->source);
				recalc.insert(ei->target);
			}
			// the source of any newly undirected edge
			for (std::set<Edge, EdgeCmp>::iterator ei = edgeLogger.addedEdges().begin();
					ei != edgeLogger.addedEdges().end(); ++ei) {
				// The newly inserted arrow is not a newly undirected one
				// Thanks to Marco Eigenmann for reported a bug here.
				if (ei->source != u_opt || ei->target != v_opt) {
					dout.level(3) << "New undirected edge: (" << ei-> source << ", " << ei->target << ")\n";
					recalcAnt.insert(ei->target);
					recalc.insert(ei->source);
				}
			}

			// Calculate anterior set of that candidate set, and add vertices that
			// have to be recalculated without the complete anterior set
			boost::dynamic_bitset<> refreshCache(getVertexCount());
			refreshCache = getAnteriorSet(recalcAnt);
			for (std::set<uint>::iterator si = recalc.begin(); si != recalc.end(); ++si)
				refreshCache.set(*si);

			// If v or u have reached the maximum degree, recalculate the optimal
			// arrow insertion for all vertices for which an insertion with new
			// parent u or v is proposed by the cache
			if (getDegree(u_opt) >= _maxVertexDegree[u_opt])
				for (uint a = 0; a < getVertexCount(); ++a)
					if (_scoreCache[a].source == u_opt)
						refreshCache.set(a);
			if (getDegree(v_opt) >= _maxVertexDegree[v_opt])
				for (uint a = 0; a < getVertexCount(); ++a)
					if (_scoreCache[a].source == v_opt)
						refreshCache.set(a);

			// Refresh cache: recalculate arrow insertions
			for (uint a = refreshCache.find_first(); a < getVertexCount(); a = refreshCache.find_next(a))
				_scoreCache[a] = getOptimalArrowInsertion(a);

			// Unregister logger
			removeLogger(&edgeLogger);
		}

		return true;
	}
	else
		return false;
}

bool EssentialGraph::greedyBackward()
{
	uint v_opt = 0;
	ArrowChange deletion, optDeletion;

	// For DEBUGGING purposes: print phase
	dout.level(1) << "== starting backward phase...\n" ;

	// Initialize optimal score gain
	optDeletion.score = _minScoreDiff;

	// TODO Allow caching for backward phase. At the moment, assumes no caching.
	for (uint v = 0; v < getVertexCount(); v++) {
		// Calculate optimal arrow insertion for given target vertex v
		deletion = getOptimalArrowDeletion(v);

		// Look for optimal score
		if (deletion.score > optDeletion.score) {
			optDeletion = deletion;
			v_opt = v;
		}
	}

	// If caching is enabled, store current phase...
	// TODO: Change that to admit actual caching also in the backward phase!!!
	if (_doCaching)
		_actualPhase = SD_BACKWARD;

	// If the score can be augmented, do it
	if (!check_interrupt() && optDeletion.score > _minScoreDiff) {
		// For DEBUGGING purposes
		dout.level(1) << "  deleting edge (" << optDeletion.source << ", " << v_opt << ") with C = "
				<< optDeletion.clique << ", S = " << optDeletion.score << "\n";
		remove(optDeletion.source, v_opt, optDeletion.clique);
		//getAdjacencyMatrix().print("A = ");
		return true;
	}
	else
		return false;
}

bool EssentialGraph::greedyTurn()
{
	uint v_opt = 0;
	ArrowChange turning, optTurning;

	// For DEBUGGING purposes: print phase
	dout.level(1) << "== starting turning phase...\n" ;

	// Initialize optimal score gain
	optTurning.score = _minScoreDiff;

	// TODO Allow caching for turning phase. At the moment, assumes no caching.
	for (uint v = 0; v < getVertexCount(); v++) {
		// Calculate optimal arrow insertion for given target vertex v
		turning = getOptimalArrowTurning(v);

		// Look for optimal score
		if (turning.score > optTurning.score) {
			optTurning = turning;
			v_opt = v;
		}
	}

	// If caching is enabled, store current phase...
	// TODO: Change that to admit actual caching also in the turning phase!!!
	if (_doCaching)
		_actualPhase = SD_TURNING;

	// Turn the highest-scoring edge
	// If the score can be augmented, do it
	if (!check_interrupt() && optTurning.score > _minScoreDiff) {
		// For DEBUGGING purposes
		dout.level(1) << "  turning edge (" << v_opt << ", " << optTurning.source << ") with C = "
				<< optTurning.clique << ", S = " << optTurning.score << "\n";
		turn(optTurning.source, v_opt, optTurning.clique);
		//getAdjacencyMatrix().print("A = ");
		return true;
	}
	else
		return false;
}

bool EssentialGraph::greedyStepDir(const step_dir direction, const ForwardAdaptiveFlag adaptive)
{
	switch (direction) {
	case SD_FORWARD:
		return greedyForward(adaptive);

	case SD_BACKWARD:
		return greedyBackward();

	case SD_TURNING:
		return greedyTurn();

	default:
		return false;
	} // SWITCH direction
}

step_dir EssentialGraph::greedyStep()
{
	uint v_opt = 0;
	step_dir optDir;
	ArrowChange change, optChange;

	// For DEBUGGING purposes: print phase
	dout.level(3) << "== looking for optimal step...\n" ;

	// Initialize optimal score gain
	optChange.score = _minScoreDiff;
	optDir = SD_NONE;

	// Look for optimal arrow insertion
	for (uint v = 0; v < getVertexCount(); v++) {
		change = getOptimalArrowInsertion(v);

		// Look for optimal score
		if (change.score > optChange.score) {
			optChange = change;
			v_opt = v;
			optDir = SD_FORWARD;
		}
	}

	// Look for optimal arrow deletion
	for (uint v = 0; v < getVertexCount(); v++) {
		change = getOptimalArrowDeletion(v);

		// Look for optimal score
		if (change.score > optChange.score) {
			optChange = change;
			v_opt = v;
			optDir = SD_BACKWARD;
		}
	}

	// Look for optimal arrow turning
	for (uint v = 0; v < getVertexCount(); v++) {
		change = getOptimalArrowTurning(v);

		// Look for optimal score
		if (change.score > optChange.score) {
			optChange = change;
			v_opt = v;
			optDir = SD_TURNING;
		}
	}

	// If caching is enabled, reset cache since it is not valid any more
	// after an arbitrary step
	if (_doCaching)
		_actualPhase = SD_NONE;

	// If the score can be augmented, perform the optimal step
	switch(optDir) {
	case SD_FORWARD :
		dout.level(3) << "  inserting edge (" << optChange.source << ", " << v_opt << ") with C = "
				<< optChange.clique << ", S = " << optChange.score << "\n";
		insert(optChange.source, v_opt, optChange.clique);
		break;
	case SD_BACKWARD :
		dout.level(1) << "  deleting edge (" << optChange.source << ", " << v_opt << ") with C = "
				<< optChange.clique << ", S = " << optChange.score << "\n";
		remove(optChange.source, v_opt, optChange.clique);
		break;
	case SD_TURNING :
		dout.level(1) << "  turning edge (" << v_opt << ", " << optChange.source << ") with C = "
				<< optChange.clique << ", S = " << optChange.score << "\n";
		turn(optChange.source, v_opt, optChange.clique);
		break;
	case SD_NONE :
		break;
	}

	return optDir;
}

bool EssentialGraph::greedyDAGForward()
{
	uint u_opt = 0, v_opt = 0;
	double diffScore, diffScore_opt;
	std::set<uint> parents, C_new;

	dout.level(2) << "= Starting forward step...\n";

	// Initialize help variables
	diffScore_opt = _minScoreDiff;
	uint p = getVertexCount();

	// Find edge that maximally increases BIC score when added
	for (uint v = 0; v < p; ++v) {
		parents = getParents(v);
		for (uint u = 0; u < p; ++u)
			if (u != v && !isAdjacent(u, v) && !gapFixed(u, v) && !existsPath(v, u)) {
				// Calculate BIC score difference for adding edge (u, v)
				C_new = parents;
				diffScore = - _score->local(v, C_new);
				C_new.insert(u);
				diffScore += _score->local(v, C_new);
				dout.level(3) << "  Score diff. for edge " << u << " --> " << v << " : " <<
						diffScore << std::endl;

				// If new score is better than previous optimum
				if (diffScore > diffScore_opt) {
					u_opt = u;
					v_opt = v;
					diffScore_opt = diffScore;
				}
			}
	}

	// Add this edge, if it increases the BIC score
	if (!check_interrupt() && diffScore_opt > _minScoreDiff) {
		dout.level(2) << "  Adding edge " << u_opt << " --> " << v_opt << std::endl;
		addEdge(u_opt, v_opt);
		return true;
	}
	else
		return false;
}

bool EssentialGraph::greedyDAGBackward()
{
	uint u_opt = 0, v_opt = 0;
	double diffScore, diffScore_opt;
	std::set<uint> parents, C_new;
	std::set<uint>::iterator ui;

	dout.level(2) << "= Starting backward step...\n";

	// Initialize help variables
	diffScore_opt = _minScoreDiff;
	uint p = getVertexCount();

	// Find edge that maximally increases BIC score when removed
	for (uint v = 0; v < p; ++v) {
		parents = getParents(v);
		for (ui = parents.begin(); ui != parents.end(); ++ui) {
			// Calculate BIC score difference when removing edge (u, v)
			C_new = parents;
			diffScore = - _score->local(v, C_new);
			C_new.erase(*ui);
			diffScore += _score->local(v, C_new);
			dout.level(3) << "  Score diff. for edge " << *ui << " --> " << v << " : " <<
					diffScore << std::endl;

			// If new score is better than previous optimum, store (u, v, C) as new optimum
			if (diffScore > diffScore_opt) {
				u_opt = *ui;
				v_opt = v;
				diffScore_opt = diffScore;
			}
		}
	}

	// Remove this edge, if this incrases the BIC score
	if (!check_interrupt() && diffScore_opt > _minScoreDiff) {
		dout.level(2) << "  Removing edge " << u_opt << " --> " << v_opt << std::endl;
		removeEdge(u_opt, v_opt);
		return true;
	}
	else
		return false;
}

bool EssentialGraph::greedyDAGTurn()
{
	uint u_opt = 0, v_opt = 0;
	double diffScore, diffScore_opt;
	std::set<uint> parents, C_new, D_new, emptyset;
	std::set<uint>::iterator ui;

	dout.level(2) << "= Starting turning step...\n";

	// Initialize help variables
	diffScore_opt = _minScoreDiff;
	uint p = getVertexCount();

	// Find edge that maximally increases BIC score when turned, i.e. when its
	// orientation is changed
	for (uint v = 0; v < p; ++v) {
		parents = getParents(v);
		for (ui = parents.begin(); ui != parents.end(); ++ui) {
			if (!existsPath(*ui, v)) {
				C_new = parents;
				D_new = getParents(*ui);
				diffScore = - _score->local(v, C_new) - _score->local(*ui, D_new);
				C_new.erase(*ui);
				D_new.insert(v);
				diffScore += _score->local(v, C_new) + _score->local(*ui, D_new);
				dout.level(3) << "  Score diff. for edge " << *ui << " --> " << v << " : " <<
						diffScore << std::endl;

				// If new score is better than previous optimum, store (u, v, C) as new optimum
				if (diffScore > diffScore_opt) {
					u_opt = *ui;
					v_opt = v;
					diffScore_opt = diffScore;
				}
			}
		}
	}

	// Turn this edge, if this incrases the BIC score
	if (!check_interrupt() && diffScore_opt > _minScoreDiff) {
		dout.level(2) << "  Turning edge " << u_opt << " --> " << v_opt << std::endl;
		removeEdge(u_opt, v_opt);
		addEdge(v_opt, u_opt);
		return true;
	}
	else
		return false;
}

bool EssentialGraph::greedyDAGStepDir(const step_dir direction)
{
	switch (direction) {
	case SD_FORWARD:
		return greedyDAGForward();

	case SD_BACKWARD:
		return greedyDAGBackward();

	case SD_TURNING:
		return greedyDAGTurn();

	default:
		return false;
	} // SWITCH direction

}

void EssentialGraph::siMySearch()
{
	// Check whether DAG is not too large (practically, the algorithm will
	// only work for substantially smaller graphs than p = 32 due to
	// memory and run time limits)
	if (getVertexCount() > 31)
		throw std::length_error("Vertex count must not exceed 31.");

	uint32_t subset, subsubset, pattern;
	bool unset;

	// TODO check for user interrupts

	// Table of best parents given a variable and a set of candidate parents
	std::vector< std::vector<uint32_t> > bestParents(getVertexCount(), std::vector<uint32_t>(1 << (getVertexCount() - 1)));
	// Table with corresponding best local scores
	std::vector< std::vector<double> > localScore(getVertexCount(), std::vector<double>(1 << (getVertexCount() - 1)));
	// Table of best sinks for all subsets of variables
	std::vector< int > bestSink(1 << getVertexCount());
	// Table of corresponding optimal scores
	std::vector< double > totalScore(1 << getVertexCount(), 0.);

	// Forward phase of DP: fill in tables of best parents and corresponding local scores
	for (uint i = 0; i < getVertexCount(); i++) {
		//std::cout << "\n" << i << ":\t";

		for (subset = 0; subset < bestParents[i].size(); subset++) {
			localScore[i][subset] = _score->local(i, _bitsToParents(i, subset));
			//std::cout << "localscore: " << localScore[i][subset] << "\n";
			bestParents[i][subset] = subset;
			for (pattern = 1; pattern < bestParents[i].size(); pattern *= 2) {
				subsubset = subset & (~pattern);
				//std::cout << "subsubset: " << subsubset << "\n";
				if (localScore[i][subsubset] > localScore[i][subset]) {
					localScore[i][subset] = localScore[i][subsubset];
					bestParents[i][subset] = bestParents[i][subsubset];
				}
			}

			//std::cout << "(" << bestParents[i][subset] << ", " << localScore[i][subset] << ")\t";
		}
	}


	// Fill in tables of best sinks and corresponding total scores
	//std::cout << "\n\nBest sinks:\n";
	for (subset = 1; subset < totalScore.size(); subset++) {
		unset = true;
		for (uint i = 0; i < getVertexCount(); i++) {
			if ((1 << i) & subset) {
				// Calculate subsubset w.r.t. candidate set for parents, not w.r.t. full set of variables
				pattern = (1 << i) - 1;
				subsubset = (subset & pattern) | ((subset & (~pattern << 1)) >> 1);
				//std::cout << "subset: " << subset << ";\ti: " << i << ";\tpattern: " << pattern << ";\tsubsubset: " << subsubset << "\n";
				if (unset || totalScore[subset] < localScore[i][subsubset] + totalScore[subset & ~(1 << i)]) {
					totalScore[subset] = localScore[i][subsubset] + totalScore[subset & ~(1 << i)];
					bestSink[subset] = i;
					unset = false;
				}
			}
		}

		//std::cout << subset << ": " << bestSink[subset] << ", " << totalScore[subset] << "\n";
	}

	// Backward phase of DP: reconstruct optimal DAG
	clear();
	std::set<uint> parents;
	std::set<uint>::iterator pi;
	subset = (1 << getVertexCount()) - 1;
	uint i;
	while (subset) {
		i = bestSink[subset];

		// Calculate subsubset w.r.t. candidate set for parents, not w.r.t. full set of variables
		pattern = (1 << i) - 1;
		subsubset = (subset & pattern) | ((subset & (~pattern << 1)) >> 1);
		//std::cout << "subset: " << subset << ";\tbest sink: " << i << ";\tsubsubset: " << subsubset << "\n";

		parents = _bitsToParents(i, bestParents[i][subsubset]);
		for (pi = parents.begin(); pi != parents.end(); ++pi)
			addEdge(*pi, i);

		subset = subset & ~(1 << i);
	}
}


std::set<uint> EssentialGraph::getOptimalTarget(uint maxSize)
{
	std::set<uint> target;

	if (maxSize == 0)
		return target;
	else if (maxSize == 1) {
		uint v = _getOptimalSingleVertexTarget();
		if (v < getVertexCount())
			target.insert(v);
		return target;
	}
	else if (maxSize == getVertexCount())
		return _getOptimalUnrestrTarget();
	else
		throw std::runtime_error("Optimal targets with size other than 1 or p are not supported.");
}

EssentialGraph castGraph(const SEXP argInEdges)
{
	Rcpp::List listInEdges(argInEdges);
	EssentialGraph result(listInEdges.size());

	for (R_len_t i = 0; i < listInEdges.size(); ++i) {
		Rcpp::IntegerVector vecParents((SEXP)(listInEdges[i]));
		// Adapt indices to C++ convention
		for (Rcpp::IntegerVector::iterator vi = vecParents.begin(); vi != vecParents.end(); ++vi)
			result.addEdge(*vi - 1, i);
	}

	return result;
}

Rcpp::List wrapGraph(const EssentialGraph& graph)
{
	Rcpp::List result;
	Rcpp::IntegerVector vecEdges;
	std::set<uint> edges;

	for (uint i = 0; i < graph.getVertexCount(); ++i) {
		edges = graph.getInEdges(i);
		Rcpp::IntegerVector vecEdges(edges.begin(), edges.end());
		// Adapt indices to R convention
		for (R_len_t i = 0; i < vecEdges.size(); ++i)
			vecEdges[i]++;
		result.push_back(vecEdges);
	}

	return result;
}


