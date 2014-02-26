/*
 * Classes for greedy estimation of causal structures
 *
 * @author Alain Hauser
 * $Id: greedy.hpp 219 2014-01-30 15:57:29Z alhauser $
 */

#ifndef GREEDY_HPP_
#define GREEDY_HPP_

#include "score.hpp"

#include <vector>
#include <set>
#include <utility>
#include <list>
#include <deque>
#include <boost/graph/adjacency_list.hpp>
#include <boost/dynamic_bitset.hpp>

enum edge_flag { NOT_PROTECTED, UNDECIDABLE, PROTECTED };

/**
 * Help functions for easier handling of set operations
 */
template <typename Key, typename Compare, typename Alloc> std::set<Key, Compare, Alloc> set_intersection(const std::set<Key, Compare, Alloc>& set1, const std::set<Key, Compare, Alloc>& set2)
{
	typename std::set<Key, Compare, Alloc> result;
	Compare comp;
	std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(result, result.end()), comp);
	return result;
}

template <typename Key, typename Compare, typename Alloc> std::set<Key, Compare, Alloc> set_union(const std::set<Key, Compare, Alloc>& set1, const std::set<Key, Compare, Alloc>& set2)
{
	typename std::set<Key, Compare, Alloc> result;
	Compare comp;
	std::set_union(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(result, result.end()), comp);
	return result;
}

template <typename Key, typename Compare, typename Alloc> std::set<Key, Compare, Alloc> set_difference(const std::set<Key, Compare, Alloc>& set1, const std::set<Key, Compare, Alloc>& set2)
{
	typename std::set<Key, Compare, Alloc> result;
	Compare comp;
	std::set_difference(set1.begin(), set1.end(), set2.begin(), set2.end(), std::inserter(result, result.end()), comp);
	return result;
}

/**
 * Classes for internal representation of graphs (specialized boost classes)
 */
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS> InternalEssentialGraph;
typedef boost::graph_traits<InternalEssentialGraph>::vertex_descriptor vertex_t;
typedef boost::graph_traits<InternalEssentialGraph>::vertex_iterator VertexIter;
typedef boost::graph_traits<InternalEssentialGraph>::adjacency_iterator AdjacencyIter;
typedef boost::graph_traits<InternalEssentialGraph>::edge_descriptor edge_t;
typedef boost::graph_traits<InternalEssentialGraph>::in_edge_iterator InEdgeIter;
typedef boost::graph_traits<InternalEssentialGraph>::out_edge_iterator OutEdgeIter;

// Edge type for internal use... the default one (edge_descriptor) is fucking complicated.
struct Edge
{
	Edge() : source(0), target(0) {};

	Edge(const uint s, const uint t) : source(s), target(t) {};

	Edge(const edge_t& e, const InternalEssentialGraph& g) : source(boost::source(e, g)), target(boost::target(e, g)) {};

	uint source, target;
};

/**
 * Comparator that yields an lexicographic ordering of edges with _inverse_
 * priority (first target, then source)
 */
struct EdgeCmp : public std::binary_function<Edge, Edge, bool>
{
	bool operator()(const Edge& first, const Edge& second) const
	{
		return first.target < second.target ||
				(first.target == second.target && first.source < second.source);
	}
};

/**
 * Helper class used as a stack for candidate cliques C \subset N
 */
class CliqueStack : public std::deque<std::set<uint> >
{
public:
	std::set<std::set<uint> > stop_sets;

	bool append(const std::set<uint>& newSet)
	{
		bool inStopSet = false;
		for (std::set<std::set<uint> >::iterator iter = stop_sets.begin(); !inStopSet && iter != stop_sets.end(); ++iter)
			inStopSet = std::includes(iter->begin(), iter->end(), newSet.begin(), newSet.end());
		//	inStopSet = (*iter == newSet);
		if (!inStopSet)
			push_back(newSet);
		return !inStopSet;
	}

	void clear_all()
	{
		clear();
		stop_sets.clear();
	}
};

/**
 * Helper classes for storing cached values
 */
struct ArrowChange
{
	uint source;
	std::set<uint> clique;
	double score;
};

struct ArrowChangeCmp : public std::binary_function<Edge, Edge, bool>
{
	bool operator()(const ArrowChange& first, const ArrowChange& second) const
	{
		return (first.score < second.score);
	}
};

enum step_dir { SD_NONE, SD_FORWARD, SD_BACKWARD, SD_TURNING };

// Forward declaration for testing
class EssentialGraphTest;
class BICScoreTest;

/**
 * Basic graph class. Support for directed and undirected edges, no loops.
 */
class EssentialGraph
{
	friend class EssentialGraphTest;
	friend class BICScoreTest;
protected:
	/**
	 * Boost graph internally storing the graph structure
	 */
	InternalEssentialGraph _graph;

	/**
	 * Fixed gaps: (undirected) graph representing edges that must not be filled
	 * in the essential graph.
	 *
	 * When only a sparse graph is allowed, _fixedGaps would be dense; then, it
	 * is more efficient to store the allowed edges than the fixed gaps. This is
	 * indicated by the flag _gapsInverted.
	 *
	 * TODO Allow for fixed edges or directed gaps; all in all, more complex
	 * restrictions
	 */
	InternalEssentialGraph _fixedGaps;
	bool _gapsInverted;

	/**
	 * Indicates whether optimal cliques and corresponding score differences
	 * should be cached during greedy search.
	 */
	bool _doCaching;

	/**
	 * Indicates whether the cache must be initialized before usage.
	 */
	step_dir _actualPhase;

	/**
	 * Map of potential edges for which cached values exist
	 */
	std::vector<ArrowChange> _scoreCache;

	/**
	 * Pointer to scoring object
	 */
	Score* _score;

	/**
	 * Constant defining minimal score difference.
	 *
	 * "Very small" score differences are often due to rounding errors and
	 * involves the danger of infinite loops
	 */
	static double _minScoreDiff;

	/**
	 * Pointer to object representing family of targets
	 */
	TargetFamily* _targets;

	/**
	 * Maximum vertex degrees, per vertex
	 */
	std::vector<uint> _maxVertexDegree;
	
	/**
	 * Vertices which are only allowed to have children, but no parents
	 * 
	 * NOTE: in order that this makes sense (i.e., is consistent with Markov
	 * equivalence classes), the corresponding vertices should also appear alone
	 * in intervention targets. However, this is not checked in the algorithm...
	 */
	boost::dynamic_bitset<> _childrenOnly; 

	/**
	 * Checks whether there is a fixed gap between two vertices.
	 */
	bool gapFixed(const uint a, const uint b) const;

	/**
	 * Checks whether there is a path from a to b in the graph that does not
	 * go through the vertices of C. The edge (a, b) is not considered, if it
	 * exists.
	 *
	 * @param	undirected	indicates whether only undirected edges shall be followed
	 */
	bool existsPath(const uint a, const uint b, const std::set<uint>& C = std::set<uint>(), const bool undirected = false);

	/**
	 * Checks whether there is a path from a to some vertex in B in the graph that
	 * does not go through the vertices in C.
	 */
	bool existsPath(const uint a, const std::set<uint>& B, const std::set<uint>& C = std::set<uint>(), const bool undirected = false) const;

	/**
	 * Checks whether there is a path from a to some vertex in B in the subgraph
	 * induced by the vertex set C
	 */
	bool existsPath(const std::set<uint>& C, const uint a, const std::set<uint>& B);

	/**
	 * Yields a LexBFS-ordering of a subset of vertices, and possibly orients
	 * the edges of the induced subgraph accordingly. Assumes
	 * that all vertices belong to the same chain component.
	 *
	 * @param	first 		first vertex of the start order
	 * @param 	last		--last: last vertex of the start order
	 * @param	orient		indicates whether edges have to be oriented according
	 * 						to the LexBFS order
	 * @param	directed	(OUT) pointer to a list of edges that become directed
	 * @return  			List (set) of oriented edges
	 */
	template <typename InputIterator> std::vector<uint> lexBFS(InputIterator first, InputIterator last, const bool orient = false, std::set<Edge, EdgeCmp>* directed = NULL)
	{
		if (directed != NULL)
			directed->clear();
		std::vector<uint> ordering;
		int length = std::distance(first, last);
		ordering.reserve(length);

		// Trivial cases: if not more than one start vertex is provided,
		// return an empty set of edges
		if (length == 1)
			ordering.push_back(*first);
		if (length <= 1)
			return ordering;

		// Create sequence of sets ("\Sigma") containing the single set
		// of all vertices in the given start order
		std::list<std::list<uint> > sets(1, std::list<uint>(first, last));
		std::list<std::list<uint> >::iterator si, newSet;
		std::list<uint>::iterator vi;

		uint a;

		while (!sets.empty()) {
			// Remove the first vertex from the first set, and remove this set
			// if it becomes empty
			a = sets.front().front();
			sets.front().pop_front();
			if (sets.front().empty())
				sets.pop_front();

			// Append a to the ordering
			ordering.push_back(a);

			// Move all neighbors of a into own sets, and orient visited edges
			// away from a
			for (si = sets.begin(); si != sets.end(); ) {
				newSet = sets.insert(si, std::list<uint>());
				for (vi = si->begin(); vi != si->end(); ) {
					if (hasEdge(a, *vi)) {
						// Orient edge to neighboring vertex, if requested, and
						// store oriented edge in return set
						if (orient)
							removeEdge(*vi, a);
						if (directed != NULL)
							directed->insert(Edge(a, *vi));

						// Move neighoring vertex
						newSet->push_back(*vi);
						vi = si->erase(vi);
					}
					else ++vi;
				}

				// If visited or newly inserted sets are empty, remove them
				// from the sequence
				if (newSet->empty())
					sets.erase(newSet);
				if (si->empty())
					si = sets.erase(si);
				else ++si;
			}
		}

		return ordering;
	}

	/**
	 * Yields a greedy coloring of a subset of vertices
	 */
	std::vector<uint> greedyColoring(std::vector<uint> vertices);

	/**
	 * Find all maximal cliques in an induced subgraph of some chain component
	 *
	 * NOTE: the function does not check whether the provided range of vertices
	 * really is a subset of some chain component!
	 */
	template <typename InputIterator> std::vector<std::set<uint> > getMaxCliques(InputIterator first, InputIterator last)
	{
		std::vector<std::set<uint> > maxCliques;

		// Trivial case: range of vertices contains at most one vertex
		if (std::distance(first, last) <= 1) {
			maxCliques.push_back(std::set<uint>(first, last));
			return maxCliques;
		}

		// For less trivial cases, first generate a LexBFS-ordering on the provided range of vertices
		std::vector<uint> ordering = lexBFS(first, last);

		// Find maximal cliques using the LexBFS-ordering
		std::set<uint> nbhdSubset(first, last);
		std::set<uint> vertices, C;
		std::vector<std::set<uint> >::iterator cliqueIter;
		bool included;
		for (int i = ordering.size() - 1; i >= 0; --i) {
			nbhdSubset.erase(ordering[i]);
			vertices = getNeighbors(ordering[i]);
			C = set_intersection(vertices, nbhdSubset);
			C.insert(ordering[i]);
			included = false;
			for (cliqueIter = maxCliques.begin(); !included && cliqueIter != maxCliques.end(); ++cliqueIter)
				included = std::includes(cliqueIter->begin(), cliqueIter->end(), C.begin(), C.end());
			if (!included)
				maxCliques.push_back(C);
		}

		return maxCliques;
	}

	/**
	 * Calculates the optimal arrow insertion for a given vertex v, that is,
	 * the best source u, clique C and the corresponding score difference.
	 */
	ArrowChange getOptimalArrowInsertion(const uint v);

	/**
	 * Calculates the optimal arrow deletion for a given vertex v, that is,
	 * the best source u, clique C and the corresponding score difference.
	 */
	ArrowChange getOptimalArrowDeletion(const uint v);

	/**
	 * Calculates the optimal arrow turning for a given vertex v, that is,
	 * the best source u, clique C and the corresponding score difference.
	 */
	ArrowChange getOptimalArrowTurning(const uint v);

	/**
	 * Yields the parent set of a node given its representation as unsigned integer.
	 *
	 * Help function for maximization of BIC via DP.
	 */
	std::set<uint> _bitsToParents(const int vertex, const uint32_t bits);

	/**
	 * Yields the "optimal" intervention target (without restriction on
	 * target size)
	 */
	std::set<uint> _getOptimalUnrestrTarget();

	/**
	 * Yields the "optimal" single-vertex intervention target
	 */
	uint _getOptimalSingleVertexTarget();
public:
	/**
	 * Constructors
	 */
	EssentialGraph() {};
	EssentialGraph(const uint vertexCount);
	// Graph(const uint vertexCount, const t_adjmat adjacency);

	/**
	 * Removes all edges from the graph
	 */
	void clear();

	/**
	 * Adds an edge to the graph. Purely technical function, does
	 * not check whether the graph is still a chain graph after the
	 * insertion.
	 */
	void addEdge(const uint a, const uint b, bool undirected = false);

	/**
	 * Removes an edge from the graph. Purely technical function, as
	 * addEdge().
	 */
	void removeEdge(const uint a, const uint b, bool bothDirections = false);

	/**
	 * Number of vertices
	 */
	uint getVertexCount() const { return boost::num_vertices(_graph); }

	/**
	 * Number of edges
	 */
	uint getEdgeCount() const { return boost::num_edges(_graph); }

	/**
	 * Checks whether some vertex is a parent, neighbor, ... of the second vertex
	 */
	bool hasEdge(const uint a, const uint b) const;
	bool isParent(const uint a, const uint b) const { return hasEdge(a, b) && !hasEdge(b, a); }
	bool isNeighbor(const uint a, const uint b) const { return hasEdge(a, b) && hasEdge(b, a); }
	bool isAdjacent(const uint a, const uint b) const { return hasEdge(a, b) || hasEdge(b, a); }

	/**
	 * Yields the parents, neighbours, ... of a certain vertex
	 */
	std::set<uint> getParents(const uint vertex) const;
	std::set<uint> getChildren(const uint vertex) const;
	std::set<uint> getNeighbors(const uint vertex) const;
	std::set<uint> getAdjacent(const uint vertex) const;
	std::set<uint> getInEdges(const uint vertex) const;

	/**
	 * Yields the degree (= number of adjacent vertices) of some vertex
	 */
	uint getDegree(const uint vertex) const;

	/**
	 * Yields the induced subgraph on some vertex subset of the current
	 * graph.
	 *
	 * The subgraph is copied, i.e., changes on it do not affect the
	 * original graph. The vertices are re-indexed, i.e., the first vertex
	 * in the selected range gets the new index 0, the second one the
	 * index 1, etc.
	 */
	template <typename InputIterator>
	EssentialGraph inducedSubgraph(InputIterator first, InputIterator last) const
	{
		EssentialGraph result(std::distance(first, last));

		InputIterator vi, vj;
		uint i, j;
		i = 0;
		for (vi = first; vi != last; ++vi, ++i) {
			j = 0;
			for (vj = first; vj != last; ++vj, ++j)
				if (hasEdge(*vi, *vj))
					result.addEdge(i, j);
		}

		return result;
	}

	/**
	 * Sets the fixed gaps
	 */
	void setFixedGaps(const EssentialGraph& fixedGaps, const bool inverted);

	/**
	 * Sets the maximum vertex degrees
	 */
	void limitVertexDegree(const std::vector<uint>& maxVertexDegree);
	void limitVertexDegree(const uint maxVertexDegree);
	void limitVertexDegree(const double maxRelativeDegree);
	
	/**
	 * Allow certain vertices to have only children, but no parents
	 */
	bool getChildrenOnly(const uint vertex) const { return _childrenOnly.test(vertex); }
	void setChildrenOnly(const uint vertex, const bool setting) { _childrenOnly.set(vertex, setting); }

	/**
	 * Yields the anterior set and posterior set of a set of vertices
	 */
	boost::dynamic_bitset<> getAnteriorSet(const std::set<uint>& A);
	boost::dynamic_bitset<> getPosteriorSet(const std::set<uint>& A);

	/**
	 * Finds the chain component of a certain vertex v.
	 */
	std::set<uint> getChainComponent(const uint v) const;

	/**
	 * Sets and gets score object
	 */
	void setScore(Score* score) { _score = score; }
	Score* getScore() const { return _score; }

	/**
	 * Sets and gets the family of targets
	 */
	void setTargets(TargetFamily* targets) { _targets = targets; }
	TargetFamily* getTargets() { return _targets; }

	/**
	 * Successively replace unprotected arrows by lines.
	 *
	 * @return	Set of previously directed edges that are now undirected
	 */
	std::set<Edge, EdgeCmp> replaceUnprotected();

	/**
	 * Getter and setter for adjacency matrix
	 */
//	arma::umat getAdjacencyMatrix() const;
//	template <typename eT> void setAdjacencyMatrix(const arma::Mat<eT> adjacency)
//	{
//		uint a, b;
//		clear();
//		for (a = 0; a < getVertexCount(); ++a)
//			for (b = 0; b < getVertexCount(); ++b)
//				if (adjacency(a, b))
//					addEdge(a, b);
//	}

	/**
	 * Yields a representative (DAG) of the equivalence class
	 *
	 * TODO: perhaps change the class of the result...
	 */
	EssentialGraph getRepresentative() const;

	/**
	 * Yields all representatives of the equivalence class
	 */
	std::vector<boost::dynamic_bitset<> > getAllRepresentatives() const;

	/**
	 * Enable caching.
	 *
	 * Sets the corresponding flag and creates a list for cached entries
	 */
	void enableCaching();

	/**
	 * Disable caching.
	 *
	 * Sets the corresponding flag and deletes the list for cached entries
	 */
	void disableCaching();

	/**
	 * Inserts a new edge according to the Ins-operator.
	 */
	void insert(const uint u, const uint v, const std::set<uint> C);

	/**
	 * Deletes an edge according to the Del-operator.
	 */
	void remove(const uint u, const uint v, const std::set<uint> C);

	/**
	 * Turns an edge according to the Turn-operator.
	 */
	void turn(const uint u, const uint v, const std::set<uint> C);

	/**
	 * Does one forward step of the greedy interventional equivalence search.
	 */
	bool greedyForward();

	/**
	 * Does one backward step of the greedy interventional equivalence search
	 */
	bool greedyBackward();

	/**
	 * Does one turning step of the greedy interventional equivalence search
	 */
	bool greedyTurn();

	/**
	 * Does one greedy step, either forward, backward, or turning, the one that
	 * yields the highest score gain.
	 */
	step_dir greedyStep();

	/**
	 * Does one forward step of the greedy DAG search (i.e., not over interventional
	 * equivalence classes).
	 */
	bool greedyDAGForward();

	/**
	 * Does one backward step of the greedy DAG search (i.e., not over interventional
	 * equivalence classes).
	 */
	bool greedyDAGBackward();

	/**
	 * Does one turning step of the greedy DAG search (i.e., not over interventional
	 * equivalence classes).
	 */
	bool greedyDAGTurn();

	/**
	 * Maximizes the BIC score by dynamic programming, as proposed by
	 * Silander and Myllymäki (2006). Only works for small graphs
	 * (technically, p <= 32; practically less due to memory and time
	 * constraints)
	 */
	void siMySearch();

	/**
	 * Yields the "optimal" intervention target of size <= q, i.e., the
	 * intervention target for which the worst case number of orientable
	 * arrows is maximal.
	 *
	 * At the moment, only q = p (no restriction on target size) and
	 * q = 1 are allowed.
	 */
	std::set<uint> getOptimalTarget(uint maxSize);
};

#endif /* GREEDY_HPP_ */
