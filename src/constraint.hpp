/**
 * Auxiliary methods for constraint-based algorithms
 *
 * @author Alain Hauser <alain.hauser@biology.unibe.ch>
 * $Id: $
 */

#include "armaLapack.hpp"
#include <vector>
#include <boost/graph/adjacency_list.hpp>

// Define BGL class for undirected graph
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> InternalUndirectedGraph;
typedef boost::graph_traits<InternalUndirectedGraph>::out_edge_iterator UndirOutEdgeIter;
typedef boost::graph_traits<InternalUndirectedGraph>::edge_iterator UndirEdgeIter;

// Define type for specification of separation sets
typedef std::vector<std::vector<arma::ivec > > SepSets;

/**
 * Virtual base class for conditional independence tests
 */
class IndepTest
{
public:
	/**
	 * Virtual test function: must be implemented in all derived classes
	 *
	 * @param	u	First variable index
	 * @param	v	Second variable index
	 * @param	S	Conditioning set
	 * @return	p-value for independence test of X_u and X_v given X_S
	 */
	virtual double test(uint u, uint v, std::vector<uint> S) const = 0;
};

/**
 * Conditional independence test based on R function
 */
class IndepTestRFunction : public IndepTest
{
protected:
	/**
	 * Sufficient statistic
	 */
	Rcpp::List* _suffStat;

	/**
	 * R function implementing the independence test
	 */
	Rcpp::Function _testFunction;

public:
	IndepTestRFunction(Rcpp::List* suffStat, Rcpp::Function testFunction) :
		_suffStat(suffStat),
		_testFunction(testFunction) {};

	virtual double test(uint u, uint v, std::vector<uint> S) const;
};

/**
 * Conditional independence test for Gaussian data
 */
class IndepTestGauss : public IndepTest
{
protected:
	/**
	 * Sufficient statistic for easier access in independence test
	 */
	uint _sampleSize;
	arma::mat _correlation;

public:
	IndepTestGauss(uint sampleSize, Rcpp::NumericMatrix& cor) :
		_sampleSize(sampleSize),
		_correlation(cor.begin(), cor.nrow(), cor.ncol(), false) {}

	virtual double test(uint u, uint v, std::vector<uint> S) const;
};

/**
 * Skeleton of a DAG model
 */
class Skeleton
{
protected:
	/**
	 * Boost graph internally storing the graph structure and the fixed edges
	 */
	InternalUndirectedGraph _graph;
	InternalUndirectedGraph _fixedEdges;

	/**
	 * Pointer to a conditional independence test
	 */
	IndepTest* _indepTest;

public:
	/**
	 * Constructor
	 */
	Skeleton(const uint vertexCount) : _graph(vertexCount), _fixedEdges(vertexCount) {};

	/**
	 * Adds an edge to the graph
	 */
	void addEdge(const uint a, const uint b) { boost::add_edge(a, b, _graph); }

	/**
	 * Adds a fixed edge to the graph.  Also adds the edge itself
	 */
	void addFixedEdge(const uint a, const uint b);

	/**
	 * Checks whether an edge is fixed
	 */
	bool isFixed(const uint a, const uint b) const;

	/**
	 * Removes an edge from the graph, if it is not fixed
	 */
	void removeEdge(const uint a, const uint b);

	/**
	 * Number of vertices
	 */
	uint getVertexCount() const { return boost::num_vertices(_graph); }

	/**
	 * Number of edges
	 */
	uint getEdgeCount() const { return boost::num_edges(_graph); }

	/**
	 * Checks whether two vertices are adjacent
	 */
	bool hasEdge(const uint a, const uint b) const;

	/**
	 * Yields the neighbours of a certain vertex
	 */
	std::set<uint> getNeighbors(const uint vertex) const;

	/**
	 * Yields the degree (= number of adjacent vertices) of some vertex
	 */
	uint getDegree(const uint vertex) const { return boost::out_degree(vertex, _graph); }

	/**
	 * Yields the adjacency matrix
	 */
	Rcpp::LogicalMatrix getAdjacencyMatrix();

	/**
	 * Sets and gets the independence test object
	 */
	void setIndepTest(IndepTest* indepTest) { _indepTest = indepTest; }
	IndepTest* getIndepTest() const { return _indepTest; }

	/**
	 * Reduces the skeleton (removes edges) until it represents a given conditional
	 * independence structure; first step of the PC algorithm
	 *
	 * TODO: return things needed by R function "skeleton"; create arguments with options
	 */
	void fitCondInd(
			const double alpha,
			Rcpp::NumericMatrix& pMax,
			SepSets& sepSet,
			std::vector<int>& edgeTests,
			int maxCondSize = -1,
			const bool NAdelete = true);
};
