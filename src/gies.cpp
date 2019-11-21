/**
 * Main file of the Greedy Interventional Equivalence Search library for R
 *
 * @author Alain Hauser
 * $Id: gies.cpp 500 2019-11-20 13:46:14Z alhauser $
 */

#include <vector>
#include <string>
#include <algorithm>
// #include <boost/lambda/lambda.hpp>
#include <boost/graph/adjacency_list.hpp>
// Experimental support for OpenMP; aim: parallelize more and more functions...
#ifdef _OPENMP
#include <omp.h>
#endif

// Define BGL class for undirected graph
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> UndirectedGraph;

#include "pcalg/constraint.hpp"
#include "pcalg/score.hpp"
#include "pcalg/greedy.hpp"
#define DEFINE_GLOBAL_DEBUG_STREAM
#include "pcalg/gies_debug.hpp"

// using namespace boost::lambda;


/**
 * Yields the local score of a vertex given its parents.
 *
 * @param	argScore		name of the score
 * @param	argPreprocData	preprocessed data; sufficient statistic and all
 * 							parameters characterizing the score to be calculated
 * @param 	argVertex		vertex index
 * @param 	argParents		vector of parents of vertex
 * @param	argOptions		additional options; at the moment: DEBUG.LEVEL
 * @return	local score value
 */
RcppExport SEXP localScore(
		SEXP argScore,
		SEXP argPreprocData,
		SEXP argVertex,
		SEXP argParents,
		SEXP argOptions)
{
	// Initialize automatic exception handling; manual one does not work any more...
	BEGIN_RCPP

	// Set debug level
	Rcpp::List options(argOptions);
	dout.setLevel(Rcpp::as<int>(options["DEBUG.LEVEL"]));
	dout.level(1) << "Calculating local score...\n";

	// Create appropriate scoring object
	Rcpp::List data(argPreprocData);
	TargetFamily targets = castTargets(data["targets"]);
	dout.level(3) << "# intervention targets: " << targets.size() << "\n";
	Score* score = createScore(Rcpp::as<std::string>(argScore), &targets, data);

	// Calculate local score and delete score object
	double result = score->local(Rcpp::as<uint>(argVertex) - 1, castVertices(argParents));
	delete score;
	return Rcpp::wrap(result);

	END_RCPP
}

/**
 * Yields the global score of a DAG.
 *
 * @param	argScore		name of the score
 * @param	argPreprocData	preprocessed data; sufficient statistic and all
 * 							parameters characterizing the score to be calculated
 * @param 	argInEdges		list of in-edges characterizing the DAG
 * @param	argOptions		additional options; at the moment: DEBUG.LEVEL
 * @return	global score value
 */
RcppExport SEXP globalScore(
		SEXP argScore,
		SEXP argPreprocData,
		SEXP argInEdges,
		SEXP argOptions)
{
	// Initialize automatic exception handling; manual one does not work any more...
	BEGIN_RCPP

	// Set debug level
	Rcpp::List options(argOptions);
	dout.setLevel(Rcpp::as<int>(options["DEBUG.LEVEL"]));

	// Create appropriate scoring object
	Rcpp::List data(argPreprocData);
	TargetFamily targets = castTargets(data["targets"]);
	Score* score = createScore(Rcpp::as<std::string>(argScore), &targets, data);

	// Calculate local score
	double result = score->global(castGraph(argInEdges));
	// TODO: check why this leads to a segfault!!!!
	delete score;
	return Rcpp::wrap(result);

	END_RCPP
}

/**
 * Yields the local MLE of a vertex given its parents.
 *
 * @param	argScore		name of the score
 * @param	argPreprocData	preprocessed data; sufficient statistic and all
 * 							parameters characterizing the score to be calculated
 * @param 	argVertex		vertex index
 * @param 	argParents		vector of parents of vertex
 * @param	argOptions		additional options; at the moment: DEBUG.LEVEL
 * @return	vector of local MLE
 */
RcppExport SEXP localMLE(
		SEXP argScore,
		SEXP argPreprocData,
		SEXP argVertex,
		SEXP argParents,
		SEXP argOptions)
{
	// Initialize automatic exception handling; manual one does not work any more...
	BEGIN_RCPP

	// Set debug level
	Rcpp::List options(argOptions);
	dout.setLevel(Rcpp::as<int>(options["DEBUG.LEVEL"]));

	// Create appropriate scoring object
	Rcpp::List data(argPreprocData);
	TargetFamily targets = castTargets(data["targets"]);
	Score* score = createScore(Rcpp::as<std::string>(argScore), &targets, data);

	// Calculate local score
	std::vector<double> result = score->localMLE(Rcpp::as<uint>(argVertex) - 1, castVertices(argParents));
	delete score;
	return Rcpp::wrap(result);

	END_RCPP
}

/**
 * Yields the global MLE of a DAG.
 *
 * @param	argScore		name of the score
 * @param	argPreprocData	preprocessed data; sufficient statistic and all
 * 							parameters characterizing the score to be calculated
 * @param 	argInEdges		list of in-edges characterizing the DAG
 * @param	argOptions		additional options; at the moment: DEBUG.LEVEL
 * @return	list of MLE vectors
 */
RcppExport SEXP globalMLE(
		SEXP argScore,
		SEXP argPreprocData,
		SEXP argInEdges,
		SEXP argOptions)
{
	// Initialize automatic exception handling; manual one does not work any more...
	BEGIN_RCPP

	// Set debug level
	Rcpp::List options(argOptions);
	dout.setLevel(Rcpp::as<int>(options["DEBUG.LEVEL"]));

	// Create appropriate scoring object
	Rcpp::List data(argPreprocData);
	TargetFamily targets = castTargets(data["targets"]);
	Score* score = createScore(Rcpp::as<std::string>(argScore), &targets, data);

	// Calculate global score
	std::vector<std::vector<double> > result = score->globalMLE(castGraph(argInEdges));
	delete score;
	return Rcpp::wrap(result);

	END_RCPP
}

/**
 * Interface to a variety of causal inference algorithms.
 *
 * @param 	argGraph			list of in-edges representing the current
 * 								essential graph
 * @param	argPreprocData		preprocessed data; sufficient statistic and all
 * 								parameters characterizing the score to be calculated
 * @param	argAlgorithm		string indicating the causal inference algorithm
 * 								to be used. Supported options: "GIES", "GDS", "DP"
 * @param	argScore			name of score object to be used. Currently supported:
 * 								"none" (= R object), "gauss.l0pen"
 * @param	argOptions			list of options specific for desired inference algorithm
 */
RcppExport SEXP causalInference(
		SEXP argGraph,
		SEXP argPreprocData,
		SEXP argAlgorithm,
		SEXP argScore,
		SEXP argOptions)
{
	// Initialize automatic exception handling; manual one does not work any more...
	BEGIN_RCPP

	// Cast debug level from options
	Rcpp::List options(argOptions);
	dout.setLevel(Rcpp::as<int>(options["DEBUG.LEVEL"]));

	// Cast graph
	dout.level(1) << "Casting graph...\n";
	EssentialGraph graph = castGraph(argGraph);
	uint p = graph.getVertexCount();

	// Cast list of targets
	dout.level(1) << "Casting options...\n";
	dout.level(2) << "  Casting list of targets...\n";
	Rcpp::List data(argPreprocData);
	TargetFamily targets = castTargets(data["targets"]);

	// Cast algorithm string
	dout.level(2) << "  Casting algorithm and options...\n";
	std::string algName = Rcpp::as<std::string>(argAlgorithm);

	// Cast score
	dout.level(2) << "  Casting score...\n";
	Score* score = createScore(Rcpp::as<std::string>(argScore), &targets, data);

	graph.setScore(score);
	graph.setTargets(&targets);

	std::vector<int> steps;
	std::vector<std::string> stepNames;
	std::stringstream ss;

	// Cast option for limits in vertex degree
	dout.level(2) << "  Casting maximum vertex degree...\n";
	Rcpp::NumericVector maxDegree((SEXP)(options["maxDegree"]));
	if (maxDegree.size() > 0) {
		if (maxDegree.size() == 1) {
			if (maxDegree[0] >= 1.) {
				uint uniformMaxDegree = static_cast<uint>(maxDegree[0]);
				graph.limitVertexDegree(uniformMaxDegree);
			}
			else {
				double maxRelativeDegree = maxDegree[0];
				graph.limitVertexDegree(maxRelativeDegree);
			}
		}
		else {
			std::vector<uint> maxDegrees = Rcpp::as< std::vector<uint> >(options["maxDegree"]);
			graph.limitVertexDegree(maxDegrees);
		}
	}

	// Cast option for required phases
	dout.level(2) << "  Casting phases...\n";
	std::vector< std::string > optPhases = Rcpp::as< std::vector<std::string> >(options["phase"]);
	std::vector< step_dir > phases(optPhases.size(), SD_FORWARD);
	for (uint i = 0; i < optPhases.size(); ++i) {
		if (optPhases[i] == "backward") {
			phases[i] = SD_BACKWARD;
		}
		else if (optPhases[i] == "turning") {
			phases[i] = SD_TURNING;
		}
	}
	dout.level(2) << "  Casting iterative...\n";
	bool doIterate = Rcpp::as<bool>(options["iterate"]);

	// Cast option for vertices which are not allowed to have parents
	// TODO: activate function in R, and check for conversion from R to C indexing convention
	std::vector<uint> childrenOnly = Rcpp::as< std::vector<uint> >(options["childrenOnly"]);
	for (std::vector<uint>::iterator vi = childrenOnly.begin(); vi != childrenOnly.end(); ++vi)
		graph.setChildrenOnly(*vi - 1, true);
	int stepLimit;

	// Cast option for fixed gaps: logical matrix, assumed to be symmetric by now
	dout.level(2) << "  Casting fixed gaps...\n";
	if (!Rf_isNull(options["fixedGaps"])) {
		Rcpp::LogicalMatrix gapsMatrix((SEXP)(options["fixedGaps"]));
		uint n_gaps = 0;
		for (uint i = 0; i < p; ++i)
			for (uint j = i + 1; j < p; ++j)
				if (gapsMatrix(i, j))
					n_gaps++;
		// Invert gaps if more than half of the possible edges are fixed gaps
		bool gapsInverted = 4*n_gaps > p*(p - 1);
		EssentialGraph fixedGaps(p);
		for (uint i = 0; i < p; ++i)
			for (uint j = i + 1; j < p; ++j)
				if (gapsMatrix(i, j) ^ gapsInverted)
					fixedGaps.addEdge(i, j, true);
		graph.setFixedGaps(fixedGaps, gapsInverted);
	}

	// Cast option for adaptive handling of fixed gaps (cf. "ARGES")
	dout.level(2) << "  Casting adaptive flag...\n";
	ForwardAdaptiveFlag adaptive(NONE);
	std::string optAdaptive = options["adaptive"];
	dout.level(2) << "Option 'adaptive': " << optAdaptive << std::endl;
	if (optAdaptive == "vstructures") {
		adaptive = VSTRUCTURES;
	}
	if (optAdaptive == "triples") {
		adaptive = TRIPLES;
	}

	// Perform inference algorithm:
	// GIES
	if (algName == "GIES") {
		dout.level(1) << "Performing GIES...\n";

		// Enable caching, if requested
		if (Rcpp::as<bool>(options["caching"]))
			graph.enableCaching();

		// Perform a greedy search with the requested phases, either iteratively or only once
		bool cont;
		int phaseCount(1);
		do {
			cont = false;
			for (uint i = 0; i < phases.size(); ++i) {
				for (steps.push_back(0);
						graph.greedyStepDir(phases[i], adaptive);
						steps.back()++) {
					cont = true;
				}
				ss.str(std::string());
				ss << optPhases[i] << phaseCount;
				stepNames.push_back(ss.str());
			}
			cont &= doIterate;
			phaseCount++;
		} while (cont);
	}

	// Single phase or step of GIES
	else if (algName == "GIES-F" || algName == "GIES-B" || algName == "GIES-T") {
		dout.level(1) << "Performing " << algName << "...\n";

		// Limit to single step if requested
		stepLimit = Rcpp::as<int>(options["maxSteps"]);
		if (stepLimit == 0)
			stepLimit = graph.getVertexCount()*graph.getVertexCount();

		// Enable caching, if requested
		if (options["caching"])
			graph.enableCaching();

		steps.push_back(0);
		if (algName == "GIES-F") {
			for (; steps.back() < stepLimit && graph.greedyForward(adaptive); steps.back()++);
			stepNames.push_back("forward1");
		}
		else if (algName == "GIES-B") {
			for (; steps.back() < stepLimit && graph.greedyBackward(); steps.back()++);
			stepNames.push_back("backward1");
		}
		else if (algName == "GIES-T") {
			for (; steps.back() < stepLimit && graph.greedyTurn(); steps.back()++);
			stepNames.push_back("turning1");
		}
	}

	// Single one or several steps of GIES into either direction
	else if (algName == "GIES-STEP") {
		dout.level(1) << "Performing " << algName << "...\n";

		// Limit to single step if requested
		stepLimit = Rcpp::as<int>(options["maxSteps"]);
		if (stepLimit == 0)
			stepLimit = graph.getVertexCount()*graph.getVertexCount();

		// Steps: 3 entries, storing number of forward, backward, and turning steps
		step_dir dir(SD_NONE), lastDir(SD_NONE);
		std::vector<int> stepCount(4);
		do {
			dir = graph.greedyStep();
			if (dir != SD_NONE) {
				if (dir != lastDir) {
					steps.push_back(1);
					stepCount[0]++;
					ss.str(std::string());
					switch(dir) {
					case SD_FORWARD:
						ss << "forward";
						break;

					case SD_BACKWARD:
						ss << "backward";
						break;

					case SD_TURNING:
						ss << "turning";
						break;
                                        
                                        default:
                                                break;
                                        }
					ss << stepCount[dir]++;
					stepNames.push_back(ss.str());
				} // IF dir
				else {
					steps.back()++;
				}
			} // IF dir
		} while (stepCount[0] < stepLimit && dir != SD_NONE);
	}

	// GDS; yields a DAG, not an equivalence class!
	else if (algName == "GDS") {
		// Perform a greedy search with the requested phases, either iteratively or only once
		bool cont;
		int phaseCount(1);
		do {
			cont = false;
			for (int i = 0; i < phases.size(); ++i) {
				for (steps.push_back(0);
						graph.greedyDAGStepDir(phases[i]);
						steps.back()++) {
					cont = true;
				}
				ss.str(std::string());
				ss << optPhases[i] << phaseCount;
				stepNames.push_back(ss.str());
			}
			cont &= doIterate;
			phaseCount++;
		} while (cont);
	}

	// DP; yields a DAG, not an equivalence class!
	else if (algName == "SiMy") {
		graph.siMySearch();
		// graph.replaceUnprotected();
	}

	// Other algorithm: throw an error
	else throw std::runtime_error(algName + ": invalid algorithm name");

	// Return new list of in-edges and steps
	delete score;
	Rcpp::IntegerVector namedSteps(steps.begin(), steps.end());
	namedSteps.names() = stepNames;

	// TODO "interrupt" zurückgeben, falls Ausführung unterbrochen wurde. Problem:
	// check_interrupt() scheint nur einmal true zurückzugeben...
	return Rcpp::List::create(
			Rcpp::Named("in.edges") = wrapGraph(graph),
			Rcpp::Named("steps") = namedSteps);

	END_RCPP
}

RcppExport SEXP representative(SEXP argGraph)
{
	// Initialize automatic exception handling; manual one does not work any more...
	BEGIN_RCPP

	// Cast graph
	EssentialGraph graph = castGraph(argGraph);

	// Get and return representative
	return wrapGraph(graph.getRepresentative());

	END_RCPP
}

RcppExport SEXP dagToEssentialGraph(SEXP argGraph, SEXP argTargets)
{
	// Initialize automatic exception handling; manual one does not work any more...
	BEGIN_RCPP

	// Cast arguments
	EssentialGraph graph = castGraph(argGraph);
	TargetFamily targets = castTargets(argTargets);

	// Calculate essential graph
	graph.setTargets(&targets);
	graph.replaceUnprotected();

	// Return essential graph
	return wrapGraph(graph);

	END_RCPP
}

/**
 * Calculate the optimal intervention target for an essential graph.
 */
RcppExport SEXP optimalTarget(SEXP argGraph, SEXP argMaxSize)
{
	// Initialize automatic exception handling; manual one does not work any more...
	BEGIN_RCPP

	// Cast arguments
	EssentialGraph graph = castGraph(argGraph);
	int maxSize = Rcpp::as<int>(argMaxSize);

	// Calculate optimal intervention target
	std::set<uint> target = graph.getOptimalTarget(maxSize);

	// Convert from C++ (0 based) to R (1 based) numbering convention.
	std::vector<uint> result(target.begin(), target.end());
	for (std::vector<uint>::iterator vi = result.begin(); vi != result.end(); ++vi)
		(*vi)++;
	return Rcpp::wrap(result);

	END_RCPP
}

/**
 * Calculate the p-value for a conditional independence test on Gaussian data
 */
RcppExport SEXP condIndTestGauss(
		SEXP argVertex1,
		SEXP argVertex2,
		SEXP argCondSet,
		SEXP argSampleSize,
		SEXP argCor)
{
	// Exception handling
	BEGIN_RCPP

	// Cast arguments; note index shift between R and C++!
	uint u = Rcpp::as<uint>(argVertex1) - 1;
	uint v = Rcpp::as<uint>(argVertex2) - 1;
	std::vector<uint> S = Rcpp::as<std::vector<uint> >(argCondSet);
	for (std::vector<uint>::iterator si = S.begin(); si != S.end(); ++si)
		(*si)--;
	uint n = Rcpp::as<uint>(argSampleSize);
	Rcpp::NumericMatrix cor(argCor);

	// Create test object and calculate p-value
	IndepTestGauss indepTest(n, cor);
	return Rcpp::wrap(indepTest.test(u, v, S));

	END_RCPP
}

/**
 * Perform undirected version of PC algorithm, i.e., estimate skeleton of DAG
 * given data
 *
 */
RcppExport SEXP estimateSkeleton(
		SEXP argAdjMatrix,
		SEXP argSuffStat,
		SEXP argIndepTest,
		SEXP argIndepTestFn,
		SEXP argAlpha,
		SEXP argFixedEdges,
		SEXP argOptions)
{
	// Exception handling
	BEGIN_RCPP

	// Cast options and set debug level
	Rcpp::List options(argOptions);
	dout.setLevel(Rcpp::as<int>(options["verbose"]));

	dout.level(1) << "Casting arguments...\n";

	// Cast sufficient statistic and significance level
	dout.level(2) << "Casting sufficient statistic...\n";
	double alpha = Rcpp::as<double>(argAlpha);
	Rcpp::List suffStat(argSuffStat);

	// Cast independence test
	dout.level(2) << "Casting independence test...\n";
	std::string indepTestName = Rcpp::as<std::string>(argIndepTest);
	IndepTest* indepTest;
	if (indepTestName == "gauss") {
		Rcpp::NumericMatrix cor((SEXP)(suffStat["C"]));
		indepTest = new IndepTestGauss(Rcpp::as<uint>(suffStat["n"]), cor);
	}
	else if (indepTestName == "rfun")
		indepTest = new IndepTestRFunction(&suffStat, Rcpp::Function(argIndepTestFn));
	// Invalid independence test name: throw error
	else throw std::runtime_error(indepTestName + ": Invalid independence test name");

	// Initialize OpenMP
	#ifdef _OPENMP
		int threads = Rcpp::as<int>(options["numCores"]);
		if (threads < 0)
			threads = omp_get_num_procs();
		omp_set_num_threads(threads);
	#endif

	// Create list of lists for separation sets
	Rcpp::LogicalMatrix adjMatrix(argAdjMatrix);
	int p = adjMatrix.nrow();
	SepSets sepSet(p, std::vector<arma::ivec>(p, arma::ivec(1)));
	for (int i = 0; i < p; ++i)
		for (int j = 0; j < p; ++j)
			sepSet[i][j].fill(-1);
	// TODO to save space, create a triangular list only

	// Cast graph and fixed edges
	dout.level(2) << "Casting graph and fixed edges...\n";
	Rcpp::LogicalMatrix fixedMatrix(argFixedEdges);
	Skeleton graph(p);
	Rcpp::NumericMatrix pMax(p, p);
	pMax.fill(-1.);
	std::vector<uint> emptySet;
	std::vector<int> edgeTests(1);
	for (int i = 0; i < p; i++) {
		#pragma omp parallel for
		for (int j = i + 1; j < p; j++) {
			if (adjMatrix(i, j) && !fixedMatrix(i, j)) {
				pMax(i, j) = indepTest->test(i, j, emptySet);
				if (pMax(i, j) >= alpha)
					sepSet[j][i].set_size(0);
				dout.level(1) << "  x = " << i << ", y = " << j << ", S = () : pval = "
						<< pMax(i, j) << std::endl;
			}
		}
	}
	for (int i = 0; i < p; i++) {
		for (int j = i + 1; j < p; j++) {
			if (fixedMatrix(i, j))
				graph.addFixedEdge(i, j);
			else if (adjMatrix(i, j)) {
				edgeTests[0]++;
				if (pMax(i, j) < alpha)
					graph.addEdge(i, j);
			}
		}
	}

	// Estimate skeleton
	graph.setIndepTest(indepTest);
	dout.level(1) << "Fitting skeleton to data...\n";
	graph.fitCondInd(alpha,
			pMax,
			sepSet,
			edgeTests,
			Rcpp::as<int>(options["m.max"]),
			Rcpp::as<bool>(options["NAdelete"]));

	// Delete test object
	delete indepTest;

	// Return list
	return Rcpp::List::create(
			Rcpp::Named("amat") = graph.getAdjacencyMatrix(),
			Rcpp::Named("pMax") = pMax,
			Rcpp::Named("sepset") = Rcpp::wrap(sepSet),
			Rcpp::Named("n.edgetests") = edgeTests);

	END_RCPP
}

