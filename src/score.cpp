/*
 * @author Alain Hauser
 * $Id: score.cpp 238 2014-02-26 16:39:25Z alhauser $
 */

#include "score.hpp"
#include "greedy.hpp"

#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <limits>
#include "gies_debug.hpp"

bool TargetFamily::protects(const uint a, const uint b) const
{
	bool aInI, bInI;
	for (std::size_t i = 0; i < size(); i++) {
		aInI = (std::find((*this)[i].begin(), (*this)[i].end(), a) != (*this)[i].end());
		bInI = (std::find((*this)[i].begin(), (*this)[i].end(), b) != (*this)[i].end());
		if (aInI ^ bInI) return true;
	}

	return false;
}

TargetFamily castTargets(const SEXP argTargets)
{
	int i;
	Rcpp::List listIventTargets(argTargets);
	std::vector<uint> vecTarget;
	std::vector<uint>::iterator vi;
	TargetFamily result;
	result.resize(listIventTargets.size());
	for (i = 0; i < listIventTargets.size(); i++) {
		vecTarget = listIventTargets[i];
		// Adapt indices to C++ convention...
		for (vi = vecTarget.begin(); vi != vecTarget.end(); ++vi)
			result[i].insert(*vi - 1);
	}
	return result;
}

std::set<uint> castVertices(SEXP argVertices)
{
	std::set<uint> result;
	std::vector<uint> vertices = Rcpp::as<std::vector<uint> >(argVertices);
	std::vector<uint>::iterator vi;

	for (vi = vertices.begin(); vi != vertices.end(); ++vi)
		result.insert(*vi - 1);

	return result;
}

double Score::global(const EssentialGraph& dag) const
{
	double result = 0.;
	uint v;

	// Standard for decomposable scores: global score is sum of local scores
	for (v = 0; v < _vertexCount; ++v)
		result += local(v, dag.getParents(v));

	return result;
}

Score* createScore(std::string name, TargetFamily* targets, Rcpp::List data)
{
	Score* result;

	dout.level(2) << "Creating score object of type '" << name << "'...\n";

	if (name == "gauss.l0pen.scatter")
		result = new ScoreGaussL0PenScatter(Rcpp::as<uint>(data["vertex.count"]), targets);
	else if (name == "none")
		result = new ScoreRFunction(Rcpp::as<uint>(data["vertex.count"]), targets);
	// Invalid score name: throw error
	else throw std::runtime_error(name + ": Invalid score name");

	// Provide preprocessed data to score object
	result->setData(data);

	return result;
}

void ScoreRFunction::setData(Rcpp::List& data)
{
	_totalDataCount = data["total.data.count"];

	// Do not really store "data", but R function object used for the
	// calculation
	dout.level(2) << "Casting R functions to calculate the score...\n";
	_rfunction.push_back(Rcpp::as<Rcpp::Function>(data["local.score"]));
	_rfunction.push_back(Rcpp::as<Rcpp::Function>(data["global.score"]));
	_rfunction.push_back(Rcpp::as<Rcpp::Function>(data["local.mle"]));
	_rfunction.push_back(Rcpp::as<Rcpp::Function>(data["global.mle"]));
}

double ScoreRFunction::local(const uint vertex, const std::set<uint>& parents) const
{
	// Adapt indices to R convention...
	std::vector<uint> shiftParents;
	shiftParents.reserve(parents.size());
	std::set<uint>::iterator si;
	for (si = parents.begin(); si != parents.end(); ++si)
		shiftParents.push_back(*si + 1);

	// Call R function for local score
	return Rcpp::as<double>((_rfunction[R_FCN_INDEX_LOCAL_SCORE])(vertex + 1, shiftParents));
}

double ScoreRFunction::global(const EssentialGraph& dag) const
{
	// Create list of in-edges to pass to R function;
	// adapt indices to R convention...
	std::vector<std::vector<uint> > inEdges(_vertexCount);
	std::set<uint> parents;
	std::set<uint>::iterator si;
	uint v;
	for (v = 0; v < _vertexCount; ++v) {
		parents = dag.getParents(v);
		inEdges[v].reserve(parents.size());
		for (si = parents.begin(); si != parents.end(); ++si)
			inEdges[v].push_back(*si + 1);
	}

	// Call R function for global score
	return Rcpp::as<double>(_rfunction[R_FCN_INDEX_GLOBAL_SCORE](inEdges));
}

std::vector<double> ScoreRFunction::localMLE(const uint vertex, const std::set<uint>& parents) const
{
	// Change parents indices to R convention
	std::vector<uint> shiftParents(parents.begin(), parents.end());
	std::vector<uint>::iterator vi;
	for (vi = shiftParents.begin(); vi != shiftParents.end(); ++vi)
		(*vi)++;

	// Return local MLE
	return Rcpp::as< std::vector<double> >(_rfunction[R_FCN_INDEX_LOCAL_MLE](vertex + 1, shiftParents));
}

std::vector< std::vector<double> > ScoreRFunction::globalMLE(const EssentialGraph& dag) const
{
	// Construct list of in-edges
	std::set<uint> parents;
	Rcpp::IntegerVector shiftParents;
	Rcpp::List inEdges(dag.getVertexCount());
	for (uint v = 0; v < _vertexCount; ++v) {
		// Get parents of vertex v and adapt their indices to the R convention
		parents = dag.getParents(v);
		shiftParents = Rcpp::IntegerVector(parents.begin(), parents.end());
		for (R_len_t i = 0; i < shiftParents.size(); ++i)
			shiftParents[i]++;

		// Add parents to list of in-edges
		inEdges[v] = shiftParents;
	}

	// Calculate (and return) global MLE
	Rcpp::List listMLE = _rfunction[R_FCN_INDEX_GLOBAL_MLE](inEdges);
	std::vector< std::vector<double> > result(listMLE.size());
	for (R_len_t i = 0; i < listMLE.size(); ++i)
		result[i] = Rcpp::as<std::vector<double> >(listMLE[i]);
	return result;
}

void ScoreGaussL0PenScatter::setData(Rcpp::List& data)
{
	std::vector<int>::iterator vi;
	//uint i;

	// Cast preprocessed data from R list
	dout.level(2) << "Casting preprocessed data...\n";
	_dataCount = Rcpp::as<std::vector<int> >(data["data.count"]);
	dout.level(3) << "# samples per vertex: " << _dataCount << "\n";
	_totalDataCount = Rcpp::as<uint>(data["total.data.count"]);
	dout.level(3) << "Total # samples: " << _totalDataCount << "\n";
	Rcpp::List scatter = data["scatter"];
	Rcpp::NumericMatrix scatterMat;
	_disjointScatterMatrices.resize(scatter.size());
	dout.level(3) << "# disjoint scatter matrices: " << scatter.size() << "\n";
	for (R_len_t i = 0; i < scatter.size(); ++i) {
		scatterMat = Rcpp::NumericMatrix((SEXP)(scatter[i]));
		_disjointScatterMatrices[i] = arma::mat(scatterMat.begin(), scatterMat.nrow(), scatterMat.ncol(), false);
	}

	// Cast index of scatter matrices, adjust R indexing convention to C++
	std::vector<int> scatterIndex = Rcpp::as<std::vector<int> >(data["scatter.index"]);
	for (std::size_t i = 0; i < scatterIndex.size(); ++i)
		_scatterMatrices[i] = &(_disjointScatterMatrices[scatterIndex[i] - 1]);

	// Cast lambda: penalty constant
	_lambda = Rcpp::as<double>(data["lambda"]);
	dout.level(3) << "Penalty parameter lambda: " << _lambda << "\n";

	// Check whether an intercept should be calculated
	_allowIntercept = Rcpp::as<bool>(data["intercept"]);
	dout.level(3) << "Include intercept: " << _allowIntercept << "\n";
}

double ScoreGaussL0PenScatter::local(const uint vertex, const std::set<uint>& parents) const
{
	double a;

	dout.level(3) << "Calculating local score...\n";

	// Cast parents set to Armadillo uvec
	arma::uvec parVec(_allowIntercept ? parents.size() + 1 : parents.size());
	std::copy(parents.begin(), parents.end(), parVec.begin());
	arma::uvec vVec(1);
	vVec[0] = vertex;

	// If intercept is allowed, add "fake parent" taking care of intercept
	if (_allowIntercept)
		parVec[parents.size()] = _vertexCount;
	dout.level(3) << "Vertex: " << vertex << "; parents (adjusted acc. to interc.): " << parVec << "\n";

	// Calculate value in the logarithm of maximum likelihood
	a = (*(_scatterMatrices[vertex]))(vertex, vertex);
	if (parVec.size()) {
		// TODO: evtl. wieder umschreiben, s.d. keine Cholesky-Zerlegung mehr
		// gebraucht wird: macht Code nämlich etwas langsamer... (wider Erwarten!)
		//arma::colvec b = arma::subvec(*(_scatterMatrices[vertex]), parents.begin(), parents.end(), vertex);
		arma::mat R;
		if (!arma::chol(R, _scatterMatrices[vertex]->submat(parVec, parVec)))
			return std::numeric_limits<double>::quiet_NaN();
		arma::colvec c = arma::solve(arma::trimatl(arma::trans(R)),
				_scatterMatrices[vertex]->submat(parVec, vVec));
		//a -= arma::as_scalar(arma::trans(b) * arma::solve(arma::submat(*(_scatterMatrices[vertex]), parents.begin(), parents.end(), parents.begin(), parents.end()), b));
		a -= arma::dot(c, c);
	}

	// Finish calculation of partial BIC score
	return -0.5*(1. + log(a/_dataCount[vertex]))*_dataCount[vertex] - _lambda*(1. + parents.size());
}

double ScoreGaussL0PenScatter::global(const EssentialGraph& dag) const
{
	double result = 0.;
	uint v;

	// L0-penalized score is decomposable => calculate sum of local scores
	for (v = 0; v < dag.getVertexCount(); ++v)
		result += local(v, dag.getParents(v));

	return result;
}

std::vector<double> ScoreGaussL0PenScatter::localMLE(const uint vertex, const std::set<uint>& parents) const
{
	std::vector<double> result(parents.size() + 2);
	std::set<uint>::iterator pi;
	arma::colvec b;
	arma::mat S_papa, S_pav;
	arma::uvec parVec;
	arma::uvec vVec(1);

	dout.level(3) << "Calculating local MLE...\n";

	// Get parents, copy them to Armadillo vector
	parVec.set_size(_allowIntercept ? parents.size() + 1 : parents.size());
	std::copy(parents.begin(), parents.end(), parVec.begin());
	if (_allowIntercept)
		parVec[parents.size()] = _vertexCount;
	vVec[0] = vertex;
	dout.level(3) << "Vertex: " << vertex << "; parents (adjusted acc. to interc.): " << parVec << "\n";

	// Initialize parameter for variance
	result[0] = _scatterMatrices[vertex]->at(vertex, vertex) / _dataCount[vertex];

	// Calculate regression coefficients
	if (parVec.size()) {
		S_pav = _scatterMatrices[vertex]->submat(parVec, vVec);
		S_papa = _scatterMatrices[vertex]->submat(parVec, parVec);
		b = arma::solve(S_papa, S_pav);

		// Correct error variance
		result[0] += (arma::dot(b, S_papa * b) - 2. * arma::dot(S_pav, b)) / _dataCount[vertex];

		// Store intercept, if requested (otherwise 0)
		result[1] = (_allowIntercept ? b(b.n_elem - 1) : 0.);

		// Copy coefficients to result vector
		std::copy(b.memptr(), b.memptr() + (_allowIntercept ? b.n_elem - 1 : b.n_elem), result.begin() + 2);
	}

	dout.level(3) << "Local MLE: " << result << "\n";

	return result;
}

std::vector< std::vector<double> > ScoreGaussL0PenScatter::globalMLE(const EssentialGraph& dag) const
{
	// Calculate local MLE for all vertices
	std::vector< std::vector<double> > result(_vertexCount);
	uint v;
	for (v = 0; v < dag.getVertexCount(); ++v)
		result[v] = localMLE(v, dag.getParents(v));

	return result;
}
