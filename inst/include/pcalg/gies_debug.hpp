/*
 * Debugging functions for GIES
 *
 * @author Alain Hauser
 * $Id: gies_debug.hpp 248 2014-03-03 11:27:22Z alhauser $
 */

#ifndef GIES_DEBUG_HPP_
#define GIES_DEBUG_HPP_

#include <iostream>
#include <Rcpp.h>

// Define default debug level
#ifndef DEBUG_OUTPUT_LEVEL
#define DEBUG_OUTPUT_LEVEL 0
#endif

namespace std {

/**
 * Sends a vector to an output stream (e.g., cout)
 */
template <typename T> ostream& operator<<(ostream& out, const vector<T>& vec) {
	out << "(";
	for (size_t i = 0; i + 1 < vec.size(); i++)
		out << vec[i] << ", ";
	if (!vec.empty())
		out << vec.back();
	out << ")";
	return out;
}

/**
 * Sends a set to an output stream (e.g., cout)
 */
template <typename T> ostream& operator<<(ostream& out, const set<T>& s) {
	typename set<T>::iterator iter;
	out << "{";
	iter = s.begin();
	if (iter != s.end()) {
		out << *iter;
		for (++iter; iter != s.end(); ++iter)
			out << ", " << *iter;
	}
	out << "}";
	return out;
}

} // NAMESPACE std

/**
 * Handling user interrupt
 */
static inline void check_interrupt_impl(void* /*dummy*/) {
	R_CheckUserInterrupt();
}

inline bool check_interrupt() {
	return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE);
}

// Macro for outputting debug messages depending on debug level
#if DEBUG_OUTPUT_LEVEL >= 1
#define DBOUT( level, message ) if ( DEBUG_OUTPUT_LEVEL >= level ) Rcpp::Rcout << message << std::endl
#else
#define DBOUT( level, message )
#endif

namespace std {

class nullbuf: public streambuf
{
	int xputc(int) { return 0; }
	streamsize xsputn(char const *, streamsize n) { return n; }
	int sync() { return 0; }
};

} // NAMESPACE std

class DebugStream
{
	int _level;
	std::nullbuf _nullbuf;
	std::ostream _nullstream;

public:
	DebugStream() : _level(0), _nullstream(&_nullbuf) {}
	DebugStream(int level) : _level(level), _nullstream(&_nullbuf) {}

	void setLevel(const int level) { _level = level;	}

	int getLevel() const { return _level; }

	std::ostream& level(const int messageLevel)
	{
		if (messageLevel <= _level)
			return Rcpp::Rcout;
		return _nullstream;
	}
};

#ifdef DEFINE_GLOBAL_DEBUG_STREAM
#define GLOBAL
#else
#define GLOBAL extern
#endif

GLOBAL DebugStream dout;

#endif /* GIES_DEBUG_HPP_ */
