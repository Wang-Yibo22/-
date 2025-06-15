#pragma once

#include <time.h>
#include <iostream>

using namespace std;  

//! Clock class for timing
class Clock
{

private:

	clock_t t0_, t1_;
	double ct_;
	bool st0_;   //!< Flag for Start method
	bool st1_;   //!< Flag for Stop method

public:

//!	Constructor
	Clock();
  
//!	Start the clock
	void Start();

//!	Stop the clock
	void Stop();
  
//!	Resume the stoped clock
	void Resume();

//!	Clear the clock
	void Clear();

//!	Return the elapsed time since the clock started
	double ElapsedTime();

};
