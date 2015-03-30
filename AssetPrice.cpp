/* 
	Trey (Carroll) Huffine
	3/15/2014
	Assignment 6 - Final Project
	AssetPrice.cpp
*/

//
//	Description
//

#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "AssetPrice.h"


using namespace std;

//******************************
//	Asset
//******************************
// Virtual Class

//******************************
//	Stock
//******************************
Stock::Stock(void) : Asset()
{
}
Stock::~Stock( void )
{
}
double Stock::getValue()
{
	return StockPrice;
}
void Stock::setValue(double price)
{
	StockPrice = price;
}

//****************************
//	Option
//****************************
Option::Option( double strike0, double spot0, double sigma, double Time, double rate ) : Asset()
{
	//  Set basic parameters
	strike = strike0;
	spot = spot0;
	volatility = sigma;
	T = Time;
	r = rate;
	nowStock.setValue( spot0 );


	//**************************************
	//	Pass stock class instead of value
	//**************************************
	d1 = calc_d1(nowStock,strike,r,volatility,T);
	d2 = calc_d2(nowStock,strike,r,volatility,T);
}
Option::~Option()
{
	//  Nothing to delete
}
double Option::calc_d1(Stock& S, double K, double r, double vix, double T)
{
	return ( log(S.getValue()/K) + (r + .5*vix*vix)*T ) / (vix*sqrt(T));
}

double Option::calc_d2(Stock& S, double K, double r, double vix, double T)
{
	return ( log(S.getValue()/K) + (r - .5*vix*vix)*T ) / (vix*sqrt(T));
}
double Option::calc_NormDens(double x)
{
	return 1/sqrt(2*3.14159265)*exp(-.5*x*x);
}
double Option::calc_CumNorm(double x)
{
	double k, N, a1, a2, a3, a4, a5, gam;

	a1 = .31938153;
	a2 = -.356563782;
	a3 = 1.781477937;
	a4 = -1.821255978;
	a5 = 1.330274429;
	gam = .2316419;

	k = 1/(1+gam*x);

	if( x >= 0.0 )
	{	N = 1 - calc_NormDens(x) * (a1*k + a2*pow(k,2) + a3*pow(k,3) + a4*pow(k,4) + a5*pow(k,5));
		return N;}
	if( x<=0 )
	{	x = -x;
		k = 1/(1+gam*x);
		N = 1 - calc_NormDens(x) * (a1*k + a2*pow(k,2) + a3*pow(k,3) + a4*pow(k,4) + a5*pow(k,5));
		return (1-N); }

	return N;
}
//
//	Value option for a range of spot prices
//
void Option::ValueSpots( double lowerSpotPrice, double upperSpotPrice, int numberOfSpots, string Name  )
{
	ofstream ValSpots(Name);

	//  Find the price increment steps for the given parameters
	double SpotStep = ( upperSpotPrice - lowerSpotPrice ) / (numberOfSpots-1);

	spot = lowerSpotPrice;

	//  Simple file output for the values at different spot prices
	for( int i = 0; i < numberOfSpots; i++ )
	{
		ValSpots << OptionType << "," << strike << "," << spot << "," << volatility << "," << T << "," << r << "," << getValue( ) << endl;
		spot += SpotStep;
	}
}
//
//	Data Access for option class
//
double Option::getStrike() { return strike; }
double Option::getSpot(){ return spot; } //  use the stock class
double Option::getVolatility(){ return volatility; }
double Option::getTime(){ return T; }
double Option::getRate(){ return r; }
string Option::getType(){ return OptionType; }


//**********************************
//	European Call
//**********************************
EuropeanCall::EuropeanCall( double strike0, double spot0, double sigma, double Time, double rate ) : Option( strike0,  spot0, sigma, Time, rate )
{
	//  set type to European call and set the BlackScholes value
	OptionType = "EC";
	BSValue = BlackScholesValue();
}
EuropeanCall::~EuropeanCall( )
{
	// nothing to delete
}
//  Find the exact value using BS
double EuropeanCall::BlackScholesValue( )
{
	return spot*calc_CumNorm(d1)-strike*exp(-r*T)*calc_CumNorm(d2);
}
//  numerically calculate the value from a binomial tree
double EuropeanCall::getValue()
{
	double currVal;

	binomialtree currOpt((*this));

	return currOpt.findval();
}
//  Find the payout if exercised early (0 if t>0)
double EuropeanCall::getPayout( double S, double t )
{
	if( t==0 )
	{
		double testP;
		testP =  S - strike;
		if (testP > 0) return testP;
		else return 0;
	}
	else
	{return 0.0;}
}
//**********************************
//	European Put
//**********************************
EuropeanPut::EuropeanPut( double strike0, double spot0, double sigma, double Time, double rate ) : Option( strike0,  spot0, sigma, Time, rate )
{
	//  Type and BS value
	OptionType = "EP";
	BSValue = BlackScholesValue();
}
EuropeanPut::~EuropeanPut( )
{
	//  Nothing to delete
}
double EuropeanPut::BlackScholesValue( )
{
	return strike*exp(-r*T)*calc_CumNorm(-d2)-spot*calc_CumNorm(-d1);
}
//  Numerical calculation using binomial tree
double EuropeanPut::getValue()
{
	double currVal;

	binomialtree currOpt((*this));

	return currOpt.findval();
}
//  Find the payout if exercised early (0 if t>0)
double EuropeanPut::getPayout(double S, double t)
{
	if( t==0 )
	{
		double testP;
		testP = strike - S;
		if (testP > 0) return testP;
		else return 0;
	}
	else
	{return 0.0;}
}
//**********************************
//	American Call
//**********************************
AmericanCall::AmericanCall( double strike0, double spot0, double sigma, double Time, double rate ) : Option( strike0,  spot0, sigma, Time, rate )
{
	// set type
	OptionType = "AC";
}
AmericanCall::~AmericanCall( )
{
	// nothing to delete
}
// NOT USED
double AmericanCall::BlackScholesValue( )
{
	return spot*calc_CumNorm(d1)-strike*exp(-r*T)*calc_CumNorm(d2);
}
double AmericanCall::getValue()
{
	//  get the value of the option numerically
	double currVal;

	binomialtree currOpt((*this));

	return currOpt.findval();
}
//  Find the payout if exercised early.  In theory this will never be used
double AmericanCall::getPayout(double S, double t)
{
		double testP;
		testP =  S - strike;
		if (testP > 0) return testP;
		else return 0;
}
//**********************************
//	American Put
//**********************************
AmericanPut::AmericanPut( double strike0, double spot0, double sigma, double Time, double rate ) : Option( strike0,  spot0, sigma, Time, rate )
{
	// set type
	OptionType = "AP";
}
AmericanPut::~AmericanPut( )
{
}
// NOT USED
double AmericanPut::BlackScholesValue( )
{
	return spot*calc_CumNorm(d1)-strike*exp(-r*T)*calc_CumNorm(d2);
}
//  Find value numerically with binomial tree
double AmericanPut::getValue()
{
	double currVal;

	binomialtree currOpt((*this));

	return currOpt.findval();
}
//  Get payout if exercised early.  Then handled in binomial tree if it applies.
double AmericanPut::getPayout(double S, double t)
{
		double testP;
		testP = strike - S;
		if (testP > 0) return testP;
		else return 0;
}

//****************************************
//	Binomial Tree
//****************************************
binomialtree::binomialtree( Option& nowOption )
{
	//  Calculate and set the binomial tree parameters
	n = 200;
	size = 0;
	P0 = nowOption.getSpot();
	for( int i = 1; i <= n+1; i++ ) size += i;
	u = exp((nowOption.getVolatility())*sqrt((nowOption.getTime())/n));
	d = 1/u;
	a = exp(nowOption.getRate()*nowOption.getTime()/n);
	q = (a - d)/(u - d);

	//  Set the matrices to the appropriate size
	stockprices.resize(size);
	optionvalues.resize(size);
	boundary.resize(n+1);

	//  initiate the boundary vector
	for( int i = 0; i <= n; i++ ) {boundary[i]=0;}

	//  set stock vector, option price vector, and delta value
	setstock();
	setoption(nowOption);
	setDelta();
	setGamma();
}
binomialtree::~binomialtree()
{
	//  Clearing pointers handled by vector class
}
//  Set the stock price based off of the inintial price
//	multiply the initial price by our variables u and d, depending on location in the tree
void binomialtree::setstock( void )
{
	int m = 0;
	for( int i = 0; i <= n; i++)
	{
		for( int j = 0; j <= i; j++ )
		{
			stockprices[m] = P0*pow(u,j)*pow(d,(i-j));
			m++;
		}
	}
}
//  Find the value of the option at each node
void binomialtree::setoption( Option& nowOption )
{
	//  Declare necessary varaibles
	double t_step = nowOption.getTime() / n;
	double t = 0;
	double currP;
	int m = size-1;
	double tracker = n;

	//  Set the values at t=0
	//  iterate from the last value in the tree backwards
	for ( int i = 0; i <= n; i++ ) 
	{
		//  Find current stock price
		currP = stockprices[m];
		//  Find payout at this stock price for t=0
		optionvalues[m] = nowOption.getPayout(currP,t);

		//  Find the location boundary for payout (basically trivial here)
		if( m < (size-1) )
		{if(optionvalues[m] > 0 && optionvalues[m+1] == 0.0)
		{ boundary[n] = currP; }}
		if( m < (size-1) )
		{if(optionvalues[m] == 0 && optionvalues[m+1] > 0.0)
		{ boundary[n] = currP;}}

		// track m to be used later.  Holds our place in the binomial tree
		m--;
	}

	//  increment time
	t = t + t_step;

	for( tracker; tracker >= 0; tracker-- )
	{
		//  Bool for finding exercise boundary
		bool boundset = true;

		for( int i = 0; i < tracker; i++ )
		{
			currP = stockprices[m];
			double exCheck;
			//  find option value - uses the length of the intervals to back out the option price at each node
			//  from the prices that were previously calculated.
			exCheck = (q*optionvalues[m+1+tracker]+(1-q)*optionvalues[m+tracker])*exp(-nowOption.getRate()*t_step);
			//  see if exercise value is more than current value
			if( nowOption.getPayout(currP,t) > exCheck )
			{ 
				//  set value to payout if it is greater than the option value
				optionvalues[m] = nowOption.getPayout(currP,t); 

				// American Put test for boundary
				if(optionvalues[m] > optionvalues[m+1] && boundset == true)
				{ boundary[tracker-1] = currP; boundset = false;}
				// American Call test for boundary
				if(optionvalues[m] > 0 && boundset == true)
				{ boundary[tracker-1] = currP;}
			}
			else
			//  set value if option is not exercised
			{ optionvalues[m] = exCheck; }
			m--;
		}
		t = t + t_step;
	}
	//for( int i = 0; i <= n; i++ ) 
	//{ cout << "boundary i " << boundary[i] << endl; }
}
// calculate the value of delta using the intial stock price and the next two values
void binomialtree::setDelta( void )
{
	delta = (optionvalues[2]-optionvalues[1])/(stockprices[2]-stockprices[1]);
}
//  return current option value
double binomialtree::findval( void )
{
	return optionvalues[0];
}
//  return vector containing the exercise boundaries
vector<double> binomialtree::ExerciseBoundary( void )
{
	return boundary;
}
// return value of delta
double binomialtree::finddelta( void )
{
	return delta;
}
void binomialtree::setGamma( void )
{
	double h = .5*(P0*u*u-P0*d*d);
	gamma = ((optionvalues[5]-optionvalues[4])/(P0*u*u-P0)-(optionvalues[4]-optionvalues[3])/(P0-P0*d*d))/(h);
}
double binomialtree::findgamma( void )
{
	return gamma;
}

//************************
//*** Helper functions ***
//************************

//  find the current token of input
bool nextToken(const string& line, int& offset, string& token)
{
	int len = line.length();
	bool rc = (offset < len);
	token = "";
	bool done = false;
	// continue through the tokens until you a "stop" value or the line is done being read.  most likely a comma since we are working is a .csv
	while (!done && (offset < len))
	{
		char c = line.at(offset);
		if (c == ',' || c == '\r' || c == '\n' || c == '\t' || c == ' ' )
		{
			done = true;
		} else {
			token += c;
		}
		offset++;
	}
	return rc;
}


void ValueOptions( string inName, string outName )
{
	ifstream inFile( inName );
	ofstream outFile( outName );

	// iterate until the file is done being read
	while( !inFile.eof() )
	{
		string AssetType, line, token;
		getline(inFile,line);
		int offset = 0;
		double thisStrike, thisSpot, thisVix, thisTime,thisRate;

		//  Parse through the lines and find the desired values
		bool rc = nextToken(line, offset, token);
		AssetType = (token.c_str());
		rc &= nextToken(line, offset, token);
		thisStrike = atof(token.c_str());
		rc &= nextToken(line, offset, token); //atof->int->string
		thisSpot = atof(token.c_str());
		rc &= nextToken(line, offset, token);
		thisVix = atof(token.c_str());
		rc &= nextToken(line, offset, token);
		thisTime = atof(token.c_str());
		rc &= nextToken(line, offset, token);
		thisRate = atof(token.c_str());

		//  Output to the .csv file based on the option type.  add price of stock at the end.
		if( AssetType == "EC" )
		{
			EuropeanCall tempOption(thisStrike,thisSpot,thisVix,thisTime,thisRate);
			outFile << AssetType << "," << thisStrike << "," << thisSpot << "," << thisVix << "," << thisTime << "," << thisRate << "," << tempOption.BlackScholesValue( ) << endl;
		}
		if( AssetType == "EP" )
		{
			EuropeanPut tempOption(thisStrike,thisSpot,thisVix,thisTime,thisRate);
			outFile << AssetType << "," << thisStrike << "," << thisSpot << "," << thisVix << "," << thisTime << "," << thisRate << "," << tempOption.BlackScholesValue( ) << endl;
		}
		if( AssetType == "AC" )
		{
			AmericanCall tempOption(thisStrike,thisSpot,thisVix,thisTime,thisRate);
			outFile << AssetType << "," << thisStrike << "," << thisSpot << "," << thisVix << "," << thisTime << "," << thisRate << "," << tempOption.getValue( ) << endl;
		}
		if( AssetType == "AP" )
		{
			AmericanPut tempOption(thisStrike,thisSpot,thisVix,thisTime,thisRate);
			outFile << AssetType << "," << thisStrike << "," << thisSpot << "," << thisVix << "," << thisTime << "," << thisRate << "," << tempOption.getValue( ) << endl;
		}
	}
}


void ValueDelta( string inName, string outName )
{
	ifstream inFile( inName );
	ofstream outFile( outName );

	// iterate until the file is done being read
	while( !inFile.eof() )
	{
		//  set vars
		string AssetType, line, token;
		getline(inFile,line);
		int offset = 0;
		double thisStrike, thisSpot, thisVix, thisTime,thisRate;


		// parse lines for the desired information
		bool rc = nextToken(line, offset, token);
		AssetType = (token.c_str());
		rc &= nextToken(line, offset, token);
		thisStrike = atof(token.c_str());
		rc &= nextToken(line, offset, token); //atof->int->string
		thisSpot = atof(token.c_str());
		rc &= nextToken(line, offset, token);
		thisVix = atof(token.c_str());
		rc &= nextToken(line, offset, token);
		thisTime = atof(token.c_str());
		rc &= nextToken(line, offset, token);
		thisRate = atof(token.c_str());

		//  Output to the .csv file based on the option type.  add price of stock and value of delta at the end.
		if( AssetType == "EC" )
		{
			EuropeanCall tempOption(thisStrike,thisSpot,thisVix,thisTime,thisRate);
			binomialtree tempTree(tempOption);
			outFile << AssetType << "," << thisStrike << "," << thisSpot << "," << thisVix << "," << thisTime << "," << thisRate << "," << tempOption.BlackScholesValue( ) << "," << tempTree.finddelta() << endl;
		}
		if( AssetType == "EP" )
		{
			EuropeanPut tempOption(thisStrike,thisSpot,thisVix,thisTime,thisRate);
			binomialtree tempTree(tempOption);
			outFile << AssetType << "," << thisStrike << "," << thisSpot << "," << thisVix << "," << thisTime << "," << thisRate << "," << tempOption.BlackScholesValue( ) << "," << tempTree.finddelta() << endl;
		}
		if( AssetType == "AC" )
		{
			AmericanCall tempOption(thisStrike,thisSpot,thisVix,thisTime,thisRate);
			binomialtree tempTree(tempOption);
			outFile << AssetType << "," << thisStrike << "," << thisSpot << "," << thisVix << "," << thisTime << "," << thisRate << "," << tempOption.getValue( ) << "," << tempTree.finddelta() << endl;
		}
		if( AssetType == "AP" )
		{
			AmericanPut tempOption(thisStrike,thisSpot,thisVix,thisTime,thisRate);
			binomialtree tempTree(tempOption);
			outFile << AssetType << "," << thisStrike << "," << thisSpot << "," << thisVix << "," << thisTime << "," << thisRate << "," << tempOption.getValue( ) << "," << tempTree.finddelta() << endl;
		}
	}
}