/* 
	Trey (Carroll) Huffine
	3/15/2014
	Assignment 6 - Final Project
	Main
*/

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

int main ()
{
	// ORDER OF INPUTS
	//	double strike0, double spot0, double sigma, double Time, double rate 

	cout << "Final Project - Option Value Calculator" << endl << endl;

	//TEST VALUES
	EuropeanCall test1(100,81,.4,1,.05);
	EuropeanPut test2(100.0,81.0,.4,1.0,.05);

	double P1 = test1.getValue();
	double P2 = test2.getValue();
	cout << "Test Values" << endl;
	cout << "********Call = " << P1 << endl;
	cout << "********Put = " << P2 << endl << endl;

	//  Value spots - 1(e)
	test2.ValueSpots(75,125,51,"ValueSpots.csv");

	//  Value options input - 1(f)
	ValueOptions( "ValueOptionsInput.csv", "ValueOptionsOutput.csv" );
	
	//  Tree depth - 2(b)
	EuropeanCall Treedepth(90,90,.25,.5,.09);

	cout << "BS Value   = " << Treedepth.BlackScholesValue() << endl;
	cout << "Tree Value = " << Treedepth.getValue() << endl << endl;

	//  European vs American comparison - 2(c)
	ofstream Put_file("PutComparison.csv");
	for( int i = 80; i <= 110; i++ )
	{
		EuropeanPut Eput(i,100,.3,2,.08);
		AmericanPut Aput(i,100,.3,2,.08);

		Put_file << i << "," << Eput.getValue() << "," << Aput.getValue() << endl;
	}

	//  Value delta input - 2(d)
	ValueDelta( "ValueDeltaInput.csv", "ValueDeltaOutput.csv" );

	//  Exercise boundary - 2(e)
	AmericanPut put1(100,90,.25,1,.05);
	AmericanPut put2(100,90,.45,1,.05);
	binomialtree testtree1(put1);
	binomialtree testtree2(put2);

	vector<double> boundtest1;
	boundtest1 = testtree1.ExerciseBoundary();

	vector<double> boundtest2;
	boundtest2 = testtree2.ExerciseBoundary();

	ofstream Bound_file("BoundComparison.csv");
	for( int i = 0; i < 100; i++ )
	{
		Bound_file << i << "," << boundtest1[i] << "," << boundtest2[i] << endl;
	}

	EuropeanPut testPut( 45, 40, .4, .5, .05 );

	cout << "HW American Put: " << testPut.getValue() << endl << endl;

	system("pause");
	return 0;
}