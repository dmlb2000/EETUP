/*
 ------------------------------------------------------------------------------- 
 PROGRAM NAME: br2en.cc
 AUTHOR:  H. Cho, Pacific Northwest National Laboratory 
 CREATION DATE:  11 October 2012 
 LAST MODIFICATION DATE: 
 DESCRIPTION:  Program to convert Bruker acquisition data file to an ENUF file.
 Link to the object library "ENUF_objects.h."  Don't forget to keep the
 version name softwareVersion updated.
 ------------------------------------------------------------------------------
 */
#include "fstream"
#include "iostream"
#include "math.h"
#include "iomanip"
#include "ENUF_objects.h"
using namespace std;
const int name_length = 200; 
const int max_chars = 1048576;
const int max_points = 524288;
const char softwareVersion[] = "br2en.101112a";
//
int main( int argc, char* argv[])
{
	char infile[name_length], outfile[name_length];
	int endian_reply, numberpts, numberFIDs, number_data_chars,
    number_complex_points, fillbytes;
    double remainder;
/*
------------------------------------------------------------------------------
 Interactively enter starting information at run time.
------------------------------------------------------------------------------
*/
	cout << "\nEnter name of Bruker acquisition file: ";
	cin.getline(infile, name_length);
	cout << "Enter name of ENUF file: ";
	cin.getline(outfile, name_length);
	cout << "Do a big endian/little endian byte swap (y/n)?: ";
    cin.getline((char *) &endian_reply, sizeof(int));
	cout << "Enter number of points (r+i) per FID: ";
	cin >> numberpts;
    number_complex_points = numberpts/2;
	cout << "Enter number of FIDs: ";
	cin >> numberFIDs;
	number_data_chars = numberpts * 4;
    remainder = fmod((double) number_data_chars, 1024.);
    if(remainder != 0.)
    {
        fillbytes = 1024 - (int) remainder;
    }
	cout << "remainder = " << remainder
    << "\nbytes to fill to 1024 = " << fillbytes << endl;
/*
------------------------------------------------------------------------------
 Connect input stream to file containing Bruker data, and terminate program
 if file not found. 
------------------------------------------------------------------------------
*/
	ifstream indata(&infile[0]);
	if(!indata)
	{
		cout << "\nError opening Bruker file.  Program aborted.\n";
		return 1;
	}
/*
------------------------------------------------------------------------------
 Define some variable names.
------------------------------------------------------------------------------ 
*/
	int I0, I1, complex_signal_int[max_points];
	float complex_signal_float[max_points];
/* 
------------------------------------------------------------------------------
 Open output file. 
------------------------------------------------------------------------------
*/
	ofstream outdata; 
	outdata.open(outfile, ios::out); 
/*
------------------------------------------------------------------------------
 Read and write data FID by FID.
------------------------------------------------------------------------------
*/
	outdata.write((char *) &softwareVersion, 12);
	outdata.write((char *) &number_complex_points, sizeof(int));
	outdata.write((char *) &numberFIDs, sizeof(int));	
	for(I1 = 0; I1 < numberFIDs; I1++)
	{
		indata.read((char *) complex_signal_int, number_data_chars);
		for(I0 = 0; I0 < numberpts; I0++)
		{
			if(endian_reply == 'y')
			{
				reverse_byte_order((char *) &complex_signal_int[I0], 4);
			}
			complex_signal_float[I0] = ((float) complex_signal_int[I0]);
		}
        outdata.write((char *) complex_signal_float, number_data_chars);
        if(remainder != 0)
        {
            indata.read((char *) &complex_signal_int[number_complex_points], fillbytes);
        }
	}
	/*
	 ------------------------------------------------------------------------------
	 Flush the output data and close the input and output files.
	 ------------------------------------------------------------------------------
	 */
	indata.close();
	outdata.flush();
    outdata.close();
}