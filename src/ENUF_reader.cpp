/*
 -------------------------------------------------------------------------------
 PROGRAM NAME: ENUF_reader.cc
 AUTHOR:  H. Cho
 CREATION DATE:  4 October 2012
 LAST MODIFICATION DATE:
 DESCRIPTION:  Program to inspect an ENUF file.
 ------------------------------------------------------------------------------
 */
#include "fstream"
#include "iostream"
#include "math.h"
using namespace std;
const int name_length = 200;
const int max_points = 1048576;
//
int main( int argc, char* argv[])
{
    char infile[name_length], outfile[name_length];
	int file_reply;
    ofstream outdata;
/*
------------------------------------------------------------------------------
 Enter starting information at run time.  Names of the input ENUF file and
 an output text file are requested.
------------------------------------------------------------------------------
*/
    cout << "Enter input ENUF filename: ";
    cin.getline(infile, name_length);
	cout << "Write text data to an output file (y/n)?: ";
    cin.getline((char *) &file_reply, sizeof(int));
	if(file_reply == 'y')
	{
		cout << "Enter output filename: ";
		cin.getline(outfile, name_length);
        outdata.open(outfile, ios::out);
	}
/*
------------------------------------------------------------------------------
 Connect input stream to file containing ENUF data, and terminate program
 if file not found.
------------------------------------------------------------------------------
*/
    ifstream indata(&infile[0]);
    if(!indata)
    {
        cout << "\nError opening ENUF file.  Program aborted.\n";
        return 1;
    }
/*
------------------------------------------------------------------------------
 Define variable names.
------------------------------------------------------------------------------
*/
    int
    I0, numberPoints, numberFIDs, first_point, last_point, index;
    float
    complex_point_float[max_points];
	char
	softwareVersion[13];
/*
------------------------------------------------------------------------------
 Read file containing data.
------------------------------------------------------------------------------
*/
	indata.read(softwareVersion, 12);
    indata.read((char *) &numberPoints, sizeof(int));
    indata.read((char *) &numberFIDs, sizeof(int));
	cout
	<< "\nENUF translator program: " << softwareVersion
    << "\nPoints per interferogram = " << numberPoints
    << "\nNumber of interferograms = " << numberFIDs
	<< endl;
/*
------------------------------------------------------------------------------
 Write data to console and to output file if specified.
------------------------------------------------------------------------------
*/
	cout << "\nEnter first complex point of range: ";
	cin >> first_point;
	cout << "\nEnter last complex point of range: ";
	cin >> last_point;
	for(I0 = 0; I0 < first_point - 1; I0++)
	{
            index = I0 * 2;
            indata.read((char *) &complex_point_float[index], 4);
            index += 1;
            indata.read((char *) &complex_point_float[index], 4);
	}
	for(I0 = first_point - 1; I0 < last_point; I0++)
	{
		index = I0 * 2;
		indata.read((char *) &complex_point_float[index], 4);
		index += 1;
		indata.read((char *) &complex_point_float[index], 4);
		cout << "\n" << I0 << "   " << complex_point_float[index-1]
		<< "   " << complex_point_float[index];
		if(file_reply == 'y')
		{
			outdata << complex_point_float[index-1]
			<< "   " << complex_point_float[index] << "\n";
		}
	}
/*
------------------------------------------------------------------------------
Close input file.
------------------------------------------------------------------------------
*/
    indata.close();
	if(file_reply == 'y')
	{
		outdata.flush();
		outdata.close();
	}
}