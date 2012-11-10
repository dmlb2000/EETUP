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
#include "rappture.h"
using namespace std;

const int max_points = 1048576;
//
int main( int argc, char* argv[])
{
    const char* infile = NULL;
    const char* fpoint_str = NULL;
    const char* lpoint_str = NULL;
    const char* path = "input.string(enuffile).current";
    const char* fpoint_path = "input.integer(firstpoint).current";
    const char* lpoint_path = "input.integer(lastpoint).current";
    const char* graph_path = "output.curve(result).component.xy";
    char output_str[1024];
    int first_point, last_point;
    RpLibrary* lib;
    int err = 0;

    if(argc < 2)
	return -1;

    lib = rpLibrary(argv[1]);
    if (lib != NULL) {
        cout << "creation of library successful" << endl;
    }
    else {
        cout << "creation of library failed" << endl;
    }
    err = rpGetString(lib,path,&infile);
    if (infile != NULL) {
        cout << "infile = " << infile << endl;
    }
    else {
        cout << "Failed to get infile." << endl;
    }

    err = rpGetString(lib,fpoint_path,&fpoint_str);
    if (fpoint_str != NULL) {
        cout << "fpoint_str = " << fpoint_str << endl;
    }
    else {
        cout << "Failed to get fpoint_str." << endl;
    }

    err = rpGetString(lib,lpoint_path,&lpoint_str);
    if (lpoint_str != NULL) {
        cout << "lpoint_str = " << lpoint_str << endl;
    }
    else {
        cout << "Failed to get lpoint_str." << endl;
    }
/*
------------------------------------------------------------------------------
 Enter starting information at run time.  Names of the input ENUF file and
 an output text file are requested.
------------------------------------------------------------------------------
*/
/*
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
*/
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
    I0, numberPoints, numberFIDs, index;
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
/*
	cout << "\nEnter first complex point of range: ";
	cin >> first_point;
	cout << "\nEnter last complex point of range: ";
	cin >> last_point;
*/
	first_point = atoi(fpoint_str);
	last_point = atoi(lpoint_str);
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
		sprintf(output_str, "%g %g\n", complex_point_float[index-1], complex_point_float[index]);
		err = rpPutString(lib,graph_path,output_str,1);
		if (err == 0) {
			cout << "rpPutString successful" << endl;
		}
		else {
			cout << "rpPutString failed" << endl;
		}
/*
		cout << "\n" << I0 << "   " << complex_point_float[index-1]
		<< "   " << complex_point_float[index];
		if(file_reply == 'y')
		{
			outdata << complex_point_float[index-1]
			<< "   " << complex_point_float[index] << "\n";
		}
*/
	}
/*
------------------------------------------------------------------------------
Close input file.
------------------------------------------------------------------------------
*/
    indata.close();
    rpFreeLibrary(&lib);
}
