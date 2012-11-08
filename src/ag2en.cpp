/*
 ------------------------------------------------------------------------------- 
 PROGRAM NAME: ag2en.cc
 AUTHOR:  H. Cho, Pacific Northwest National Laboratory 
 CREATION DATE:  20 September 2012 
 LAST MODIFICATION DATE: 
 DESCRIPTION:  Program to convert VNMR data file to an ENUF file.
 Link to the object library "ENUF_objects.h."  Don't forget to keep the
 version name softwareVersion updated.
 ------------------------------------------------------------------------------
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "fstream"
#include "iostream"
#include "iomanip"
#include "math.h"
#include "ENUF_objects.h"
using namespace std;
const int name_length = 200; 
const int max_chars = 1048576;
const int max_points = 524288;
const char softwareVersion[] = "ag2en.100412a";

/*
 * Return the usage to stdout C++ style.
 */
void usage(char *progname)
{
	cout << progname << " -i FILE -o FILE [-b]" << endl << endl;
	cout << "\t-i FILE - Input VNMR file name." << endl;
	cout << "\t-o FILE - Output ENUF file name." << endl;
	cout << "\t-b      - Byte Swap the data." << endl;
}

int main( int argc, char* argv[])
{
	char infile[name_length], outfile[name_length];
	bool endian_reply = false;
	int c;
	while ((c = getopt (argc, argv, "i:o:bh")) != -1)
	{
		switch (c)
		{
			case 'i':
			strcpy(optarg, infile);
			break;
			case 'o':
			strcpy(optarg, outfile);
			break;
			case 'b':
			endian_reply = true;
			break;
			case 'h':
			usage(argv[0]);
			return 0;
			break;
			case '?':
			if(optopt == 'i')
				cerr << "Option 'i' requires an argument" << endl;
			else if(optopt == 'o')
				cerr << "Option 'o' requires an argument" << endl;
			else if (isprint (optopt))
				cerr << "Unknown Option -" << optopt << endl;
			else
				cerr << "Unknown option character `\\x" << hex << optopt << "'." << endl;
			return 1;
			break;
			default:
			usage(argv[0]);
			return 1;
			break;
		}
	}
	
/*
------------------------------------------------------------------------------
 Interactively enter starting information at run time.
------------------------------------------------------------------------------
*/
/*
	cout << "\nEnter name of VNMR file: ";
	cin.getline(infile, name_length);
	cout << "Enter name of ENUF file: ";
	cin.getline(outfile, name_length);
	cout << "Do a big endian/little endian byte swap (y/n)?: ";
    cin.getline((char *) &endian_reply, sizeof(int));
*/
/*
------------------------------------------------------------------------------
 Connect input stream to file containing VNMR data, and terminate program
 if file not found. 
------------------------------------------------------------------------------
*/
	ifstream indata(&infile[0]);
	if(!indata)
	{
		cout << "\nError opening VNMR file.  Program aborted.\n";
		return 1;
	}
/*
------------------------------------------------------------------------------
 Define some variable names.
------------------------------------------------------------------------------ 
*/
	datafilehead filehead;
	datablockhead blockhead;
	hypercmplxbhead hcblockhead;
	int
	I0, I1, I2, numpts, fileheadsize, blockheadsize, number_data_chars,
	numberFIDs, complex_signal_int[max_points];
	fileheadsize = (7*sizeof(int)) + (2*sizeof(short));
	blockheadsize = (4*sizeof(short)) + (4*sizeof(float)) + sizeof(int);
	char first_status_char[2];
	char* hdr_ptr;
	float scale_factor = 1.0, complex_signal_float[max_points];
	short complex_signal_short[max_points];
/*
------------------------------------------------------------------------------
 Read vnmr fileheader and reverse the byte order of individual data if
 necessary. 
------------------------------------------------------------------------------
*/
	indata.read((char *) &filehead, fileheadsize);
	if(endian_reply)
	{
		hdr_ptr = (char*) &filehead.nblocks;
		for(I0 = 0; I0 < 6; I0 += 1)
		{ 
			reverse_byte_order(hdr_ptr, sizeof(int));
			hdr_ptr += sizeof(int);
		}
		for(I0 = 0; I0 < 2; I0 += 1)
		{ 
		reverse_byte_order(hdr_ptr, sizeof(short));
		hdr_ptr += sizeof(short);
		}
		reverse_byte_order(hdr_ptr, sizeof(int));
		hdr_ptr += sizeof(int);
	}
//
	numpts = filehead.np / 2; 
	number_data_chars = filehead.np * sizeof(float); 
	numberFIDs = filehead.nblocks*filehead.ntraces;
	display_VNMR_parameters(&fileheadsize, &numpts, &filehead, &blockhead,
							&hcblockhead);
	memcpy((void *) &first_status_char, (void *) &filehead.status, 2); 
	if ((first_status_char[0] & 0x8) == 0x8)
	{ 
		cout << "\nData are stored as floating point." << endl; 
	}
/* 
------------------------------------------------------------------------------
 Open output file. 
------------------------------------------------------------------------------
*/
	ofstream outdata; 
	outdata.open(outfile, ios::out); 
/*
------------------------------------------------------------------------------
 Read the vnmr block header and data, and then write the data to the 
 output file.  This is repeated block-by-block.
 Sometimes VNMR files have a second header block for hypercomplex data sets.
 When this header is present, filehead.nbheaders is equal to two, and this
 second header block must be read before the actual NMR data.
 Check the filehead.status parameter to determine if the data are single
 precision floating point or integer, and adjust the type conversion
 appropriately.
 Check the filehead.ebytes parameter for two-byte or four-byte integer
 data, and convert input data to floating point.
 VNMR data use a left handed coordinate system to define signal phases..
 To correct for this, the imaginary part of the complex data points (i.e., 
 every second word in the fid data stream) must be negated.
------------------------------------------------------------------------------
*/
	outdata.write((char *) &softwareVersion, 12);
	outdata.write((char *) &numpts, sizeof(int));
	outdata.write((char *) &numberFIDs, sizeof(int));	
	for(I1 = 0; I1 < filehead.nblocks; I1++)
	{
		indata.read((char *) &blockhead, blockheadsize);
		if(endian_reply)
		{
			hdr_ptr = (char*) &blockhead.scale;
			for(I0 = 0; I0 < 4; I0 += 1)
			{ 
				reverse_byte_order(hdr_ptr, sizeof(short));
				hdr_ptr += sizeof(short);
			}
			reverse_byte_order(hdr_ptr, sizeof(int));
			hdr_ptr += sizeof(int);
			for(I0 = 0; I0 < 4; I0 += 1)
			{ 
				reverse_byte_order(hdr_ptr, sizeof(float));
				hdr_ptr += sizeof(float);
			}
		}
/* 
------------------------------------------------------------------------------
 Read the hypercomplex block if present. The code does not yet implement the
 optional byte reversal for the hypercomplex block.
------------------------------------------------------------------------------
*/
		if (filehead.nbheaders > 1)
		{
			indata.read((char *) &hcblockhead, blockheadsize);
		}
//
		for(I0 = 0; I0 < filehead.ntraces; I0++) 
		{
            if ((first_status_char[0] & 0x8) == 0x8)
			{ 
                indata.read((char *) complex_signal_float, filehead.tbytes); 
                for(I2 = 0; I2 < filehead.np; I2++)
				{
					if(endian_reply)
					{
						reverse_byte_order((char *) &complex_signal_float[I2],
										   sizeof(float));
					}
					complex_signal_float[I2] *= (float) pow(-1., (double) I2);  
				} 
			}
            else  
			{
				scale_factor = (float) pow(2, (double) blockhead.scale);
                if(filehead.ebytes == 2) 
				{ 
					indata.read((char *) complex_signal_short, filehead.tbytes); 
					for(I2 = 0; I2 < filehead.np; I2++) 
					{
						if(endian_reply)
						{
							reverse_byte_order((char *) &complex_signal_short[I2],
											   sizeof(short));
						}
						complex_signal_float[I2] = ((float) complex_signal_short[I2])
						* scale_factor * ((float) pow(-1., (double) I2));
					} 
				}
                else
				{
					indata.read((char *) complex_signal_int, filehead.tbytes); 
					for(I2 = 0; I2 < filehead.np; I2++)
					{
						if(endian_reply)
						{
							reverse_byte_order((char *) &complex_signal_int[I2],
											   sizeof(int));
						}
						complex_signal_float[I2] = ((float) complex_signal_int[I2])
						* scale_factor * ((float) pow(-1., (double) I2)); 
					}
				}
			}
            outdata.write((char *) complex_signal_float, number_data_chars);
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
