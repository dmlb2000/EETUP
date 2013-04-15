/*
 -------------------------------------------------------------------------------
 PROGRAM NAME: nt2en.cc
 AUTHOR:  Herman Cho, Pacific Northwest National Laboratory
 CREATION DATE:  20 September 2012
 LAST MODIFICATION DATE:
 DESCRIPTION:  Program to convert NTNMR data file to an ENUF file.
 Link to the object library "ENUF_objects.h."  Don't forget to keep the
 version name softwareVersion updated.
 ------------------------------------------------------------------------------
 */
#include "fstream"
#include "iostream"
#include "iomanip"
#include "math.h"
#include "ENUF_objects.h"
using namespace std;
const int name_length = 200;
const int max_points = 1048576;
const char softwareVersion[] = "nt2en.100912a";
//
int main( int argc, char* argv[])
{
    char infile[name_length], outfile[name_length], metadatafile[name_length],
	data_source[name_length], JCAMP_origin[name_length],
	JCAMP_owner[name_length], JCAMP_sample[name_length];
    int endian_reply, file_reply;
/*
------------------------------------------------------------------------------
 Interactively enter starting information at run time.
------------------------------------------------------------------------------
*/
	cout << "Enter file with NTNMR-to-ENUF translation inputs: ";
    cin.getline(data_source, name_length);
    get_file_data(data_source, infile, outfile, JCAMP_origin, JCAMP_owner,
				  JCAMP_sample, &name_length);
	cout << "Do a big endian/little endian byte swap (y/n)?: ";
    cin.getline((char *) &endian_reply, sizeof(int));
    cout << "Write a JCAMP metadata file (y/n)?: ";
    cin.getline((char *) &file_reply, sizeof(int));
	if(file_reply == 'y')
	{
		cout << "Enter metadata filename: ";
		cin.getline(metadatafile, name_length);
	}
/*
------------------------------------------------------------------------------
 Connect input stream to file containing NTNMR data, and terminate program
 if file not found.
------------------------------------------------------------------------------
*/
    cout << infile << endl;
    ifstream indata(infile);
    if(!indata)
    {
        cout << "\nError opening NTNMR file.  Program aborted.\n";
        return 1;
    }
/*
------------------------------------------------------------------------------
 Define variable names.
------------------------------------------------------------------------------
*/
    char  temp_array[16384];
    int  I0, I1, index, numberFIDs;
    float  complex_point_float[max_points];
	TECMAG_structure_header structure_header;
    TECMAG_structure structure;
    TECMAG_data_header data_header;
	TECMAG2_structure_header structure2_header;
	TECMAG2_structure structure2;
	pulse_sequence pulse_program;
/*
------------------------------------------------------------------------------
 Read 20-byte structure header from input file and store in structure_header.
------------------------------------------------------------------------------
*/
    indata.read((char *) structure_header.version_ID, 20);
	if(endian_reply == 'y')
	{
		reverse_byte_order((char *) &structure_header.structure_length, 4);
	}
/*
------------------------------------------------------------------------------
 Read 1024-byte NTNMR structure header.  The structure is read in sections
 since the computer may not store them in a contiguous block in memory.
 The byte orders of individual integers or floats are reversed if necessary.
------------------------------------------------------------------------------
*/
    indata.read(temp_array, 1024);
	if(endian_reply == 'y')
	{
		for(I0 = 0; I0 < 72; I0 += 4)
		{ 
			reverse_byte_order(&temp_array[I0], 4);
		}
		for(I0 = 76; I0 < 120; I0 += 8)
		{
			reverse_byte_order(&temp_array[I0], 8);
		}
		for(I0 = 240; I0 < 336; I0 += 8)
		{
			reverse_byte_order(&temp_array[I0], 8);
		}
		for(I0 = 336; I0 < 344; I0 += 2)
		{
			reverse_byte_order(&temp_array[I0], 2);
		}
		for(I0 = 344; I0 < 348; I0 += 4)
		{
			reverse_byte_order(&temp_array[I0], 4);
		}
		for(I0 = 368; I0 < 372; I0 += 2)
		{
			reverse_byte_order(&temp_array[I0], 2);
		}
		for(I0 = 388; I0 < 400; I0 += 2)
		{
			reverse_byte_order(&temp_array[I0], 2);
		}
		for(I0 = 400; I0 < 424; I0 += 8)
		{
			reverse_byte_order(&temp_array[I0], 8);
		}
		for(I0 = 440; I0 < 464; I0 += 8)
		{
			reverse_byte_order(&temp_array[I0], 8);
		}
		for(I0 = 464; I0 < 536; I0 += 2)
		{
			reverse_byte_order(&temp_array[I0], 2);
		}
		reverse_byte_order(&temp_array[536], 8);
		for(I0 = 544; I0 < 564; I0 += 2)
		{
			reverse_byte_order(&temp_array[I0], 2);
		}		
	}
	memcpy((void *) structure.npts, (void *) &temp_array[0], 76);
	memcpy((void *) &structure.magnet_field, (void *) &temp_array[76], 164);
	memcpy((void *) structure.sw, (void *) &temp_array[240], 96);
	memcpy((void *) &structure.spectrum_direction, (void *) &temp_array[336], 8);
	memcpy((void *) &structure.bDigRec, (void *) &temp_array[344], 24);
	memcpy((void *) &structure.transmitter_gain, (void *) &temp_array[368], 20);
	memcpy((void *) &structure.set_spin_rate, (void *) &temp_array[388], 12);
	memcpy((void *) &structure.lock_freq_mhz, (void *) &temp_array[400], 40);
	memcpy((void *) &structure.set_temperature, (void *) &temp_array[440], 24);
	memcpy((void *) structure.shims, (void *) &temp_array[464], 72);
	memcpy((void *) &structure.shim_FWHM, (void *) &temp_array[536], 8);
	memcpy((void *) &structure.HH_dcpl_attn, (void *) &temp_array[544], 480);
/*
------------------------------------------------------------------------------
 Read data header.
------------------------------------------------------------------------------
*/
	indata.read((char *) data_header.DATA_TAG, 12);
	if(endian_reply == 'y')
	{
		reverse_byte_order((char *) &data_header.DATA_length, 4);
	}
/*
------------------------------------------------------------------------------
 Open ENUF output file.
------------------------------------------------------------------------------
*/
    ofstream outdata;
    outdata.open(outfile, ios::out);
/*
------------------------------------------------------------------------------
 Read time-domain data from NTNMR file.  If the NTNMR file holds a
 multidimensional experiment, the traces are read in one at a time.
------------------------------------------------------------------------------
*/
	numberFIDs = (int) (structure.actual_npts[1]*structure.actual_npts[2]
                        *structure.actual_npts[3]);
	outdata.write((char *) &softwareVersion, 12);
	outdata.write((char *) &structure.actual_npts[0], 4);
	outdata.write((char *) &numberFIDs, sizeof(int));
    cout << indata.tellg() << endl;
    cout << numberFIDs << endl;
    cout << structure.npts[0] << endl;
    for(I0 = 0; I0 < numberFIDs; I0++)
    {
        for(I1 = 0; I1 < structure.npts[0]; I1++)
        {
            index = I1 * 2;
            indata.read((char *) &complex_point_float[index], 4);
			index += 1;
            indata.read((char *) &complex_point_float[index], 4);
        }
        index += 1;
		if(endian_reply == 'y')
		{
			for(I1 = 0; I1 < 2*structure.npts[0]; I1++)
			{
				reverse_byte_order((char *) &complex_point_float[I1], 4);
			}
		}
/*
------------------------------------------------------------------------------
 NTNMR uses a left handed coordinate system to define signal phases.  To
 correct for this, the imaginary part of the complex data points (i.e., every
 second word in the fid data stream) must be negated.
------------------------------------------------------------------------------
*/
        for(I1 = 0; I1 < structure.npts[0]; I1++)
        {
            complex_point_float[(2*I1)+1] *= -1.;
        }
/*
------------------------------------------------------------------------------
 Write the output file, and loop back to read the next t2 trace (if file is
 multi-dimensional).
------------------------------------------------------------------------------
*/
        outdata.write((char *) complex_point_float, 4*index);
    }
/*
------------------------------------------------------------------------------
 Read end of data file and display metadata to user.  This version of nt2en
 does not read the TECMAG2 structure and pulse sequence parameters in full
 at the moment.
------------------------------------------------------------------------------
*/
    indata.read((char *) structure2_header.TMAG2_tag, 12);
	if(endian_reply == 'y')
	{
		reverse_byte_order((char *) &structure2_header.structure2_length, 4);
    }
    indata.read((char *) &structure2.placeholder[0], 2048);
    indata.read((char *) pulse_program.PSEQ_tag, 16);
	cout << "\nENUF translator program: " << softwareVersion;
	display_NTNMR_parameters(&structure_header, &structure, &data_header,
							 &pulse_program);
    if(file_reply == 'y')
	{
        write_metadata(&structure_header, &structure, &data_header,
                       &pulse_program, metadatafile, infile, softwareVersion,
					   JCAMP_origin, JCAMP_owner, JCAMP_sample);
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
