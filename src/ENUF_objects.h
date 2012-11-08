/*
 -------------------------------------------------------------------------------
 PROGRAM NAME: ENUF_objects.h
 AUTHOR:  H. Cho, Pacific Northwest National Laboratory
 CREATION DATE: 20 September 2012
 LAST MODIFICATION DATE:
 DESCRIPTION:  Code for NMR data translation software.
 ------------------------------------------------------------------------------
 */
/*
 * Include Statements
 */
#include <cstring>


using namespace std;
/*
 -------------------------------------------------------------------------------
 The structures datafilehead, datablockhead, and hypercmplxbhead describe the
 format of a VNMR fid file header.
 ------------------------------------------------------------------------------
 */
struct datafilehead
{
    int nblocks;
    int ntraces;
    int np;
    int ebytes;
    int tbytes;
    int bbytes;
    short vers_id;
    short status;
    int nbheaders;
};
struct datablockhead
{
    short scale;
    short status;
    short index;
    short mode;
    int ctcount;
    float lpval;
    float rpval;
    float lvl;
    float tlt;
};
struct hypercmplxbhead
{
    short s_apare1;
    short status;
    short s_spare2;
    short s_spare3;
    int l_spare1;
    float lpval1;
    float rpval1;
    float f_spare1;
    float f_spare2;
};
/*
------------------------------------------------------------------------------
 AUTHOR:  Herman Cho
 CREATION DATE:  4 October 2012
 LAST MODIFICATION DATE:
 This function displays VNMR parameters.
------------------------------------------------------------------------------
*/
void display_VNMR_parameters(int *fileheadsize, int *numpts,
                             datafilehead *filehead,
                             datablockhead *blockhead,
                             hypercmplxbhead *hcblockhead)
{
    cout
    << "\nFile header size (bytes) = " << *fileheadsize
	<< "\nNumber of blocks = " << (*filehead).nblocks
	<< "\nNumber of traces per block = " << (*filehead).ntraces
	<< "\nNumber of elements per trace = " << (*filehead).np
	<< "\nNumber of points, direct dimension  = " << *numpts
	<< "\nNumber of bytes per element = " << (*filehead).ebytes
	<< "\nNumber of bytes per trace = " << (*filehead).tbytes
	<< "\nNumber of block headers per block = "
	<< (*filehead).nbheaders
	<< "\nFileheader status word = " << (*filehead).status
	<< endl;
};
/*
------------------------------------------------------------------------------
 The structures below are in NTNMR .tnt files.
------------------------------------------------------------------------------
*/
struct TECMAG_structure_header
{
    char version_ID[8];
    char TMAG_tag[4];
    bool bool_value1[4];
    int structure_length;
};
struct TECMAG_structure
{
    int npts[4];
    int actual_npts[4];
    int acq_points;
    int npts_start[4];
    int scans;
    int actual_scans;
    int dummy_scans;
    int repeat_times;
    int sadimension;
    char space1[4];
    
    double magnet_field;
    double ob_freq[4];
    double base_freq[4];
    double offset_freq[4];
    double ref_freq;
    double NMR_frequency;
    char space2[44];
    
    double sw[4];
    double dwell[4];
    double filter;
    double experiment_time;
    double acq_time;
    double last_delay;
    short spectrum_direction;
    short hardware_sideband;
    short Taps;
    short Type;
    bool bDigRec;
    int nDigitalCenter;
    char space3[16];
    short transmitter_gain;
    short receiver_gain;
    char space4[16];
    
    unsigned short set_spin_rate;
    unsigned short actual_spin_rate;
    
    short lock_field;
    short lock_power;
    short lock_gain;
    short lock_phase;
    double lock_freq_mhz;
    double lock_ppm;
    double H2O_freq_ref;
    char space5[16];
    
    double set_temperature;
    double actual_temperature;
    
    double shim_units;
    short shims[36];
    double shim_FWHM;
    
    short HH_dcpl_attn;
    short DF_DN;
    short F1_tran_mode[7];
    short dec_BW;
    char grd_orientation[4];
    char space6[296];
    
    char date[32];
    char nucleus[16];
    char nucleus_2D[16];
    char nucleus_3D[16];
    char nucleus_4D[16];
    char sequence[32];
    char lock_solvent[16];
    char lock_nucleus[16];
};
struct TECMAG_data_header
{
    char DATA_TAG[4];
    bool bool_value2[4];
    int DATA_length;
};
struct TECMAG2_structure_header
{
    char TMAG2_tag[4];
    bool bool_value3[4];
    int structure2_length;
};
struct TECMAG2_structure
{
	char placeholder[2048];
};
struct pulse_sequence
{
	char PSEQ_tag[4];
	bool bool_value4[4];
	char Sequence_ID[8];
};
/*
------------------------------------------------------------------------------
 AUTHOR:  Herman Cho
 CREATION DATE:  22 July 2004
 LAST MODIFICATION DATE:
 This function takes a character array N_chars chars in length, and reverses
 the order of the chars.  This operation is necessary for conversions of
 unformatted data from PC to Unix-based systems.
------------------------------------------------------------------------------
 */
void reverse_byte_order( char *byte_array, int N_chars)
{
	char
	dummy[N_chars];
	int
	I0;
/*
------------------------------------------------------------------------------
 Create a char array called dummy with the reversed byte order of
 byte_array.  The re-ordering of the first N_bytes-1 chars is performed
 in the loop, and the last outside the loop in order to avoid
 incrementing the byte_array pointer past the end of the array.
------------------------------------------------------------------------------
*/
	for(I0 = N_chars - 1; I0 >= 1; I0--)
	{
		dummy[I0] = *byte_array;
		byte_array++;
	}
	dummy[0] = *byte_array;
/*
------------------------------------------------------------------------------
 Reset the byte_array pointer to the start of the array, and overwrite
 the objects pointed to by byte_array with dummy.
------------------------------------------------------------------------------
 */
	byte_array -= N_chars - 1;
	for(I0 = 0; I0 <= N_chars - 2; I0++)
	{
		*byte_array = dummy[I0];
		byte_array++;
	}
	*byte_array = dummy[N_chars-1];
}
/*
------------------------------------------------------------------------------
 AUTHOR:  Herman Cho
 CREATION DATE:  27 September 2012
 LAST MODIFICATION DATE:
 This function displays parameters stored in NTNMR .tnt files.
------------------------------------------------------------------------------
*/
void display_NTNMR_parameters(TECMAG_structure_header *structure_header_start,
                              TECMAG_structure *structure_start,
                              TECMAG_data_header *data_header_start,
                              pulse_sequence *pulse_program)
{
	int I0;
//
	cout << "\nVersion ID = ";
	for(I0 = 0; I0 < 8; I0++)
	{
        if (int ((*structure_header_start).version_ID[I0]) == 0) break;
		cout << (*structure_header_start).version_ID[I0];
	}
//
	cout << "\nSequence ID: ";
	for(I0 = 0; I0 < 8; I0++)
	{
        if (int ((*pulse_program).Sequence_ID[I0]) == 0) break;
		cout << (*pulse_program).Sequence_ID[I0];
	}
//
	cout << "\nExperiment date: ";
	for(I0 = 0; I0 < 32; I0++)
	{
        if (int ((*structure_start).date[I0]) == 0) break;
		cout << (*structure_start).date[I0];
	}
//
    cout << "\nObserve nucleus: ";
	for(I0 = 0; I0 < 16; I0++)
	{
        if (int ((*structure_start).nucleus[I0]) == 0) break;
		cout << (*structure_start).nucleus[I0];
	}
//
	cout
    << "\nStructure length (bytes) = "
	<< (*structure_header_start).structure_length
	<< "\nTotal number of data bytes = "
	<< (*data_header_start).DATA_length
    << "\nMagnetic field = " << (*structure_start).magnet_field
	<< endl;
}
/*
 ------------------------------------------------------------------------------
 AUTHOR:  Herman Cho
 CREATION DATE:  17 October 2012
 LAST MODIFICATION DATE:
 This function writes metadata from an NTNMR .tnt file into a JCAMP-DX file.
 ------------------------------------------------------------------------------
 */
void write_metadata(TECMAG_structure_header *structure_header_start,
                    TECMAG_structure *structure_start,
                    TECMAG_data_header *data_header_start,
                    pulse_sequence *pulse_program,
                    char *metadatafile,
					char *infile,
                    const char *softwareVersion,
					char *JCAMP_origin,
					char *JCAMP_owner,
					char *JCAMP_sample)
{
	int I0;
    ofstream metadata;
    metadata.open(metadatafile, ios::out);
    //
	metadata
	<< "##TITLE= NMR parameters from NTNMR"
	<< "\n##JCAMP-DX= 5.00   $$SOFTWARE SOURCE= " << softwareVersion
	<< "\n##DATA TYPE= NMR parameters"
	<< "\n##ORIGIN= " << JCAMP_origin
	<< "\n##OWNER= " << JCAMP_owner
    << "\n##SPECTROMETER/DATA SYSTEM= Tecmag Discovery 300"
	<< "\n##SAMPLE DESCRIPTION= " << JCAMP_sample
	<< "\n##DATE= ";
	for(I0 = 0; I0 < 10; I0++)
	{
        if (int ((*structure_start).date[I0]) == 0) break;
		metadata << (*structure_start).date[I0];
	}
    metadata << "\n##TIME= ";
	for(I0 = 12; I0 < 20; I0++)
	{
        if (int ((*structure_start).date[I0]) == 0) break;
		metadata << (*structure_start).date[I0];
	}
    metadata.setf(ios_base::fixed,ios_base::floatfield);
	metadata
	<< "\n##SOURCE REFERENCE= " << infile
    << "\n##.OBSERVE FREQUENCY (MHz)= "
    << (*structure_start).ob_freq[0]
    << "\n##$OBSERVE CHANNEL (MHz)= CH1"
    << "\n##.OBSERVE NUCLEUS= ";
    for(I0 = 0; I0 < 16; I0++)
	{
        if (int ((*structure_start).nucleus[I0]) == 0) break;
		metadata << (*structure_start).nucleus[I0];
	}
    metadata
    << setprecision(6)
    << "\n##$FREQUENCY CH1 (MHz)= "
    << (*structure_start).ob_freq[0]
    << "\n##$FREQUENCY CH2 (MHz)= "
    << (*structure_start).ob_freq[1]
    << "\n##$FREQUENCY CH3 (MHz)= "
    << (*structure_start).ob_freq[2]
    << "\n##$FREQUENCY CH4 (MHz)= "
    << (*structure_start).ob_freq[3]
    << setprecision(1)
	<< "\n##$SPECTRAL WIDTH D1 (Hz)= "
    << (*structure_start).sw[0]
    << "\n##$SPECTRAL WIDTH D2 (Hz)= "
    << (*structure_start).sw[1]
    << "\n##$SPECTRAL WIDTH D3 (Hz)= "
    << (*structure_start).sw[2]
    << "\n##$SPECTRAL WIDTH D4 (Hz)= "
    << (*structure_start).sw[3]
    << "\n##$POINTS ACQUIRED D1= "
	<< (*structure_start).actual_npts[0]
    << "\n##$POINTS ACQUIRED D2= "
	<< (*structure_start).actual_npts[1]
    << "\n##$POINTS ACQUIRED D3= "
	<< (*structure_start).actual_npts[2]
    << "\n##$POINTS ACQUIRED D4= "
	<< (*structure_start).actual_npts[3]
	<< "\n##$TEMPERATURE SET POINT (K)= "
    << (*structure_start).set_temperature
	<< "\n##$TEMPERATURE ACTUAL (K)= "
    << (*structure_start).actual_temperature
	<< "\n##END= " << endl;
    metadata.flush();
    metadata.close();
}
/*
 ------------------------------------------------------------------------------
 AUTHOR:  Herman Cho
 CREATION DATE:  22 October 2012
 LAST MODIFICATION DATE:
 This function reads a descriptor for an NTNMR file.
 ------------------------------------------------------------------------------
 */
void get_file_data(char *data_source,
				   char *infile,
				   char *outfile,
                   char *JCAMP_origin,
				   char *JCAMP_owner,
				   char *JCAMP_sample,
				   const int *name_length)
{
	int string_length;
    char dummy_string[*name_length];
    const char null_character = '\0';
    ifstream source_inputs(data_source);
//
    source_inputs.getline(dummy_string, *name_length);
    string_length = (int) source_inputs.gcount() - 2;
    memcpy((void *) infile, (void *) dummy_string, string_length);
    memcpy((void *) &infile[string_length], (void *) &null_character, 1);
//
    source_inputs.getline(dummy_string, *name_length);
    string_length = (int) source_inputs.gcount() - 2;
    memcpy((void *) outfile, (void *) dummy_string, string_length);
    memcpy((void *) &outfile[string_length], (void *) &null_character, 1);
//
	source_inputs.getline(dummy_string, *name_length);
    string_length = (int) source_inputs.gcount() - 2;
    memcpy((void *) JCAMP_origin, (void *) dummy_string, string_length);
    memcpy((void *) &JCAMP_origin[string_length], (void *) &null_character, 1);
//
	source_inputs.getline(dummy_string, *name_length);
    string_length = (int) source_inputs.gcount() - 2;
    memcpy((void *) JCAMP_owner, (void *) dummy_string, string_length);
    memcpy((void *) &JCAMP_owner[string_length], (void *) &null_character, 1);
//
	source_inputs.getline(dummy_string, *name_length);
    string_length = (int) source_inputs.gcount() - 2;
    memcpy((void *) JCAMP_sample, (void *) dummy_string, string_length);
    memcpy((void *) &JCAMP_sample[string_length], (void *) &null_character, 1);
//		
    source_inputs.close();
}
