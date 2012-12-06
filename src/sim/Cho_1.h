/*
 ------------------------------------------------------------------------------
 FILE NAME: Cho_1.h 
 AUTHOR: Herman Cho, Pacific Northwest National Laboratory
 MODIFICATION DATE: 17 August 2007
 Previous versions of Cho_1.h contained a lot of non-standard C++ code that
 used to be tolerated by the Gnu C++ compiler, but no more.  The non
 standard code has been fixed or removed in this version.
 ------------------------------------------------------------------------------
 */
#include <cstring>

using namespace std;
/*
 ------------------------------------------------------------------------------
 AUTHOR:  H. Cho 
 CREATION DATE:  4 June 1999  
 LAST MODIFICATION DATE:  14 June 1999 
 This function creates a high field quadrupolar Hamiltonian evaluated to
 first order in coherent averaging theory.  The Hamiltonian is evaluated 
 in a coordinate system that is rotated with respect to its principal axis 
 system by the spherical polar angles theta and phi.  Note that the 
 Hamiltonian does not possess cylindrical symmetry, and therefore for 
 complete orientational averaging, the Hamiltonian should be evaluated for 
 orientations differing by a rotation over a third angle corresponding to a 
 rotation about the static field direction.  
 ------------------------------------------------------------------------------
 */ 
inline gen_op Hint_Q( double theta, double phi, double eta,  
					 double ang_mom, double AQ, double omega_L, 
					 spin_sys ID, int numID)
{
	double a20, two_theta, sin_2theta, sin_theta, cos_theta,
	two_phi, sin_2phi, cos_2phi, real_part;           
	complex a21, a21cc, a22, a22cc, imag_part, i(0,1);
	const double threehalf = 1.5; 
	const double sqrt_threehalf = sqrt(threehalf);  
	gen_op H_Q; 
	spin_op IzIz(Iz(ID,numID)*Iz(ID,numID)); 
	//
	two_theta = theta + theta;
	two_phi = phi + phi;
	sin_2theta = sin(two_theta);
	sin_theta = sin(theta);
	cos_theta = cos(theta);
	sin_2phi = sin(two_phi);
	cos_2phi = cos(two_phi);
	// 
	a20 = sqrt_threehalf * ((3.0*cos_theta*cos_theta) - 1.0  
							+ (eta*cos_2phi*sin_theta*sin_theta)); 
	// 
	real_part = (threehalf*sin_2theta)
	- (eta*sin_theta*cos_theta*cos_2phi); 
	imag_part = -i*eta*sin_theta*sin_2phi;    
	a21 = imag_part + real_part;  
	a21cc = -imag_part + real_part; 
	//    
	real_part = (threehalf*sin_theta*sin_theta) 
	+ (0.5*eta*cos_2phi*(1+(cos_theta*cos_theta)));
	imag_part = i*eta*sin_2phi*cos_theta; 
	a22 = imag_part + real_part; 
	a22cc = -imag_part + real_part; 
	// 
	H_Q = (AQ * a20 / (sqrt(6.0))) 
	* ((3.0*IzIz) - (ang_mom*(ang_mom+1.0)*Ie(ID,numID)));  
	H_Q -= (0.5*AQ*AQ/omega_L)  
	* (
	   (0.5*sqrt_threehalf*a21*a20  
		*((4.0*IzIz) - (4.0*Iz(ID,numID)) + Ie(ID,numID))
		*Ip(ID,numID))  
	   + (0.5*sqrt_threehalf*a21cc*a20
		  *((4.0*IzIz) + (4.0*Iz(ID,numID)) + Ie(ID,numID)) 
		  *Im(ID,numID))
	   + (sqrt_threehalf*a22*a20
		  *(Ie(ID,numID) - Iz(ID,numID))
		  *Ip(ID,numID)*Ip(ID,numID)) 
	   - (sqrt_threehalf*a22cc*a20
		  *(Ie(ID,numID) + Iz(ID,numID))
		  *Im(ID,numID)*Im(ID,numID)) 
	   + (a21*a21cc
		  *(
			(((4.0*ang_mom*(ang_mom+1))-1.0)*Ie(ID,numID)) 
			- (8.0*IzIz) 
			)*Iz(ID,numID) 
		  ) 
	   - (a22*a22cc
		  *( 
			(((2.0*ang_mom*(ang_mom+1))-1.0)*Ie(ID,numID)) 
			- (2.0*IzIz) 
			)*Iz(ID,numID) 
		  )  
	   ); 
	// 
	return H_Q; 
}  
/*
 ------------------------------------------------------------------------------
 AUTHOR:  Herman Cho
 CREATION DATE:  4 June 1999 
 LAST MODIFICATION DATE:  14 June 1999 
 This function constructs the anisotropic part of the chemical shift 
 Hamiltonian.  The Hamiltonian is rotated from the chemical shift 
 principal axis system (CS-PAS) into an intermediate coordinate system 
 related by the three Euler angles alp, bet, and gam, and then rotated 
 from this intermediate system into the rotating frame (with the static 
 field aligned along the z-axis) via the rotations of angle theta and 
 phi.  
 The chemical shift Hamiltonian is given in units of radians/sec 
 (angular frequency).  
 ------------------------------------------------------------------------------
 */
inline gen_op Hint_cs( double theta, double phi, double eta, 
					  double delta, double alp, double bet, 
					  double gam, double cos_bet, double sin_bet, 
					  double cos_2gam, double sin_2bet, 
					  spin_sys ID, int numID)
{ 
	double cos_theta, sin_theta, two_theta, sin_2theta; 
	const double pi = 3.14159265358979; 
	gen_op H_cs; 
	// 
	two_theta = theta + theta; 
	sin_2theta = sin(two_theta);  
	sin_theta = sin(theta); 
	cos_theta = cos(theta); 
	// 
	H_cs = 0.5*pi*delta*Iz(ID,numID)*
	(
	 ( 
	  ((3.0*cos_theta*cos_theta)-1.0) 
	  *(((3.0*cos_bet*cos_bet)-1.0)
		-(eta*sin_bet*sin_bet*cos_2gam))
	  )
	 +(
	   sin_2theta*((3.0*cos(phi-alp)*sin_2bet)
				   +(eta*sin_bet*
					 (((1.0+cos_bet)*cos(phi-alp-gam-gam))
					  -((1.0-cos_bet)*cos(phi-alp+gam+gam)))))
	   )
	 +(
	   (sin_theta*sin_theta)*
	   ((3.0*sin_bet*sin_bet*cos(2.0*(phi-alp)))
		-(0.5*eta*
		  (((1.0+cos_bet)*(1.0+cos_bet)*
			cos(2.0*(phi-alp-gam)))+
		   ((1.0-cos_bet)*(1.0-cos_bet)*
			cos(2.0*(phi-alp+gam))))))
	   )
	 );
	// 
	return H_cs; 
} 
/*
 ------------------------------------------------------------------------------
 AUTHOR:  Herman Cho
 CREATION DATE:  26 June 1999
 LAST MODIFICATION DATE:  29 April 2002.
 This function opens and writes a Felix file containing complex time domain 
 data in the old Felix ASCII format. 
 ------------------------------------------------------------------------------
 */
void felix_time_out_ASCII( int numpts, double dwell, double omega_L,  
						  char *outfile, double twopi, double *real_part, 
						  double *imag_part)
{  
	const float rzero = 0.0; 
	const int one1 = 1, 
	zero0= 0; 
	// 
	double
	v_L, swidth;
	int
	I0;
	ofstream outdata; 
	//
	swidth = 1.0 / dwell; 
	v_L = omega_L / twopi;
	//
	outdata.open(outfile, ios::out); 
	outdata.setf( ios::scientific, ios::floatfield); 
	outdata.precision(8); 
	//
	outdata << "params      16" << endl; 
	// 
	outdata << setw(16) << numpts << ',' 
	<< setw(16) << swidth << endl 
	<< setw(16) << one1 << ',' 
	<< setw(16) << v_L << endl 
	<< setw(16) << zero0 << ',' 
	<< setw(16) << rzero << endl 
	<< setw(16) << one1 << ',' 
	<< setw(16) << rzero << endl; 
	// 
	I0 = 5; 
	do {
		outdata << setw(16) << zero0 << ',' 
		<< setw(16) << rzero << endl; 
		I0 += 1; 
	}
	while(I0 < 17); 
	//
	outdata << "data" << setw(10) << numpts;  
	//
	I0 = 0; 
begin:; 
	if (I0 == numpts) 
		goto end; 
	outdata << endl; 
	outdata << setw(16) << *real_part 
	<< setw(15) << *imag_part; 
	I0 += 1; 
	real_part++; 
	imag_part++; 
	//
	if (I0 == numpts) 
		goto end; 
	outdata << setw(15) << *real_part
	<< setw(15) << *imag_part;
	I0 += 1; 
	real_part++; 
	imag_part++; 
	goto begin; 
end:; 
	// 
	outdata << endl; 
	outdata.close(); 
} 
/*
 ------------------------------------------------------------------------------
 AUTHOR:  H. Cho
 CREATION DATE:  29 June 1999
 LAST MODIFICATION DATE:  29 June 1999
 This function creates a high field quadrupolar Hamiltonian evaluated to
 first order in coherent averaging theory, with non-secular terms discarded.
 The Hamiltonian is evaluated in a coordinate system that is rotated with 
 respect to its principal axis system by the spherical polar angles theta 
 and phi.  
 The eigenvalues of this Hamiltonian are the energies calculated in the 
 laboratory frame using second order time-independent perturbation theory. 
 ------------------------------------------------------------------------------
 */
inline gen_op Hint_Q_0( double theta, double phi, double eta,
					   double ang_mom, double AQ, double omega_L,
					   spin_sys ID, int numID)
{
	double a20, two_theta, sin_2theta, sin_theta, cos_theta,
	two_phi, sin_2phi, cos_2phi, real_part;          
	complex a21, a21cc, a22, a22cc, imag_part, i(0,1);
	const double threehalf = 1.5;
	const double sqrt_threehalf = sqrt(threehalf);
	gen_op H_Q; 
	spin_op IzIz(Iz(ID,numID)*Iz(ID,numID));
	//
	two_theta = theta + theta;
	two_phi = phi + phi;
	sin_2theta = sin(two_theta);
	sin_theta = sin(theta);
	cos_theta = cos(theta);
	sin_2phi = sin(two_phi);
	cos_2phi = cos(two_phi);
	//
	a20 = sqrt_threehalf * ((3.0*cos_theta*cos_theta) - 1.0
							+ (eta*cos_2phi*sin_theta*sin_theta));
	//
	real_part = (threehalf*sin_2theta)
	- (eta*sin_theta*cos_theta*cos_2phi);
	imag_part = -i*eta*sin_theta*sin_2phi;
	a21 = imag_part + real_part;
	a21cc = -imag_part + real_part;
	//
	real_part = (threehalf*sin_theta*sin_theta)
	+ (0.5*eta*cos_2phi*(1+(cos_theta*cos_theta)));
	imag_part = i*eta*sin_2phi*cos_theta;
	a22 = imag_part + real_part;
	a22cc = -imag_part + real_part;
	//
	//
	H_Q = (AQ * a20 / (sqrt(6.0)))
	* ((3.0*IzIz) - (ang_mom*(ang_mom+1.0)*Ie(ID,numID)));
	H_Q -= (0.5*AQ*AQ/omega_L)
	* (
	   (a21*a21cc
		*(
		  (((4.0*ang_mom*(ang_mom+1))-1.0)*Ie(ID,numID))
		  - (8.0*IzIz)
		  )*Iz(ID,numID)
		)
	   - (a22*a22cc
		  *(
			(((2.0*ang_mom*(ang_mom+1))-1.0)*Ie(ID,numID))
			- (2.0*IzIz)
			)*Iz(ID,numID)
		  )
	   );
	//
	return H_Q;
}
/*
 ------------------------------------------------------------------------------
 AUTHOR:  Herman Cho
 CREATION DATE:  8 June 2004
 LAST MODIFICATION DATE:
 This function creates a character string in the format of a new Felix 
 preheader, header, and frame.
 The preheader is 4 bytes in length, and the header is 256 4-byte words in 
 length.
 The spectrometer frequencies must be passed in MHz and the dwell times 
 in units of seconds.
 ------------------------------------------------------------------------------
 */
void format_felix_header( int numpts, 
						 float dwell_time_f2, 
						 float direct_frequency, 
						 float dwell_time_f1, 
						 float indirect_frequency,
						 int number_of_frames, 
						 char *header_pointer) 
{
	//
	const int
	first_word = 0x04030201, header_size = 256, version_number = 200,
	frame_size = 32, IEEE_float = 0, unused_one = 1, time_domain = 0, 
	complex_data_type = 1, point_axis_type = 1;
	float
	spectral_width_f2, spectral_width_f1;
	//
	spectral_width_f2 = 1.0 / dwell_time_f2;
	/*
	 ------------------------------------------------------------------------------
	 Write the byte key that precedes all data and the header size word.
	 Afterwards, write header data.
	 ------------------------------------------------------------------------------
	 */
	memcpy(header_pointer, (char *)&first_word, sizeof(int)); 
	header_pointer += 4;
	memcpy(header_pointer, (char *)&header_size, sizeof(int)); 
	header_pointer += 4;
	memcpy(header_pointer, (char *)&number_of_frames, sizeof(int)); 
	header_pointer += 4;
	memcpy(header_pointer, (char *)&IEEE_float, sizeof(int)); 
	header_pointer += 4;
	memcpy(header_pointer, (char *)&unused_one, sizeof(int)); 
	header_pointer += 4;
	memcpy(header_pointer, (char *)&frame_size, sizeof(int));
	header_pointer += 4;
	memcpy(header_pointer, (char *)&version_number, sizeof(int));
	/*
	 ------------------------------------------------------------------------------
	 Write the first frame. 
	 ------------------------------------------------------------------------------
	 */
	header_pointer += 360; 
	memcpy(header_pointer, (char *)&numpts, sizeof(int)); 
	header_pointer += 4;
	memcpy(header_pointer, (char *)&complex_data_type, sizeof(int)); 
	header_pointer += 4;
	memcpy(header_pointer, (char *)&time_domain, sizeof(int));
	header_pointer += 4;
	memcpy(header_pointer, (char *)&point_axis_type, sizeof(int)); 
	header_pointer += 52; 
	memcpy(header_pointer, (char *)&spectral_width_f2, sizeof(float)); 
	header_pointer += 4; 
	memcpy(header_pointer, (char *)&direct_frequency, sizeof(float)); 
	/*
	 ------------------------------------------------------------------------------
	 Write the second frame.
	 ------------------------------------------------------------------------------
	 */
	header_pointer += 60;
	memcpy(header_pointer, (char *)&numpts, sizeof(int));
	header_pointer += 4;
	memcpy(header_pointer, (char *)&complex_data_type, sizeof(int));
	header_pointer += 4;
	memcpy(header_pointer, (char *)&time_domain, sizeof(int));
	header_pointer += 4;
	memcpy(header_pointer, (char *)&point_axis_type, sizeof(int));
	header_pointer += 52;
	memcpy(header_pointer, (char *)&spectral_width_f1, sizeof(float));
	header_pointer += 4;
	memcpy(header_pointer, (char *)&indirect_frequency, sizeof(float));
}
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
 CREATION DATE:  29 April 2002
 LAST MODIFICATION DATE: 20 August 2007
 This function opens and writes a Felix file containing a complex 1-D time
 domain signal in the new Felix format.  The header is specified to be 256
 4-byte words long.  The Larmor frequency omega_L must be passed in units of
 radians/sec.
 The mods on 20 Aug 2007 correct some non-standard C++ code in the prior
 version.
 ------------------------------------------------------------------------------
 */
void felix_time_out_new( int numpts, double dwell, double omega_L,
						char *outfile, double twopi, double *real_part,
						double *imag_part)
{
	//
	const int MaxPoints = 524288, NumberFrames = 1;
	float
	dwell_time_f1, dwell_time_f2, DirectFrequency,
	IndirectFrequency, ComplexPointFloat[MaxPoints];
	int
	I0, index, NumDataWords, NumFelixChars;
	char
	HeaderChars[1032];
	//
	ofstream outdata;
	//
	DirectFrequency = (float) (omega_L / twopi);
	IndirectFrequency = DirectFrequency;
	dwell_time_f2 = (float) dwell;
	dwell_time_f1 = dwell_time_f2;
	//
	outdata.open(outfile, ios::out);
	/*
	 ------------------------------------------------------------------------------
	 Write the Felix preheader, header, and frames.  The function
	 "format_felix_header" is defined further below.
	 ------------------------------------------------------------------------------
	 */
	format_felix_header(numpts, dwell_time_f2, DirectFrequency,
						dwell_time_f1, IndirectFrequency,
						NumberFrames, HeaderChars);
	outdata.write(HeaderChars, 1032);
	/*
	 ------------------------------------------------------------------------------
	 Put the real and imaginary data in a single array, write the array to the
	 Felix file, and then close.
	 ------------------------------------------------------------------------------
	 */
	for(I0 = 0; I0 < numpts; I0++)
	{
		index = I0 * 2;
		ComplexPointFloat[index] = (float) *real_part;
		index += 1;
		ComplexPointFloat[index] = (float) *imag_part;
		real_part++;
		imag_part++;
	}
	NumDataWords = 2 * numpts;
	NumFelixChars = NumDataWords * sizeof(float);
	outdata.write((char *) &NumDataWords, sizeof(int));
	outdata.write((char *) ComplexPointFloat, NumFelixChars);
	outdata.flush();
	outdata.close();
}
