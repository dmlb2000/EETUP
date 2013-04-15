/*
 ------------------------------------------------------------------------------- 
 PROGRAM NAME: Quad_NMR_Simulator.cc
 AUTHOR:  Herman Cho 
 CREATION DATE:  30 October 2012 
 LAST MODIFICATION DATE:  
 DESCRIPTION:  Program to calculate high-field spectra of quadrupolar 
 nuclides with chemical shift, quadrupolar, and heteronuclear dipolar
 interactions present.  
 Main inputs to this program are the spin system, pulse, and phase cycle 
 parameters.
 The user can specify a subset of the angle pairs in the computation of
 the averaged spectrum.  This makes it possible to parallelize a powder
 average calculation by splitting up the powder average among separate 
 processors.
 Gaussian (cgs) units are used throughout. 
 The simulated interferogram is written to a Felix format file "outfile." 
 Program must be linked to GAMMA object library.  
 Source code references only the mathematical shell of the GAMMA library.  
 When porting this program from hyde to king, some arrays needed to be
 reduced in size, or the executable code would not work properly.
 ------------------------------------------------------------------------------
 */
#include "iomanip" 
#include "gamma.h" 
#include "fftw3.h"
#include "Cho_1.h"
using namespace std;
complex i(0,1);
const double pi = 3.14159265358979, 
hbar = 1.054588757E-27, 
e_charge = 4.80324238E-10,  
twopi = 2.0*pi; 
const int WordLength = 160; 
//
int main( int argc, char* argv[])
{
	double 
	dwell, IAngMom, omega_L, center_ppm, sw_Hz, GammaIbar, GammaSbar;
	int 
	I0, I1, I2, numpts, NumberSSpins, progressnumber = 10000;
	char
	Infile[WordLength], ISNetworkFile[WordLength], outfile[WordLength],
	OrientationFile[WordLength];
	char 
	*outfile_ptr;
/*
 ------------------------------------------------------------------------------
 Interactively enter starting information at run time.  The following
 data are requested:
 - Infile: name of file containing spin system data.  For the observe spin,
   this includes the Larmor frequency and EFG and chemical shift data.
   For the coupled S-spin, this contains the number of S spins, all of the
   I-S internuclear distances, and the two angles relating the I-S
   internuclear vectors to the I-spin's EFG PAS.
 - ISNetworkFile: name of file containing spin system isotope information
 - outfile: name of Felix format file where simulated interferogram is
   saved
 - omega_L: Larmor frequency, in Hz.  This quantity is immediately converted
   to radians/sec (units of angular frequency)
 ------------------------------------------------------------------------------
*/
    cout << "Enter acquisition dwell time, in seconds: ";
    cin >> dwell;
    sw_Hz = 1.0 / dwell;
//
    cout << "Enter number of points to acquire: ";
    cin >> numpts;
//
    cout << "Enter frequency at 0 ppm, in Hz: ";
	cin >> omega_L;
	omega_L *= twopi;
//
    cout << "Enter spectrum center position, in ppm: ";
	cin >> center_ppm;
//
	cout << "Enter name of file containing network isotope information: ";
	cin.get(); 
	cin.getline(ISNetworkFile, WordLength);
//
	cout << "Enter name of file containing EFG and chemical shift data: ";
	cin.getline(Infile, WordLength);
//
	cout << "Enter output filename: ";
	cin.getline(outfile, WordLength);
	outfile_ptr = &outfile[0]; 
/*
 ------------------------------------------------------------------------------
 Connect input streams to files containing spin system data, isotope data,
 orientation data, and return error message and terminate program if files
 are not found.
 The isotopes and number of spins are specified in the gamma spin_sys file
 ISNetworkFile.  The observe isotope must be the last nucleus listed in
 ISNetworkFile, and it must be quadrupolar.  All nuclei before the observe
 nucleus (I) are assumed to be of another isotope (S).  The data in
 ISNetworkFile are used to define the spin_sys object called ISNetwork.
 Input stream for spin system data is "indata", and for orientation data
 "orientationdata."
 Number of phase cycle steps is read from phase cycle file and number of
 S spins is read from spin system file.
 Note that gamma stores gyromagnetic ratios in SI units (rad/sec/Tesla),
 whereas this program uses cgs units.  To convert to the cgs unit of
 magnetic induction (Gauss), the gyromagnetic ratios must be divided by 1e04.
 ------------------------------------------------------------------------------
*/
	spin_sys ISNetwork;  
	ISNetwork.read(ISNetworkFile); 
	NumberSSpins = ISNetwork.spins() - 1;
	IAngMom = ISNetwork.qn(NumberSSpins); 
	GammaIbar = ISNetwork.gamma(NumberSSpins) / 10000.;
	GammaSbar = ISNetwork.gamma(0) / 10000.;
	cout << " GammaI = " << GammaIbar/twopi << "   GammaS = " << 
	GammaSbar/twopi << "\n";
	// 
	ifstream indata(&Infile[0]);
	if(!indata)
	{
        cout << "\nError opening spin system file.  Program aborted.\n";
        return 1;
	}
/*
 ------------------------------------------------------------------------------
 Define variables for calculations.
 ------------------------------------------------------------------------------
*/
	const int 
	max_data = 524288, maxangles = 65536;
	double  
	theta[maxangles], phi[maxangles], eta, AQ, Vxx, Vyy, Vzz, Q_moment, 
	Q_SpectralWidth, CS_SpectralWidth, Total_SpectralWidth, delta_xx,
    delta_yy, delta_zz, alp, bet, gam, cos_bet, sin_bet, sin_2bet, cos_2gam,
    iso_cs, delta_cs, eta_cs,
	ISprefactor[NumberSSpins], rho[NumberSSpins],
	cos2psi[NumberSSpins], sin2psi[NumberSSpins],
	sin_2psi[NumberSSpins], psi, ISdistance, FirstSphericalHarmonic; 
	int 
	orientations, numberchars, reply, UpperBound, LowerBound; 
	gen_op 
	Hint, H_csiso, rho0, rho1, 
	Ud1, i_td1_Hint, Ud2, i_td2_Hint, i_dwell_Hint, Udwell, detectop;
	complex 
	i_dwell, spur;

	fftw_complex *data;
	fftw_plan fp;
/*
 ------------------------------------------------------------------------------
 Enter orientation data as (theta,phi) pairs, either from a file or manually.
 The program allows tracking of the progress of the calculation by printing
 a short message every time the number of orientations evaluated reaches a
 multiple of the parameter progressnumber.
 The parameters LowerBound and UpperBound specify the range of orientations
 within the (theta,phi) array to be used in the calculation.
 ------------------------------------------------------------------------------
*/
	cout << 
	"Will orientation data be entered now (0) or read from a file (1)?: ";
	cin >> reply;
	//
	if(reply == 1) 
	{
        cout << "Enter name of file with orientation data: ";
        cin.get();
        cin.getline(OrientationFile, WordLength);
        cout << "Enter progress tracking number: ";
        cin >> progressnumber;
        ifstream orientationdata(&OrientationFile[0]);
        if(!orientationdata)
		{
			cout << "\nError opening orientation data file.  Program aborted.\n";
			return 1;
		}
        orientationdata.read((char *) &orientations, sizeof(int));
        numberchars = sizeof(double) * orientations; 
        cout << "There are " << orientations
		<< " orientations in the file." << endl;
        cout << "Enter starting orientation number (0 or higher): ";
        cin >> LowerBound;
        cout << "Enter final orientation number (" 
		<< orientations-1 << " or lower): ";
        cin >> UpperBound;
        UpperBound += 1;
        orientationdata.read((char *) theta, numberchars); 
        orientationdata.read((char *) phi, numberchars); 
        orientationdata.close(); 
	}
	else
	{
        cout << "Enter number of orientations: ";
        cin >> orientations;
        LowerBound = 0;
        UpperBound = orientations;
	double h, kp;
	int n = orientations;
	double p = 1.0/2.0;
	double a = 1 - 2*p/(n-3);
	double b = p*(n+1)/(n-3);
	double r[maxangles];
	r[0]     = 0.0;
	theta[0] = PI;
	phi[0]  = 0;
        for(I0 = 1; I0 < n-1; I0++)
        {
		kp = a*I0 + b;
		h = -1 + 2*(kp-1)/(n-1);
		r[I0]     = sqrt(1-(h*h));
		theta[I0] = acos(h);
		phi[I0] = fmod((phi[I0-1] + 3.6/sqrt(n)*2/(r[I0-1]+r[I0])), (2*pi));
        }
	theta[n-1] = 0;
	phi[n-1]   = 0;
#if 0

        for(I0 = 0; I0 < orientations; I0++)
		{
			cout << endl << "Enter theta (in deg) for orientation number " 
			<< I0 << ": "; 
			cin >> theta[I0]; 
			cout << "Enter phi (in deg) for orientation number " << I0 << ": ";
			cin >> phi[I0];
			theta[I0] *= pi / 180.;
			phi[I0] *= pi / 180.;
		}
#endif
	}
/*
 ------------------------------------------------------------------------------
 Read file containing observe spin system data, starting with chemical shift
 information.
 Principal components of the chemical shift tensor are used to compute
 asymmetry parameters, isotropic shifts, and isotropic Hamiltonian, H_csiso.
 Data must be ordered in the following way:
 - xx component of CS tensor, in ppm
 - yy component of CS tensor, in ppm
 - zz component of CS tensor, in ppm
 - Euler angles alpha, beta, and gamma relating CS-PAS to the quadrupolar PAS
   (in degrees).
 After reading chemical shift data, the offset is subtracted and the shifts
 are converted to units of Hz.
 The isotropic shift Hamiltonian is created in units of radians/sec (angular
 frequency units).
------------------------------------------------------------------------------
*/
	indata >> delta_xx >> delta_yy >> delta_zz;
	indata >> alp >> bet >> gam;
    delta_xx += -center_ppm;
    delta_yy += -center_ppm;
    delta_zz += -center_ppm;
    delta_xx *= omega_L / 1.0e6;
    delta_yy *= omega_L / 1.0e6;
    delta_zz *= omega_L / 1.0e6;
	alp *= pi / 180.0;
	bet *= pi / 180.0;
	gam *= pi / 180.0;
	cos_bet = cos(bet);
	sin_bet = sin(bet);
	sin_2bet = sin(bet+bet);
	cos_2gam = cos(gam+gam);
	iso_cs = (delta_xx+delta_yy+delta_zz)/3.0;
	if(delta_zz == delta_yy)
	{
        delta_cs = 0.0;
        eta_cs = 0.0;
	}
	else
	{
        delta_cs = delta_zz-iso_cs;
        eta_cs = (delta_yy-delta_xx)/delta_cs;
	}
	H_csiso = twopi * iso_cs * Iz(ISNetwork,0); 
/*
 ------------------------------------------------------------------------------
 Continue reading spin system data from input stream indata.
 Quadrupolar data are next.  Values must be arrayed in the following order:
 - Q_moment: quadrupolar moment of nuclide, in cm^2
 - Vxx: xx-component of quadrupolar EFG in the quadrupolar PAS, in
   statvolts/cm^2
 - Vyy: yy-component of quadrupolar EFG in the quadrupolar PAS, in
   statvolts/cm^2
 - Vzz: zz-component of quadrupolar EFG in the quadrupolar PAS, in
   statvolts/cm^2
   The fundamental charge e_charge is defined above in units of esu, and
   the quantity AQ is calculated in radians/sec (units of angular frequency).
 ------------------------------------------------------------------------------
*/
	indata >> Q_moment >> Vxx >> Vyy >> Vzz; 
	AQ = (e_charge * Vzz * Q_moment ) 
	/ (4.0 * hbar * IAngMom * ((2.0*IAngMom)-1.0)); 
	eta = (Vxx-Vyy) / Vzz;
/*
 ------------------------------------------------------------------------------
 Now read heteronuclear dipolar coupling data, and close input file upon
 completion.
 Data for individual S-spins are arrayed in the following way:
 - I-S internuclear separation in cm (ISdistance)
 - Longitudinal angle psi specifying angle between the I spin's EFG principal
   axis and the IS internuclear vector (in degrees)
 - Azimuthal angle rho (in degrees)
 ------------------------------------------------------------------------------
*/
	for(I0 = 0; I0 < NumberSSpins; I0++)
	{
        indata >> ISdistance >> psi >> rho[I0];
        psi *= pi / 180.0;
        rho[I0] *= pi / 180.0;
        cos2psi[I0] = cos(psi) * cos(psi);
        sin2psi[I0] = sin(psi) * sin(psi); 
        sin_2psi[I0] = sin(psi+psi);
        ISprefactor[I0] = -0.5*GammaIbar*GammaSbar*hbar
		/(ISdistance*ISdistance*ISdistance);
	}
	indata.close();
/*
 ------------------------------------------------------------------------------
 Compute a spectral width that will contain the entire line, including tbe
 first order quadrupolar lineshape.  The maximum inhomogeneous linewdth is
 estimated to be the sum of the width of the quadrupolar powder line and the
 width of the CSA pattern.  The desired spectral width is two times the 
 maximum inhomogeneous linewidth, and is specified by halving the dwell time
 until it is large enough.  The frequency resolution is kept constant by
 doubling the number of acquisition points every time the spectral width is 
 halved.
 ------------------------------------------------------------------------------
*/    
    Q_SpectralWidth = abs(3. * AQ * ((2.0*IAngMom)-1.0) / pi);
    CS_SpectralWidth  = abs(delta_zz - delta_xx) / (2. * pi);
    Total_SpectralWidth = 2. * (Q_SpectralWidth + CS_SpectralWidth);
    while (1./dwell < Total_SpectralWidth)
	{
		dwell *= 0.5;
		numpts *= 2;
	}
/*
 ------------------------------------------------------------------------------
 Initialize variables, and calculate expressions that appear frequently
 to avoid repetition.
 ------------------------------------------------------------------------------
*/
	i_dwell = i * dwell; 
	rho0 = Ix(ISNetwork,NumberSSpins);
	detectop = Ix(ISNetwork,NumberSSpins) + (i*Iy(ISNetwork,NumberSSpins));
/*
 ------------------------------------------------------------------------------
 Null the one-dimensional vector that will be used to store the simulated FID.
 ------------------------------------------------------------------------------
*/
	data = (fftw_complex *)malloc(sizeof(fftw_complex)*numpts);
	for(I0 = 0; I0 < numpts; I0++)
	{
       	data[I0][0] = 0.0;
        data[I0][1] = 0.0;
	}
/*
 ------------------------------------------------------------------------------
 Initialize some parameters, and start the master I0 for loop, which iterates
 through the different orientations.
 ------------------------------------------------------------------------------
*/
	for(I0 = LowerBound; I0 < UpperBound; I0++)
	{
        Hint = H_csiso; 
/*
 ------------------------------------------------------------------------------
 Construct anisotropic part of the chemical shift interaction and add it to
 the internal Hamiltonian, Hint.  Note that this Hamiltonian has been
 initialized in the line above by setting it equal to the isotropic shift
 Hamiltonian.
 ------------------------------------------------------------------------------
*/
        Hint += Hint_cs(theta[I0],phi[I0],eta_cs,delta_cs,alp,bet,gam, 
                        cos_bet,sin_bet,cos_2gam,sin_2bet,ISNetwork,
                        NumberSSpins);  
/*
 ------------------------------------------------------------------------------
 Define quadrupolar Hamiltonian (evaluated to first order) and add it to
 Hint.
 ------------------------------------------------------------------------------
*/
        Hint += Hint_Q_0(theta[I0],phi[I0],eta,IAngMom,AQ,omega_L,
                         ISNetwork,NumberSSpins);  
/*
 ------------------------------------------------------------------------------
 Construct heteronuclear dipolar interaction Hamiltonian (if needed) term by
 term and add it to Hint.
 ------------------------------------------------------------------------------
*/
        if (NumberSSpins > 0)
		{
            FirstSphericalHarmonic = (3.0*cos(theta[I0])*cos(theta[I0]))-1.0;
            for(I2 = 0; I2 < NumberSSpins; I2++)
			{
				Hint += ISprefactor[I2]
				* (Iz(ISNetwork,I2)*Iz(ISNetwork,NumberSSpins))
				* (
				   (FirstSphericalHarmonic*((3.0*cos2psi[I2])-1.0))
				   +(3.0*sin_2psi[I2]*sin(theta[I0]+theta[I0])
					 *cos(phi[I0]-rho[I2]))
				   +(3.0*sin2psi[I2]*sin(theta[I0])*sin(theta[I0])
					 *cos(rho[I2]+rho[I2]-phi[I0]-phi[I0]))
				   );
			} 
		}
/*
 ------------------------------------------------------------------------------
 Define free evolution propagator.
 ------------------------------------------------------------------------------
*/
        i_dwell_Hint = i_dwell * Hint;
        Udwell = exp(-i_dwell_Hint);
        rho1 = rho0;
        spur = trace(rho1,detectop);
        data[0][0] += Re(spur);
        data[0][1] += Im(spur);
/*
 ------------------------------------------------------------------------------
 Each I1 loop is a single point acquisition of the FID.
 ------------------------------------------------------------------------------
*/
        for(I1 = 1; I1 < numpts; I1++)
			{
				rho1.sim_trans_ip(Udwell);
				spur = trace(rho1,detectop);
				data[I1][0] += Re(spur); 
				data[I1][1] += Im(spur); 
			}
        if(fmod((double) I0, (double) progressnumber) == 0)
		{
            cout << "lower bound = " << LowerBound 
			<< "; upper bound = " << UpperBound 
			<< "; count = " << I0 << endl;
		}
	}
/*
 ------------------------------------------------------------------------------
 Write the felix output file.
 ------------------------------------------------------------------------------
*/
	/*felix_time_out_new(numpts, dwell, omega_L, outfile_ptr, twopi, 
					   real_part, imag_part); 
	*/
	//for(I0 = 0; I0 < numpts; I0++)
	//{
	//	if(I0&1)
	//	{
	//		data[I0][0] = -data[I0][0];
	//		data[I0][1] = -data[I0][1];
	//	}
	//}
	/*fp = fftw_plan_dft_1d(numpts, data, data, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(fp); */
	/* print to gnuplot stuff */
	ofstream output(outfile_ptr);
	//output << "Real,Imaginary" << endl;
	cout << numpts << endl;
	float tmp_float;
	for(I0 = 0; I0 < numpts; I0++)
	{
		tmp_float = (float)data[I0][0];
		output.write((char *)&(tmp_float), sizeof(float));
	}
	for(I0 = 0; I0 < numpts; I0++)
	{
		tmp_float = (float)data[I0][1];
		output.write((char *)&(tmp_float), sizeof(float));
	}
	//	output << data[I0][0] << "," << data[I0][1] << endl;
	output.close();
	//fftw_free(data);
	//fftw_destroy_plan(fp);
}
