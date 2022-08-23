/* dsm.c
 * triton_dsm simulator
 * Eric Fogleman, Ian Galton -- 1997-2000
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

/* windows platform -- use rngs.c for drand48()  */
#ifdef _WIN32
#include "rngs.h"
#define srand48(x) (PlantSeeds((x)))
#define drand48()  (Random())
#endif

#ifdef _WIN32
#define CP_CMD_STRING "copy "
#else
#define CP_CMD_STRING "cp "
#endif

#define            SUCCESS                    0
#define            FAIL                      -1
#define            NOT_APPLICABLE             0
#define            MAX_TOKEN_SIZE             240
#define            YES                        1
#define            NO                         0
#define            NUM_WARM_UP_SAMPLES        100
#define            PI                         3.1415926

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define ABS(x) (((x) < 0) ? -(x) : (x))
#define SGN(x) (((x) > 0) ? 1 : -1)
#define SQR(x) ((x)*(x))

typedef struct {
    double reference;     /* Reference voltage */
    double offset;        /* Comparator offset voltage */
} FADC_COMPARATOR;        /* State of each comparator in the Flash ADC */

typedef struct {
    int num_comparators;
    double *ref_ladder_resistors;
    double step_size;
    FADC_COMPARATOR *comparators;
} FADC_STATE;             /* Overall state of the Flash ADC */

typedef struct {
  int num_comparators;
  double *ref_ladder_resistors;
  double step_size;
  FADC_COMPARATOR *p_comparators;
  FADC_COMPARATOR *n_comparators;
} DFADC_STATE;  /* state of the differential flash adc */

typedef struct {
    int *input;           /* Pointer to output of previous Switching Block in tree */

    int top_output;       /* Value of top output of Switching Block */
    int bottom_output;    /* Value of bottom output of Switching Block */
    int s_kr;             /* Switching sequence value */

    int integrator1_accum;  /* to implement 2nd order encoder */
    int integrator2_accum;
    long int integrator2_overload_count; 
} DAC_ENCODER_SB;             /* Switching Block state  */


typedef struct {
    int num_layers;
    int num_one_bit_dacs;
    int num_switching_blocks;
    int use_second_order_mmdac;
    int skr_int2_max;
    int skr_int2_min;
    DAC_ENCODER_SB **sb;  /* 2-d array of DAC_ENCODER_SBs */
} DAC_ENCODER_STATE;

typedef struct {
    int state;            /* State of dac: 1 == high, 0 == low */
    double dac_error_hi;  /* High state DAC error */
    double dac_error_lo;  /* Low state DAC error */
} DAC_ELEMENT;   /* State of each one bit DAC element*/

typedef struct {
    int num_one_bit_dacs;

    double nominal_one_bit_dac_hi;
    double nominal_one_bit_dac_lo;
    DAC_ELEMENT *one_bit_dacs;
} DAC_BANK_STATE;

typedef struct {
    int s_kr;
    int tri_state_accum;
} REQUANTIZER;


typedef struct {
    DAC_BANK_STATE dac;
    double integrator_gain;
    double integrator_output;
} DSM_STAGE;                      /* State of each delta-sigma modulator stage */


typedef struct {
    int num_stages;
    int num_requantizers;
    double adc_step_size;
    DSM_STAGE *stages;
    DFADC_STATE adc;
    DAC_ENCODER_STATE dac_encoder;
    REQUANTIZER *requantizers; 
} DSM_STATE;           /* Overall state of the pipelined ADC */


typedef struct {
    long int n;                           /* Master sample index */
    int num_periodogram_bins;
    int num_periodogram_averages;  
    int per_no;                           /* Periodogram number */
    int bin_no;                           /* Bin number within periodogram */
    double input_freq;

    /* Pointers to buffers for ADC output and psd estimates
     */
    double *dsm_output;
    double *dsm_noise;
    double *dsm_output_psd;
    double *dsm_noise_psd;

    /* Pointers to buffers for psd estimation process
     */
    double *real, *imag;
    double *psd_window, psd_U;
    double *dc_term_fft_real, *dc_term_fft_imag;
    FILE *output_fp;
} SIMULATION_STATE;     /* Simulation state information */



typedef struct {
    int num_stages;

    struct parms_stage_info     {
        double integrator_gain;
        double integrator_gain_std_dev;
        double dac_step_size;
        double dac_step_size_std_dev;
    } *stage_info;                        /* The number of stages is user-specified */

    int use_second_order_mmdac;  /* 1 ==> Jared's mmdac, 0 ==> first order */
    int num_skr_int2_bits;  /* number of bits in second order mmdac second integrator */
    int num_requantizers;                 /* Number 1-bit requantizers following ADC */
    
    int num_adc_comparators;
    double max_adc_input;
    double adc_unit_res_std_dev;
    double adc_comp_offset_std_dev;

    double input_amplitude;                
    double input_freq;                    /* sample-rate = 1 */
    double input_offset;
    double cm_amplitude;                
    double cm_freq;                    /* sample-rate = 1 */
    double cm_offset;
    int num_periodogram_bins;
    int num_periodogram_averages;
    int oversampling_ratio;               /* For SINAD calculation */
    int random_number_seed;

    int argc;                             /* Command line ...        */
    char **argv ;                         /*                 ...info */
} PARMS;                                  /* User Specified Parameters */

extern double gaussian_rv();
extern double drand48();
/* extern double log(); */
/* extern double log10(); */

void copy();
void fft();
void zero_fill();


main(int argc, char *argv[]) {
    PARMS parms;
    SIMULATION_STATE sim;
    DSM_STATE dsm;

    char system_command_buffer[256];
    int ret = SUCCESS;
    int i, j, pow_of_two; 

    double input, cm, in_band_noise, in_band_signal_and_noise;
    int output;

    /* Open the parameter file and get the simulation parameters.
    */
    ret |= get_parameters(argc, argv, &parms);

    /* Allocate memory and initialize things.
    */
    ret |= sim_initialize(&sim, &parms);
    ret |= dsm_initialize(&dsm, &parms);
            
    if (ret != SUCCESS)
	exit(-1);

    /* Go through enough "warm up" samples so that all state variables are initialized
     */
    for (i= 0; i < NUM_WARM_UP_SAMPLES; i++)     {
	    ret |= next_cm_val(&sim, &parms, &cm);  /* advance cm noise signal */
        ret |= next_input_val(&sim, &parms, &input);  /* advance input signal and index */
	    ret |= dsm_step(&dsm, input, cm, &output);
    }
 
    /* Generate the simulation output data
     */
    for (sim.per_no = 0; sim.per_no < sim.num_periodogram_averages; sim.per_no++)     {	
	    printf("Computing data for periodogram no.: %d\n", sim.per_no + 1);
	
	    for (sim.bin_no = 0; sim.bin_no < sim.num_periodogram_bins; sim.bin_no++) 	{
	        /* Step the delta-sigma modulator.
	         */
            ret |= next_cm_val(&sim, &parms, &cm);  /* advance cm noise signal */
	        ret |= next_input_val(&sim, &parms, &input);  /* advance input signal and index */
	        ret |= dsm_step(&dsm, input, cm, &output);

	        /* Update the periodogram arrays
	         */
	        sim.dsm_output[sim.bin_no] = output;
	    }

	    /* Process the new periodograms and add the results to the PSD estimation arrays
	     */
	    advance_psd_estimation(&sim, sim.dsm_output, YES, NO, sim.dsm_output_psd);
	    advance_psd_estimation(&sim, sim.dsm_output, YES, YES, sim.dsm_noise_psd);
    }

    /* Normalize the psd estimates
     */
    normalize_psd_estimate(&sim, sim.dsm_output_psd);
    normalize_psd_estimate(&sim, sim.dsm_noise_psd);

    /* Calculate and print SINAD
     */
    calculate_in_band_power(sim.dsm_noise_psd, sim.num_periodogram_bins, parms.oversampling_ratio, &in_band_noise);
    calculate_in_band_power(sim.dsm_output_psd, sim.num_periodogram_bins, parms.oversampling_ratio, &in_band_signal_and_noise);
    printf("Total in-band noise variance = %e v^2\n", in_band_noise);
    printf("SINAD = %f dB\n", 10*log10((in_band_signal_and_noise - in_band_noise) / in_band_noise));

    /* print the overload information */
    printf( "\nSwitching block overload information follows:\n" );
    for (i = 0, pow_of_two = 1; i < dsm.dac_encoder.num_layers; i++, pow_of_two *= 2) {
        printf( "\nLayer number %d\n", i );
        for ( j = 0; j < pow_of_two; j++ ) {
            printf( "s_[%d][%d] ==> %d overloads\n", i, j, dsm.dac_encoder.sb[i][j].integrator2_overload_count );
        }
    }

    /* Only save data for frequencies 0 - NyquistRate/2
    */
    for (i = 0; i < sim.num_periodogram_bins / 2; i++) {
        fprintf(sim.output_fp, " %e", (double)i / sim.num_periodogram_bins);
        fprintf(sim.output_fp, " %e", sim.dsm_output_psd[i]);
        fprintf(sim.output_fp, " %e", sim.dsm_noise_psd[i]);
        fprintf(sim.output_fp, "\n");
    }
     
    /* Copy the input file to the .log file.
    */
    strcpy(system_command_buffer, CP_CMD_STRING);   /* 'copy' on PC, 'cp' for *nix */
    strcat(system_command_buffer, argv[1]);
    strcat(system_command_buffer, " ");
    strcat(system_command_buffer, argv[1]);
    strcat(system_command_buffer, ".log");
    strcat(system_command_buffer, "\n");
    system(system_command_buffer);
    printf("done...\n");
}


int sim_initialize(SIMULATION_STATE *sim, PARMS *parms) {
    char output_file_name[256];

    srand48((long) parms->random_number_seed);
    sim->n = 0;
    sim->num_periodogram_bins = parms->num_periodogram_bins;
    sim->num_periodogram_averages = parms->num_periodogram_averages;
    sim->input_freq = parms->input_freq;

    /* Open the output file 
     */
    strcpy(output_file_name, parms->argv[1]);
    strcat(output_file_name, ".dat");
    sim->output_fp = fopen(output_file_name, "w");
    if (sim->output_fp == NULL) {
	    printf("Could not open output file.\n");
	    return(FAIL);
    }
    
    /* Allocate buffer for delta-sigma modulator output, and noise
     */
    sim->dsm_output = (double*)calloc((size_t)sim->num_periodogram_bins, sizeof(double));
    sim->dsm_noise = (double*)calloc((size_t)sim->num_periodogram_bins, sizeof(double));
    
    /* Allocate buffer for PSD of delta-sigma modulator output, and noise
     */
    sim->dsm_output_psd = (double*)calloc((size_t)sim->num_periodogram_bins, sizeof(double));
    sim->dsm_noise_psd = (double*)calloc((size_t)sim->num_periodogram_bins, sizeof(double));

    /* Allocate buffer for zero imaginary part of sample data 
     */
    sim->real = (double*)calloc((size_t)sim->num_periodogram_bins, sizeof(double));
    sim->imag = (double*)calloc((size_t)sim->num_periodogram_bins, sizeof(double));

    /* Allocate for hanning window and psd of hanning window to use when subtracting DC
     */
    sim->psd_window = (double*)calloc((size_t)sim->num_periodogram_bins, sizeof(double));
    sim->dc_term_fft_real = (double*)calloc((size_t)sim->num_periodogram_bins, sizeof(double));
    sim->dc_term_fft_imag = (double*)calloc((size_t)sim->num_periodogram_bins, sizeof(double));

    initialize_psd_estimation(sim, parms);

    return(SUCCESS);
}



int next_input_val(SIMULATION_STATE *sim, PARMS *parms, double *input) {
    /* Generate the sinusoidal plus offset test signal.
     */
    *input = parms->input_amplitude*sin(2*PI*(sim->n)*(parms->input_freq)) + parms->input_offset;
    (sim->n) ++;     
    return(SUCCESS);
}

int next_cm_val(SIMULATION_STATE *sim, PARMS *parms, double *cm) {
    /* Generate the sinusoidal plus offset cm noise
     */
    *cm = parms->cm_amplitude*sin(2*PI*(sim->n)*(parms->cm_freq)) + parms->cm_offset;
    return(SUCCESS);
}

int dsm_initialize(DSM_STATE *dsm, PARMS *parms) {
    int ret = SUCCESS, num_dac_encoder_layers, num_one_bit_dacs, i;

    dsm->num_stages = parms->num_stages;
    dsm->num_requantizers = parms->num_requantizers;
    num_dac_encoder_layers = (int)( log((double)parms->num_adc_comparators) / log(2.0) );
    num_one_bit_dacs = parms->num_adc_comparators / (int) pow(2.0, (double) parms->num_requantizers);

    /* Allocate memory for the delta-sigma stage state information
     */
    dsm->stages = (DSM_STAGE*)calloc((size_t)parms->num_stages, sizeof(DSM_STAGE));

    /* Initialize the delta-sigma stages
     */
    for (i = 0; i < dsm->num_stages; i++) {
        dsm->stages[i].integrator_output = 0;
        dsm->stages[i].integrator_gain = parms->stage_info[i].integrator_gain + 
            parms->stage_info[i].integrator_gain_std_dev*gaussian_rv();

	    ret |= dac_bank_initialize(&(dsm->stages[i].dac), num_one_bit_dacs, 
            parms->stage_info[i].dac_step_size, parms->stage_info[i].dac_step_size_std_dev);
    }

    /* Initialize the ADC 
    */
    ret |= dfadc_initialize(&(dsm->adc), -parms->max_adc_input, parms->max_adc_input, 
        parms->num_adc_comparators, parms->adc_unit_res_std_dev, parms->adc_comp_offset_std_dev);
     
    /* Initialize the DAC digital encoder
    */
    ret |= dac_encoder_initialize(&(dsm->dac_encoder), num_dac_encoder_layers - parms->num_requantizers, 
        parms->use_second_order_mmdac, parms->num_skr_int2_bits );

    /* Allocate memory for and initialize the requantizers (if applicable)
     */
    if (parms->num_requantizers > 0) {
        dsm->requantizers = (REQUANTIZER*)calloc(parms->num_requantizers, sizeof(REQUANTIZER));
        for (i = 0; i < parms->num_requantizers; i++)
            dsm->requantizers[i].s_kr = dsm->requantizers[i].tri_state_accum = 0;
    }
    else dsm->requantizers = NULL;

    return(SUCCESS);      
}


int dsm_step(DSM_STATE *dsm, double input, double cm, int *output) {
    int i;
    int adc_output, requantized_output;
    double dac_bank_output;

    /* Step the flash ADC, the requantizers (if applicable), and the digital encoder for the DAC banks
     */
    dfadc_step(&(dsm->adc), dsm->stages[dsm->num_stages-1].integrator_output, cm, &adc_output);

    requantized_output = adc_output;
    for (i = 0; i < dsm->num_requantizers; i++)
        requantizer_step(&(dsm->requantizers[i]), requantized_output, &requantized_output);
    dac_encoder_step(&(dsm->dac_encoder), requantized_output);

    /* Step all but the first stage of the delta-sigma modulator
     */
    for (i = dsm->num_stages - 1; i > 0; i--) {
        /* Step the DAC bank and the delaying integrator for this stage
        */
        dac_bank_step(&(dsm->stages[i].dac), &(dsm->dac_encoder), &dac_bank_output);
        dsm->stages[i].integrator_output += 
            dsm->stages[i].integrator_gain * (dsm->stages[i-1].integrator_output - dac_bank_output);
    }

    /* Step the DAC bank and the delaying integrator for the first stage
     */
    dac_bank_step(&(dsm->stages[0].dac), &(dsm->dac_encoder), &dac_bank_output);
    dsm->stages[0].integrator_output += 
        dsm->stages[0].integrator_gain * (input  - dac_bank_output);

    *output = requantized_output - dsm->stages[0].dac.num_one_bit_dacs/2;  /* subtract dc offset in dfadc output */
    return(SUCCESS);      
}

int requantizer_step(REQUANTIZER *requantizer, int input, int *output) {
    if (floor(input/2.0) == input/2.0) /* if input is even */
        requantizer->s_kr = 0;
    else {
        if (requantizer->tri_state_accum == 0) {
	        requantizer->s_kr = SGN(drand48() - .5);
	        requantizer->tri_state_accum = -requantizer->s_kr;
	    }
	    else {
	        requantizer->s_kr = requantizer->tri_state_accum;
	        requantizer->tri_state_accum = 0;
	    }
    }
    *output = (input + requantizer->s_kr)/2;
    return(SUCCESS);      
}
	
int dac_encoder_initialize(DAC_ENCODER_STATE *encoder, int num_layers, int use_second_order_mmdac, int num_skr_int2_bits ) {
    int i, j;

    encoder->num_layers = num_layers;
    encoder->num_one_bit_dacs = pow(2.0, (double) num_layers);
    encoder->num_switching_blocks = encoder->num_one_bit_dacs - 1;
    encoder->use_second_order_mmdac = use_second_order_mmdac;
    encoder->skr_int2_max = pow( 2, (double)( num_skr_int2_bits - 1 ) ) - 1;  
    encoder->skr_int2_min = -pow(2, (double)(num_skr_int2_bits - 1));  

    /* Allocate memory for and inialize the switching block control logic.
    * The structure is sb[i][j], where i is the layer, and j is the element in the i-th layer\
    * For a 32-element DAC
    * the encoder has log2(32)=5 layers, 
    * the layers from 0 to 4 have 2^0, 2^1, 2^2, 2^4 elements (31 total)
    */
    encoder->sb = (DAC_ENCODER_SB**)calloc((size_t)num_layers, sizeof(DAC_ENCODER_SB *));
    for (i = 0; i < num_layers; i++) {
	    /* at i-th layer, j = 0 to 2^i */
        encoder->sb[i] = (DAC_ENCODER_SB*)calloc((size_t)pow(2.0, (double) i), sizeof(DAC_ENCODER_SB));

	    for (j = 0; j < pow(2.0, (double)i); j++) {
	        encoder->sb[i][j].s_kr = -1;
	        encoder->sb[i][j].integrator1_accum = 0;
	        encoder->sb[i][j].integrator2_accum = 0;
	        encoder->sb[i][j].integrator2_overload_count = 0;
	    
	        if (i == 0) continue;  /* layer 0 is the encoder input, so no hookup not needed */
	
            /* connect inputs for layers i > 0 */
	        if (floor(j/2.0) == j/2.0) {/* if j/2 == an even number */
		        encoder->sb[i][j].input = &(encoder->sb[i-1][j/2].top_output);
	        } 
	        else {
		        encoder->sb[i][j].input= &(encoder->sb[i-1][j/2].bottom_output);
	        }
	    }
    }
    return(SUCCESS);      
}


int dac_encoder_step(DAC_ENCODER_STATE *encoder, int input) {
    int dummy, i, j, k, pow_of_two, r_i;
    double sb_input;
    int odd_input;
    double s_kr;

    /* The input is applied to the first sb 
     */
    encoder->sb[0][0].input = &input;

    /* Step the switching blocks.
     */
    for (i = k = 0, pow_of_two = 1; i < encoder->num_layers; i++, pow_of_two *= 2)
    {
        r_i = SGN(drand48() - .5);              /* This layer's independent random bit */

	    /* step switching blocks in the layer */
        for (j = 0; j < pow_of_two; j++, k++) 	{
	        sb_input = *(encoder->sb[i][j].input);
	        odd_input = ( floor( ( sb_input )/2.0 ) != ( sb_input )/2.0 );

	        /* update second order s_kr state
            */
	        if ( encoder->use_second_order_mmdac ) { 
                                
                /* compute s_kr */
                if ( encoder->sb[i][j].integrator1_accum == 0 ) {
                    encoder->sb[i][j].s_kr = ( encoder->sb[i][j].integrator2_accum == 0 ) ? -r_i : SGN( encoder->sb[i][j].integrator2_accum ); 
                }
                else {
                    encoder->sb[i][j].s_kr = SGN( encoder->sb[i][j].integrator1_accum );
                }    
                
                /* update integrator 2 */
                encoder->sb[i][j].integrator2_accum += encoder->sb[i][j].integrator1_accum;
                if ( encoder->sb[i][j].integrator2_accum > encoder->skr_int2_max ) {  /* saturate +fs */     
                    encoder->sb[i][j].integrator2_accum = encoder->skr_int2_max;
                    encoder->sb[i][j].integrator2_overload_count += 1;
                }
                else if ( encoder->sb[i][j].integrator2_accum < encoder->skr_int2_min ) {  /* saturate -fs */
                    encoder->sb[i][j].integrator2_accum = encoder->skr_int2_min;
                    encoder->sb[i][j].integrator2_overload_count += 1;
                }
                
                /* update integrator 1 */
                encoder->sb[i][j].integrator1_accum -= ( odd_input ) ? encoder->sb[i][j].s_kr : 0;
                if ( abs( encoder->sb[i][j].integrator1_accum ) > 1 ) {
                    printf( "sb[%d][%d].integrator1_accum = %d\n", i, j, encoder->sb[i][j].integrator1_accum );
                }
            }
            
            /* first order s_kr state
            */	     
            else {
                if ( odd_input ) { 
                    if (encoder->sb[i][j].integrator1_accum == 0) {
                        encoder->sb[i][j].s_kr = r_i;
                        encoder->sb[i][j].integrator1_accum = -encoder->sb[i][j].s_kr;
                    }
                    else {
                        encoder->sb[i][j].s_kr = encoder->sb[i][j].integrator1_accum;
                        encoder->sb[i][j].integrator1_accum = 0;
                    }
                }
            }   
            /* update s_kr out of block based on input parity and s_kr state
            */
            s_kr = ( odd_input ) ? encoder->sb[i][j].s_kr : 0;
            encoder->sb[i][j].top_output = (sb_input + s_kr)/2;
            encoder->sb[i][j].bottom_output = (sb_input - s_kr)/2;
        }
    }
    return(SUCCESS);
}

int dac_bank_initialize(DAC_BANK_STATE *dac_bank, int num_one_bit_dacs, double step_size, double step_size_std_dev) {
    int i;

    dac_bank->num_one_bit_dacs = num_one_bit_dacs;
    dac_bank->nominal_one_bit_dac_hi = step_size/2.0;
    dac_bank->nominal_one_bit_dac_lo = -step_size/2.0;

    /* Allocate memory for and set up the one bit DACs.
    */
    dac_bank->one_bit_dacs = (DAC_ELEMENT*)calloc((size_t)dac_bank->num_one_bit_dacs, sizeof(DAC_ELEMENT));
	
    for (i = 0; i < dac_bank->num_one_bit_dacs; i++)     {
	    dac_bank->one_bit_dacs[i].dac_error_hi = step_size * step_size_std_dev * gaussian_rv() / 2.0;
	    dac_bank->one_bit_dacs[i].dac_error_lo = step_size * step_size_std_dev * gaussian_rv() / 2.0;
    }
    return(SUCCESS);      
}




int dac_bank_step(DAC_BANK_STATE *dac_bank, DAC_ENCODER_STATE *encoder, double *output) {
    int i, j;
    DAC_ELEMENT *one_bit_dac;

    /* Transfer the final layer switching block outputs to the DAC element inputs.
     */
    for (j = 0, i = encoder->num_layers - 1; j < dac_bank->num_one_bit_dacs/2; j++) {
        dac_bank->one_bit_dacs[2*j].state = encoder->sb[i][j].top_output;
        dac_bank->one_bit_dacs[2*j + 1].state = encoder->sb[i][j].bottom_output;
    }
	
    /* Operate the bank of 1-bit DACs to obtain the overall DAC output.
     */
    for (i = 0, *output = 0; i < dac_bank->num_one_bit_dacs; i++) {
	    one_bit_dac = &(dac_bank->one_bit_dacs[i]);
	    if (one_bit_dac->state == 1)
            *output += dac_bank->nominal_one_bit_dac_hi + one_bit_dac->dac_error_hi;
	    else if (one_bit_dac->state == 0)
            *output += dac_bank->nominal_one_bit_dac_lo + one_bit_dac->dac_error_lo;
	    else
            printf("Number conservation rule violated.\n");
    }
    
    return(SUCCESS);
}



int fadc_initialize(FADC_STATE *adc, double min_input, double max_input, int num_comparators, 
    double ladder_unit_res_std_dev, double offset_std_dev) {
    int i;
    double total_ladder_resistance, ladder_current, last_reference;

    adc->num_comparators = num_comparators;
    adc->comparators = (FADC_COMPARATOR*)calloc((size_t)num_comparators, sizeof(FADC_COMPARATOR));
    adc->ref_ladder_resistors = (double*)calloc((size_t)(num_comparators + 1), sizeof(double));  

    /* Choose resistor values for the reference voltage ladder.  The first (last) reference voltage is
       nominally a half step-size above (below) the minimum (maximum) input voltage.
    */
    for (i = 1, total_ladder_resistance = 0; i < num_comparators; i++) {
        adc->ref_ladder_resistors[i] = 1 + ladder_unit_res_std_dev * gaussian_rv();
        total_ladder_resistance += adc->ref_ladder_resistors[i];
    }
    adc->ref_ladder_resistors[0] = 0.5 + ladder_unit_res_std_dev * gaussian_rv() / sqrt(2.0);
    total_ladder_resistance += adc->ref_ladder_resistors[0];
    adc->ref_ladder_resistors[num_comparators] = 0.5 + ladder_unit_res_std_dev * gaussian_rv() / sqrt(2.0);
    total_ladder_resistance += adc->ref_ladder_resistors[num_comparators];

    /* Using the resistor values just generated, calculate the reference levels such that the top and
       bottom of the resistor ladder have the appropriate voltages (i.e., max_input and min_input).
       Also, generate the comparator offset values.
    */
    ladder_current = (max_input - min_input)/total_ladder_resistance;

    for (i = 0, last_reference = min_input; i < num_comparators; i++) {
        adc->comparators[i].reference = last_reference + ladder_current*adc->ref_ladder_resistors[i];
        last_reference = adc->comparators[i].reference;
        adc->comparators[i].offset = offset_std_dev * gaussian_rv();
    }

    /* Calculate the nominal quantization step size of the ADC 
     */
    adc->step_size = (max_input - min_input)/num_comparators;
    return(SUCCESS);      
}


int dfadc_initialize(DFADC_STATE *adc, double min_input, double max_input, int num_comparators, 
    double ladder_unit_res_std_dev, double offset_std_dev) {

    int i;
    double total_ladder_resistance, ladder_current, last_reference;
    
    adc->num_comparators = num_comparators;
    adc->p_comparators = (FADC_COMPARATOR*)calloc((size_t)num_comparators/2, sizeof(FADC_COMPARATOR));
    adc->n_comparators = (FADC_COMPARATOR*)calloc((size_t)num_comparators/2, sizeof(FADC_COMPARATOR));
    adc->ref_ladder_resistors = (double*)calloc((size_t)(num_comparators/2 + 1), sizeof(double));  

    /* Reference ladder shared between the two sets of comparators; there are num_comps/2 taps.
       The first (last) ref voltage is nominally 1/2 (3/2) steps from the min (max) input voltage 
    */
    for (i = 1, total_ladder_resistance = 0; i < num_comparators/2; i++) {
        adc->ref_ladder_resistors[i] = 2 + ladder_unit_res_std_dev * gaussian_rv() * sqrt(2.0);
	    total_ladder_resistance += adc->ref_ladder_resistors[i];
    }
    adc->ref_ladder_resistors[0] = 0.5 + ladder_unit_res_std_dev * gaussian_rv() * sqrt(0.5);
    total_ladder_resistance += adc->ref_ladder_resistors[0];
    adc->ref_ladder_resistors[num_comparators/2] = 1.5 + ladder_unit_res_std_dev * gaussian_rv() * sqrt(1.5);
    total_ladder_resistance += adc->ref_ladder_resistors[num_comparators/2];

    /* Using the resistor values just generated, calculate the reference levels such that the top and
       bottom of the resistor ladder have the appropriate voltages (i.e., max_input and min_input).
       Also, generate the comparator offset values.
    */
    ladder_current = (max_input - min_input)/total_ladder_resistance;

    for (i = 0, last_reference = min_input; i < num_comparators/2; i++) {
        adc->p_comparators[i].reference = last_reference + ladder_current*adc->ref_ladder_resistors[i];
        adc->n_comparators[i].reference = adc->p_comparators[i].reference;
	    last_reference = adc->p_comparators[i].reference;
	    adc->p_comparators[i].offset = offset_std_dev * gaussian_rv();
	    adc->n_comparators[i].offset = offset_std_dev * gaussian_rv();
    }

    /* Calculate the nominal quantization step size of the ADC 
     */
    adc->step_size = (max_input - min_input)/num_comparators;
    
    return(SUCCESS);      
}


int fadc_step(FADC_STATE *adc, double input, int *output) {
    int i;
    /* Note: the number of possible output values is one plus the number
       of comparators.
    */
    for (i = 0; i < adc->num_comparators; i++)     {
        if (input <  adc->comparators[i].reference + adc->comparators[i].offset)
	    break;
    }
    *output = i;
    return(SUCCESS);      
}


int dfadc_step(DFADC_STATE *adc, double input, double cm, int *output) {
    int ip, in;

    /* positive half of differential signal */
    for (ip = 0; ip < (adc->num_comparators)/2; ip++) {
        if (input+cm <  adc->p_comparators[ip].reference + adc->p_comparators[ip].offset)
	    break;
    }

    /* negative half of differential signal */
    for (in = 0; in < (adc->num_comparators)/2; in++) {
        if (-input+cm <  adc->n_comparators[in].reference + adc->n_comparators[in].offset)
	    break;
    }  
    /* total number of output levels = num_comparators+1 */
    /* add dc offset to make compatible w/ dac encoder input */
    *output = ip-in+(adc->num_comparators)/2;  
    return(SUCCESS);      
}

/* gaussian_rv() is a version of the numerical recipies function modified to use
   drand48() instead of ran1()
*/
double gaussian_rv() {
	static int iset=0;
	static double gset;
	double fac,r,v1,v2;
	double log(), sqrt();

	if  (iset == 0) {
	    do {
		    v1=2.0*drand48()-1.0;
		    v2=2.0*drand48()-1.0;
		    r=v1*v1+v2*v2;
	    } while (r >= 1.0 || r == 0.0);

	    fac=sqrt(-2.0*log(r)/r);
	    gset=v1*fac;
	    iset=1;
	    return(v2*fac);
	} 
	else {
	    iset=0;
	    return(gset);
	}
}

int get_parameters(int argc, char *argv[], PARMS *parms) {
    int i, ret = SUCCESS, status = SUCCESS;
    char token[MAX_TOKEN_SIZE];
    FILE *fp;
    
    /* Open the input file
     */
    fp = fopen(argv[1], "r");

    if (fp == NULL || argc < 2) {
	    printf("Could not open input file.\n");
	    return(FAIL);
    }

    /* Store access to the command line information
     */
    parms->argv = argv;
    parms->argc = argc;

    /* Parse the input file
     */
    status |= get_token(fp, token);
    parms->num_stages = atof(token);

    parms->stage_info = (struct parms_stage_info*)calloc((size_t) parms->num_stages, sizeof(struct parms_stage_info));
    for (i = 0; i < parms->num_stages; i++) {
        status |= get_token(fp, token);
	    parms->stage_info[i].integrator_gain = atof(token);

	    status |= get_token(fp, token);
	    parms->stage_info[i].integrator_gain_std_dev = atof(token);

	    status |= get_token(fp, token);
	    parms->stage_info[i].dac_step_size = atof(token);

	    status |= get_token(fp, token);
	    parms->stage_info[i].dac_step_size_std_dev = atof(token);
    }
    status |= get_token(fp, token);
    parms->use_second_order_mmdac = atof(token);
	    	  
    status |= get_token(fp, token);
    parms->num_skr_int2_bits = atof(token);
	    	    
    status |= get_token(fp, token);
    parms->num_requantizers = atof(token);
	    
    status |= get_token(fp, token);
    parms->num_adc_comparators = atof(token);
	    
    status |= get_token(fp, token);
    parms->max_adc_input = atof(token);

    status |= get_token(fp, token);
    parms->adc_unit_res_std_dev = atof(token);

    status |= get_token(fp, token);
    parms->adc_comp_offset_std_dev = atof(token);

    status |= get_token(fp, token);
    parms->input_amplitude = atof(token);

    status |= get_token(fp, token);
    parms->input_freq = atof(token);

    status |= get_token(fp, token);
    parms->input_offset = atof(token);

    status |= get_token(fp, token);
    parms->cm_amplitude = atof(token);

    status |= get_token(fp, token);
    parms->cm_freq = atof(token);

    status |= get_token(fp, token);
    parms->cm_offset = atof(token);

    status |= get_token(fp, token);
    parms->num_periodogram_bins = atof(token);

    status |= get_token(fp, token);
    parms->num_periodogram_averages = atof(token);

    status |= get_token(fp, token);
    parms->oversampling_ratio = atof(token);

    status |= get_token(fp, token);
    parms->random_number_seed = atof(token);

    fclose(fp);

    if (status != SUCCESS) {
	    printf("Incomplete parameter file.\n");
	    ret = FAIL;
    }
    return(ret);
}    

int get_token(FILE *fp, char *token) {
    int c;

    while (((c = getc(fp)) == ' ' || c == '\t' || c == '\n') && c != EOF);

    /* Allow for comments following '%' and continuing until '\n' or EOF.
    */
    while (c == '%') {
        while((c = getc(fp)) != '\n' && c!= EOF);
        if (c == EOF)
	    break;
        while (((c = getc(fp)) == ' ' || c == '\t' || c == '\n') && c != EOF);
    }
        
    if (c != EOF) {
        *token = c;
        token++;
        while ((c = getc(fp)) != ' ' && c != '\t' && c != '\n' && c != ',' && c != EOF) {
            *token = c;
            token++;
        }
    }
        
    *token = '\0';

    if (c != EOF)
	c = SUCCESS;
    
    return(c);
}


int calculate_in_band_power(double *psd, int num_points, int osr, double *in_band_power) {
    int k;
    for (*in_band_power = 0, k = 0; k < (num_points/2)/osr; k++)
        (*in_band_power) += psd[k]/(num_points/2);
    return(SUCCESS);
}


void hanning_window(double *window, int window_length) {
    int j;
    
    for (j = 1; j <= window_length; j++)
	window[j - 1] = 0.5 * (1 - cos((2.0 * PI * j) / (window_length + 1)));  
}

int initialize_psd_estimation(SIMULATION_STATE *sim, PARMS *parms) {
    int k;

    hanning_window(sim->psd_window, parms->num_periodogram_bins);
    for (k = 0, sim->psd_U = (double) 0.0; k < parms->num_periodogram_bins; k++)
	    sim->psd_U += (sim->psd_window[k] * sim->psd_window[k]);    
    sim->psd_U /= parms->num_periodogram_bins;

    /* Compute the fft of the window to obtain a dc term fft for dc removal
     */
    copy(sim->dc_term_fft_real, sim->psd_window, parms->num_periodogram_bins);
    zero_fill(sim->dc_term_fft_imag, parms->num_periodogram_bins);
    fft(sim->dc_term_fft_real, sim->dc_term_fft_imag, parms->num_periodogram_bins, 0);

    return(SUCCESS);
}


int remove_sinusoidal_component(double* input, int num_points, double freq, double* window, double* output) {
    int k;
    double sin_square_sum, cos_square_sum, sin_cos_sum;
    double cos_signal_sum, sin_signal_sum;
    double denominator, cos_scale, sin_scale;

    /* Calculate the constants required to compute the cosine and sine scale factors.
     */
    sin_square_sum = cos_square_sum = sin_cos_sum = cos_signal_sum = sin_signal_sum = 0;
    for (k = 0; k < num_points; k++) {
        sin_square_sum += SQR(window[k] * sin(2 * PI * k * freq));
        cos_square_sum += SQR(window[k] * cos(2 * PI * k * freq));
        sin_cos_sum += SQR(window[k]) * sin(2 * PI * k * freq) * cos(2 * PI * k * freq);

        cos_signal_sum += SQR(window[k]) * input[k] * cos(2 * PI * k * freq);
        sin_signal_sum += SQR(window[k]) * input[k] * sin(2 * PI * k * freq);
    }

    /* Calculate the scale factors.
     */
    denominator = sin_square_sum * cos_square_sum - SQR(sin_cos_sum);
    cos_scale = (sin_square_sum * cos_signal_sum - sin_cos_sum * sin_signal_sum) / denominator;
    sin_scale = (cos_square_sum * sin_signal_sum - sin_cos_sum * cos_signal_sum) / denominator;

    /* Subtract the sinusoidal term from the sequence and copy the result to the output array.
     */
    for (k = 0; k < num_points; k++)
        output[k] = input[k] - cos_scale * cos(2 * PI * k * freq) - sin_scale * sin(2 * PI * k * freq);

    return(SUCCESS);
}


int advance_psd_estimation(SIMULATION_STATE *sim, double *time_sequence, int subtract_DC, 
    int subtract_sinusoid, double *psd) {

    int j;
    double dc_scale;

    /* If supposed to remove the sinusoidal signal, do it.  In either case, copy the time sequence
       into the real part of the fft array.  Zero-fill the imaginary part.
     */
    if (subtract_sinusoid) {
        remove_sinusoidal_component(time_sequence,
				    sim->num_periodogram_bins,
				    sim->input_freq,
				    sim->psd_window,
				    sim->real);
    }
    else {
        for (j = 0; j < sim->num_periodogram_bins; j++)
	    sim->real[j] = time_sequence[j];
    }      

    zero_fill(sim->imag, sim->num_periodogram_bins);
    
    /* window the input sequence
     */
    for (j = 0; j < sim->num_periodogram_bins; j++)
        sim->real[j] *= sim->psd_window[j];
	
    fft(sim->real, sim->imag, sim->num_periodogram_bins, 0);    
       
    /* If supposed to remove dc offset, subtract scaled version of the fft of the dc term
     */
    if (subtract_DC) {	
        dc_scale = sim->real[0]/sim->dc_term_fft_real[0];
        for (j = 0; j <sim->num_periodogram_bins; j++) {
            sim->real[j] -= dc_scale * sim->dc_term_fft_real[j];
            sim->imag[j] -= dc_scale * sim->dc_term_fft_imag[j];
        } 
    }   

    /* finally, complete the computation of the periodogram
     */
    for (j = 0; j < sim->num_periodogram_bins; j++)
        psd[j] += (sim->real[j] * sim->real[j]) + (sim->imag[j] * sim->imag[j]);
    
    return(SUCCESS);    
}

int normalize_psd_estimate(SIMULATION_STATE *sim, double *psd) {
    int j;
    
    for (j = 0; j < sim->num_periodogram_bins; j++)
        psd[j] /= (sim->psd_U * sim->num_periodogram_bins * sim->num_periodogram_averages);
    return(SUCCESS);    
}



/*
Introduced into UseNet By: aeusemrs@csun.uucp (Mike Stump) Sun Jul 26 11:26:31 PDT 1987
HEADER:
TITLE:		Fast Fourier Transform;
DATE:		05/18/1985;
DESCRIPTION:	"Performs fast fourier transform using method described
		by E. O. Brigham.  For details of the method, refer
		to Brigham's book. THE FAST FOURIER TRANSFORM";
KEYWORDS: 	Fourier, transform;
FILENAME:	FFT.C;
WARNINGS:
  "This program is self-contained, all that is needed is a manner of getting
  the data into the array real_data (& imag_data, if applicable).  The
  transformed data will reside in these two arrays upon return with the
  original data being destroyed."
AUTHORS:	Jim Pisano;
COMPILERS:	DeSmet;
REFERENCES:	AUTHORS:	E. O. Brigham;
		TITLE:		"THE FAST FOURIER TRANSFORM";
		CITATION:	"";
	ENDREF
*/

/*	file name fft.c
*	program name fft() ... Fast Fourier Transform
*
*	Perform fast fourier transform using method described by E. O. Brigham.
*  For details of the method, refer to Brigham's book
*
*	Translated to C from FORTRAN by
*
*		Jim Pisano
*		P.O. Box 3134
*		University Station
*		Charlottesville, VA 22903
*
*  This program is in the public domain & may be used by anyone for commercial
*  or non-commercial purposes.
*
*  real_data ... ptr. to real part of data to be transformed
*  imag_data ... ptr. to imag  "   "   "   "  "      "
*  inv ..... Switch to flag normal or inverse transform
*  n_pts ... Number of real data points
*  nu ...... logarithm in base 2 of n_pts e.g. nu = 5 if n_pts = 32.
*
*  This program is self-contained, all that is needed is a manner of getting
*  the data into the array real_data (& imag_data, if applicable).  The
*  transformed data will reside in these two arrays upon return with the
*  original data being destroyed.
*/

int bit_swap(int i, int nu)
{
	int ib, i1, i2;

	ib = 0;

	for( i1 = 0; i1 < nu; i1++ )
	{
		i2  = i / 2;
		ib = ib * 2 + i - 2 * i2;
		i   = i2;
	}
	return( ib );
}
/*
* Simple exchange routine where *x1 & *x2 are swapped
*/
void dswap(double *x1,  double *x2)
{
	double temp_x;

	temp_x = *x1;
	*x1 = *x2;
	*x2 = temp_x;
}

void copy(double *to, double *from, int length)
{
    while (length)
    {
	length--;
        to[length] = from[length];
    }
}

void zero_fill(double *seq, int seq_length)
{
    int j;
    
    for (j = 0; j < seq_length; j++)
	seq[j] = (double) 0.0;
}

void zero_pad(double *seq, int start_index, int seq_length)
{
    int j;
    
    for (j = start_index; j < seq_length; j++)
	seq[j] = 0.0;
}

void time_domain_alias(double *output_seq, 
		      int output_length, 
		      double *input_seq, 
		      int input_length)
{
    int i;

    for (i = 0; i < output_length; i++)
	output_seq[i] = 0.0;
    
    for (i = 0; i < input_length; i++)   
 	output_seq[(i % output_length)] += input_seq[i];
}

void fft(double *real_data, double *imag_data, int n_pts, int inv)
{
	int n2, nu, n_pts_copy, j, l, i, ib, k, k1, k2;
	int sgn;
	double tr, ti, arg, nu1;	/* intermediate values in calcs. */
	double c, s;	       /* cosine & sine components of Fourier trans. */

	for (nu = 0, n_pts_copy = 1; n_pts_copy < n_pts; nu++, n_pts_copy *= 2) ;

	n2 = n_pts / 2;
	nu1 = nu - 1.0;
	k = 0;
/*
* sign change for inverse transform
*/
	sgn = inv ? -1 : 1;	

/*
* Calculate the componets of the Fourier series of the function
*/
	for( l = 0; l != nu; l++ )
	{
		do
		{
			for( i = 0; i != n2; i++ )
			{
				j = k / ((int)pow(2.0, (double)nu1));
				ib = bit_swap( j, nu );
				arg = 2.0 * PI * ib / n_pts;
				c = cos( arg );
				s = sgn * sin( arg );
				k1 = k;
				k2 = k1 + n2;
				tr = *(real_data+k2) * c + *(imag_data+k2) * s;
				ti = *(imag_data+k2) * c - *(real_data+k2) * s;
				*(real_data+k2) = *(real_data+k1) - tr;
				*(imag_data+k2) = *(imag_data+k1) - ti;
				*(real_data+k1) = *(real_data+k1) + tr;
				*(imag_data+k1) = *(imag_data+k1) + ti;
				k++;
			}
			k +=  n2;
		} while( k < n_pts - 1);
		k = 0;
		nu1 -= 1.0;
		n2 /= 2;
	}

	for( k = 0; k < n_pts; k++ )
	{
		
		ib = bit_swap(k, nu);
	
		if( ib > k)
		{
			dswap( (real_data+k), (real_data+ib) );
			dswap( (imag_data+k), (imag_data+ib) );
		}
	}
/*
* If calculating the inverse transform, must divide the data by the number of
* data points.
*/
	if( inv )
		for( k = 0; k != n_pts; k++)
		{
			*(real_data+k) /= n_pts;
			*(imag_data+k) /= n_pts;
		}
}
