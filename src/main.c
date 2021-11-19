#include "main.h"
#include <stdarg.h>
#include <string.h>
#include <stdio.h>

#include "arm_math.h"
#include "arm_const_structs.h"
#include "arm_common_tables.h"

#define DEBUG_SAMPLING_FREQUENCY 1
#define FFT_USE_WINDOWING 1
#define FFT_USE_FILTER_HIGH 1
#define FFT_USE_FILTER_LOW 1

#define SAMPLING_FREQUENCY (4000) //44000) /* Hz */

#define SAMPLE_SIZE 1024UL
#define FULLBUFFERSIZE (SAMPLE_SIZE*2)

#define AUDIO_APP_LOG_ENABLED               1
#define FFT_LOG_ENABLED                     1
#define GRAPH_WINDOW_HEIGHT                 40                              //!< Graph window height used in draw function.
#define BLOCKS_TO_TRANSFER                  1 //2
#define FFT_COMPLEX_INPUT                   (SAMPLE_SIZE * 2UL)
#define FFT_OUTPUT_SIZE                     SAMPLE_SIZE //(I2S_DATA_BLOCK_WORDS / 2UL)

#define FFT_TEST_SAMPLE_FREQ_HZ             SAMPLING_FREQUENCY //(MCLK_FREQ_HZ/128.0f) 
#define FFT_TEST_SAMPLE_FREQ_FACTOR         (3.7878) /* for 4kHz sampling frequency */   // todo create calibration mecanisme get n point and compute factor           
#define FFT_TEST_SAMPLE_RES_HZ              ((FFT_TEST_SAMPLE_FREQ_HZ / FFT_COMPLEX_INPUT) * FFT_TEST_SAMPLE_FREQ_FACTOR)

#define FFT_TEST_HIGH_CUTT_OFF_FREQ         50 /* Hz */
#define FFT_TEST_LOW_CUTT_OFF_FREQ          50 /* Hz */

ADC_HandleTypeDef hadc1;
DMA_HandleTypeDef hdma_adc1;
DMA_HandleTypeDef hdma_dac1;
TIM_HandleTypeDef htim2;
UART_HandleTypeDef huart1;

uint32_t adcData[FULLBUFFERSIZE];

static volatile uint32_t* inbufPtr = NULL;
static volatile uint32_t* outbufPtr;

char sz_traceBuffer[10000];
size_t sz_traceBufferSize;

const static arm_cfft_instance_f32 *S;

uint16_t blocksOffset = BLOCKS_TO_TRANSFER;
uint16_t volatile blocksTransferred;
bool loop = FALSE;

static volatile bool fftIsBusy = FALSE;

static float m_fft_output_f32[FFT_OUTPUT_SIZE];             //!< FFT output data. Frequency domain.
static float m_fft_input_f32[FFT_COMPLEX_INPUT] = {0};     //!< FFT input array for complex numbers. Time domain.
static uint32_t     fftSamples[SAMPLE_SIZE] = {0};

BEEP_protocol_s     fft_result;
BEEP_protocol_s     settings;


/* Private function prototypes -----------------------------------------------*/
void SystemClock_Config(void);
static void MX_GPIO_Init(void);
static void MX_DMA_Init(void);
static void MX_ADC1_Init(void);
static void MX_TIM2_Init(void);

/* Private user code ---------------------------------------------------------*/

/**
 * @brief Macro to be used in a formatted string to a pass float number to the log.
 *
 * Use this macro in a formatted string instead of the %f specifier together with
 * @ref NRF_LOG_FLOAT macro.
 * Example: NRF_LOG_INFO("My float number" TRACE_FLOAT_MARKER "\r\n", TRACE_FLOAT(f)))
 */
//#define TRACE_FLOAT_MARKER "%s%d.%02d"

#define TRACE_FLOAT_MARKER "%s%d.%d"

/**
 * @brief Macro for dissecting a float number into two numbers (integer and residuum).
 */

#define EXP2(a, b) a ## b
#define EXP(a, b) EXP2(a ## e,b)
#define POW(C,x) EXP(C, x)
#define TRACE_FLOAT_PRECISION 4

#define TRACE_FLOAT(val) (uint32_t)(((val) < 0 && (val) > -1.0) ? "-" : ""),   \
                           (int32_t)(val),                                       \
                           (int32_t)((((val) > 0) ? (val) - (int32_t)(val)       \
                                                : (int32_t)(val) - (val))* 10000)

#define TRACE(fmt, ...) trace(fmt, ##__VA_ARGS__)

void trace(const char * format, ... ) {
  va_list argptr;

  va_start(argptr, format);
  sz_traceBufferSize = vsnprintf(sz_traceBuffer, sizeof(sz_traceBuffer), format, argptr);
  va_end(argptr);

//  if(HAL_UART_Transmit_IT(&huart1, (uint8_t*)sz_traceBuffer, sz_traceBufferSize)!= HAL_OK) {
//    Error_Handler();
//  }

  if(HAL_UART_Transmit(&huart1, (uint8_t*)sz_traceBuffer, sz_traceBufferSize, 1000)!= HAL_OK) {
    Error_Handler();
  }
}

#if (FFT_USE_WINDOWING == 1)
#if 0
/* Windowing type */
#define FFT_WIN_TYP_RECTANGLE 0x00 /* rectangle (Box car) */
#define FFT_WIN_TYP_HAMMING 0x01 /* hamming */
#define FFT_WIN_TYP_HANN 0x02 /* hann */
#define FFT_WIN_TYP_TRIANGLE 0x03 /* triangle (Bartlett) */
#define FFT_WIN_TYP_NUTTALL 0x04 /* nuttall */
#define FFT_WIN_TYP_BLACKMAN 0x05 /* blackman */
#define FFT_WIN_TYP_BLACKMAN_NUTTALL 0x06 /* blackman nuttall */
#define FFT_WIN_TYP_BLACKMAN_HARRIS 0x07 /* blackman harris*/
#define FFT_WIN_TYP_FLT_TOP 0x08 /* flat top */
#define FFT_WIN_TYP_WELCH 0x09 /* welch */

#define twoPi 6.28318531
#define fourPi 12.56637061
#define sixPi 18.84955593
#define sq(x) ((x)*(x))

#define FFT_FORWARD 0x01
#define FFT_REVERSE 0x00

void Windowing(double *vData, uint16_t samples, uint8_t windowType, uint8_t dir)
{// Weighing factors are computed once before multiple use of FFT
// The weighing function is symetric; half the weighs are recorded

	double samplesMinusOne = (double) samples - 1.0;
	for (uint16_t i = 0; i < (samples >> 1); i++) {
		double indexMinusOne = (double) i;
		double ratio = (indexMinusOne / samplesMinusOne);
		double weighingFactor = 1.0;
		// Compute and record weighting factor
		switch (windowType) {
		case FFT_WIN_TYP_RECTANGLE: // rectangle (box car)
			weighingFactor = 1.0;
			break;
		case FFT_WIN_TYP_HAMMING: // hamming
			weighingFactor = 0.54 - (0.46 * cos(twoPi * ratio));
			break;
		case FFT_WIN_TYP_HANN: // hann
			weighingFactor = 0.54 * (1.0 - cos(twoPi * ratio));
			break;
		case FFT_WIN_TYP_TRIANGLE: // triangle (Bartlett)
			weighingFactor = 1.0 - ((2.0 * abs(indexMinusOne - (samplesMinusOne / 2.0))) / samplesMinusOne);
			break;
		case FFT_WIN_TYP_NUTTALL: // nuttall
			weighingFactor = 0.355768 - (0.487396 * (cos(twoPi * ratio))) + (0.144232 * (cos(fourPi * ratio))) - (0.012604 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN: // blackman
			weighingFactor = 0.42323 - (0.49755 * (cos(twoPi * ratio))) + (0.07922 * (cos(fourPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN_NUTTALL: // blackman nuttall
			weighingFactor = 0.3635819 - (0.4891775 * (cos(twoPi * ratio))) + (0.1365995 * (cos(fourPi * ratio))) - (0.0106411 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_BLACKMAN_HARRIS: // blackman harris
			weighingFactor = 0.35875 - (0.48829 * (cos(twoPi * ratio))) + (0.14128 * (cos(fourPi * ratio))) - (0.01168 * (cos(sixPi * ratio)));
			break;
		case FFT_WIN_TYP_FLT_TOP: // flat top
			weighingFactor = 0.2810639 - (0.5208972 * cos(twoPi * ratio)) + (0.1980399 * cos(fourPi * ratio));
			break;
		case FFT_WIN_TYP_WELCH: // welch
			weighingFactor = 1.0 - sq((indexMinusOne - samplesMinusOne / 2.0) / (samplesMinusOne / 2.0));
			break;
		}
		if (dir == FFT_FORWARD) {
			vData[i] *= weighingFactor;
			vData[samples - (i + 1)] *= weighingFactor;
		}
		else {
			vData[i] /= weighingFactor;
			vData[samples - (i + 1)] /= weighingFactor;
		}
	}
}
#else
//https://github.com/sidneycadot/WindowFunctions
//https://github.com/kichiki/WaoN/blob/master/fft.c
//https://github.com/kfrlib/kfr

/* Reference: "Numerical Recipes in C" 2nd Ed.
 * by W.H.Press, S.A.Teukolsky, W.T.Vetterling, B.P.Flannery
 * (1992) Cambridge University Press.
 * ISBN 0-521-43108-5
 * Sec.13.4 - Data Windowing
 */
double
parzen (int i, int nn)
{
  return (1.0 - fabs (((double)i-0.5*(double)(nn-1))
		      /(0.5*(double)(nn+1))));
}

double
welch (int i, int nn)
{
  return (1.0-(((double)i-0.5*(double)(nn-1))
	       /(0.5*(double)(nn+1)))
	  *(((double)i-0.5*(double)(nn-1))
	    /(0.5*(double)(nn+1))));
}

double
hanning (int i, int nn)
{
  return ( 0.5 * (1.0 - cos (2.0*M_PI*(double)i/(double)(nn-1))) );
}

/* Reference: "Digital Filters and Signal Processing" 2nd Ed.
 * by L. B. Jackson. (1989) Kluwer Academic Publishers.
 * ISBN 0-89838-276-9
 * Sec.7.3 - Windows in Spectrum Analysis
 */
double
hamming (int i, int nn)
{
  return ( 0.54 - 0.46 * cos (2.0*M_PI*(double)i/(double)(nn-1)) );
}

double
blackman (int i, int nn)
{
  return ( 0.42 - 0.5 * cos (2.0*M_PI*(double)i/(double)(nn-1))
	  + 0.08 * cos (4.0*M_PI*(double)i/(double)(nn-1)) );
}

double
steeper (int i, int nn)
{
  return ( 0.375
	  - 0.5 * cos (2.0*M_PI*(double)i/(double)(nn-1))
	  + 0.125 * cos (4.0*M_PI*(double)i/(double)(nn-1)) );
}

/* apply window function to data[]
 * INPUT
 *  flag_window : 0 : no-window (default -- that is, other than 1 ~ 6)
 *                1 : parzen window
 *                2 : welch window
 *                3 : hanning window
 *                4 : hamming window
 *                5 : blackman window
 *                6 : steeper 30-dB/octave rolloff window
 */
void
windowing (int n, const double *data, int flag_window, double scale, double *out)
{
  int i;
  for (i = 0; i < n; i ++) {
    switch (flag_window) {
	case 1: // parzen window
	  out [i] = data [i] * parzen (i, n) / scale;
	  break;

	case 2: // welch window
	  out [i] = data [i] * welch (i, n) / scale;
	  break;

	case 3: // hanning window
	  out [i] = data [i] * hanning (i, n) / scale;
	  break;

	case 4: // hamming window
	  out [i] = data [i] * hamming (i, n) / scale;
	  break;

	case 5: // blackman window
	  out [i] = data [i] * blackman (i, n) / scale;
	  break;

	case 6: // steeper 30-dB/octave rolloff window
	  out [i] = data [i] * steeper (i, n) / scale;
	  break;

	default:
	  //fprintf (stderr, "invalid flag_window\n");
	case 0: // square (no window)
	  out [i] = data [i] / scale;
	  break;
	  }
  }
}

#endif
#endif

void HAL_ADC_ConvHalfCpltCallback(ADC_HandleTypeDef* hadc)
{
  inbufPtr = &adcData[0];

  // Check if the FFT is ready to accept new data
  //if(!fftIsBusy)
  //{
  //  memcpy(fftSamples, &adcData[0], SAMPLE_SIZE);
  //  fftIsBusy = TRUE;
  //}

  HAL_GPIO_TogglePin(GPIOA, GPIO_PIN_8);
}

void HAL_ADC_ConvCpltCallback(ADC_HandleTypeDef* hadc)
{
	inbufPtr = &adcData[SAMPLE_SIZE];

  // Check if the FFT is ready to accept new data
  if(!fftIsBusy)
  {
    memcpy(fftSamples, &adcData[SAMPLE_SIZE], SAMPLE_SIZE);
    fftIsBusy = TRUE;
  }

  HAL_GPIO_TogglePin(GPIOA, GPIO_PIN_8);
}

#if (DEBUG_SAMPLING_FREQUENCY == 1)
void HAL_TIM_OC_DelayElapsedCallback(TIM_HandleTypeDef* htim)
{
  //HAL_GPIO_TogglePin(GPIOA, GPIO_PIN_8);
}
#endif

bool processDSP(void) {
  // Check whether a block has been released by the I2S peripheral.
  if(fftIsBusy) {
    blocksTransferred++;

    // Only process the I2S data of the I2S measurement result we'll actually use. When loop==true each sample is processed.
    if((blocksTransferred < blocksOffset) && (blocksOffset != 0) && !loop) {
      fftIsBusy = FALSE;
      return;
    }

#if 1
    //print_plotter(fftSamples, SAMPLE_SIZE);
#endif

#if (FFT_USE_WINDOWING == 1)   
    double f64_fftSamples[SAMPLE_SIZE];
    float32_t f32_fftSamples[SAMPLE_SIZE];
    uint32_t i;

    for (i=0;i<SAMPLE_SIZE;i++) {
      //f64_fftSamples[i] = (double) (__LL_ADC_CALC_DATA_TO_VOLTAGE(3300, fftSamples[i], LL_ADC_RESOLUTION_12B) / 1000.0) - 1.25;
      f64_fftSamples[i] = (double) (__LL_ADC_CALC_DATA_TO_VOLTAGE(3300, fftSamples[i], LL_ADC_RESOLUTION_12B) / 1000.0);
      //p_input[i] = (float32_t) __LL_ADC_CALC_DATA_TO_VOLTAGE(3300, samples[i], LL_ADC_RESOLUTION_12B);
      //p_input[i] = (float32_t) samples[i];
    }

#if 0
#if 0
    #define FFT_WIN_TYP_RECTANGLE 0x00 /* rectangle (Box car) */
    #define FFT_WIN_TYP_HAMMING 0x01 /* hamming */
    #define FFT_WIN_TYP_HANN 0x02 /* hann */
    #define FFT_WIN_TYP_TRIANGLE 0x03 /* triangle (Bartlett) */
    #define FFT_WIN_TYP_NUTTALL 0x04 /* nuttall */
    #define FFT_WIN_TYP_BLACKMAN 0x05 /* blackman */
    #define FFT_WIN_TYP_BLACKMAN_NUTTALL 0x06 /* blackman nuttall */
    #define FFT_WIN_TYP_BLACKMAN_HARRIS 0x07 /* blackman harris*/
    #define FFT_WIN_TYP_FLT_TOP 0x08 /* flat top */
    #define FFT_WIN_TYP_WELCH 0x09 /* welch */
#endif
    Windowing(f64_fftSamples, SAMPLE_SIZE, FFT_WIN_TYP_HANN, FFT_FORWARD);

#else
/*  flag_window : 0 : no-window (default -- that is, other than 1 ~ 6)
 *                1 : parzen window
 *                2 : welch window
 *                3 : hanning window
 *                4 : hamming window
 *                5 : blackman window
 *                6 : steeper 30-dB/octave rolloff window
 */
    windowing (SAMPLE_SIZE, f64_fftSamples, 4 /* hamming */, 1.0, f64_fftSamples);
#endif

    for (i=0;i<SAMPLE_SIZE;i++)
      f32_fftSamples[i] = (float32_t) f64_fftSamples[i];

    /* Convert the uint32_t array containing the samples of the Audio ADC to an float array with complex numbers. The real part will be placed on the even 
     * indexes and the imaginary part will be set to 0 on all uneven indexes. This means that the complex input array is twice the size of the number of
     * samples.
     */
    fft_generate_complexNumbersArray(m_fft_input_f32, FFT_COMPLEX_INPUT, f32_fftSamples, SAMPLE_SIZE);
#else 
    /* Convert the uint32_t array containing the samples of the Audio ADC to an float array with complex numbers. The real part will be placed on the even 
     * indexes and the imaginary part will be set to 0 on all uneven indexes. This means that the complex input array is twice the size of the number of
     * samples.
     */
    fft_generate_complexNumbersArray(m_fft_input_f32, FFT_COMPLEX_INPUT, fftSamples, SAMPLE_SIZE);
#endif
    
    //print_plotter(m_fft_input_f32, FFT_COMPLEX_INPUT);

    /* Use CFFT module to process the data.
     * Ensure that the used module is equal to the number of complex number of samples.
     */
  switch (SAMPLE_SIZE) {
    case 16:
      S = &arm_cfft_sR_f32_len16;
      break;
    case 32:
      S = &arm_cfft_sR_f32_len32;
      break;
    case 64:
      S = &arm_cfft_sR_f32_len64;
      break;
    case 128:
      S = &arm_cfft_sR_f32_len128;
      break;
    case 256:
      S = &arm_cfft_sR_f32_len256;
      break;
    case 512:
      S = &arm_cfft_sR_f32_len512;
      break;
    case 1024:
      S = &arm_cfft_sR_f32_len1024;
      break;
    case 2048:
      S = &arm_cfft_sR_f32_len2048;
      break;
    case 4096:
      S = &arm_cfft_sR_f32_len4096;
      break;
  }
  
  arm_cfft_f32(S, m_fft_input_f32, 0, 1);

  double Fmax = SAMPLING_FREQUENCY; //maximum frequency
  double df = Fmax / FFT_COMPLEX_INPUT; //delta frequency

#if (FFT_USE_FILTER_HIGH == 1)
  //high-pass filter
  double FcutHigh = FFT_TEST_HIGH_CUTT_OFF_FREQ / df; //set cutoff frequency as FFT_TEST_HIGH_CUTT_OFF_FREQ Hz
  
  for(uint32_t i = 0; i< FFT_COMPLEX_INPUT; i+=2) { //set frequencies <FFT_TEST_HIGH_CUTT_OFF_FREQ Hz to zero
    if ((i<FcutHigh*2)||(i>(FFT_COMPLEX_INPUT-FcutHigh*2))) {
      m_fft_input_f32[i] = 0; // Real part.
      m_fft_input_f32[i+1] = 0; // Img part.
    }
  }
#endif

#if (FFT_USE_FILTER_LOW == 1)
  //low-pass filter
  double FcutLow = FFT_TEST_LOW_CUTT_OFF_FREQ / df; //set cutoff frequency as FFT_TEST_HIGH_CUTT_OFF_FREQ Hz

  for(uint32_t i = 0; i< FFT_COMPLEX_INPUT; i+=2) { //set frequencies >FFT_TEST_HIGH_CUTT_OFF_FREQ Hz to zero
    if (((i>FcutLow*2)&&(i<FFT_OUTPUT_SIZE))||((i>FFT_OUTPUT_SIZE) &&(i<(FFT_COMPLEX_INPUT-FcutLow*2)))) {
      m_fft_input_f32[i] = 0; // Real part.
      m_fft_input_f32[i+1] = 0; // Img part.
    }
  }  
#endif

  // Calculate the magnitude at each bin using Complex Magnitude Module function.
  arm_cmplx_mag_f32(m_fft_input_f32, m_fft_output_f32, FFT_OUTPUT_SIZE);

  //m_fft_output_f32[0] = 0.0; // Remove the DC offset

#if FFT_LOG_ENABLED
    //print_plotter(m_fft_output_f32, FFT_OUTPUT_SIZE); // full 
    print_plotter(m_fft_output_f32, FFT_OUTPUT_SIZE / 2); // N/2 (nyquist-theorm)
#endif

    // Compress the data to the beep format for storage and transmission.
    fft_create_result(&fft_result,
                      m_fft_output_f32, 
                      settings.param.audio_config.fft_count,
                      settings.param.audio_config.fft_start * 2,
                      settings.param.audio_config.fft_stop  * 2,
                      FFT_OUTPUT_SIZE / 2);

    /*if(audio.callback != NULL) {
      fft_result.param.meas_result.source = audio.source;
      audio.callback(&fft_result.param.meas_result);
    }*/

#if FFT_LOG_ENABLED
    //TRACE("Received block: %u\r\n", blocksTransferred);
#endif
    fftIsBusy = FALSE;
  }

  // Sample the TLV with the I2S interface until enough samples have been taken or sample for as long as the loop boolean is set.
  if(blocksTransferred >= blocksOffset && !loop) {
    //I2S_stop();
    //I2C_uninit();
    //audioFFT_deinit();
    //TLV_reset();
    //audio_app_nextState(AUDIO_PROCESS);

    //todo stop audio sampling
    //uinit adc / dma / timer / gpio / ...

    //return TRUE;
  }

  return FALSE;
}

void print_plotter(float const * p_data, uint16_t size) {
  uint16_t i;
  for (i = 0; i<size; i++)
    TRACE(TRACE_FLOAT_MARKER"\r\n", TRACE_FLOAT(p_data[i]));
}

void print_data(uint32_t const * p_block, uint16_t size) {
  uint16_t i;
  int16_t val;

  TRACE("ADC = [");
  
  for (i = 0; i<size; i++) {
    //val = (uint16_t)(p_block[i] >> 8);
    val = (uint16_t)(p_block[i]);
		if(i < (size - 1))
			TRACE("%i, ", val);
		else
			TRACE("%i ];\r\n", val);
	}	
}

/**
 * @brief Function for generating sine wave samples for FFT calculations.
 *
 * This function fill up output array with generated sine wave data with proper sampling frequency.
 * Must be executed before fft_process function.
 *
 * @param[in] p_input     Input array to fill with sine wave generated data.
 * @param[in] size        Input array size.
 * @param[in] sample_freq Sine wave sampling frequency.
 * @param[in] sine_freq   Sine wave frequency.
 * @param[in] add_noise   Flag for enable or disble adding noise for generated data.
 */
#if (FFT_USE_WINDOWING == 1)
void fft_generate_complexNumbersArray(float  *p_input, uint16_t p_inputSize, float *samples, uint16_t samplesSize) {
  uint32_t i, y;

  for (i=0;i<p_inputSize;i+=2) {
    y = (i) ? i / 2 : i;

    // Real part.
    p_input[i] = samples[y];

    // Img part.
    p_input[i+1] = 0;
  }
}
#else
void fft_generate_complexNumbersArray(float  *p_input, uint16_t p_inputSize, uint32_t *samples, uint16_t samplesSize)
{
  uint32_t i, y;

  for (i=0;i<p_inputSize;i+=2) {
    y = (i) ? i / 2 : i;

    // Real part.
    //p_input[complex_idx] = (float32_t) (__LL_ADC_CALC_DATA_TO_VOLTAGE(3300, samples[y], LL_ADC_RESOLUTION_12B) / 1000.0) - 1.25;
    p_input[i] = (float32_t) (__LL_ADC_CALC_DATA_TO_VOLTAGE(3300, samples[y], LL_ADC_RESOLUTION_12B) / 1000.0);

    // Img part.
    p_input[i+1] = 0;
  }
}
#endif

void fft_create_result(BEEP_protocol_s * ret, float * fft, uint16_t count, uint16_t start, uint16_t stop, uint16_t fftSize) {
  uint16_t binOffset;
  uint16_t i, j, sumNbins, diff;
  float binSum;
  FFT_RESULTS * result;

  float32_t max_value = 0;
  uint32_t  max_val_index = 0;
  
  // Search FFT max value in input array with an offset of 10 to ignore the DC part.
  uint32_t offset = 0; //10
  arm_max_f32(&fft[offset], fftSize - offset, &max_value, &max_val_index);
  max_val_index += offset;
  
  TRACE("fSample: %u Hz, res:%u Hz\r\n", (uint32_t)SAMPLING_FREQUENCY, (uint32_t) FFT_TEST_SAMPLE_RES_HZ); /* todo compute factor for frequency accuracy */
  TRACE("Max value: "TRACE_FLOAT_MARKER" at [%u] = %u Hz\r\n\r\n", TRACE_FLOAT(max_value), max_val_index, (uint32_t)(((float) max_val_index * (float) FFT_TEST_SAMPLE_RES_HZ) /* * 3.7878  factor frequency */));

//  #define DEFAULT_FFT_DIVIDER         10 //count
//  #define DEFAULT_FFT_START           0 //start
//  #define DEFAULT_FFT_STOP            UINT8_MAX //stop

  memset(ret, 0, sizeof(BEEP_protocol_s));
  result = &ret->param.meas_result.result.fft;
  result->bins = 0;
  result->start = start / 2; //0 / 2 = 0
  result->stop = stop  / 2; //255 / 2 = 127
  ret->param.meas_result.type = AUDIO_ADC;

  if(fftSize < count || fftSize == 0 || count > FFT_MAX_BINS || count == 0 || start >= stop)
    return;

  diff = stop - start; // 127
  sumNbins = diff / count; // 127 / 10 = 12

  /*
   * Round up the number of bins to sum when the number of bins times the size is smaller than the fft size
   * start=0, stop=255
   * 256 / 20 = 13.8
   * 13 * 20 = 240. 240 < 256 -> sumNbins = 14
   */
  if((count * sumNbins) < diff) //(10 * 12) < 127 then 121 
    {
        sumNbins++;
    }
    binOffset = start; 
     
    for(i=0; i<count; i++) //0 to 10
    {
        binSum = 0.0;

        for(j=0; j<sumNbins; j++) // 0 to 121
        {
            if(binOffset >= stop) 
            {
                break;
            }
            binSum += fft[binOffset];
            binOffset++;
        }

        #if 0
            // Average bin value
            float binAverage = binSum / (float) j;
            result->values[i] = (binAverage > UINT16_MAX) ? UINT16_MAX : (uint16_t) binAverage;
        #else
            // Total sum value per bin
            result->values[i] = (binSum > UINT16_MAX) ? UINT16_MAX : (uint16_t) binSum;
        #endif
        result->bins++;
    }

#if FFT_LOG_ENABLED
    TRACE("FFT result start %u/%u Hz, stop: %u/%u Hz, bin count: %u, samples/bin %u ", start, start * FFT_TEST_SAMPLE_RES_HZ,  stop, stop * FFT_TEST_SAMPLE_RES_HZ, count, sumNbins);

    for(i=0; i<count; i++)
    {
        // Total sum value per bin
        TRACE("[%u] = %u, ", i, result->values[i]);
    }

    TRACE("\r\n");
#endif
}

/**
 * @brief Function for drawing line and processed data informations.
 * @param[in] input_sine_freq Input sine wave frequency.
 * @param[in] is_noisy        Flag if data is noisy.
 * @param[in] chart_width     Drawing chart height.
 */
static void draw_fft_header(float32_t input_sine_freq)
{
    TRACE("fSample: %u Hz, res:%u Hz\r\n", (uint32_t)input_sine_freq, (uint32_t) FFT_TEST_SAMPLE_RES_HZ);
}

void draw_fft_max(float * m_fft_output_f32, uint16_t data_size)
{
  float32_t max_value = 0;
  uint32_t  max_val_index = 0;
  
  // Search FFT max value in input array with an offset of 10 to ignore the DC part.
  uint32_t offset = 0;//10;
  arm_max_f32(&m_fft_output_f32[offset], data_size-offset, &max_value, &max_val_index);
  max_val_index += offset;
  
  TRACE("fSample: %u Hz, res:%u Hz\r\n", (uint32_t)SAMPLING_FREQUENCY, (uint32_t) FFT_TEST_SAMPLE_RES_HZ);
  TRACE("Max value: "TRACE_FLOAT_MARKER" at [%u] = %u Hz\r\n\r\n", TRACE_FLOAT(max_value), max_val_index, (max_val_index * FFT_TEST_SAMPLE_RES_HZ));
}

int main(void)
{
  HAL_Init();

  SystemClock_Config(); /* 16 mHz */

  MX_GPIO_Init();
  MX_UART1_Init();
  MX_DMA_Init();
  MX_ADC1_Init();
  MX_TIM2_Init();
	
  memset(adcData, 0x00, sizeof(adcData));
  
  HAL_TIM_Base_Start(&htim2);

#if (DEBUG_SAMPLING_FREQUENCY == 1)
  HAL_GPIO_WritePin(GPIOA, GPIO_PIN_8, GPIO_PIN_SET);
  HAL_TIM_OC_Start(&htim2, TIM_CHANNEL_1); //debug timer
#endif

	HAL_ADC_Start_DMA(&hadc1, (uint32_t *) adcData, sizeof(adcData));
  
	//settings.param.audio_config.channel;
  //settings.param.audio_config.gain;
  //settings.param.audio_config.volume;
  settings.param.audio_config.fft_count = DEFAULT_FFT_DIVIDER;
  settings.param.audio_config.fft_start = DEFAULT_FFT_START;
  settings.param.audio_config.fft_stop = DEFAULT_FFT_STOP;    
  //settings.param.audio_config.min6dB;   

  blocksTransferred = 0;

  //TRACE("Start Applications\r\n");

  while (TRUE) {
		if(processDSP() == TRUE)
      break;
  }
}

void SystemClock_Config(void) { /* 16 mHz */
  RCC_OscInitTypeDef RCC_OscInitStruct = {};
  RCC_ClkInitTypeDef RCC_ClkInitStruct = {};

  /** Configure the main internal regulator output voltage
  */
  __HAL_RCC_PWR_CLK_ENABLE();
  __HAL_PWR_VOLTAGESCALING_CONFIG(PWR_REGULATOR_VOLTAGE_SCALE3 /*PWR_REGULATOR_VOLTAGE_SCALE1*/);
  /** Initializes the RCC Oscillators according to the specified parameters
  * in the RCC_OscInitTypeDef structure.
  */
  RCC_OscInitStruct.OscillatorType = RCC_OSCILLATORTYPE_HSI;
  RCC_OscInitStruct.HSIState = RCC_HSI_ON;
  RCC_OscInitStruct.HSICalibrationValue = RCC_HSICALIBRATION_DEFAULT;
  RCC_OscInitStruct.PLL.PLLState = RCC_PLL_OFF;
  //RCC_OscInitStruct.PLL.PLLSource = RCC_PLLSOURCE_HSI;
  //RCC_OscInitStruct.PLL.PLLM = 14;
  //RCC_OscInitStruct.PLL.PLLN = 336;
  //RCC_OscInitStruct.PLL.PLLP = RCC_PLLP_DIV4;
  //RCC_OscInitStruct.PLL.PLLQ = 8;
  if (HAL_RCC_OscConfig(&RCC_OscInitStruct) != HAL_OK) {
    Error_Handler();
  }
  /** Initializes the CPU, AHB and APB buses clocks
  */
  RCC_ClkInitStruct.ClockType = RCC_CLOCKTYPE_HCLK | RCC_CLOCKTYPE_SYSCLK | RCC_CLOCKTYPE_PCLK1 | RCC_CLOCKTYPE_PCLK2;
  RCC_ClkInitStruct.SYSCLKSource = RCC_SYSCLKSOURCE_HSI; //RCC_SYSCLKSOURCE_PLLCLK;
  RCC_ClkInitStruct.AHBCLKDivider = RCC_SYSCLK_DIV1;
  RCC_ClkInitStruct.APB1CLKDivider = RCC_HCLK_DIV1;
  RCC_ClkInitStruct.APB2CLKDivider = RCC_HCLK_DIV1;

  if (HAL_RCC_ClockConfig(&RCC_ClkInitStruct, FLASH_LATENCY_3) != HAL_OK) {
    Error_Handler();
  }
}

static void MX_ADC1_Init(void)
{
  ADC_ChannelConfTypeDef sConfig = {0};

  hadc1.Instance = ADC1;
  hadc1.Init.ClockPrescaler = ADC_CLOCK_SYNC_PCLK_DIV4;
  hadc1.Init.Resolution = ADC_RESOLUTION_12B;
  hadc1.Init.ScanConvMode = DISABLE;
  hadc1.Init.ContinuousConvMode = DISABLE;
  hadc1.Init.DiscontinuousConvMode = DISABLE;
  hadc1.Init.ExternalTrigConvEdge = ADC_EXTERNALTRIGCONVEDGE_RISING;
  hadc1.Init.ExternalTrigConv = ADC_EXTERNALTRIGCONV_T2_TRGO;
  hadc1.Init.DataAlign = ADC_DATAALIGN_RIGHT;
  hadc1.Init.NbrOfConversion = 1;
  hadc1.Init.DMAContinuousRequests = ENABLE;
  hadc1.Init.EOCSelection = ADC_EOC_SINGLE_CONV;
  if (HAL_ADC_Init(&hadc1) != HAL_OK)
  {
    Error_Handler();
  }

  sConfig.Channel = ADC_CHANNEL_2;
  sConfig.Rank = 1;
  sConfig.SamplingTime = ADC_SAMPLETIME_3CYCLES;
  if (HAL_ADC_ConfigChannel(&hadc1, &sConfig) != HAL_OK)
  {
    Error_Handler();
  }
}

static void MX_TIM2_Init(void)
{
  TIM_ClockConfigTypeDef sClockSourceConfig = {0};
  TIM_MasterConfigTypeDef sMasterConfig = {0};
  TIM_OC_InitTypeDef sConfigOC = {0};

  htim2.Instance = TIM2;
  htim2.Init.Prescaler = 0;
  htim2.Init.CounterMode = TIM_COUNTERMODE_UP;
  htim2.Init.Period = (HAL_RCC_GetHCLKFreq() / SAMPLING_FREQUENCY) - 1;
  htim2.Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  htim2.Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_Base_Init(&htim2) != HAL_OK)
  {
    Error_Handler();
  }
  sClockSourceConfig.ClockSource = TIM_CLOCKSOURCE_INTERNAL;
  if (HAL_TIM_ConfigClockSource(&htim2, &sClockSourceConfig) != HAL_OK)
  {
    Error_Handler();
  }
  sMasterConfig.MasterOutputTrigger = TIM_TRGO_UPDATE;
  sMasterConfig.MasterSlaveMode = TIM_MASTERSLAVEMODE_DISABLE;
  if (HAL_TIMEx_MasterConfigSynchronization(&htim2, &sMasterConfig) != HAL_OK)
  {
    Error_Handler();
  }

#if (DEBUG_SAMPLING_FREQUENCY == 1)
  sConfigOC.OCMode = TIM_OCMODE_TOGGLE;
  sConfigOC.Pulse = 0;
  if (HAL_TIM_OC_ConfigChannel(&htim2, &sConfigOC, TIM_CHANNEL_1) != HAL_OK)
  {
    Error_Handler();
  }

  HAL_NVIC_ClearPendingIRQ(TIM2_IRQn);  // make sure that any pending interrupt is cleared
  HAL_NVIC_EnableIRQ(TIM2_IRQn);  

  __HAL_TIM_ENABLE_IT(&htim2, TIM_IT_CC1);
#endif
}

static void MX_DMA_Init(void) 
{
  __HAL_RCC_DMA2_CLK_ENABLE();
  __HAL_RCC_DMA1_CLK_ENABLE();

  HAL_NVIC_SetPriority(DMA1_Stream5_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(DMA1_Stream5_IRQn);

  HAL_NVIC_SetPriority(DMA2_Stream0_IRQn, 0, 0);
  HAL_NVIC_EnableIRQ(DMA2_Stream0_IRQn);
}

static void MX_GPIO_Init(void)
{
  GPIO_InitTypeDef GPIO_InitStruct = {0};

  __HAL_RCC_GPIOA_CLK_ENABLE();
  //__HAL_RCC_GPIOB_CLK_ENABLE();

#if (DEBUG_SAMPLING_FREQUENCY == 1)
  /*Configure GPIO pin : PA8 */
  GPIO_InitStruct.Pin = GPIO_PIN_8;
  GPIO_InitStruct.Mode = GPIO_MODE_OUTPUT_PP;
  GPIO_InitStruct.Pull = GPIO_NOPULL;
  GPIO_InitStruct.Speed = GPIO_SPEED_FREQ_VERY_HIGH; //GPIO_SPEED_FREQ_LOW;
  HAL_GPIO_Init(GPIOA, &GPIO_InitStruct);
#endif  
}

void MX_UART1_Init(void)
{
  huart1.Instance = USART1;        
  huart1.Init.BaudRate = 115200;
  huart1.Init.WordLength = UART_WORDLENGTH_8B;
  huart1.Init.StopBits = UART_STOPBITS_1;
  huart1.Init.Parity = UART_PARITY_NONE;
  huart1.Init.Mode = UART_MODE_TX_RX;
  huart1.Init.HwFlowCtl = UART_HWCONTROL_NONE;
  if (HAL_UART_Init(&huart1) != HAL_OK) {
    Error_Handler();
  }
}

void Error_Handler(void)
{

}

#ifdef  USE_FULL_ASSERT
/**
  * @brief  Reports the name of the source file and the source line number
  *         where the assert_param error has occurred.
  * @param  file: pointer to the source file name
  * @param  line: assert_param error line source number
  * @retval None
  */
void assert_failed(uint8_t *file, uint32_t line)
{ 
  /* USER CODE BEGIN 6 */
  /* User can add his own implementation to report the file name and line number,
     tex: printf("Wrong parameters value: file %s on line %d\r\n", file, line) */
  /* USER CODE END 6 */
}
#endif 