#include "main.h"
#include <stdarg.h>
#include <string.h>
#include <stdio.h>

#include "arm_math.h"
#include "arm_const_structs.h"
#include "arm_common_tables.h"

#define DEBUG_SAMPLING_FREQUENCY      1
#define FFT_USE_WINDOWING             1
#define FFT_USE_FILTER_HIGH           0
#define FFT_USE_FILTER_LOW            0
#define FTT_USE_FREQ_FACTOR           1
#define FFT_LOG_ENABLED               1

#define FFT_SAMPLE_FREQ_HZ            (4000) /* only : 2000 / 4000 / 8000 / 16000 Hz */
#define FFT_OUTPUT_SIZE               1024UL
#define FTT_COUNT                     2 
#define FFT_ADC_SAMPLING_SIZE         (FFT_OUTPUT_SIZE * 2)
#define FFT_ADC_BUFFER_SIZE           (FFT_OUTPUT_SIZE * (FTT_COUNT + (FTT_COUNT % 2)))
#define FFT_COMPLEX_INPUT             (FFT_OUTPUT_SIZE * 2UL)

#if (FTT_USE_FREQ_FACTOR == 1)
  #if (FFT_SAMPLE_FREQ_HZ == 2000)
    #define FFT_SAMPLE_FREQ_FACTOR      (1.9526) /* For sampling 2000Hz, res 1.91Hz */
  #elif (FFT_SAMPLE_FREQ_HZ == 4000)
    #define FFT_SAMPLE_FREQ_FACTOR      (1.94) /* For sampling 4000Hz, res 3.79Hz */
  #elif  (FFT_SAMPLE_FREQ_HZ == 8000)
    #define FFT_SAMPLE_FREQ_FACTOR      (1.9466) /* For sampling 8000Hz, res 7.6Hz */
  #elif  (FFT_SAMPLE_FREQ_HZ == 16000)
    #define FFT_SAMPLE_FREQ_FACTOR      (1.9451) /* For sampling 16000Hz, res 15.20Hz */      
  #else
    //try sample 4100Hz for beep with ratio 1.95
    #define FFT_SAMPLE_FREQ_FACTOR      (1.950) /* For sampling 16000Hz */      
  #endif
#else
  /*
    Apply in input 100Hz, 200Hz, 500Hz, 1000Hz, 1500Hz, 2000Hz, etc ....
    And read the frequency in trace (for exemple 54.69HZ): Max value: 15.6448 at [7] = 54.69 Hz" 
    And divide each target frequency by each measure : ex 1000 / 54.69 HZ = Ratio for 1000 Hz target
    Finish by compute average of each ratio (FFT_SAMPLE_FREQ_FACTOR = average of each ratio)
  */
  #define FFT_SAMPLE_FREQ_FACTOR      (1.0)
#endif  
#define FFT_SAMPLE_RES_HZ             ((float32_t)(((float32_t)FFT_SAMPLE_FREQ_HZ / (float32_t)FFT_COMPLEX_INPUT) * (float32_t)FFT_SAMPLE_FREQ_FACTOR))
#define FFT_HIGH_CUTT_OFF_FREQ        50 /* Hz */
#define FFT_LOW_CUTT_OFF_FREQ         50 /* Hz */

ADC_HandleTypeDef hadc1;
DMA_HandleTypeDef hdma_adc1;
DMA_HandleTypeDef hdma_dac1;
TIM_HandleTypeDef htim2;
UART_HandleTypeDef huart1;

uint32_t adcData[FFT_ADC_BUFFER_SIZE];

static uint8_t adcDataCount = 0;
static uint8_t fftPerformed = 0;
static bool fftStarted = FALSE;

FFT_RESULTS fftResult[FTT_COUNT];

char sz_traceBuffer[10000];
size_t sz_traceBufferSize;

const static arm_cfft_instance_f32 *S;

static float m_fft_output_f32[FFT_OUTPUT_SIZE];             //!< FFT output data. Frequency domain.
static float m_fft_input_f32[FFT_COMPLEX_INPUT] = {0};     //!< FFT input array for complex numbers. Time domain.

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

void HAL_ADC_ConvHalfCpltCallback(ADC_HandleTypeDef* hadc) {
  if(++adcDataCount >= FTT_COUNT)
    HAL_ADC_Stop_DMA(&hadc1);
  
  if(!fftStarted)
    fftStarted = TRUE;

  HAL_GPIO_TogglePin(GPIOA, GPIO_PIN_8);
}

void HAL_ADC_ConvCpltCallback(ADC_HandleTypeDef* hadc) {
  if(++adcDataCount < FTT_COUNT) 
    HAL_ADC_Start_DMA(&hadc1, (uint32_t*) &adcData[/*pwr of 2*/adcDataCount*FFT_OUTPUT_SIZE], FFT_ADC_BUFFER_SIZE);
  
  if(!fftStarted) 
    fftStarted = TRUE;
  
  HAL_GPIO_TogglePin(GPIOA, GPIO_PIN_8);
}

#if (DEBUG_SAMPLING_FREQUENCY == 1)
void HAL_TIM_OC_DelayElapsedCallback(TIM_HandleTypeDef* htim)
{
  //HAL_GPIO_TogglePin(GPIOA, GPIO_PIN_8);
}
#endif

bool processDSP(void) {
  uint16_t adcIndex = fftPerformed*FFT_OUTPUT_SIZE;

  if(!fftStarted)
    return FALSE;
  
  if(fftPerformed > FTT_COUNT)
    return FALSE;

  TRACE("New fft %d:\r\n", fftPerformed);  

#if (FFT_USE_WINDOWING == 1)   
    //double f64_fftSamples[FFT_OUTPUT_SIZE];
    float32_t f32_fftSamples[FFT_OUTPUT_SIZE];
    uint32_t i;

#if 0
    /* Trace input sample for current fft */
    TRACE("ADC samples:\r\n");
    for(i=adcIndex;i<adcIndex + FFT_OUTPUT_SIZE;i++) {
     TRACE("[ %d | %u ] , ", i, adcData[i]);
    }
    TRACE("\r\n");
#endif

  float32_t dcComponent = 0.0;

  for(i=adcIndex;i<adcIndex + FFT_OUTPUT_SIZE;i++) {
      f32_fftSamples[i - adcIndex] = (float32_t) (__LL_ADC_CALC_DATA_TO_VOLTAGE(3300, adcData[i], LL_ADC_RESOLUTION_12B) / 1000.0);

      dcComponent += f32_fftSamples[i - adcIndex];
  }

  dcComponent = dcComponent / (float32_t) FFT_OUTPUT_SIZE;

  TRACE("FFT DC component %f :\r\n", dcComponent);

  for(i=0;i<FFT_OUTPUT_SIZE;i++)
      f32_fftSamples[i] -= dcComponent;

/*  flag_window : 0 : no-window (default -- that is, other than 1 ~ 6)
 *                1 : parzen window
 *                2 : welch window
 *                3 : hanning window
 *                4 : hamming window
 *                5 : blackman window
 *                6 : steeper 30-dB/octave rolloff window
 */
    //windowing (FFT_OUTPUT_SIZE, f64_fftSamples, 4 /* hamming */, 1.0, f64_fftSamples);

    /* Convert the uint32_t array containing the samples of the Audio ADC to an float array with complex numbers. The real part will be placed on the even 
     * indexes and the imaginary part will be set to 0 on all uneven indexes. This means that the complex input array is twice the size of the number of
     * samples.
     */
    fft_generate_complexNumbersArray(m_fft_input_f32, FFT_COMPLEX_INPUT, f32_fftSamples, FFT_OUTPUT_SIZE);
#else 
    /* Convert the uint32_t array containing the samples of the Audio ADC to an float array with complex numbers. The real part will be placed on the even 
     * indexes and the imaginary part will be set to 0 on all uneven indexes. This means that the complex input array is twice the size of the number of
     * samples.
     */
    fft_generate_complexNumbersArray(m_fft_input_f32, FFT_COMPLEX_INPUT, fftSamples, FFT_OUTPUT_SIZE);
#endif
    
    //print_plotter(m_fft_input_f32, FFT_COMPLEX_INPUT);

    /* Use CFFT module to process the data.
     * Ensure that the used module is equal to the number of complex number of samples.
     */
  switch (FFT_OUTPUT_SIZE) {
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

#if (FFT_USE_FILTER_HIGH == 1)
  //high-pass filter
  float32_t FcutHigh = FFT_HIGH_CUTT_OFF_FREQ / FFT_SAMPLE_RES_HZ; 
  
  for(uint32_t i=0;i<FFT_COMPLEX_INPUT;i+=2) { //set frequencies <FFT_HIGH_CUTT_OFF_FREQ Hz to zero
    if(((float32_t)i < (FcutHigh * 2)) || ((float32_t)i > ((float32_t) FFT_COMPLEX_INPUT - (FcutHigh * 2)))) {
      m_fft_input_f32[i] = 0; // Real part.
      m_fft_input_f32[i+1] = 0; // Img part.
    }
  }
#endif

#if (FFT_USE_FILTER_LOW == 1)
  //low-pass filter
  float32_t FcutLow = FFT_LOW_CUTT_OFF_FREQ / FFT_SAMPLE_RES_HZ;

  for(uint32_t i = 0; i< FFT_COMPLEX_INPUT; i+=2) { //set frequencies >FFT_HIGH_CUTT_OFF_FREQ Hz to zero
    if ((((float32_t)i > (FcutLow*2)) && (i < FFT_OUTPUT_SIZE)) || ((i > FFT_OUTPUT_SIZE) &&((float32_t)i < ((float32_t)FFT_COMPLEX_INPUT - (FcutLow * 2))))) {
      m_fft_input_f32[i] = 0; // Real part.
      m_fft_input_f32[i+1] = 0; // Img part.
    }
  }  
#endif

  // Calculate the magnitude at each bin using Complex Magnitude Module function.
  arm_cmplx_mag_f32(m_fft_input_f32, m_fft_output_f32, FFT_OUTPUT_SIZE);

  m_fft_output_f32[0] = 0.0; // Remove first bin correspond to the DC component

#if FFT_LOG_ENABLED
    //print_plotter(m_fft_output_f32, FFT_OUTPUT_SIZE); // full 
    //print_plotter(m_fft_output_f32, FFT_OUTPUT_SIZE / 2); // N/2 (nyquist-theorm)

    print_fft_max(m_fft_output_f32, FFT_OUTPUT_SIZE / 2); // N/2 (nyquist-theorm)
#endif

    // Compress the data to the beep format for storage and transmission.
#if 1
    fft_create_result(&fftResult[fftPerformed],
                      m_fft_output_f32, 
                      20 /*settings.param.audio_config.fft_count*/,
                      0 /*settings.param.audio_config.fft_start * 2*/,
                      13 /*settings.param.audio_config.fft_stop  * 2*/,
                      FFT_OUTPUT_SIZE / 2);
#else
    fft_create_result(&fft_result,
                      m_fft_output_f32, 
                      settings.param.audio_config.fft_count,
                      settings.param.audio_config.fft_start * 2,
                      settings.param.audio_config.fft_stop  * 2,
                      FFT_OUTPUT_SIZE / 2);
#endif

  // Sample the TLV with the I2S interface until enough samples have been taken or sample for as long as the loop boolean is set.
  //if(blocksTransferred >= blocksOffset && !loop) {
    //I2S_stop();
    //I2C_uninit();
    //audioFFT_deinit();
    //TLV_reset();
    //audio_app_nextState(AUDIO_PROCESS);

    //todo stop audio sampling
    //uinit adc / dma / timer / gpio / ...

  //  return FALSE;
  //}

  if(++fftPerformed >= adcDataCount)
    fftStarted = FALSE;

  return FALSE;
}

void print_plotter(float const * p_data, uint16_t size) {
  uint16_t i;
  for (i = 0; i<size; i++)
    TRACE(TRACE_FLOAT_MARKER"\r\n", TRACE_FLOAT(p_data[i]));
}

void print_fft_max(float * m_fft_output_f32, uint16_t data_size) {
  float32_t max_value = 0;
  uint32_t  max_val_index = 0;
  
  // Search FFT max value in input array with an offset of 10 to ignore the DC part.
  uint32_t offset = 0;//10;
  arm_max_f32(&m_fft_output_f32[offset], data_size-offset, &max_value, &max_val_index);
  max_val_index += offset;
  
  TRACE("fSample: %u Hz, res:%0.2f Hz\r\n", (uint32_t) FFT_SAMPLE_FREQ_HZ, (float32_t) FFT_SAMPLE_RES_HZ);
  TRACE("Max value: "TRACE_FLOAT_MARKER" at [%u] = %0.2f Hz\r\n\r\n", TRACE_FLOAT(max_value), max_val_index, (float32_t) max_val_index * FFT_SAMPLE_RES_HZ);
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
void fft_generate_complexNumbersArray(float32_t  *p_input, uint16_t p_inputSize, float32_t *samples, uint16_t samplesSize) {
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
void fft_generate_complexNumbersArray(float32_t  *p_input, uint16_t p_inputSize, uint32_t *samples, uint16_t samplesSize)
{
  uint32_t i, y;

  double sum = 0;

  for (i=0;i<p_inputSize;i+=2) {
    y = (i) ? i / 2 : i;

    // Real part.
    p_input[i] = (float32_t) (__LL_ADC_CALC_DATA_TO_VOLTAGE(3300, samples[y], LL_ADC_RESOLUTION_12B) / 1000.0);

    sum += p_input[i];

    TRACE("FFT sum value %f, value inserted %f\r\n", sum, p_input[i]);

    // Img part.
    p_input[i+1] = 0;
  }

  sum = sum / (double)((p_inputSize / 2));

  TRACE("FFT sum value %f :\r\n", (float32_t) sum);

}
#endif

#if 1
void fft_create_result(FFT_RESULTS *p_fftResult, float *p_fftValue, uint16_t binOutputCount, uint16_t binInputOffset, uint16_t binCountInOneOutput, uint16_t fftSize) {
  uint16_t binsOutput, binInputCurrent, i;
  float32_t binOutputSum;
  bool end = FALSE;

  if((p_fftResult == NULL) || (p_fftValue == NULL) || !binOutputCount || (binOutputCount > FFT_MAX_BINS) || !fftSize)
    return;

  p_fftResult->start = binInputOffset; // start input bin
  p_fftResult->stop = binCountInOneOutput; // n input bin for one output bin

  //TRACE("FFT compute ouput:\r\n");

  for(binsOutput=0;binsOutput<binOutputCount;binsOutput++) {
    binOutputSum = 0.0;

    //TRACE("%02d -> ", binsOutput);

    for(i=0;i<binCountInOneOutput;i++) {
      binInputCurrent = binInputOffset + ((binsOutput * binCountInOneOutput) + i);

      if(binInputCurrent >= fftSize) {
        end = TRUE;
        break;
      }  

      binOutputSum += p_fftValue[binInputCurrent];
      //TRACE("%03d: %0.3f   ", binInputCurrent, p_fftValue[binInputCurrent]);
    }

    //TRACE("\r\n");

    p_fftResult->values[binsOutput] = (binOutputSum > UINT16_MAX) ? UINT16_MAX : (uint16_t) binOutputSum; // check diff on 32 bits !!!

    if(end)
      break;
  }

  p_fftResult->bins = binsOutput; //n output bin computed

#if 1 //FFT_LOG_ENABLED
  float32_t startFrequency = binInputOffset * FFT_SAMPLE_RES_HZ;
  float32_t binFrequency = (float32_t) binCountInOneOutput * FFT_SAMPLE_RES_HZ;
  uint16_t binEnd = binInputOffset + (p_fftResult->bins * binCountInOneOutput);

  if(binEnd >= fftSize)
    binEnd = fftSize;

  TRACE("FFT result start: %u / %0.2f Hz, end: %u / %0.2f Hz, out bin freq: %0.2f Hz, out bin count: %u, in bin by out bin: %u\r\n", binInputOffset, startFrequency,  binEnd, (float32_t) binEnd * FFT_SAMPLE_RES_HZ, binFrequency, p_fftResult->bins, binCountInOneOutput);

  for(binsOutput=0;binsOutput<p_fftResult->bins;binsOutput++) {
    TRACE("%02d -> brut s_bin%f_%fHz,\t\ts_bin%04d_%04dHz = %u\r\n", binsOutput,
                                               (float32_t) (startFrequency + ((float32_t) binsOutput * binFrequency)), //todo round
                                               (float32_t) (startFrequency + ((float32_t) binsOutput * binFrequency) + binFrequency), //todo round 
                                               (uint16_t) (startFrequency + ((float32_t) binsOutput * binFrequency)), //todo round
                                               (uint16_t) (startFrequency + ((float32_t) binsOutput * binFrequency) + binFrequency), //todo round
                                               p_fftResult->values[binsOutput]);


  }
#endif
}
#else
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
  
  TRACE("fSample: %u Hz, res:%0.2f Hz\r\n", (uint32_t)FFT_SAMPLE_FREQ_HZ, (float32_t) FFT_SAMPLE_RES_HZ); /* todo compute factor for frequency accuracy */
  TRACE("Max value: "TRACE_FLOAT_MARKER" at [%u] = %0.2f Hz\r\n\r\n", TRACE_FLOAT(max_value), max_val_index, (float32_t) max_val_index * FFT_SAMPLE_RES_HZ);

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
  if((count * sumNbins) < diff) { //(10 * 12) < 127 then 121 
    sumNbins++;
  }
  binOffset = start; 
     
  for(i=0; i<count; i++) { //0 to 10
    binSum = 0.0;

    for(j=0; j<sumNbins; j++) { // 0 to 121
      if(binOffset >= stop) 
        break;
    
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
  TRACE("FFT result start %u/%0.2f Hz, stop: %u/%0.2f Hz, bin count: %u, samples/bin %u ", start, (float32_t) start * FFT_SAMPLE_RES_HZ,  stop, (float32_t) stop * FFT_SAMPLE_RES_HZ, count, sumNbins);

  for(i=0;i<count;i++) {
    // Total sum value per bin
    TRACE("[%u] = %u, ", i, result->values[i]);
  }

  TRACE("\r\n");

  //print_plotter(result->values, count);
#endif
}
#endif

int main(void) {
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

	HAL_ADC_Start_DMA(&hadc1, (uint32_t *) adcData, FFT_ADC_BUFFER_SIZE);

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
  hadc1.Init.DMAContinuousRequests = DISABLE; //ENABLE; /* circular or not */
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
  htim2.Init.Period = (HAL_RCC_GetHCLKFreq() / FFT_SAMPLE_FREQ_HZ) - 1;
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
