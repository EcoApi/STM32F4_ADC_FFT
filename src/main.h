/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.h
  * @brief          : Header for main.c file.
  *                   This file contains the common defines of the application.
  ******************************************************************************
  * @attention
  *
  * <h2><center>&copy; Copyright (c) 2019 STMicroelectronics.
  * All rights reserved.</center></h2>
  *
  * This software component is licensed by ST under BSD 3-Clause license,
  * the "License"; You may not use this file except in compliance with the
  * License. You may obtain a copy of the License at:
  *                        opensource.org/licenses/BSD-3-Clause
  *
  ******************************************************************************
  */
/* USER CODE END Header */

/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __MAIN_H
#define __MAIN_H

#ifdef __cplusplus
extern "C" {
#endif

/* Includes ------------------------------------------------------------------*/
#include "stm32f4xx_hal.h"

#define TRUE 1
#define FALSE 0

typedef uint8_t bool;
typedef float float32_t;

#define FFT_MAX_BINS 12

typedef struct
{
    uint8_t     bins;
    uint8_t     start;
    uint8_t     stop;
    uint16_t    values[FFT_MAX_BINS];
}FFT_RESULTS;

typedef enum
{
	DS18B20		= 0,
	BME280		= 1,
	HX711		= 2,
	AUDIO_ADC	= 3,
    nRF_ADC		= 4,
    SQ_MIN		= 5,
    ATECC		= 6,
	BUZZER		= 7,
	LORAWAN		= 8,
	MX_FLASH	= 9,
    nRF_FLASH	= 10,
	APPLICATIE	= 11,
}SENSOR_TYPE;

//#define DEFAULT_AUDIO_CHANNEL       AIN_IN2RP
#define DEFAULT_AUDIO_COARSE_GAIN   (20*2)
#define DEFAULT_AUDIO_VOLUME        0
#define DEFAULT_AUDIO_INPUT_DIV     false
#define DEFAULT_FFT_DIVIDER         10
#define DEFAULT_FFT_START           0
#define DEFAULT_FFT_STOP            UINT8_MAX

typedef struct
	{
		uint8_t		channel;
        uint8_t     gain;
        int8_t      volume;
        uint8_t     fft_count;
        uint8_t     fft_start;
        uint8_t     fft_stop;    
        bool        min6dB;    
	}AUDIO_CONFIG_s;

    // The communication interface origin.
	typedef enum
	{
		INTERNAL_SOURCE,
		BLE_SOURCE,
		LORAWAN_SOURCE,
		UNKNOWN_SOURCE,
	}CONTROL_SOURCE;

typedef struct
{
	SENSOR_TYPE		type;
  CONTROL_SOURCE	source;

    union
    {
	    //DS18B20_RESULTS_s		ds18B20;
        //BME280_RESULT_s			bme280;
        //ADC_s					saadc;	
        //HX711_CONV_s			hx711;
        FFT_RESULTS             fft;
        //SQ_ORIENTATION_s		sq;
	}result;
}MEASUREMENT_RESULT_s;

typedef struct
{
	//BEEP_CID        command;
    //CONTROL_SOURCE  source;
		
		// Union to generalize data access.
    union
	{
		//RESPONSE_s				error;
        //STATUS_s				status;
        //FH_VERSION_s			version;
        //DS18B20_STATE_s			ds_state;
        //DS18B20_CONFIG_s		ds_config;
        AUDIO_CONFIG_s          audio_config;
        MEASUREMENT_RESULT_s	meas_result;	// Measurement storage container for all sensors.
        //ATTECC_s				atecc_id;
        //BUZZER_s				buzz;
        //SQ_MIN_s				sq;
        //LORAWAN_s				lorawan_key;
        //nRF_FLASH_s				flash;
        //uint32_t                size;
        //ALARM_CONFIG_s          alarm;
        //BME_CONFIG_s            bme_config;
	}param;
}BEEP_protocol_s;



/* Exported functions prototypes ---------------------------------------------*/
void Error_Handler(void);

#ifdef __cplusplus
}
#endif

#endif /* __MAIN_H */

/************************ (C) COPYRIGHT STMicroelectronics *****END OF FILE****/
