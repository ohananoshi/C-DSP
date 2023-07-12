/* Header: my_DSP.h

    A set of mathematical functions useful for Digital Signal Processing.

    Author: Guilherme Arruda

    GitHub: https://github.com/ohananoshi/C_Projects/tree/main/Digital_Signal_Processing

    Created on: 23/02/2023

    Last updated: 03/03/2023
*/

#pragma once

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define MAX_ARRAY_SIZE 65535

typedef double (*f_x)(double);

double point_ret_to_polar(double point_x, double point_y, bool axis){

    if(axis){
        return sqrt(pow(point_x, 2.0) + pow(point_y, 2.0));
    }else return atan2(point_y, point_x);
}

double* signal_generate(f_x func,
                        double interval_start,
                        double interval_end,
                        double signal_source_length){

    if(signal_source_length > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    double* generated_signal = (double*)malloc((unsigned int)(signal_source_length+1)*sizeof(double));
    double step = (interval_end - interval_start)/signal_source_length;

    for(unsigned int i = 0; i <= signal_source_length; i++){
       generated_signal[i] = func(interval_start + ((double)i)*step);
    }

    return generated_signal;
}

double signal_mean(double* signal_source, unsigned int signal_source_length){

    if(signal_source_length > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    double mean = 0;

    for(unsigned int i = 0; i < signal_source_length; i++){
        mean += signal_source[i];
    }

    return mean/signal_source_length;
}

double signal_variance(double* signal_source, double signal_mean, unsigned int signal_source_length){

    if(signal_source_length > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    double variance;

    for(unsigned int i = 0; i < signal_source_length; i++){
        variance += pow(signal_source[i] - signal_mean, 2.0);
    }

    return variance;
}

double signal_std_deviation(double signal_variance){
    return sqrt(signal_variance);
}

double* signal_convolution(double* signal_source, unsigned int signal_source_length, double* impulse_response, unsigned int impulse_response_length){

    if(signal_source_length > MAX_ARRAY_SIZE  || impulse_response_length > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    unsigned int conv_out_length = signal_source_length+impulse_response_length-1;
    double *conv_signal = (double*)malloc(conv_out_length*sizeof(double));

    for(int i = 0; i < conv_out_length; i++){
        conv_signal[i] = 0;
    }

    for(int i = 0; i < signal_source_length; i++){
        for(int j = 0; j < impulse_response_length; j++){
            conv_signal[i+j] += signal_source[i]*impulse_response[j];
        }
    }

    return conv_signal;
}

double* signal_first_difference(double* signal_source, unsigned int signal_source_length){

    if(signal_source_length > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    double* out_signal = (double*)malloc(signal_source_length*sizeof(double));

    for(int i = 1; i < signal_source_length; i++){
        out_signal[i] = signal_source[i] - signal_source[i-1];
    }

    return out_signal;
}

double* signal_running_sum(double* signal_source, unsigned int signal_source_length){

    if(signal_source_length > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    double* out_signal = (double*)malloc(signal_source_length*sizeof(double));

    out_signal[0] = signal_source[0];

    for(int i = 1; i < signal_source_length; i++){
        out_signal[i] = signal_source[i] + out_signal[i-1];
    }

    return out_signal;
}

double* signal_DFT(double* signal_source, unsigned int signal_source_length, bool mode){

    if(signal_source_length > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    double* out_signal = (double*)malloc((signal_source_length/2)*sizeof(double));

    for(unsigned int k = 0; k < (signal_source_length/2); k++){
        out_signal[k] = 0;
    }

    if(mode){
        for(unsigned int i = 0; i < (signal_source_length/2); i++){
            for(unsigned int j = 0; j < signal_source_length; j++){
                out_signal[i] += signal_source[j]*cos(2*M_PI*i*j/signal_source_length);
            }
        }
    }else{
        for(unsigned int i = 0; i < (signal_source_length/2); i++){
            for(unsigned int j = 0; j < signal_source_length; j++){
                out_signal[i] += -signal_source[j]*sin(2*M_PI*i*j/signal_source_length);
            }
        }
    }

    return out_signal;
}

double* signal_complex_DFT(double* real_signal_source, double* imaginary_signal_source, unsigned int signal_source_length, bool mode){

    if((signal_source_length*2) > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    double* out_signal = (double*)malloc(signal_source_length*sizeof(double));

    for(unsigned int k = 0; k < signal_source_length; k++){
        out_signal[k] = 0;
    }

    if(mode){
        for(unsigned int i = 0; i < signal_source_length; i++){
            for(unsigned int j = 0; j < signal_source_length; j++){
                out_signal[i] += real_signal_source[i]*cos(2*M_PI*i*j/signal_source_length) + imaginary_signal_source[i]*sin(2*M_PI*i*j/signal_source_length);
            }
        }
    }else{
        for(unsigned int i = 0; i < signal_source_length; i++){
            for(unsigned int j = 0; j < signal_source_length; j++){
                out_signal[i] += -imaginary_signal_source[i]*(cos(2*M_PI*i*j/signal_source_length) + sin(2*M_PI*i*j/signal_source_length));
            }
        }
    }

    return out_signal;
}

double* signal_inverse_DFT(double* real_signal_source, double* imaginary_signal_source, unsigned int signal_source_length){

    if((signal_source_length*2) > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Output array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    double* out_signal = (double*)malloc((2*signal_source_length)*sizeof(double));
    double* temp_real_signal = (double*)malloc(signal_source_length*sizeof(double));
    double* temp_imaginary_signal = (double*)malloc(signal_source_length*sizeof(double));

    for(unsigned int i = 0; i < signal_source_length; i++){
        temp_real_signal[i] = real_signal_source[i]/signal_source_length;
        temp_imaginary_signal[i] = -imaginary_signal_source[i]/signal_source_length;

    }

    for(unsigned int i = 0; i < 2*signal_source_length; i++){
        out_signal[i] = 0;
    }

    temp_real_signal[0] = real_signal_source[0]/2;
    temp_imaginary_signal[0] = -imaginary_signal_source[0]/2;

    for(int i = 0; i < signal_source_length; i++){
        for(int j = 0; j < (2*signal_source_length); j++){
            out_signal[j] += temp_real_signal[i]*cos(M_PI*j*i/signal_source_length) + temp_imaginary_signal[i]*sin(M_PI*i*j/signal_source_length);
        }
    }

    free(temp_real_signal);
    free(temp_imaginary_signal);

    return out_signal;
}

double* signal_amplitude(double* real_signal_source, double* imaginary_signal_source, unsigned int signal_source_length){

    if(signal_source_length > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    double* amplitude = (double*)malloc(signal_source_length*sizeof(double));

    for(int i = 0; i < signal_source_length; i++){
        amplitude[i] = sqrt(pow(real_signal_source[i],2.0) + pow(imaginary_signal_source[i], 2.0));
    }

    return amplitude;
}

double* signal_phase(double* real_signal_source, double* imaginary_signal_source, unsigned int signal_source_length){

    if(signal_source_length > MAX_ARRAY_SIZE){
        fprintf(stderr, "ERROR. Array size exceeds %d elements.", MAX_ARRAY_SIZE);
        system("pause");
        exit(2);
    }

    double* phase = (double*)malloc(signal_source_length*sizeof(double));

    for(int i = 0; i < signal_source_length; i++){
        phase[i] = atan2(imaginary_signal_source[i], real_signal_source[i]);
    }

    return phase;
}
