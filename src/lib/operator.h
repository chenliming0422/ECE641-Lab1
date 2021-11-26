#ifndef _OPERATOR_H_
#define _OPERATOR_H_

void get_conv_data()

double circ_conv_size3(double data[3][3], double filter[3][3]);
double circ_conv_size5(double data[3][3], double filter[5][5]);

double pair_wise_error_filter3(double data[3][3], double filter[3][3]);
double pair_wise_error_filter5(double data[3][3], double filter[5][5]);

#endif