#include "operator.h"

double circ_conv_size3(double data[3][3], double filter[3][3])
{
    double out_val = 0.0;
    int i, j;

    for (i = 0; i < filter_size; i++)
    {
        for(j = 0; j < filter_size; j++)
        {
            out_val += filter[i][j] * data[i][j];
        } 
    }

    return out_val;
}

double circ_conv_size5(double data[5][5], double filter[5][5])
{
    double out_val = 0.0;
    int i, j;

    for (i = 0; i < filter_size; i++)
    {
        for(j = 0; j < filter_size; j++)
        {
            out_val += filter[i][j] * data[i][j];
        } 
    }

    return out_val;
}