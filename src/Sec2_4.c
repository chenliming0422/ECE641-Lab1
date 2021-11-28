#include <math.h>
#include "lib\tiff.h"
#include "lib\allocate.h"
#include "lib\randlib.h"
#include "lib\typeutil.h"

void error(char *name);
double prior_model(double **img, int M, int N, int u, int v, double g[3][3]);
double pair_wise_gaussian(double **img, int M, int N, int u, int v, double g[3][3]);
double circ_conv2d(double **img, int M, int N, int u, int v, double h[5][5]);
void update_error(double **e, int M, int N, int u, int v, double num, double h[5][5]);
double cost_function(double **y, double **x, int M, int N, double sigma_w, double sigma_x, double h[5][5], double g[3][3]); 

int main (int argc, char **argv) 
{
    FILE *fp;
    struct TIFF_img input_img, output_img;
    double **img1, **img2, **y, **x, **e;
    int32_t pixel;
    double sigma_x, sigma_w;
    double cost[20] = {0};
    int32_t i,j, iter;
    double v, theta_1, theta_2;
    int M,N;
    double filter[5][5] = {{1/81.0, 2/81.0, 3/81.0, 2/81.0, 1/81.0},
                           {2/81.0, 4/81.0, 6/81.0, 4/81.0, 2/81.0},
                           {3/81.0, 6/81.0, 9/81.0, 6/81.0, 3/81.0},
                           {2/81.0, 4/81.0, 6/81.0, 4/81.0, 2/81.0},
                           {1/81.0, 2/81.0, 3/81.0, 2/81.0, 1/81.0}};
    double prediction[3][3] = {{1/12.0, 1/6.0, 1/12.0}, {1/6.0, 0.0, 1/6.0}, {1/12.0, 1/6.0, 1/12.0}};

    if ( argc != 2 ) error( argv[0] );

      /* open image file */
    if ( ( fp = fopen ( argv[1], "rb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file %s\n", argv[1] );
        exit(1);
    }

    /* read image */
    if ( read_TIFF ( fp, &input_img ) ) {
        fprintf ( stderr, "error reading file %s\n", argv[1] );
        exit(1);
    }

    /* close image file */
     fclose ( fp );

    /* check the type of image data */
    if ( input_img.TIFF_type != 'g' ) {
        fprintf ( stderr, "error:  image must be gray scale\n" );
        exit ( 1 );
    }

    /* Allocate image of double precision floats */
    img1 = (double **)get_img(input_img.width,input_img.height,sizeof(double));
    img2 = (double **)get_img(input_img.width,input_img.height,sizeof(double));
    y = (double **)get_img(input_img.width,input_img.height,sizeof(double));
    x = (double **)get_img(input_img.width,input_img.height,sizeof(double));
    e = (double **)get_img(input_img.width,input_img.height,sizeof(double));


    M = input_img.height;
    N = input_img.width;

    /* copy image pixels to double array */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        img1[i][j] = input_img.mono[i][j];
    }

    /* blurring fitler */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        img2[i][j] = circ_conv2d(img1, M, N, i, j, filter);
    }

    /* Add noise */
    srandom2(1);
    
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        y[i][j] = img2[i][j] + 4 * normal();
    }

    get_TIFF ( &output_img, input_img.height, input_img.width, 'g' );
    
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)y[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;        
        output_img.mono[i][j] = pixel;
        y[i][j] = pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/Y.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file Y.tif\n");
        exit ( 1 );
    }
    
    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }

    /* MAP estimation using ICD*/
   sigma_w = 4.0;
   sigma_x = 16.95;

    /* Initialize x <- y */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        x[i][j] = y[i][j];
    }

    /* Initialize e <- y - Hx*/
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        e[i][j] = y[i][j] - circ_conv2d(x, M, N, i, j, filter);
    }

    theta_2 = pow((9/81.0), 2) + 8 * pow((2/81.0), 2) + 4 * (pow((6/81.0), 2) + pow((4/81.0), 2) + pow((3/81.0), 2) + pow((1/81.0), 2));
    theta_2 = theta_2 / (sigma_w * sigma_w);

    for (iter = 0; iter < 20; iter++)
    {
        for ( i = 0; i < input_img.height; i++ )
        for ( j = 0; j < input_img.width; j++ ) {
            v = x[i][j];
            theta_1 = circ_conv2d(e, M, N, i, j, filter);
            theta_1 = (-1.0) * theta_1 / (sigma_w * sigma_w);

            x[i][j] = prior_model(x, M, N, i, j, prediction);
            x[i][j] = (theta_2 * v - theta_1 + (1.0/(sigma_x * sigma_x)) * x[i][j]) / (theta_2 + (1.0/(sigma_x * sigma_x)));
  
            if(x[i][j] < 0) x[i][j] = 0.0;
            
            update_error(e, M, N, i, j, x[i][j] - v, filter);
        }

        cost[iter] = cost_function(y, x, M, N, sigma_w, sigma_x, filter, prediction);
        printf("iter %2d : cost = %.15f\n", iter, cost[iter]);
    }

    
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)x[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;
        output_img.mono[i][j] = (int32_t)pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/X.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file X.tif\n");
        exit ( 1 );
    }
    
    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }

    /* close output image file */
    fclose ( fp );

    /* de-allocate space which was used for the images */   
    free_TIFF ( &(input_img) );
    free_TIFF ( &(output_img) );

    free_img( (void**)img1 );
    free_img( (void**)img2 );

    return(0);

}

void error(char *name)
{
    printf("usage:  %s  image.tiff \n\n",name);
    printf("this program reads in a 24-bit color TIFF image.\n");
    printf("It then horizontally filters the green component, adds noise,\n");
    printf("and writes out the result as an 8-bit image\n");
    printf("with the name 'green.tiff'.\n");
    printf("It also generates an 8-bit color image,\n");
    printf("that swaps red and green components from the input image");
    exit(1);
}

double prior_model(double **img, int M, int N, int u, int v, double g[3][3])
{
    int i, j;
    double val = 0.0;
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
        val += (img[(u+i-1+M)%M][(v+j-1+N)%N]) * g[i][j];
    }
    return val;
}

double pair_wise_gaussian(double **img, int M, int N, int u, int v, double g[3][3])
{
    int i, j;
    double val = 0.0;
    for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
        val += pow((img[u][v] - img[(u+i-1+M)%M][(v+j-1+N)%N]), 2) * g[i][j];
    }
    return val;
}


double circ_conv2d(double **img, int M, int N, int u, int v, double h[5][5])
{
    int i, j;
    double val = 0.0;
    for (i = 0; i < 5; i++)
    for (j = 0; j < 5; j++){
        val += (img[(u+i-2+M)%M][(v+j-2+N)%N]) * h[i][j];
    }
    return val;
}

void update_error(double **e, int M, int N, int u, int v, double num, double h[5][5])
{
    int i, j;
    for (i = 0; i < 5; i++)
    for (j = 0; j < 5; j++) {
        e[(u+i-2+M)%M][(v+j-2+N)%N] -= h[i][j] * num;
    }
}

double cost_function(double **y, double **x, int M, int N, double sigma_w, double sigma_x, double h[5][5], double g[3][3])
{
    int i, j;
    double val = 0.0;
    double tmp = 0.0;
    for (i = 0; i < M; i++)
    for (j = 0; j < N; j++){
        tmp = y[i][j] - circ_conv2d(x, M, N, i, j, h);
        val += (1.0/(2 * sigma_w * sigma_w)) * pow(tmp, 2);
        val += (1.0/(2 * sigma_x * sigma_x)) * pair_wise_gaussian(x, M, N, i, j, g);
    }

    return val;
}

