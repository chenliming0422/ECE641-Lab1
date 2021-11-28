#include <math.h>
#include "lib\tiff.h"
#include "lib\allocate.h"
#include "lib\randlib.h"
#include "lib\typeutil.h"

void error(char *name);
double prior_model(double **img, int M, int N, int u, int v, double g[3][3]);
double pair_wise_gaussian(double **img, int M, int N, int u, int v, double g[3][3]);

int main (int argc, char **argv) 
{
    FILE *fp;
    struct TIFF_img input_img, output_img;
    double **img1, **x, **y;
    int pixel;
    double sigma_x, sigma_w, sigma_x_2;
    double cost[20] = {0.0};
    double cost_1, cost_2;
    int32_t i,j, iter;
    int M,N;
    double decrease_val = 0;
    double temp;
    double sigma_ratio;

    double prediction[3][3] = {{1/12.0, 1/6.0, 1/12.0}, 
                               {1/6.0, 0.0, 1/6.0}, 
                               {1/12.0, 1/6.0, 1/12.0}};

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
    y = (double **)get_img(input_img.width,input_img.height,sizeof(double));
    x = (double **)get_img(input_img.width,input_img.height,sizeof(double));

    /* copy image pixels to double array */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        img1[i][j] = input_img.mono[i][j];
    }

    /* Set seed for random noise generator */
    srandom2(1);

    /* Add noise to image */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        img1[i][j] += 16.0*normal();
    }

    /* set up structure for output achromatic image */
    /* to allocate a full color image use type 'c' */
    get_TIFF ( &output_img, input_img.height, input_img.width, 'g' );

    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)img1[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;
        output_img.mono[i][j] = pixel;
        y[i][j] = pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/noise.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file noise.tif\n");
        exit ( 1 );
    }

    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }

    /* close output image file */
    fclose ( fp );

    /* compute MAP estimate using 20 iterations of ICD optimization */
    /* sigma_x = 16.95 from section 2.1 */
    sigma_x = 16.95;
    sigma_w = 16.0;

    sigma_ratio = pow(sigma_w, 2.0) / pow(sigma_x, 2.0);

    /* Initialize with ML estimate */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        x[i][j] = y[i][j];
    }

    /* iteration */
    M = input_img.height;
    N = input_img.width;
    for (iter = 0; iter < 20; iter++){
        for ( i = 0; i < input_img.height; i++ )
        for ( j = 0; j < input_img.width; j++ ) {
            temp = prior_model(x, M, N, i, j, prediction);
            temp = (y[i][j] + sigma_ratio * temp) / (1.0 + sigma_ratio);

            if(temp > 0) {
                x[i][j] = temp;
            }
            else {
                x[i][j] = 0.0;
            } 
        }

        cost_1 = 0.0;
        cost_2 = 0.0;
        for ( i = 0; i < input_img.height; i++ )
        for ( j = 0; j < input_img.width; j++ ) {
            cost_1 += pow((y[i][j] - x[i][j]), 2.0);
            cost_2 += pair_wise_gaussian(x, M, N, i, j, prediction) / 2.0;
        }
        cost[iter] = cost_1 / (2.0 * sigma_w * sigma_w) + cost_2 / (2.0 * sigma_x * sigma_x);


        if(iter>0) decrease_val = cost[iter-1] - cost[iter]; 
        //printf("iteration %02d: cost = %.15f, decrease: %.15f\n", iter, cost[iter], decrease_val);
        printf("%.15f\n", cost[iter]);
    }

    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)x[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;
        output_img.mono[i][j] = pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/MAPestimate1.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file MAPestimate1.tif\n");
        exit ( 1 );
    }

    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }

    /* close output image file */
    fclose ( fp );

    /* compute MAP estimate using 20 iterations of ICD optimization: 2 */
    sigma_x_2 = 5.0 * sigma_x * sigma_x;
    sigma_ratio = pow(sigma_w, 2.0) / sigma_x_2;
    
    /* Initialize with ML estimate */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        x[i][j] = y[i][j];
    }

    for (iter = 0; iter < 20; iter++){
        for ( i = 0; i < input_img.height; i++ )
        for ( j = 0; j < input_img.width; j++ ) {
            temp = prior_model(x, M, N, i, j, prediction);
            temp = (y[i][j] + sigma_ratio * temp) / (1.0 + sigma_ratio);

            if(temp > 0) {
                x[i][j] = temp;
            }
            else {
                x[i][j] = 0.0;
            } 
        }
    }

    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)x[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;
        output_img.mono[i][j] = pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/MAPestimate2.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file MAPestimate2.tif\n");
        exit ( 1 );
    }

    /* write output image */
    if ( write_TIFF ( fp, &output_img ) ) {
        fprintf ( stderr, "error writing TIFF file" );
        exit ( 1 );
    }

    /* close output image file */
    fclose ( fp );

    /* compute MAP estimate using 20 iterations of ICD optimization: 3 */
    sigma_x_2 = (1/5.0) * sigma_x * sigma_x;
    sigma_ratio = pow(sigma_w, 2.0) / sigma_x_2;

    /* Initialize with ML estimate */
    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        x[i][j] = y[i][j];
    }

    for (iter = 0; iter < 20; iter++){
        for ( i = 0; i < input_img.height; i++ )
        for ( j = 0; j < input_img.width; j++ ) {
            temp = prior_model(x, M, N, i, j, prediction);
            temp = (y[i][j] + sigma_ratio * temp) / (1.0 + sigma_ratio);

            if(temp > 0) {
                x[i][j] = temp;
            }
            else {
                x[i][j] = 0.0;
            } 
        }
    }

    for ( i = 0; i < input_img.height; i++ )
    for ( j = 0; j < input_img.width; j++ ) {
        pixel = (int32_t)x[i][j];
        if(pixel > 255) pixel = 255;
        if(pixel < 0) pixel = 0;
        output_img.mono[i][j] = pixel;
    }

    /* open output image file */
    if ( ( fp = fopen ( "../output/MAPestimate3.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file MAPestimate3.tif\n");
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
    free_img( (void**)y);
    free_img( (void**)x );
    
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
        val += pow((img[u][v] - img[(u+i-1+M)%M][(v+j-1+N)%N]), 2.0) * g[i][j];
    }
    return val;
}
