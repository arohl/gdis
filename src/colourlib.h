/**********************************************************/
/***			colorLib.h 			***/
/**********************************************************/

#define ALPHA_MASK 0xff000000
void    hsv2rgb(float,float,float,float*,float*,float*);
void 	rgb2hsv(float,float,float,float *,float *,float *);

void    hsv2lrgb(float,float,float,unsigned long *);


