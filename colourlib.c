/**************************************************/
/***			colourLib.c						***/
/*** Somehow from Xrgb							***/
/**************************************************/

#include <glib.h>
#include "colourlib.h"



void hsv2lrgb(float hue,float saturation,float value,unsigned long *color)
{
	unsigned long ired,igreen,iblue,ialpha;	
	float red,green,blue;

    ialpha = (ALPHA_MASK & *color); /*** Keep whatever value of alpha already exists ***/
    hsv2rgb(hue, saturation, value, &red, &green, &blue);
    ired =   (unsigned long)(red * 255.0);
    igreen = (unsigned long)(green * 255.0);
    iblue =  (unsigned long)(blue * 255.0);
    iblue =  iblue << 16;
    igreen = igreen << 8;
    *color = (iblue | igreen | ired | ialpha);
}
	
		
void hsv2rgb(float hue, float saturation, float value, float *red, float *green, float *blue)
{
    int i;
    float f, m, n, k;

    hue = (hue > 360.0) ? 360.0 : (hue < 0.0) ? 0.0 : hue;
    saturation = (saturation > 1.0) ? 1.0 : (saturation < 0.0) ? 0.0 : saturation;
    value = (value > 1.0) ? 1.0 : (value < 0.0) ? 0.0 : value;

    if ( saturation == 0.0 ) { 
        *red = *green = *blue = value;
    } else { 
	if ( hue == 360.0 ) {
            hue = 0.0;
        } else {
            hue /= 60.0;
            i = hue; /* truncate */
            f = hue - (float)i;
            m = value * (1.0 - saturation);
            n = value * (1.0 - (saturation * f));
            k = value * (1.0 - (saturation * (1.0 - f)));

            switch (i) {
            case 0:
                *red = value;
                *green = k;
                *blue = m;
                break;
            case 1:
                *red = n;
                *green = value;
                *blue = m;
                break;
            case 2:
                *red = m;
                *green = value;
                *blue = k;
                break;
            case 3:
                *red = m;
                *green = n;
                *blue = value;
                break;
            case 4:
                *red = k;
                *green = m;
                *blue = value;
                break;
            case 5:
            default:
                *red = value;
                *green = m;
                *blue = n;
                break;
            }
        }
    }

    *red = (*red > 1.0) ? 1.0 : (*red < 0.0) ? 0.0 : *red;
    *green = (*green > 1.0) ? 1.0 : (*green < 0.0) ? 0.0 : *green;
    *blue = (*blue > 1.0) ? 1.0 : (*blue < 0.0) ? 0.0 : *blue;

    return;
}

/************************************************************/
/*	RGB to HSV converter									*/
/* From an algorithm in:									*/
/*  Foley, J. D., van Dam, A., et al. (1990). 				*/
/*	Computer Graphics - Principles and Practice, 			*/
/*	Addison-Wesley.											*/
/************************************************************/
void rgb2hsv(float red, float green, float blue, float *hue, float *saturation, float *value )
{

  float myMin,myMax, delta;
  
  myMin = 1.0;
  myMax = 0.0;
	
  if (red   < 0.0) red   = 0.0;
  if (green < 0.0) green = 0.0;
  if (blue  < 0.0) blue  = 0.0;
  
	myMax = MAX(red, green);
	myMax = MAX(myMax, blue);
	
	myMin = MIN(red, green);
	myMin = MIN(myMin, blue);
	
	*value = myMax;
	
	delta = myMax - myMin;
	
	if ( myMax > 0.0 ){
		*saturation = (myMax - myMin)/myMax;
	} else {
		*saturation = 0.0;
	}
	
	if (*saturation == 0.0){
		*hue = 0.0; 		/* supposed to be "undefined", but that doesn't really work for me */
	} else {
		if (red == myMax){
			*hue = (green - blue)/delta;			/* resulting colour is between yellow and magenta */
		} else if (green == myMax){
			*hue = 2 + (blue - red)/delta;		/* resulting colour is between cyan and yellow */
		} else if (blue == myMax){
			*hue = 4 + (red - green)/delta;		/* resulting colour is between magenta and cyan */
		}
		*hue = *hue * 60;							/* convert hue to degrees */
		
		if ( *hue < 0.0) *hue+= 360;				/* Make sure hue is non-negative */
		if ( *hue < 0.0) *hue = 0.0;
		*hue = *hue / 360.0;						/* finally convert to a value in the range [0.0 - 1.0] */
	}
}

