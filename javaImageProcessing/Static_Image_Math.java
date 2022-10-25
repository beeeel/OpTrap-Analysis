public class Static_Image_Math{

// Centroid by weighted variance
public static double[] getCentroid(short[] pixels, int width, int height, short thresh) {
    // Finds centroid for object based on diff ^2 from mean
    // Good for bright objects with gray BG and no other objects in frame
    double mean;
    double[] CenterOfMass = {0, 0};
    int i, count = 0;
    double v, v2, sum = 0.0, sum2 = 0.0, x2sum=0.0, y2sum=0.0;
    
    for (int y = 0; y < height; y++) {
        i = y * width;
        for (int x = 0; x < width; x++) {
            v = (pixels[i] > thresh) ? pixels[i] : 0;
            sum += v;
            count ++;
            i ++;
        }
    }
    mean = sum / count;
    
    for (int y=0; y<(height); y++) {
        i = y*width;
        for (int x=0; x<(width); x++) {
            v = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE - mean : 0;
	        v2 = v*v;
            sum2 += v2;
            x2sum += x*v2;
            y2sum += y*v2;
            i++;
        }
    }
    CenterOfMass[0] = x2sum/sum2+0.5;
    CenterOfMass[1] = y2sum/sum2+0.5;
    return CenterOfMass;
}

// Centroid by brightness squared
public static double[] getCentroidBasicSq(short[] pixels, int width, int height, short thresh) {
    // Finds centroid for object based on brightness squared
    // Good for bright objects with gray BG and no other objects in frame
    double mean;
    double[] CenterOfMass = {0, 0};
    int i, count = 0;
    double v, v2, sum = 0.0, sum2 = 0.0, x2sum=0.0, y2sum=0.0;
    
    for (int y=0; y<(height); y++) {
        i = y*width;
        for (int x=0; x<(width); x++) {
            v = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE : 0;
            v2 = v*v;
            sum2 += v2;
            x2sum += x*v2;
            y2sum += y*v2;
            i++;
        }
    }
    CenterOfMass[0] = x2sum/sum2+0.5;
    CenterOfMass[1] = y2sum/sum2+0.5;
    return CenterOfMass;
}

// Centroid by simple brightness
public static double[] getCentroidBasic(short[] pixels, int width, int height, short thresh) {
    // Finds centroid for object based on brightness
    // Good for bright objects with gray BG and no other objects in frame
    double mean;
    double[] CenterOfMass = {0, 0};
    int i, count = 0;
    double v, sum = 0.0, xsum=0.0, ysum=0.0;
    
    for (int y=0; y<(height); y++) {
        i = y*width;
        for (int x=0; x<(width); x++) {
            v = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE : 0;
            sum += v;
            xsum += x*v;
            ysum += y*v;
            i++;
        }
    }
    CenterOfMass[0] = xsum/sum+0.5;
    CenterOfMass[1] = ysum/sum+0.5;
    return CenterOfMass;
}

// Centroid by variance method to find dark object
public static double[] getDarkCentroid(short[] pixels, int width, int height) {
    // Find centre of object which is darker than the mean of the image
    // Ignores all pixels where v >= mean
    double mean;
    double[] CenterOfMass = {0, 0};
    int i, count = 0;
    double v, v2, sum = 0.0, sum2 = 0.0, x2sum=0.0, y2sum=0.0;
    
    for (int y = 0; y < height; y++) {
        i = y * width;
        for (int x = 0; x < width; x++) {
             v = pixels[i];
             sum += v;
             i ++;
             count ++;
        }
    }
    mean = sum / count;
    
    for (int y=0; y<(height); y++) {
        i = y*width;
        for (int x=0; x<(width); x++) {
                v = pixels[i]+Double.MIN_VALUE - mean;
                v2 = v*v;
                if (v < 0){
                    sum2 += v2;
                    x2sum += x*v2;
                    y2sum += y*v2;
                }
                i++;
        }
    }
    CenterOfMass[0] = x2sum/sum2+0.5;
    CenterOfMass[1] = y2sum/sum2+0.5;
    return CenterOfMass;
}


public static double getMean(short[] pixels, int width, int height, int YSkip) {
    // Finds mean of image
    double mean;
    int i, count = 0;
    double v, sum = 0.0;
    
    for (int y = YSkip; y < height-YSkip; y++) {
        i = y * width;
        for (int x = 0; x < width; x++) {
             v = pixels[i];
             sum += v;
             i ++;
             count ++;
        }
    }
    mean = sum / count;
    return mean;
}
    
// Calculate centroids by weighted variance measure
public static double[] getTwoCentroids(short[] pixels, int width, int height, short thresh) {
    // Finds centroid for object based on diff ^2 from mean
    // Good for bright objects with gray BG and no other objects in frame
    // Get one centroid for left half of image and one for right half
    double mean;
    double[] CenterOfMass = {0, 0, 0, 0};
    int i, countl = 0, countr = 0;
    double vl, v2l, suml = 0.0, sum2l = 0.0, x2suml=0.0, y2suml=0.0;
    double vr, v2r,             sum2r = 0.0, x2sumr=0.0, y2sumr=0.0;
    
    for (int y = 0; y < height; y++) {
	// Index of first pixel on row
        i = y * width;
        for (int x = 0; x < width; x++) {
            vl = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE : 0;
            suml += vl;
            i ++;
            countl ++;
        }
    }
    mean = suml / countl;
    
    for (int y=0; y<(height); y++) {
	    i = y*width;
    	for (int x=0; x<(width/2); x++) {
            vl = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE - mean : 0;
		    v2l = vl*vl;
    		sum2l += v2l;
			x2suml += x*v2l;
	    	y2suml += y*v2l;
		    i++;
	    }
	    for (int x=width/2; x<(width); x++) {
            vr = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE - mean : 0;
		    v2r = vr*vr;
    		sum2r += v2r;
			x2sumr += x*v2r;
	    	y2sumr += y*v2r;
		    i++;
	}
	}                                   
	
	CenterOfMass[0] = x2suml/sum2l+0.5;
    CenterOfMass[1] = y2suml/sum2l+0.5;
    CenterOfMass[2] = x2sumr/sum2r+0.5;
    CenterOfMass[3] = y2sumr/sum2r+0.5;
    return CenterOfMass;
    }

// Calculate centroids by centre of brightness squared measure
public static double[] getTwoCentroidsBasicSq(short[] pixels, int width, int height, short thresh) {
    // Finds centroid for object based brightness
    // Get one centroid for left half of image and one for right half
    double mean;
    double[] CenterOfMass = {0, 0, 0, 0};
    int i, countl = 0, countr = 0;
    double vl, v2l, suml = 0.0, sum2l = 0.0, x2suml=0.0, y2suml=0.0;
    double vr, v2r,             sum2r = 0.0, x2sumr=0.0, y2sumr=0.0;
    
    for (int y=0; y<(height); y++) {
		i = y*width;
		for (int x=0; x<(width/2); x++) {
            
            vl = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE : 0;
			v2l = vl*vl;
			sum2l += v2l;
			x2suml += x*v2l;
			y2suml += y*v2l;
			
			i++;
		}
		for (int x=width/2; x<(width); x++) {
			
            vr = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE : 0;
			v2r = vr*vr;
			sum2r += v2r;
			x2sumr += x*v2r;
			y2sumr += y*v2r;

			i++;
		}
	}                                   
	
    CenterOfMass[0] = x2suml/sum2l+0.5;
    CenterOfMass[1] = y2suml/sum2l+0.5;
    CenterOfMass[2] = x2sumr/sum2r+0.5;
    CenterOfMass[3] = y2sumr/sum2r+0.5;
    return CenterOfMass;
}

// Calculate centroids by basic centre of brightness measure
public static double[] getTwoCentroidsBasic(short[] pixels, int width, int height, short thresh) {
    // Finds centroid for object based brightness
    // Get one centroid for left half of image and one for right half
    double mean;
    double[] CenterOfMass = {0, 0, 0, 0};
    int i, countl = 0, countr = 0;
    double vl, suml = 0.0, xsuml=0.0, ysuml=0.0;
    double vr, sumr = 0.0, xsumr=0.0, ysumr=0.0;
    
    for (int y=0; y<(height); y++) {
		i = y*width;
		for (int x=0; x<(width/2); x++) {
			
            vl = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE : 0;
			suml += vl;
			xsuml += x*vl;
			ysuml += y*vl;
			
			i++;
		}
		for (int x=width/2; x<(width); x++) {
			
            vr = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE : 0;
			sumr += vr;
			xsumr += x*vr;
			ysumr += y*vr;

			i++;
		}
	}                                   
	
    CenterOfMass[0] = xsuml/suml+0.5;
    CenterOfMass[1] = ysuml/suml+0.5;
    CenterOfMass[2] = xsumr/sumr+0.5;
    CenterOfMass[3] = ysumr/sumr+0.5;
    return CenterOfMass;
}

// Calculate centroids by centre of brightness squared measure
public static double[] getTwoCentroidsBasicSqROI(short[] pixels, int width, int height, int subWidth, short thresh) {
    // Finds centroid for object based brightness
    // Get one centroid for left half of image and one for right half
    double mean;
    double[] CenterOfMass = {0, 0, 0, 0};
    int i, countl = 0, countr = 0;
    double vl, v2l, suml = 0.0, sum2l = 0.0, x2suml=0.0, y2suml=0.0;
    double vr, v2r,             sum2r = 0.0, x2sumr=0.0, y2sumr=0.0;
    
    mean = getMean(pixels, width, height, 0);

    for (int y=0; y<(height); y++) {
		i = y*width;
		for (int x=0; x<(subWidth); x++) {
            
            vl = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE - mean : 0;
			v2l = vl*vl;
			sum2l += v2l;
			x2suml += x*v2l;
			y2suml += y*v2l;
			
			i++;
		}
    	for (int x=(width-subWidth); x<(width); x++) {
			
            vr = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE - mean : 0;
			v2r = vr*vr;
			sum2r += v2r;
			x2sumr += x*v2r;
			y2sumr += y*v2r;

			i++;
		}
	}                                   
	
    CenterOfMass[0] = x2suml/sum2l+0.5;
    CenterOfMass[1] = y2suml/sum2l+0.5;
    CenterOfMass[2] = x2sumr/sum2r+0.5;
    CenterOfMass[3] = y2sumr/sum2r+0.5;
    return CenterOfMass;
}

// Calculate centroids for part of image by basic centre of brightness measure
public static double[] getTwoCentroidsBasicROI(short[] pixels, int width, int height, int subWidth, short thresh) {
    // Finds centroid for object based brightness
    // Get one centroid for left half of image and one for right half
    double mean;
    double[] CenterOfMass = {0, 0, 0, 0};
    int i, countl = 0, countr = 0;
    double vl, suml = 0.0, xsuml=0.0, ysuml=0.0;
    double vr, sumr = 0.0, xsumr=0.0, ysumr=0.0;
    
    for (int y=0; y<(height); y++) {
	    i = y*width;
    	for (int x=0; x<(subWidth); x++) {
		
            vl = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE : 0;
		    suml += vl;
    		xsuml += x*vl;
	    	ysuml += y*vl;
		
		    i++;
    	}
	    i = y*width + (width-subWidth);
    	for (int x=(width-subWidth); x<(width); x++) {
		
            vr = (pixels[i] > thresh) ? pixels[i] + Double.MIN_VALUE : 0;
		    sumr += vr;
    		xsumr += x*vr;
	    	ysumr += y*vr;

		    i++;
    	}
	}                                   
	
    CenterOfMass[0] = xsuml/suml+0.5;
    CenterOfMass[1] = ysuml/suml+0.5;
    CenterOfMass[2] = xsumr/sumr+0.5;
    CenterOfMass[3] = ysumr/sumr+0.5;
    return CenterOfMass;
}

// END OF CLASS
}

//        for (int y=0; y<(height); y++) {
//        i = y*width;
//        for (int x=0; x<(width); x++) {
//        vl = (x < (width/2)) & (pixels[i]+Double.MIN_VALUE - mean);
//        v2l = vl*vl;
//        sum2l += v2l;
//        x2suml += x*v2l;
//        y2suml += y*v2l;
//        
//        vr = (x >= (width/2)) & (pixels[i]+Double.MIN_VALUE - mean);
//        v2r = vr*vr;
//        sum2r += v2r;
//        x2sumr += x*v2r;
//        y2sumr += y*v2r;
//        
//        i++;
//        
//        }
//        }

//        for (int y=0; y<(height); y++) {
//        i = y*width;
//        for (int x=0; x<(width); x++) {
//        if (x < (width/2))
//        		    {
//        		    vl = pixels[i]+Double.MIN_VALUE - mean;
//        		    v2l = vl*vl;
//        		    sum2l += v2l;
//        		    x2suml += x*v2l;
//        		    y2suml += y*v2l;
//        		    i++;
//        		    }
//        		    else
//        		    {
//        		    vr = pixels[i]+Double.MIN_VALUE - mean;
//        		    v2r = vr*vr;
//        		    sum2r += v2r;
//        		    x2sumr += x*v2r;
//        		    y2sumr += y*v2r;
//        		    i++;
//        		    }
//        		    }
//        		    }                                       
