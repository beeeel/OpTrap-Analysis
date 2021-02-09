public class Static_Image_Math{

    // Basically steal ImageJ's centroid calculation
    public static double[] getCentroid(short[] pixels, int width, int height) {
        // Finds centroid for object based on diff ^2 from mean
        // Good for bright objects with gray BG and no other objects in frame
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
                    if (v < mean){
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


public static double getMean(short[] pixels, int width, int height) {
        // Finds mean of image
        double mean;
        int i, count = 0;
        double v, sum = 0.0;
        
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
        return mean;
        }
    
// Basically steal ImageJ's centroid calculation
public static double[] getTwoCentroids(short[] pixels, int width, int height) {
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
                 vl = pixels[i];
                 suml += vl;
                 i ++;
                 countl ++;
            }
        }
        mean = suml / countl;
      
        for (int y=0; y<(height); y++) {
		i = y*width;
		for (int x=0; x<(width/2); x++) {
			
			vl = pixels[i]+Double.MIN_VALUE - mean;
			v2l = vl*vl;
			sum2l += v2l;
			x2suml += x*v2l;
			y2suml += y*v2l;
			
			i++;
		}
		for (int x=width/2; x<(width); x++) {
			
			vr = pixels[i]+Double.MIN_VALUE - mean;
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
