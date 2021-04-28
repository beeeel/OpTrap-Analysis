public class oneCentroid{

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
    // Kinda crap tbh
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
    // Good for bright objects with no BG and no other objects in frame
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

public static double[] getCentroidQuart(short[] pixels, int width, int height, short thresh) {
    // Finds centroid for object based on diff ^4 from mean
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
	        v2 = v*v*v*v;
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

// END OF CLASS
}

