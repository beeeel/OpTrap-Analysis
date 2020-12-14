package org.will.hardiman;

public class Short_Math {

    // Basically steal ImageJ's centroid calculation
    static double[] getShortCentroid(short pixels, int width, int height) {
        // do your calculations with pixels here
        double mean;
        double variance;
        double[] CenterOfMass = {0, 0};
        int i, count = 0;
        double v, v2, sum = 0.0, sum2 = 0.0, x2sum=0.0, y2sum=0.0;
        
        for (int y = 0; y < height; y++) {
            i = y * width;
            for (int x = 0; x < width; x++) {
                 v = pixels;
                 sum += v;
                 i ++;
                 count ++;
            }
        }
        mean = sum / count;
        
        for (int y=0; y<(height); y++) {
            i = y*width;
            for (int x=0; x<(width); x++) {
                    v = pixels+Double.MIN_VALUE - mean;
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
}