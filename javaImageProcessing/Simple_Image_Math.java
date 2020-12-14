package org.will.hardiman;

public class Simple_Image_Math {

public double mean;
public double variance;
public double xCenterOfMass;
public double yCenterOfMass;
public int width;
public int height;
public float[] pixels;

    public void SimpleImageMath(int width, int height) {
        this.width = width;
        this.height = height;
        }
        
    public void setPixels(Object pixels){
        this.pixels = (float[])pixels;
        getStats();
        getCentroid();
        }
        
    // Basically steal ImageJ's centroid calculation
    private void getCentroid() {
        // do your calculations with pixels here
        int i;
        double v, v2, sum2=0.0, x2sum=0.0, y2sum=0.0;
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
        xCenterOfMass = x2sum/sum2+0.5;
        yCenterOfMass = y2sum/sum2+0.5;
    }
    
    private void getStats() {
        // Calculate the mean of the image
        int i, count = 0;
        double v, sum = 0.0, sum2 = 0.0;
        for (int y = 0; y < height; y++) {
            i = y * width;
            for (int x = 0; x < width; x++) {
                 v = pixels[i];
                 sum += v;
                 sum2 += v * v;
                 i ++;
                 count ++;
            }
        }
        mean = sum / count;
        double mean2 = mean * mean;                
    }
}
