public class Static_Image_Math{

    // Basically steal ImageJ's centroid calculation
    public static double[] getCentroid(short[] pixels, int width, int height) {
        // Finds centroid for object based on diff ^2 from mean
        // Good for bright objects with gray BG and no other objects in frame
        double mean;
        double variance;
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
        double variance;
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
}
