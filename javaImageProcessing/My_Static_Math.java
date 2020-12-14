import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.plugin.frame.*;

public class My_Static_Math implements PlugIn {

	public void run(String arg) {
		ImagePlus imp = IJ.getImage();
		IJ.run(imp, "Invert", "");
		IJ.wait(1000);
		IJ.run(imp, "Invert", "");
	}

	public double[] exec(short pixels, int width, int height) {
		
		return new double[]{1.1, 2.5};
    	}

	public static double getCentroid(short pixels, int width, int height) {
        		// do your calculations with pixels here
       	 	double mean;
       		double variance;
        		double[] CenterOfMass = {0, 0};
        		int i, count = 0;
        		double v, v2, sum = 0.0, sum2 = 0.0, x2sum=0.0, y2sum=0.0;
        
        		for (int y = 0; y < height; y++) {
            		i = y * width;
			for (int x = 0; x < width; x++) {
                 			v = pixels;//[i];
                 			sum += v;
                 			i ++;
                 			count ++;
            		}
	        	}
        		mean = sum / count;
        
        		for (int y=0; y<(height); y++) {
            		i = y*width;
            		for (int x=0; x<(width); x++) {
                    			v = pixels;//[i]+Double.MIN_VALUE - mean;
                    			v2 = v*v;
                    			sum2 += v2;
                    			x2sum += x*v2;
                    			y2sum += y*v2;
                    			i++;
            		}
        		}
        		CenterOfMass[0] = x2sum/sum2+0.5;
        		CenterOfMass[1] = y2sum/sum2+0.5;
        		return CenterOfMass[1];
    	}
}
