import java.lang.Math;

// Get one centroid for each area of FOV, cut at columns specified in subWidth
public class nCentroids{

// Calculate centroids by basic centre of brightness measure
public static double[] getNCentroidsBasic(short[] pixels, int width, int height, int[] subWidth, short thresh) 
{
    // Finds centroid for each object based on brightness
    int i, j, count, nB = length(subWidth);// Number of beads is number of subwidths - last one being the fullwidth
    double CenterOfMass[], v[], sum[], xsum[], ysum[];
    
    count = new int[ nB ];
    CentreOfMass = new double[ 2*nB ]; // 2D CoM per particle
    sum  = new double[ nB ];
    xsum = new double[ nB ];
    ysum = new double[ nB ];

    for (int y=0; y<(height); y++) // Iterator over the rows
    {
        // Start pixel count at the beginning of the row
		i = y*width;
        // Persistent counter along the row
        j = 0;
        for (int b=0; b<nB; b++) // Iterator based on bead number
        {
            for (int x=0; x<(subWidth[b]); x++) // Iterator along the row (within the subROI)
            {     
            // Apply theshold
            v = (pixels[i] > thresh) ? pixels[i] - thresh + Double.MIN_VALUE : 0;

            // Running totals
            sum[b]  += v;
			xsum[b] += x*v;
			ysum[b] += y*v;
			
            // Increment pixel count
			i++;
            }
		}
	}                                   
	
    for (int b=0; b < nB; b++) // Calculate CoM from running totals
    {
        CenterOfMass[2*b]   = xsum[b]/sum[b]+0.5;
        CenterOfMass[2*b+1] = ysum[b]/sum[b]+0.5;
    }
    return CenterOfMass;
}

// END OF CLASS
}

