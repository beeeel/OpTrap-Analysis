import java.io.File;
import java.io.FileOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataOutputStream;
import java.io.IOException;

public class myTools{

	// Writes short array to file as byte array. Needs to be given the length of the array
	public static void writeShortArrayToFile(short[][] array, int len, int inner, String path)
	{ 
		File file;
		FileOutputStream fos = null;
		byte[] bytesArray;

		try
		{
		file = new File(path);
		fos = new FileOutputStream(file);
		if (!file.exists()){
			file.createNewFile();
		}
		ByteArrayOutputStream bos = new ByteArrayOutputStream();
		DataOutputStream dos = new DataOutputStream(bos);
		for (int i = 0; i < len; i++)
		{
			for (int idx = 0; idx < inner; idx++)
			{
		 		dos.writeShort(array[i][idx]);
		 	}
		 	dos.flush();
			bytesArray = bos.toByteArray();
			fos.write(bytesArray);
			bos.reset();
		}
		fos.flush();
		fos.close();
		}
		catch (IOException ioe) {
			ioe.printStackTrace();
		}
		finally {
	  try {
	     if (fos != null)
	     {
		 fos.close();
	     }
          }
	  catch (IOException ioe) {
	     System.out.println("Error in closing the Stream");
	  }
       }
	}



	// Writes double array to file as byte array. Needs to be given the length of the array
	public static void writeDoubleArrayToFile(double[] array, int len, String path)
	{ 
		File file;
		FileOutputStream fos = null;
		byte[] bytesArray;

		try
		{
		file = new File(path);
		fos = new FileOutputStream(file);
		if (!file.exists()){
			file.createNewFile();
		}
		ByteArrayOutputStream bos = new ByteArrayOutputStream();
		DataOutputStream dos = new DataOutputStream(bos);
		for (int idx = 0; idx < len; idx++)
		{
			dos.writeDouble(array[idx]);
		}
		dos.flush();
		bytesArray = bos.toByteArray();
		fos.write(bytesArray);
		fos.flush();
		fos.close();
		}
		catch (IOException ioe) {
			ioe.printStackTrace();
		}
		finally {
	  try {
	     if (fos != null)
	     {
		 fos.close();
	     }
          }
	  catch (IOException ioe) {
	     System.out.println("Error in closing the Stream");
	  }
       }
	}

}
	

