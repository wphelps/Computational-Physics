import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class printToFile {
	public static void main(String args[]) {
		try {
			PrintWriter dhiraj = new PrintWriter(new FileWriter("/Users/wphelps/Desktop/test.txt", false));
			dhiraj.print("Blah blah test, blah. newline \n second line... hopefully");
			dhiraj.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.out.println("Could not write to file!");
		}
	}
}
