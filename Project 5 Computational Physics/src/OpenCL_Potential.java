import static org.jocl.CL.CL_CONTEXT_PLATFORM;
import static org.jocl.CL.CL_DEVICE_TYPE_ALL;
import static org.jocl.CL.CL_MEM_COPY_HOST_PTR;
import static org.jocl.CL.CL_MEM_READ_ONLY;
import static org.jocl.CL.CL_MEM_READ_WRITE;
import static org.jocl.CL.CL_TRUE;
import static org.jocl.CL.clBuildProgram;
import static org.jocl.CL.clCreateBuffer;
import static org.jocl.CL.clCreateCommandQueue;
import static org.jocl.CL.clCreateContext;
import static org.jocl.CL.clCreateKernel;
import static org.jocl.CL.clCreateProgramWithSource;
import static org.jocl.CL.clEnqueueNDRangeKernel;
import static org.jocl.CL.clEnqueueReadBuffer;
import static org.jocl.CL.clGetDeviceIDs;
import static org.jocl.CL.clGetPlatformIDs;
import static org.jocl.CL.clReleaseCommandQueue;
import static org.jocl.CL.clReleaseContext;
import static org.jocl.CL.clReleaseKernel;
import static org.jocl.CL.clReleaseMemObject;
import static org.jocl.CL.clReleaseProgram;
import static org.jocl.CL.clSetKernelArg;

import org.jocl.CL;
import org.jocl.Pointer;
import org.jocl.Sizeof;
import org.jocl.cl_command_queue;
import org.jocl.cl_context;
import org.jocl.cl_context_properties;
import org.jocl.cl_device_id;
import org.jocl.cl_kernel;
import org.jocl.cl_mem;
import org.jocl.cl_platform_id;
import org.jocl.cl_program;

public class OpenCL_Potential {
	 private static String programSource = "__kernel void " +
	 "sampleKernel(__global const float *x,"
	 + "__global const float *y," + " __global const float *z," +
	 "__global const float *x2,"
	 + "__global const float *y2," + " __global const float *z2," +
	 "  __global float *c)" + "{"
	 + "    int gid = get_global_id(0);" +
	 "    float xdiff = x[gid] - x2[gid];"
	 + "    float ydiff = y[gid] - y2[gid];" +
	 "    float zdiff = z[gid] - z2[gid];"
	 + "    float radius = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);"
	 + "c[gid] = .001/pow(radius,12) - .001 / pow(radius,6);" + "}";
	
	public int N;
	public float[] x_input;
	public float[] y_input;
	public float[] z_input;
	Pointer x;
	Pointer y;
	Pointer z;

	public float[] x2_input;
	public float[] y2_input;
	public float[] z2_input;
	Pointer x2;
	Pointer y2;
	Pointer z2;
	public double result;
	public float[] potential;
	Pointer out;

	public void setN(int N) {
		this.N = N;
	}

	public OpenCL_Potential() {
		// The platform, device type and device number
		// that will be used

	}
	public void Run() {
		// The platform, device type and device number
		// that will be used
		final int platformIndex = 0;
		final long deviceType = CL_DEVICE_TYPE_ALL;
		final int deviceIndex = 0;

		// Enable exceptions and subsequently omit error checks in this sample
		CL.setExceptionsEnabled(true);

		// Obtain the number of platforms
		int numPlatformsArray[] = new int[1];
		clGetPlatformIDs(0, null, numPlatformsArray);
		int numPlatforms = numPlatformsArray[0];

		// Obtain a platform ID
		cl_platform_id platforms[] = new cl_platform_id[numPlatforms];
		clGetPlatformIDs(platforms.length, platforms, null);
		cl_platform_id platform = platforms[platformIndex];

		// Initialize the context properties
		cl_context_properties contextProperties = new cl_context_properties();
		contextProperties.addProperty(CL_CONTEXT_PLATFORM, platform);

		// Obtain the number of devices for the platform
		int numDevicesArray[] = new int[1];
		clGetDeviceIDs(platform, deviceType, 0, null, numDevicesArray);
		int numDevices = numDevicesArray[0];

		// Obtain a device ID
		cl_device_id devices[] = new cl_device_id[numDevices];
		clGetDeviceIDs(platform, deviceType, numDevices, devices, null);
		cl_device_id device = devices[deviceIndex];

		// Create a context for the selected device
		cl_context context = clCreateContext(contextProperties, 1, new cl_device_id[]{device}, null, null, null);

		// Create a command-queue for the selected device
		cl_command_queue commandQueue = clCreateCommandQueue(context, device, 0, null);

		x = Pointer.to(x_input);
		y = Pointer.to(y_input);
		z = Pointer.to(z_input);

		x2 = Pointer.to(x2_input);
		y2 = Pointer.to(y2_input);
		z2 = Pointer.to(z2_input);
		potential = new float[N];
		out = Pointer.to(potential);
		cl_mem memObjects[] = new cl_mem[7];
		memObjects[0] = clCreateBuffer(context, // x
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_float * x_input.length, x, null);
		memObjects[1] = clCreateBuffer(context, // x
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_float * y_input.length, y, null);
		memObjects[2] = clCreateBuffer(context, // x
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_float * z_input.length, z, null);
		memObjects[3] = clCreateBuffer(context, // x
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_float * x2_input.length, x2, null);
		memObjects[4] = clCreateBuffer(context, // x
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_float * y2_input.length, y2, null);
		memObjects[5] = clCreateBuffer(context, // x
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, Sizeof.cl_float * z2_input.length, z2, null);
		memObjects[6] = clCreateBuffer(context, CL_MEM_READ_WRITE, Sizeof.cl_float * N, null, null);
		cl_program program = clCreateProgramWithSource(context, 1, new String[]{programSource}, null, null);

		// Build the program
		clBuildProgram(program, 0, null, null, null, null);

		// Create the kernel
		cl_kernel kernel = clCreateKernel(program, "sampleKernel", null);

		// Set the arguments for the kernel
		for (int i = 0; i < memObjects.length; i++) {
			clSetKernelArg(kernel, i, Sizeof.cl_mem, Pointer.to(memObjects[i]));
		}

		// Set the work-item dimensions
		long global_work_size[] = new long[]{N};
		long local_work_size[] = new long[]{1};

		// Execute the kernel
		clEnqueueNDRangeKernel(commandQueue, kernel, 1, null, global_work_size, local_work_size, 0, null, null);

		// Read the output data
		clEnqueueReadBuffer(commandQueue, memObjects[6], CL_TRUE, 0, N * Sizeof.cl_float, out, 0, null, null);

		for (int i = 0; i < memObjects.length; i++) {
			clReleaseMemObject(memObjects[i]);
		}

		clReleaseKernel(kernel);
		clReleaseProgram(program);
		clReleaseCommandQueue(commandQueue);
		clReleaseContext(context);

		result = 0.0;
		for (int i = 0; i < potential.length; i++) {
			result += (double) potential[i];
			//System.out.println(result);
		}

	}

}
