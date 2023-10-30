#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "arrayhelp.h"
#include "filehelp.h"
#include <hdf5.h>
#include <hdf5_hl.h>

/* linspace - generates a linearly spaced array of numbers between two limits

   Inputs:
	x0 (double)    - the number to start with
	x1 (double)    - the number to end with
	n (int)        - the number of points to make
   Output:
	(double *)     - a pointer to an array with the values (length: n)

   Notes: The output pointer must be checked to ensure that it is not null and
 	must later be deallocated.

   Example: 
*/

double * linspace(double x0, double x1, int n) {
	int i;
	double *linspacearray = malloc(sizeof(double)*n);
	double steps = 1;

	if (n>1) steps = n-1;

	if (linspacearray != NULL) {
		for (i=0; i<steps; i++) {
			linspacearray[i] = x0 + i*(x1-x0)/((double)steps);
		}
		/* make the last one exactly right, no round-off */
		if (steps>1) linspacearray[n-1] = x1;
	}
	return(linspacearray);
}

/* logspace - generates a logarithmically spaced array of numbers between limits

   Inputs:
	x0    - the number to start with: base^x0
	x1    - the number to end with: base^x1
	n     - the number of points to make
	base  - the base of the logarithm

   Output:
	(double *)     - a pointer to an array with the values (length: n)

   Notes: The output pointer must be checked to ensure that it is not null and
 	must later be deallocated.

   Examples:

   y=logspace(-3,2,100,10) // 100 points between 0.001 and 100 
   y=logspace(log10(0.001),log10(100),100,10) // same 
*/ 

double * logspace(double x0, double x1, int n, double base) {
	int i;
	double *logspacearray = malloc(sizeof(double)*n);
	double stepsize = 0;
	double stepvalue = x0;

	if (n>0) stepsize = (x1-x0)/(n-1);

	if (logspacearray != NULL) {
		for (i=0; i<n-1; i++) {
			logspacearray[i] = pow(base,stepvalue);
			stepvalue += stepsize;
		}
		/* make the last one exactly right, no round-off */
		if (n>1) logspacearray[n-1] = pow(base,x1);
	}
	return(logspacearray);
}

/* printdoublearray - print a double array to standard out

   Inputs:
	double printarray[]  - An array of doubles
	int length           - The length of the array

   Output:
	None (void)

   Notes: The output pointer must be checked to ensure that it is not null and
 	must later be deallocated.

   Examples:

   y=logspace(-3,2,10,10); // 10 points between 0.001 and 100 
   printdoublearray(y,10); // print the array
   free(y); // free the memory
*/ 


void printdoublearray(double printarray[], int length) {
	int i;

	for (i=0;i<length;i++) {
		printf("%f ",printarray[i]);
	}
}

named_array named_array_create(char name[], double array_values[], int array_length) {
	named_array new_array;
	int i;

	strcpy(new_array.name,name);
	new_array.array_length = array_length;
	new_array.array_values = (double*)malloc(array_length*sizeof(double));
	if (new_array.array_values != NULL) {
		for (i=0;i<array_length;i++) {
			new_array.array_values[i] = array_values[i];
		}
	} else {
		printf("Error making new array.\n");
		exit(EXIT_FAILURE);
	}
	return(new_array);
}

void named_array_delete(named_array *todelete_array) {
	if ((*todelete_array).array_values != NULL) {
		free((*todelete_array).array_values);
		(*todelete_array).array_values = NULL;
		(*todelete_array).array_length = 0;
		strcpy((*todelete_array).name,"\n");
	}
}

int named_array_fwrite(FILE *stream, named_array *to_write) {
	int output;
	size_t countsize;

	output = fprintf(stream, "%s\n%d\n", (*to_write).name, (*to_write).array_length);

	if (output==0) {
		lendian_fwrite(((*to_write).array_values), sizeof(double), (*to_write).array_length, stream);
	}

	return(output);
}


int named_array_fread(FILE *stream, named_array *to_read) {
	char name[255];
	int array_length;
	double *data;
	int output;
	size_t countsize;

	fscanf(stream, "%s\n%d\n", name, &array_length);

	data = malloc(sizeof(double)*array_length);

	lendian_fread(data, sizeof(double), array_length, stream);

	strcpy((*to_read).name,name);
	(*to_read).array_length = array_length;
	(*to_read).array_values = data;

	return(0);
}

/* named_array_h5read - read an HDF5 array into a named array

   Inputs:
	const char *filename - The HDF5 file name
	const char *datapath - The data path to be read
	named_array *to_read - The new named_array that is read

   Output:
	status (always 0)

   This function opens a dataset in an HDF5 file. The dataset must be a 1xN or Mx1
   array of doubles. A named array with the name of the array (the dataset path is trimmed off)
   is created.

   Examples:

   // if you have a dataset called "/myDataset1/arrayName" in file "myfile.h5"
   named_array na2;
   int status;
   status = named_array_h5read("myfile.h5","/myDataset1/arrayName", &na2);
*/ 


int named_array_h5read(const char *filename, const char *datapath, named_array *to_read) {
	int status = 0;
	char name[255];
	char h5pathsep = '/';
	int i, j, lastSep;

	hid_t       file, dataset;         /* handles */
	hid_t       datatype, filespace;   
	hid_t       memspace; 
	H5T_class_t class;                 /* datatype class */
	H5T_order_t order;                 /* data order */
	size_t      size;                  
	hsize_t elements_bufsize;
	double *elements;

	hsize_t     dims[2];                     /* dataset and chunk dimensions*/ 
	herr_t      status_n, out_status;
	int rank, ctr;

	for (i=0;i<strlen(datapath);i++) { 
		if (h5pathsep==datapath[i]) {
			lastSep = i;
		}
	}
	for (i=lastSep+1;i<strlen(datapath);i++) {
		name[i - (lastSep+1)] = datapath[i];
	}
	for (i=strlen(datapath) - (lastSep+1);i<255;i++) {
		name[i] = (char)0;
	}

	//printf("name is %s.\n", name);

	//printf("About to open file %s with dataset %s\n", filename, datapath);

	file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset = H5Dopen1(file, datapath);

	datatype  = H5Dget_type(dataset);     // datatype handle
	class     = H5Tget_class(datatype);
	order     = H5Tget_order(datatype);
	//if (order == H5T_ORDER_LE) printf("Little endian order \n");
	H5Tclose(datatype);

	elements_bufsize = H5Dget_storage_size(dataset);
	//printf("Datasize: %llu\n", elements_bufsize);
	elements = (double *) malloc(elements_bufsize);

	filespace = H5Dget_space(dataset);    /* Get filespace handle first. */
	rank      = H5Sget_simple_extent_ndims(filespace);
	status_n  = H5Sget_simple_extent_dims(filespace, dims, NULL);
	/*printf("dataset rank %d, dimensions %lu x %lu\n",
	   rank, (unsigned long)(dims[0]), (unsigned long)(dims[1]));*/

	memspace = H5Screate_simple(rank,dims,NULL);
 
	/*
	* Read dataset back and display.
	*/
	out_status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace,
		     H5P_DEFAULT, elements);

	H5Sclose(filespace);
	H5Sclose(memspace);
	H5Dclose(dataset);
	H5Fclose(file);

        strcpy((*to_read).name,name);
        (*to_read).array_length = dims[0]*dims[1];
        (*to_read).array_values = elements;

	/* for testing
	ctr = 0;
	for (i=0;i<dims[0];i++) {
		for (j=0;j<dims[1];j++) {
			printf("%f ", elements[ctr]);
			ctr = ctr + 1;
		}
		printf("\n");
	}
	*/

	return(status);
}

int named_array_h5write(const char *filename, const char *datapath, named_array *to_read) {
	int status = 0;








	return(status);
}



 /* thanks to Bard for this one, well mostly Bard */

int write_double_array_to_hdf5(const double *array, size_t n_rows, size_t n_cols, const char *filename, const char *dataset_name) {
  // Open the HDF5 file for writing.
  hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file_id == H5I_INVALID_HID) {
    fprintf(stderr, "Error creating HDF5 file: %s\n", filename);
    return(1);
  }

  // Create the HDF5 dataset.
  hsize_t dims[] = {n_rows, n_cols};
  hid_t dataset_id = H5Dcreate(file_id, dataset_name, H5T_NATIVE_DOUBLE, H5Screate_simple(2, dims, NULL), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  if (dataset_id == H5I_INVALID_HID) {
    fprintf(stderr, "Error creating HDF5 dataset: %s\n", dataset_name);
    H5Fclose(file_id);
    return(1);
  }

  // Write the double array to the HDF5 dataset.
  herr_t status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, array);
  if (status < 0) {
    fprintf(stderr, "Error writing double array to HDF5 dataset: %s\n", dataset_name);
    H5Dclose(dataset_id);
    H5Fclose(file_id);
    return(1);
  }

  // Close the HDF5 dataset and file.
  H5Dclose(dataset_id);
  H5Fclose(file_id);
  return(0);
}


