// Predict transmembrane lipid exposure using hhblits pssms and support vector machines
// (C) Tim Nugent May 2012

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "svm.h"

#define MAXSEQLEN 5000

using namespace std;

int main(int argc, const char* argv[]){

	float nullfreqs[20];
	float freqs[MAXSEQLEN][20];
	int pos = 0, i = 1;
	int window = 3;
	int precision = 5;
	int getfreqs = 0;
	int found_null_freqs = 0;
	char buf[1024];
	string model_file = "model/lipex.model";
	string hhmake_file;

	if (argc < 2){
		printf("Usage : %s [-m <svm model file>] <hhm output>\n", argv[0]);
		exit(1);
	}
	while(i < argc){
  		if( argv[i][0] == '-'){
    			i++;
    			switch(argv[i-1][1]){
      				case 'm' : {model_file=argv[i]; break;}
   			}   
   		}
    		i++;
  	}
	hhmake_file = argv[i-1];

	// Parse hhmake output
	FILE *fin = fopen (hhmake_file.c_str(), "r");
	if (fin != NULL){
		//cout << "#pos\tA\tC\tD\tE\tF\tG\tH\tI\tK\tL\tM\tN\tP\tQ\tR\tS\tT\tV\tW\tY" << endl;    
		while (fgets (buf, 1024, fin)) {
			string line = buf;
			if(line.compare(0,40,"Amino acid frequencies WITH pseudocounts") == 0){	
				getfreqs = 1;
			}
			if(getfreqs&& line.compare(0,14,"Writing HMM to") == 0){
				getfreqs = 0;
			}
			if (sscanf(buf, "NULL%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f",&nullfreqs[0],&nullfreqs[1],&nullfreqs[2],&nullfreqs[3],&nullfreqs[4],&nullfreqs[5],&nullfreqs[6],&nullfreqs[7],&nullfreqs[8],&nullfreqs[9],&nullfreqs[10],&nullfreqs[11],&nullfreqs[12],&nullfreqs[13],&nullfreqs[14],&nullfreqs[15],&nullfreqs[16],&nullfreqs[17],&nullfreqs[18],&nullfreqs[19]) == 20){
				for(int i = 0;i < 20;i++){
					// Inverse log, *100
					nullfreqs[i] = 100*pow(2,(float)nullfreqs[i]/-1000);
					found_null_freqs = 1;
					//cout << nullfreqs[i] << endl;
				}
			}
			if (getfreqs && 20 == sscanf(buf, "%*d%*c%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*c[4]%*f",&freqs[pos][0],&freqs[pos][1],&freqs[pos][2],&freqs[pos][3],&freqs[pos][4],&freqs[pos][5],&freqs[pos][6],&freqs[pos][7],&freqs[pos][8],&freqs[pos][9],&freqs[pos][10],&freqs[pos][11],&freqs[pos][12],&freqs[pos][13],&freqs[pos][14],&freqs[pos][15],&freqs[pos][16],&freqs[pos][17],&freqs[pos][18],&freqs[pos][19])){
				pos++;
			}
		
		}
		fclose (fin);
	}else{
		cout << "Couldn't open file " << hhmake_file << endl;
		exit(1);
	}
	if(!found_null_freqs || !pos){
		cout << "Failed to parse hhmake output" << endl;
		exit(1);
	}

	// Load libsvm model file
    	struct svm_model* model;
    	if((model = svm_load_model(model_file.c_str())) == 0){
		cout << "Couldn't load model file " << model_file << endl << "Try passing it via -m flag" << endl;
        	exit(1);
	}

	// Libsvm containers
	int max_nr_attr = 60;
	int nr_class=svm_get_nr_class(model);
	struct svm_node *feature_vector = (struct svm_node *) malloc(max_nr_attr*sizeof(struct svm_node));
	double *prob_estimates = (double *) malloc(nr_class*sizeof(double));

	// Generate sliding window
	for(int i = 0;i < pos;i++){		
		int feature = 1;
		for(int w = (int)-window/2;w < 1+(int)window/2;w++){
			for(int j = 0;j < 20;j++){
				if(i+w >= 0 && i+w < pos){
					// Log of frequency ratio
					double p = log(freqs[i][j]/nullfreqs[j]);
					// Sigmoidal transformation
					p = 1/(1+exp(-p));
					// Fill feature vector
					feature_vector[feature-1].index = feature;
					feature_vector[feature-1].value = p;
				}
				feature++;
			}		
		}
		// Mark end of feature vector with index -1
		feature_vector[feature-1].index = -1;
		svm_predict_probability(model,feature_vector,prob_estimates);	
		cout << setprecision(precision) << prob_estimates[0] << endl;
	}
	svm_free_and_destroy_model(&model);
	return(1);
}
