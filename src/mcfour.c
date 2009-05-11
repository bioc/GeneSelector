void mcfour(int *R, int *U, int *lx, int *ll, int *mm, int *maxrank, double *L) {
 int i, j, m, temp1, temp2;
 for(m = 0; m < *mm; m++) 
     for(i=0; i < *ll; i++){
	 for(j=(i+1); j < *ll; j++){
             temp1 = R[U[i] + (*lx)* m - 1];
             temp2 = R[U[j] + (*lx)* m - 1];
	     if(temp1 <= *maxrank && temp2 <= *maxrank){  
		 if(temp2 < temp1) L[i + (*ll)*j] = L[i + (*ll)*j] + 1;
                 else L[j + (*ll)*i] = L[j + (*ll)*i] + 1;
	     } 
	 }
     }
 for(i=0; i < *ll; i++){
	 for(j=(i+1); j < *ll; j++){
	     if(L[i + (*ll)*j] > L[j + (*ll)*i]){ 
		L[j + (*ll)*i] = 0; 
	        L[i + (*ll)*j] = 1;
	     }
    else{
	if(L[i + (*ll)*j] == L[j + (*ll)*i]){  
	L[i + (*ll)*j]    = 0.5;
	L[j + (*ll)*i]    = 0.5;
       }
    
	else{ 
          L[j + (*ll)*i] = 1; 
	  L[i + (*ll)*j] = 0;
        } 
    }
  }
  } 
} 
