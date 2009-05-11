void kendall(int *r, int *rprime, int *lx, double *w, double *res) {
 int i, j;
 double tau;
 for(i=0; i < *lx; i++)
     for(j=(i+1); j < *lx; j++)
         tau = tau + w[i] * w[j] * (((r[i] - r[j]) * (rprime[i] - rprime[j])) < 0); 
  res[0] = tau;
}



