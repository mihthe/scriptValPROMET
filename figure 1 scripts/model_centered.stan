data{
    array[5000] int m_surg;
     vector[5000] Su;
    array[5000] int y;
    array[5000] int n;
     vector[5000] m_bio;
     vector[5000] Si;
    array[5000] int L;
}
parameters{
     vector[5] b;
     vector[5] e;
     real a;
     real e_bar;
     real d;
     real g;
     real b_bar;
     real<lower=0> sigmae;
     real<lower=0> sigmab;
}
model{
     vector[5000] lambda;
    sigmab ~ exponential( 2 );
    sigmae ~ exponential( 2 );
    b_bar ~ normal( -1 , 0.2 );
    g ~ normal( -1 , 0.2 );
    d ~ normal( -1 , 0.2 );
    e_bar ~ normal( -1 , 0.2 );
    a ~ normal( -1 , 0.2 );
    e ~ normal( e_bar , sigmae );
    b ~ normal( b_bar , sigmab );
    for ( i in 1:5000 ) {
        lambda[i] = a + b[L[i]] * Si[i] + g * m_bio[i] * n[i] + d * m_bio[i] * y[i] + e[L[i]] * Su[i];
        lambda[i] = exp(lambda[i]);
    }
    m_surg ~ poisson( lambda );
}
