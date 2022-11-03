void get_amp_and_frq(const float *floatbuf,float *amp,float *freq,int clength);
void put_amp_and_frq(float *floatbuf,const float *amp, const float *freq,int clength);
void do_spectral_shiftp(float *amp, float *freq,double pitch,int clength);
