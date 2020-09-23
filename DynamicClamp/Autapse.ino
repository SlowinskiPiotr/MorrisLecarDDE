// A self-delayed simple, state-dependent autapse model
// For speed, the function is pre-computed
// and put in lookup tables stored as global variables.

// Declare the lookup table variables
float autapse[1501] = {0.0};                     // Pre-calculate autapse function
int V_hist[buffer_size] = {350};                 // buffer for delayed voltage


// Generate the lookup tables
void GenerateAutapseLUT() {
  float v;
  for (int x=0; x<1501; x++) {
    v = (float)x/10 - 100.0;                     // We use membrane potentials between -100 mV and +50 mV
  autapse[x] = 0.5*(1.0+tanh( (v+20.0)/0.5));
  }
}

// At every time step, calculate delay and autapse value
float AutapseCurrent(float v) {

  while (delay_sum > tau) {
    ind_delay++;
    if (ind_delay >= buffer_size) {
      ind_delay = 0;
    }
    delay_sum -= t_hist[ind_delay];
  }

  int vIdx = (int)v*10.0 + 1000;
  vIdx = constrain(vIdx,0,1500);
  V_hist[ind_curr] = vIdx;
  t_hist[ind_curr] = dt;

  ind_curr++;
  if (ind_curr >= buffer_size) {
    ind_curr = 0;
  }

  delay_sum += dt;

  float current = kappa*autapse[V_hist[ind_delay]];
  return current;
}
