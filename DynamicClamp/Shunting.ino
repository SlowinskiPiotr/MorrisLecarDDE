// A constant conductance with a reversal potential of -80 mV
float Shunting(float v) {
  float current = -gLeak * (v + 80);
  return current;
}
