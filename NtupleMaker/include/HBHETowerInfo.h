struct HBHETowerInfo {
  float energy;
  float energyEM;
  float energyHAD;
  int nCaloHits;
  float energy_25ns;
  float energyEM_25ns;
  float energyHAD_25ns;
  int nCaloHits_25ns;

  HBHETowerInfo(){
    energy=0.f; energyHAD=0.f; energyEM=0.f; nCaloHits=0;
    energy_25ns=0.f; energyHAD_25ns=0.f; energyEM_25ns=0.f; nCaloHits_25ns=0;
  }
};

