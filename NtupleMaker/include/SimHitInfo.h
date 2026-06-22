struct SimHitInfo {
  float energy;
  float energyEM;
  float energyHAD;
  int nCaloHits;

  float energy_25ns;
  float energyEM_25ns;
  float energyHAD_25ns;
  int nCaloHits_25ns;

  float energy_50ns;
  float energyEM_50ns;
  float energyHAD_50ns;
  int nCaloHits_50ns;

  float energy_75ns;
  float energyEM_75ns;
  float energyHAD_75ns;
  int nCaloHits_75ns;

  float energy_100ns;
  float energyEM_100ns;
  float energyHAD_100ns;
  int nCaloHits_100ns;

  float energy_200ns;
  float energyEM_200ns;
  float energyHAD_200ns;
  int nCaloHits_200ns;

  SimHitInfo(){
    energy=0; energyHAD=0; energyEM=0; nCaloHits=0;
    energy_25ns=0; energyHAD_25ns=0; energyEM_25ns=0; nCaloHits_25ns=0;
    energy_50ns=0; energyHAD_50ns=0; energyEM_50ns=0; nCaloHits_50ns=0;
    energy_75ns=0; energyHAD_75ns=0; energyEM_75ns=0; nCaloHits_75ns=0;
    energy_100ns=0; energyHAD_100ns=0; energyEM_100ns=0; nCaloHits_100ns=0;
    energy_200ns=0; energyHAD_200ns=0; energyEM_200ns=0; nCaloHits_200ns=0;
  }
};

