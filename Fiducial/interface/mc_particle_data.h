#ifndef MC_PARTICLE_DATA_H
#define MC_PARTICLE_DATA_H

#include "TLorentzVector.h"
#include <iostream>

class MCParticleData {

 private:
  TLorentzVector four_vector_;
  int pid_;
  int mom_pid_;
  int status_;
  int mc_parentage_;
    
 public:
  void SetFourVector(TLorentzVector four_vector) {four_vector_ = four_vector;};
  void SetPID(int pid) {pid_ = pid;};
  void SetMomPID(int mom_pid){ mom_pid_ = mom_pid;};
  void SetStatus(int status) {status_ = status; };
  void SetMCParentage(int mc_parentage){mc_parentage_ = mc_parentage; };

  TLorentzVector GetFourVector() { return four_vector_; };
  int GetPID() {return pid_;} ;
  int GetMomPID() {return mom_pid_; } ;
  int GetStatus() {return status_; } ;
  int GetMCParentage() { return mc_parentage_; };
};
  
#endif
