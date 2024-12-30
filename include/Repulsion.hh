#ifndef REPULSION_HH
#define REPULSION_HH

#include <cmath>
#include <limits>
#include <iostream>
#include <vector>
#include <memory>

#include "G4Fragment.hh"
#include "G4SystemOfUnits.hh"
#include "Nucleon.hh"

namespace RepulsionStage {

constexpr double fm = 1e-15 * CLHEP::m;
constexpr double theta = 0.3;

constexpr double totalTime = 200 * fm / CLHEP::c_light;
constexpr double iterations = 100;
constexpr double max_adaptive_delta = std::numeric_limits<double>::max();

G4FragmentVector CalculateRepulsion(G4FragmentVector frags, aamcc::NucleonVector nucleons, const std::vector<int>& maps);

class BHNode {
 public:
  int totalA; // total proton count
  G4ThreeVector cr; // mean coordinates of the charges in box
  G4ThreeVector ctr; // coordinates of the box center
  std::vector<std::shared_ptr<BHNode>> children; // child nodes
  int index; // -1 if > 1 particles, index in nucleons vector otherwise
  double size; // size of the box

  BHNode() = default;
  BHNode(double size, G4ThreeVector ctr) : size(size), ctr(ctr), totalA(0), cr({0.0, 0.0, 0.0}), index(-1) {};
  ~BHNode() = default;
  void Divide();
  inline G4ThreeVector GetCr() const { return cr;}
  inline double GetTotalA() const { return totalA;}
};

class BHTree {
 public:
  explicit BHTree(const aamcc::NucleonVector* nucleons, G4FragmentVector* frags, const std::vector<int>* maps);

  std::vector<G4ThreeVector> Iterate(double time_delta);

  double GetAdaptiveTimeDelta() const;

 private:
  std::shared_ptr<BHNode> rootnode_;
  G4FragmentVector* frags_;
  const std::vector<int>* maps_;
  std::vector<G4ThreeVector> fs_;

  std::shared_ptr<BHNode> BuildBHTree(const aamcc::NucleonVector* nucleons);

  std::shared_ptr<BHNode> InitializeRoot(const aamcc::NucleonVector* nucleons);

  void GetForces(const std::shared_ptr<BHNode>& node);

  G4ThreeVector Force(const std::shared_ptr<BHNode>& rootnode, const std::shared_ptr<BHNode>& node) const;

  G4ThreeVector DuoForce(G4ThreeVector vfrom, G4ThreeVector target, double from_totalA) const;

  void InsertNucleon(std::shared_ptr<BHNode>& node, const G4ThreeVector cords, int pIndex);
};

}

#endif