#ifndef REPULSION_HH
#define REPULSION_HH

#include <cmath>
#include <limits>
#include <iomanip>
#include <iostream>
#include <array>
#include <memory>
#include <vector>
#include <cassert>
#include <memory>

#include "G4Fragment.hh"
#include "G4SystemOfUnits.hh"
#include "Nucleon.hh"

namespace RepulsionStage {

G4FragmentVector CalculateRepulsion(G4FragmentVector frags, aamcc::NucleonVector nucleons, const std::vector<int>& maps);

class Vector3d {
 public:
  double x;
  double y;
  double z;
  Vector3d() {
    x = 0;
    y = 0;
    z = 0;
  }

  Vector3d(double x1, double y1, double z1 = 0) {
    x = x1;
    y = y1;
    z = z1;
  }
  Vector3d(const Vector3d &vec);
  Vector3d operator+(const Vector3d &vec);
  Vector3d &operator+=(const Vector3d &vec);
  Vector3d operator-(const Vector3d &vec);
  Vector3d &operator-=(const Vector3d &vec);
  Vector3d operator*(double value);
  Vector3d &operator*=(double value);
  Vector3d operator/(double value);
  Vector3d &operator/=(double value);
  Vector3d &operator=(const Vector3d &vec);
  double dot_product(const Vector3d &vec);
  Vector3d cross_product(const Vector3d &vec);
  const double magnitude() const;
  Vector3d normalization();
  double square() const;

  double distance(const Vector3d &vec);
  double show_X();
  double show_Y();
  double show_Z();
  void disp();
};

class BHNode {
 public:
  int totalA; // total proton count
  Vector3d cr; // mean coordinates of the charges in box
  Vector3d ctr; // coordinates of the box center
  std::vector<std::shared_ptr<BHNode>> children; // child nodes
  int index; // -1 if > 1 particles, index in nucleons vector otherwise
  double size; // size of the box

 public:
  BHNode() = default;
  BHNode(double size, Vector3d ctr) : size(size), ctr(ctr), totalA(0), cr({0.0, 0.0, 0.0}), index(-1) {};
  ~BHNode() = default;
  void Divide();
  inline Vector3d GetCr() const { return cr;}
  inline double GetTotalA() const { return totalA;}
};

std::shared_ptr<BHNode> BuildBhtree(const aamcc::NucleonVector* nucleons);

std::shared_ptr<BHNode> InitializeRoot(const aamcc::NucleonVector* nucleons);

Vector3d Force(const std::shared_ptr<BHNode>& rootnode, const std::shared_ptr<BHNode>& node, const std::vector<int>& maps);

Vector3d DuoForce(Vector3d vfrom, Vector3d target, double from_totalA);

std::vector<double> GetFragmentMasses(G4FragmentVector& fragments);

void GetForces(const std::shared_ptr<BHNode>& rootnode, const std::shared_ptr<BHNode>& node, const std::vector<int>& maps, std::vector<Vector3d>& fs);

void LazyGetForces(const aamcc::NucleonVector* nucleons, const std::vector<int>& maps, std::vector<Vector3d>& fs);

std::vector<Vector3d> Iterate(G4FragmentVector& frags, const std::vector<Vector3d>& fs, double time_delta, const std::vector<double>& frag_masses);

inline double GetGamma(const G4LorentzVector& momentum, double mass);

void InsertNucleon(std::shared_ptr<BHNode>& node, const Vector3d cords, int pIndex);

double GetAdaptiveTimeDelta(const G4FragmentVector& frags, const std::vector<Vector3d>& fs);

}

#endif