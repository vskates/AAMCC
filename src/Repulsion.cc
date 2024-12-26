#include "Repulsion.hh"

namespace RepulsionStage {

const double fm = 1e-15;
const double MeV = 1.602176634e-13;
const double theta = 0.3;
const double MeVc = 1.602176634e-13 / 2.99792458e8;
const double k_const = 8.9875517873681764e9;
const double e_charge = 1.602176634 * 1e-19;
const double n_mass = 1.67492749804 * 1e-27;
const double p_mass = 1.67262192369e-27;
const double speed_light = 2.99792458e8;

const double totalTime = 200 * fm / speed_light;
const double iterations = 100;

std::shared_ptr<BHNode> InitializeRoot(const aamcc::NucleonVector* nucleons) {
  double minX = (*nucleons)[0].GetX(), maxX = minX;
  double minY = (*nucleons)[0].GetY(), maxY = minY;
  double minZ = (*nucleons)[0].GetZ(), maxZ = minZ;

  for (const auto& n : *nucleons) {
    minX = std::min(minX, n.GetX());
    maxX = std::max(maxX, n.GetX());
    minY = std::min(minY, n.GetY());
    maxY = std::max(maxY, n.GetY());
    minZ = std::min(minZ, n.GetZ());
    maxZ = std::max(maxZ, n.GetZ());
  }

  Vector3d cr = {0.0, 0.0, 0.0};
  size_t nucleons_sz = nucleons->size();
  for (const auto& nuc : *nucleons) {
    cr.x += nuc.GetX();
    cr.y += nuc.GetY();
    cr.z += nuc.GetZ();
  }
  cr.x /= static_cast<double>(nucleons_sz);
  cr.y /= static_cast<double>(nucleons_sz);
  cr.z /= static_cast<double>(nucleons_sz);

  double maxRange = std::max({maxX - minX, maxY - minY, maxZ - minZ});
  auto rootnode = std::make_shared<BHNode>(maxRange, Vector3d((minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2));
  return rootnode;
}

std::shared_ptr<BHNode> BuildBhtree(const aamcc::NucleonVector* nucleons) {
  if (nucleons->empty()) {
    return nullptr;
  }

  std::shared_ptr<BHNode> rootnode = InitializeRoot(nucleons);

  for (size_t i = 0; i < nucleons->size(); i++) {
    InsertNucleon(rootnode, {(*nucleons)[i].GetX(), (*nucleons)[i].GetY(), (*nucleons)[i].GetZ()}, i);
  }
  return rootnode;
}

void BHNode::Divide() {
  for (size_t i = 0; i < 8; i++) {
    Vector3d offset(
      (i & 1 ? 1 : -1) * size / 4.0,
      (i & 2 ? 1 : -1) * size / 4.0,
      (i & 4 ? 1 : -1) * size / 4.0
    );
    children.push_back(std::make_shared<BHNode>(size / 2.0, ctr + offset));
  }
}

void InsertNucleon(std::shared_ptr<BHNode>& node, const Vector3d cords, int pIndex) {
  if (node->totalA == 0) {
    node->totalA = 1;
    node->cr = cords;
    node->index = pIndex;
    return;
  }
  if (node->totalA == 1) {
    node->Divide();

    int index = 0;
    if (node->cr.x > node->ctr.x) index |= 1;
    if (node->cr.y > node->ctr.y) index |= 2;
    if (node->cr.z > node->ctr.z) index |= 4;
    std::shared_ptr<BHNode> targetNode = node->children[index];
    InsertNucleon(targetNode, node->cr, node->index);

    index = 0;
    if (cords.x > node->ctr.x) index |= 1;
    if (cords.y > node->ctr.y) index |= 2;
    if (cords.z > node->ctr.z) index |= 4;
    targetNode = node->children[index];
    InsertNucleon(targetNode, cords, pIndex);

    node->cr = (node->cr + cords) / 2;
    node->totalA += 1;
    node->index = -1;
    return;
  }
  node->cr = (node->cr * node->totalA + cords) / (node->totalA + 1);
  node->totalA += 1;

  int index = 0;
  if (cords.x > node->ctr.x) index |= 1;
  if (cords.y > node->ctr.y) index |= 2;
  if (cords.z > node->ctr.z) index |= 4;
  std::shared_ptr<BHNode> targetNode = node->children[index];
  InsertNucleon(targetNode, cords, pIndex);
  return;
}

inline double GetGamma(const G4LorentzVector& momentum, double mass) {
  double p2 = momentum.x() * momentum.x() + momentum.y() * momentum.y() + momentum.z() * momentum.z();
  return std::sqrt(1.0 + p2 / (mass * mass * speed_light * speed_light));
}

double GetAdaptiveTimeDelta(const G4FragmentVector& frags, const std::vector<Vector3d>& fs) {
  double min_time = 100;
  for (size_t i = 0; i < frags.size(); i++) {
    min_time = std::min(min_time, 0.05 * frags[i]->GetMomentum().vect().mag() / fs[i].magnitude());
  }
  return min_time;
}

G4FragmentVector CalculateRepulsion(G4FragmentVector frags, aamcc::NucleonVector nucleons, const std::vector<int>& maps) {
  for (auto& n : nucleons) {
    n.x *= fm;
    n.y *= fm;
    n.z *= fm;
  }

  std::vector<double> frag_masses = GetFragmentMasses(frags);

  double time = 0.0;
  double delta_time = totalTime / iterations;

  while (time < totalTime) {
    std::shared_ptr<BHNode> rootnode = BuildBhtree(&nucleons);

    if (rootnode == nullptr) {
      return frags;
    }

    std::vector<Vector3d> fs(frags.size(), {0.0, 0.0, 0.0});

    GetForces(rootnode, rootnode, maps, fs);

    double temp_timedelta = delta_time;
    auto adaptive_tm = GetAdaptiveTimeDelta(frags, fs);
    if ((adaptive_tm < temp_timedelta) && (totalTime / adaptive_tm < (2 * iterations))) {
      temp_timedelta = adaptive_tm;
    }

    auto delta_x = Iterate(frags, fs, temp_timedelta, frag_masses);

    for (int i = 0; i < nucleons.size(); i++) {
      nucleons[i].x += delta_x[maps[i]].x;
      nucleons[i].y += delta_x[maps[i]].y;
      nucleons[i].z += delta_x[maps[i]].z;
    }
    time += temp_timedelta;
  }

  for (auto& frag : frags) {
    frag->SetMomentum(frag->GetMomentum() / MeVc);
  }
  return frags;
}

std::vector<Vector3d> Iterate(G4FragmentVector& frags, const std::vector<Vector3d>& fs, double time_delta, const std::vector<double>& frag_masses) {
  std::vector<Vector3d> delta_x(frags.size());

  for (size_t i = 0; i < frags.size(); ++i) {
    if (frags[i]->GetZ_asInt() == 0) {
      continue;
    }

    G4LorentzVector currentMomentum = frags[i]->GetMomentum();
    G4ThreeVector spatialMomentum = currentMomentum.vect();

    G4ThreeVector halfstepMomentumIncrement(
        0.5 * time_delta * (fs[i].x),
        0.5 * time_delta * (fs[i].y),
        0.5 * time_delta * (fs[i].z)
    );

    G4ThreeVector newMomentum = spatialMomentum + halfstepMomentumIncrement;
    double gamma = GetGamma(currentMomentum, frag_masses[i]);
    G4ThreeVector velocity = newMomentum / (frag_masses[i] * gamma);

    delta_x[i] = Vector3d(time_delta * velocity.x(),
                            time_delta * velocity.y(),
                            time_delta * velocity.z());

    newMomentum = newMomentum + halfstepMomentumIncrement;
    G4LorentzVector finalMomentum(std::sqrt(newMomentum.mag() * newMomentum.mag() + std::pow(frag_masses[i] * speed_light, 2)), newMomentum);

    frags[i]->SetMomentum(finalMomentum);
  }
  return delta_x;
}

void LazyGetForces(const aamcc::NucleonVector* nucleons, const std::vector<int>& maps, std::vector<Vector3d>& fs) {
  for (size_t i = 0; i < nucleons->size(); i++) {
    for (size_t j = i + 1; j < nucleons->size(); j++) {
      auto nuc = (*nucleons)[i];
      auto nuc2 = (*nucleons)[j];
      if (nuc.isospin != 1 || nuc2.isospin != 1) {
        continue;
      }
      Vector3d r2 = {nuc2.GetX(), nuc2.GetY(), nuc2.GetZ()};
      Vector3d r = {nuc.GetX(), nuc.GetY(), nuc.GetZ()};
      Vector3d tempf = DuoForce(r2, r, 1);
      fs[maps[i]] += tempf;
      fs[maps[j]] -= tempf;
    }
  }
}

std::vector<double> GetFragmentMasses(G4FragmentVector& fragments) {
  std::vector<double> ms;
  for (auto& frag : fragments) {
    double mass = (frag->GetZ_asInt() * p_mass) + (frag->GetA_asInt() - frag->GetZ_asInt()) * n_mass;
    frag->SetMomentum(MeVc * frag->GetMomentum());
    ms.push_back(mass);
  }
  return ms;
}

void GetForces(const std::shared_ptr<BHNode>& rootnode, const std::shared_ptr<BHNode>& node, const std::vector<int>& maps, std::vector<Vector3d>& fs) {
  if (node->totalA == 0) {
    return;
  }
  if (node->totalA == 1) {
    fs[maps[node->index]] += Force(rootnode, node, maps);
    return;
  }
  for (const auto& child : node->children) {
    GetForces(rootnode, child, maps, fs);
  }
  return;
}

Vector3d Force(const std::shared_ptr<BHNode>& rootnode, const std::shared_ptr<BHNode>& node, const std::vector<int>& maps) {
  if ((rootnode->totalA == 0) || ((rootnode->index != -1) && maps[rootnode->index] == maps[node->index])) {
    return {0.0, 0.0, 0.0};
  }
  if (rootnode->totalA == 1) {
    return DuoForce(rootnode->GetCr(), node->GetCr(), rootnode->GetTotalA());
  }
  if ((rootnode->size / (node->GetCr() - rootnode->GetCr()).magnitude()) < theta) {
    return DuoForce(rootnode->GetCr(), node->GetCr(), rootnode->GetTotalA());
  }
  Vector3d totalForce;
  for (const auto& child : rootnode->children) {
    totalForce += Force(child, node, maps);
  }
  return totalForce;
}

Vector3d DuoForce(Vector3d vfrom, Vector3d target, double from_totalA) {
  Vector3d vec = target - vfrom;
  Vector3d fos = vec * k_const * e_charge * (e_charge * from_totalA) / std::pow(vec.magnitude(), 3);
  return fos;
}


Vector3d::Vector3d(const Vector3d &vec) {
  x = vec.x;
  y = vec.y;
  z = vec.z;
}

Vector3d Vector3d::operator+(const Vector3d &vec) { return Vector3d(x + vec.x, y + vec.y, z + vec.z); }

Vector3d &Vector3d::operator+=(const Vector3d &vec) {
  x += vec.x;
  y += vec.y;
  z += vec.z;
  return *this;
}

Vector3d Vector3d::operator-(const Vector3d &vec) { return Vector3d(x - vec.x, y - vec.y, z - vec.z); }

Vector3d &Vector3d::operator-=(const Vector3d &vec) {
  x -= vec.x;
  y -= vec.y;
  z -= vec.z;
  return *this;
}

Vector3d Vector3d::operator*(double value) { return Vector3d(x * value, y * value, z * value); }

Vector3d &Vector3d::operator*=(double value) {
  x *= value;
  y *= value;
  z *= value;
  return *this;
}

Vector3d Vector3d::operator/(double value) {
  assert(value != 0);
  return Vector3d(x / value, y / value, z / value);
}

Vector3d &Vector3d::operator/=(double value) {
  assert(value != 0);
  x /= value;
  y /= value;
  z /= value;
  return *this;
}

Vector3d &Vector3d::operator=(const Vector3d &vec) {
  x = vec.x;
  y = vec.y;
  z = vec.z;
  return *this;
}

double Vector3d::dot_product(const Vector3d &vec) { return x * vec.x + vec.y * y + vec.z * z; }

Vector3d Vector3d ::cross_product(const Vector3d &vec) {
  double ns = y * vec.z - z * vec.y;
  double nj = z * vec.x - x * vec.z;
  double nk = x * vec.y - y * vec.x;
  return Vector3d(ns, nj, nk);
}

const double Vector3d::magnitude() const { return sqrt(square()); }

double Vector3d::square() const { return x * x + y * y + z * z; }

Vector3d Vector3d::normalization() {
  assert(magnitude() != 0);
  *this /= magnitude();
  return *this;
}

double Vector3d::distance(const Vector3d &vec) {
  Vector3d dist = *this - vec;
  return dist.magnitude();
}

double Vector3d::show_X() { return x; }

double Vector3d::show_Y() { return y; }

double Vector3d::show_Z() { return z; }

void Vector3d::disp() { std::cout << x << " " << y << " " << z << std::endl; }

}