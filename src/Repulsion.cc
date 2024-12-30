#include "Repulsion.hh"

namespace RepulsionStage {

G4FragmentVector CalculateRepulsion(G4FragmentVector frags, aamcc::NucleonVector nucleons, const std::vector<int>& maps) {
  if (nucleons.empty()) {
    return frags;
  }

  for (auto& n : nucleons) {
    n.x *= fm;
    n.y *= fm;
    n.z *= fm;
  }

  double time = 0.0;
  double delta_time = totalTime / static_cast<double>(iterations);

  while (time < totalTime) {
    BHTree bhtree(&nucleons, &frags, &maps);

    double temp_timedelta = std::min(delta_time, bhtree.GetAdaptiveTimeDelta());

    auto r_delta = bhtree.Iterate(temp_timedelta);

    for (int i = 0; i < nucleons.size(); i++) {
      nucleons[i].x += r_delta[maps[i]].x();
      nucleons[i].y += r_delta[maps[i]].y();
      nucleons[i].z += r_delta[maps[i]].z();
    }
    time += temp_timedelta;
  }

  for (auto& n : nucleons) {
    n.x /= fm;
    n.y /= fm;
    n.z /= fm;
  }
  return frags;
}

BHTree::BHTree(const aamcc::NucleonVector* nucleons, G4FragmentVector* frags, const std::vector<int>* maps) : frags_(frags), maps_(maps) {
  fs_.assign(frags->size(), {0.0, 0.0, 0.0});
  rootnode_ = BuildBHTree(nucleons);
  GetForces(rootnode_);
}

std::shared_ptr<BHNode> BHTree::InitializeRoot(const aamcc::NucleonVector* nucleons) {
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

  G4ThreeVector cr = {0.0, 0.0, 0.0};
  size_t nucleons_sz = nucleons->size();
  for (const auto& nuc : *nucleons) {
    cr.setX(cr.x() + nuc.GetX());
    cr.setY(cr.y() + nuc.GetY());
    cr.setZ(cr.z() + nuc.GetZ());
  }
  cr.setX(cr.x() / static_cast<double>(nucleons_sz));
  cr.setY(cr.y() / static_cast<double>(nucleons_sz));
  cr.setZ(cr.z() / static_cast<double>(nucleons_sz));

  double maxRange = std::max({maxX - minX, maxY - minY, maxZ - minZ});
  return std::make_shared<BHNode>(maxRange, G4ThreeVector((minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2));
}

std::shared_ptr<BHNode> BHTree::BuildBHTree(const aamcc::NucleonVector* nucleons) {
  rootnode_ = InitializeRoot(nucleons);

  for (size_t i = 0; i < nucleons->size(); i++) {
    InsertNucleon(rootnode_, {(*nucleons)[i].GetX(), (*nucleons)[i].GetY(), (*nucleons)[i].GetZ()}, i);
  }
  return rootnode_;
}

void BHTree::InsertNucleon(std::shared_ptr<BHNode>& node, const G4ThreeVector cords, int pIndex) {
  if (node->totalA == 0) {
    node->totalA = 1;
    node->cr = cords;
    node->index = pIndex;
    return;
  }
  if (node->totalA == 1) {
    node->Divide();

    int index = 0;
    if (node->cr.x() > node->ctr.x()) index |= 1;
    if (node->cr.y() > node->ctr.y()) index |= 2;
    if (node->cr.z() > node->ctr.z()) index |= 4;
    std::shared_ptr<BHNode> targetNode = node->children[index];
    InsertNucleon(targetNode, node->cr, node->index);

    index = 0;
    if (cords.x() > node->ctr.x()) index |= 1;
    if (cords.y() > node->ctr.y()) index |= 2;
    if (cords.z() > node->ctr.z()) index |= 4;
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
  if (cords.x() > node->ctr.x()) index |= 1;
  if (cords.y() > node->ctr.y()) index |= 2;
  if (cords.z() > node->ctr.z()) index |= 4;
  std::shared_ptr<BHNode> targetNode = node->children[index];
  InsertNucleon(targetNode, cords, pIndex);
  return;
}

double BHTree::GetAdaptiveTimeDelta() const {
  double min_time = max_adaptive_delta;
  for (size_t i = 0; i < frags_->size(); i++) {
    if (!frags_->at(i)->GetZ_asInt() || !fs_[i].mag() || !frags_->at(i)->GetMomentum().vect().mag()) {
      continue;
    }
    min_time = std::min(min_time, 0.05 * frags_->at(i)->GetMomentum().vect().mag() / fs_[i].mag());
  }
  return min_time;
}

std::vector<G4ThreeVector> BHTree::Iterate(double time_delta) {
  std::vector<G4ThreeVector> r_delta(frags_->size());

  for (size_t i = 0; i < frags_->size(); ++i) {
    if (frags_->at(i)->GetZ_asInt() == 0) {
      continue;
    }
    G4LorentzVector p = frags_->at(i)->GetMomentum();

    G4ThreeVector half_dp = fs_[i] * time_delta * 0.5;

    p = G4LorentzVector((p.vect() + half_dp), std::sqrt((p.vect() + half_dp).mag2() + p.m2()));
    G4ThreeVector mid_v = p.vect() / p.e();
    r_delta[i] = G4ThreeVector(time_delta * mid_v.x(), time_delta * mid_v.y(), time_delta * mid_v.z());

    p = G4LorentzVector((p.vect() + half_dp), std::sqrt((p.vect() + half_dp).mag2() + p.m2()));

    frags_->at(i)->SetMomentum(p);
  }
  return r_delta;
}

void BHTree::GetForces(const std::shared_ptr<BHNode>& node) {
  if (node->totalA == 0) {
    return;
  }
  if (node->totalA == 1) {
    fs_[maps_->at(node->index)] += Force(rootnode_, node);
    return;
  }
  for (const auto& child : node->children) {
    GetForces(child);
  }
  return;
}

G4ThreeVector BHTree::Force(const std::shared_ptr<BHNode>& rootnode, const std::shared_ptr<BHNode>& node) const {
  if ((rootnode->totalA == 0) || ((rootnode->index != -1) && maps_->at(rootnode->index) == maps_->at(node->index))) {
    return {0.0, 0.0, 0.0};
  }
  if (rootnode->totalA == 1) {
    return DuoForce(rootnode->GetCr(), node->GetCr(), rootnode->GetTotalA());
  }
  if ((rootnode->size / (node->GetCr() - rootnode->GetCr()).mag()) < theta) {
    return DuoForce(rootnode->GetCr(), node->GetCr(), rootnode->GetTotalA());
  }
  G4ThreeVector totalForce = {0.0, 0.0, 0.0};
  for (const auto& child : rootnode->children) {
    totalForce += Force(child, node);
  }
  return totalForce;
}

G4ThreeVector BHTree::DuoForce(G4ThreeVector vfrom, G4ThreeVector target, double from_totalA) const {
  G4ThreeVector vec = target - vfrom;
  G4ThreeVector fos = vec * CLHEP::elm_coupling * from_totalA / std::pow(vec.mag(), 3);
  return fos;
}

void BHNode::Divide() {
  for (size_t i = 0; i < 8; i++) {
    G4ThreeVector offset(
      (i & 1 ? 1 : -1) * size / 4.0,
      (i & 2 ? 1 : -1) * size / 4.0,
      (i & 4 ? 1 : -1) * size / 4.0
    );
    children.push_back(std::make_shared<BHNode>(size / 2.0, ctr + offset));
  }
}

}