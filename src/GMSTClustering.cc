#include "GMSTClustering.hh"
#include "Repulsion.hh"

GMSTCluster::GMSTCluster(G4int Z_in, G4int A_in) : A(A_in), Z(Z_in){};

GMSTCluster::~GMSTCluster(){};

GMSTClustering::GMSTClustering() { CritDist = 100; };

GMSTClustering::GMSTClustering(G4double CD_in, G4double Aa_in, G4double Ab_in) {
  CritDist = CD_in;
  SpecAa = Aa_in;
  SpecAb = Ab_in;
};

GMSTClustering::~GMSTClustering(){};

Graph GMSTClustering::ClusterToGraph(aamcc::NucleonVector *nucleons, G4double A) {
  Graph g(A, A * (A - 1) / 2);
  //  making full graph of nucleons
  for (G4int iArray = 0; iArray < nucleons->size(); iArray++) {
    aamcc::Nucleon *nucleon = &(nucleons->at(iArray));
    for (G4int iArray_pairs = iArray + 1; iArray_pairs < nucleons->size(); iArray_pairs++) {
      aamcc::Nucleon *nucleon_pair = &(nucleons->at(iArray_pairs));
      g.addEdge(iArray, iArray_pairs, std::sqrt(pow(nucleon->GetX() - nucleon_pair->GetX(), 2) + pow(nucleon->GetY() - nucleon_pair->GetY(), 2) + pow(nucleon->GetZ() - nucleon_pair->GetZ(), 2)));
    }
  }
  return g;
};

void GMSTClustering::SetCDExEn(G4double Ex, G4int A) {
  if (Ex / G4double(A) < 2.17 * MeV) {
    CritDist = d0;
  } else {
    G4double ex = Ex / G4double(A);
    G4double dep = std::exp(-std::pow((ex / b_opt), a_opt)) * c_opt + d_opt;
    CritDist = d0 * std::pow(dep, 1. / 3.);
  }
  if (alphaPow > 10) {
    CritDist = d0;
  }
}

std::vector<G4FragmentVector> GMSTClustering::GetClusters(aamcc::NucleonVector *nucleons_in, G4double ExA, G4double ExB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB) {
  G4int A = 0;
  G4int Z = 0;
  G4int Ab = 0;
  G4int Zb = 0;
  nucleonVector = nucleons_in;
  SpecAa = nucleonVector->GetA("A");
  SpecAb = nucleonVector->GetA("B");

  auto *nucleons = new aamcc::NucleonVector();    // nucleons from Side A
  auto *nucleons_B = new aamcc::NucleonVector();  // nucleons from Side B

  for (G4int iArray = 0; iArray < nucleons_in->size(); iArray++) {
    aamcc::Nucleon *nucleon = &(nucleons_in->at(iArray));
    if (nucleon->isParticipant == 0 && nucleon->Nucl == "A") {
      A += 1;
      nucleons->push_back(nucleons_in->at(iArray));
    }
    if (nucleon->isParticipant == 0 && nucleon->Nucl == "A" && nucleon->isospin == 1) {
      Z += 1;
    }
    if (nucleon->isParticipant == 0 && nucleon->Nucl == "B") {
      Ab += 1;
      nucleons_B->push_back(nucleons_in->at(iArray));
    }
    if (nucleon->isParticipant == 0 && nucleon->Nucl == "B" && nucleon->isospin == 1) {
      Zb += 1;
    }
  }
  // Creating a graph of nucleons (side A)
  Graph g = this->ClusterToGraph(nucleons, A);

  // Creating a graph of nucleons (side B)
  Graph g_B = this->ClusterToGraph(nucleons_B, Ab);

  // Applying MST + critical distance cut + DFS clustering (side A)
  vector<vector<G4int> > clusters;
  if (ExA > 0 && ExB > 0) this->SetCDExEn(ExA, A);
  // G4cout<<"CritDist = "<<CritDist<<"\n";
  CritDistA = this->CritDist;
  clusters = g.AdvancedKruskalMST(this->CritDist);

  // Applying MST + critical distance cut + DFS clustering (side B)
  vector<vector<G4int> > clusters_B;
  if (ExA > 0 && ExB > 0) this->SetCDExEn(ExB, Ab);
  clusters_B = g_B.AdvancedKruskalMST(this->CritDist);

  std::vector<G4FragmentVector> outClusters;
  G4FragmentVector output_vector_A;
  G4FragmentVector output_vector_B;

  std::vector<int> rmapsA;  // repulsion maps with fragments for size A
  std::vector<int> rmapsB;  // repulsion maps with fragments for size B
  int rcountA = 0;
  int rcountB = 0;
  aamcc::NucleonVector rnucsA;  // protons for size A
  aamcc::NucleonVector rnucsB;  // protons for size B

  // Filling a Clister Vector (side A)
  for (G4int i = 0; i < clusters.size(); ++i) {
    G4int Z_clust = 0;
    G4int A_clust = 0;
    for (G4int j = 0; j < clusters[i].size(); ++j) {
      aamcc::Nucleon *nucleon = &(nucleons->at((clusters[i])[j]));

      if (nucleon->isospin == 1) {
        rnucsA.push_back(*nucleon);
        rmapsA.push_back(rcountA);

        Z_clust += 1;
      }
      A_clust += 1;
    }

    G4Fragment *frag = new G4Fragment();
    frag->SetA(A_clust);
    frag->SetZ(Z_clust);
    frag->SetMomentum({0, 0, 0, 0});
    output_vector_A.push_back(frag);

    ++rcountA;
  }

  // Filling a Clister Vector (side B)
  for (G4int i = 0; i < clusters_B.size(); ++i) {
    G4int Z_clust = 0;
    G4int A_clust = 0;
    for (G4int j = 0; j < clusters_B[i].size(); ++j) {
      aamcc::Nucleon *nucleon = &(nucleons_B->at((clusters_B[i])[j]));

      if (nucleon->isospin == 1) {
        rnucsB.push_back(*nucleon);
        rmapsB.push_back(rcountB);

        Z_clust += 1;
      }
      A_clust += 1;
    }

    G4Fragment *frag = new G4Fragment();
    frag->SetA(A_clust);
    frag->SetZ(Z_clust);
    frag->SetMomentum({0, 0, 0, 0});
    output_vector_B.push_back(frag);

    ++rcountB;
  }
  outClusters.push_back(output_vector_A);
  outClusters.push_back(output_vector_B);

  delete nucleons;
  delete nucleons_B;

  return CalculateMomentum(outClusters, ExA, ExB, boostA, boostB, rnucsA, rnucsB, rmapsA, rmapsB);
};

std::vector<G4FragmentVector> GMSTClustering::CalculateMomentum(std::vector<G4FragmentVector> noMomClusters, G4double ExEnA, G4double ExEnB, CLHEP::Hep3Vector boostA, CLHEP::Hep3Vector boostB, aamcc::NucleonVector rnucsA, aamcc::NucleonVector rnucsB, std::vector<int> rmapsA, std::vector<int> rmapsB) {
  std::vector<G4FragmentVector> momClusters = noMomClusters;
  std::vector<G4double> MstMassVector_A;
  MstMassVector_A.reserve(noMomClusters.at(0).size());

  std::vector<G4double> MstMassVector_B;
  MstMassVector_B.reserve(noMomClusters.at(1).size());

  G4double SumMassMst = 0;
  G4double SumMassMstEx = 0;

  if (!noMomClusters.at(0).empty()) {
    for (G4int i = 0; i < noMomClusters.at(0).size(); ++i) {
      G4int clfrag_A = (noMomClusters.at(0)[i])->GetA();
      G4int clfrag_Z = (noMomClusters.at(0)[i])->GetZ();
      // Excitation energy computing
      G4double energy = 0;
      if (!((clfrag_Z == 0 && (clfrag_A == 1)) || (clfrag_Z == 1 && (clfrag_A == 1)))) {
        energy = ExEnA * G4double(clfrag_A) / G4double(SpecAa);
      }

      G4double NuclearMass = G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z) + energy;

      SumMassMst += G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z);

      SumMassMstEx += NuclearMass;

      MstMassVector_A.push_back(NuclearMass);
    }

    G4double PrefragmentMass_A = SumMassMst + ExEnA;

    std::vector<G4LorentzVector *>* momentumVectorA;
    // if is commented as a part of a long history
    if (PrefragmentMass_A < (SumMassMstEx + 1e-5 * MeV)) {
      PrefragmentMass_A += 1e-5 * MeV;
    }
    momentumVectorA = phaseSpaceDecay.Decay(PrefragmentMass_A, MstMassVector_A);
    //}

    // ----------- Momentums A
    for (int I = 0; I < momClusters.at(0).size(); ++I) {
      momClusters.at(0).at(I)->SetMomentum(*momentumVectorA->at(I));
    }

    momClusters.at(0) = RepulsionStage::CalculateRepulsion(momClusters.at(0), rnucsA, rmapsA);
    for (int I = 0; I < momClusters.at(0).size(); ++I) {
      G4LorentzVector momentum = momClusters.at(0).at(I)->GetMomentum();
      momentum.boost(boostA);
      momClusters.at(0).at(I)->SetMomentum(momentum);
    }
    // ----------

    momentumVectorA->clear();
  }
  // side B
  SumMassMst = 0;
  SumMassMstEx = 0;
  if (!noMomClusters.at(1).empty()) {
    for (G4int i = 0; i < noMomClusters.at(1).size(); ++i) {
      G4int clfrag_A = (noMomClusters.at(1)[i])->GetA();
      G4int clfrag_Z = (noMomClusters.at(1)[i])->GetZ();
      // Excitation energy computing
      G4double energy = 0;
      if (!((clfrag_Z == 0 && (clfrag_A == 1)) || (clfrag_Z == 1 && (clfrag_A == 1)))) {
        energy = ExEnB * G4double(clfrag_A) / G4double(SpecAb);
      }

      G4double NuclearMass = G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z) + energy;

      SumMassMst += G4NucleiProperties::GetNuclearMass(clfrag_A, clfrag_Z);

      SumMassMstEx += NuclearMass;

      MstMassVector_B.push_back(NuclearMass);
    }

    G4double PrefragmentMass_B = SumMassMst + ExEnB;
    std::vector<G4LorentzVector *> *momentumVectorB;
    // if is commented as a part of a history
    if (PrefragmentMass_B < (SumMassMstEx + 1e-5 * MeV)) {
      PrefragmentMass_B += 1e-5 * MeV;
    }
    momentumVectorB = phaseSpaceDecay.Decay(PrefragmentMass_B, MstMassVector_B);
    //}

    // ----------- Momentums B
    for (int I = 0; I < momClusters.at(1).size(); ++I) {
      momClusters.at(1).at(I)->SetMomentum(*momentumVectorB->at(I));
    }

    momClusters.at(1) = RepulsionStage::CalculateRepulsion(momClusters.at(1), rnucsB, rmapsB);
    for (int I = 0; I < momClusters.at(1).size(); ++I) {
      G4LorentzVector momentum = momClusters.at(1).at(I)->GetMomentum();
      momentum.boost(boostB);
      momClusters.at(1).at(I)->SetMomentum(momentum);
    }
    // ----------

    momentumVectorB->clear();
  }

  MstMassVector_A.clear();
  MstMassVector_B.clear();
  noMomClusters.clear();

  return momClusters;
};

Graph::Graph(G4int V, G4int E) {
  this->V = V;
  this->E = E;
  adj = new list<G4int>[V];
}

Graph::Graph() {
  this->V = 0;
  this->E = 0;
}

Graph::~Graph() { delete[] adj; }

void Graph::addEdge(G4int u, G4int v, G4double w) { edges.push_back({w, {u, v}}); }

vector<vector<G4int> > Graph::AdvancedKruskalMST(G4double CD_in) {
  // Sort edges in increasing order on basis of cost
  sort(edges.begin(), edges.end());

  // Create disjoint sets
  DisjointSets ds(V);

  // Iterate through all sorted edges
  vector<pair<G4double, iPair> >::iterator it;
  for (it = edges.begin(); it != edges.end(); it++) {
    G4int u = it->second.first;
    G4int v = it->second.second;

    G4int set_u = ds.find(u);
    G4int set_v = ds.find(v);

    // Check if the selected edge is creating
    // a cycle or not (Cycle is created if u
    // and v belong to same set)
    if (set_u != set_v) {
      // Current edge will be in the MST
      if (it->first < CD_in) {
        this->addConn(u, v);
      }
      // Merge two sets
      ds.merge(set_u, set_v);
    }
  }

  return this->connectedComponents();
}

void Graph::addConn(G4int v, G4int w) {
  adj[v].push_back(w);
  adj[w].push_back(v);
}

vector<vector<G4int> > Graph::connectedComponents() {
  vector<vector<G4int> > clusters;
  // Mark all the vertices as not visited
  G4bool *visited = new G4bool[V];
  for (G4int v = 0; v < V; v++) visited[v] = false;
  for (G4int v = 0; v < V; v++) {
    if (visited[v] == false) {
      // vector to fill with nucleon's numbers from one cluster
      vector<G4int> clust_inside;
      // print all reachable vertices from v
      DFSUtil(v, visited, &clust_inside);
      // push cluster to vector of clusters
      clusters.push_back(clust_inside);
    }
  }
  delete[] visited;
  return clusters;
}

void Graph::DFSUtil(G4int v, G4bool visited[], vector<G4int> *clust_inside) {
  // Mark the current node as visited and print it
  visited[v] = true;
  (*clust_inside).push_back(v);

  // Recur for all the vertices adjacent to this vertex
  list<G4int>::iterator i;
  for (i = adj[v].begin(); i != adj[v].end(); ++i)
    if (!visited[*i]) DFSUtil(*i, visited, &(*clust_inside));
}

DisjointSets::DisjointSets(G4int n) {
  // Allocate memory
  this->n = n;
  parent = new G4int[n + 1];
  rnk = new G4int[n + 1];

  // Initially, all vertices are in different sets and have rank 0.
  for (G4int i = 0; i <= n; i++) {
    rnk[i] = 0;
    // every element is parent of itself
    parent[i] = i;
  }
}

DisjointSets::~DisjointSets() {
  delete[] parent;
  delete[] rnk;
}

G4int DisjointSets::find(G4int u) {
  /* Make the parent of the nodes in the path
     from u--> parent[u] point to parent[u] */
  if (u != parent[u]) parent[u] = find(parent[u]);
  return parent[u];
}

void DisjointSets::merge(G4int x, G4int y) {
  x = find(x), y = find(y);

  /* Make tree with smaller height
     a subtree of the other tree  */
  if (rnk[x] > rnk[y])
    parent[y] = x;
  else  // If rnk[x] <= rnk[y]
    parent[x] = y;

  if (rnk[x] == rnk[y]) rnk[y]++;
}