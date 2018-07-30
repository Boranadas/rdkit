#ifndef RD_ATOMPAIRGEN_H_2018_06
#define RD_ATOMPAIRGEN_H_2018_06

#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <GraphMol/Fingerprints/FingerprintUtil.h>

namespace RDKit {
namespace AtomPair {
using namespace AtomPairs;

class AtomPairAtomInvGenerator : public AtomInvariantsGenerator {
  const bool df_includeChirality;

 public:
  AtomPairAtomInvGenerator(bool includeChirality = false);

  std::vector<std::uint32_t> *getAtomInvariants(const ROMol &mol) const;

  std::string infoString() const;
};

/*!
  /brief class that holds atom-pair fingerprint specific arguments

 */
template <typename OutputType>
class AtomPairArguments : public FingerprintArguments<OutputType> {
 public:
  const bool df_includeChirality;
  const bool df_use2D;
  const unsigned int d_minDistance;
  const unsigned int d_maxDistance;

  OutputType getResultSize() const;

  std::string infoString() const;

  /*!
    /brief construct a new AtomPairArguments object

    /param countSimulation  if set, use count simulation while generating the
    fingerprint
    /param includeChirality if set, chirality will be used in the atom
    invariants, this is ignored if atomInvariantsGenerator is present for
    the /c FingerprintGenerator that uses this
    /param use2D            if set, the 2D (topological) distance matrix will be
    used
    /param minDistance      minimum distance between atoms to be considered in a
    pair, default is 1 bond
    /param maxDistance      maximum distance between atoms to be considered in a
    pair, default is maxPathLen-1 bonds
    /param countBounds      boundries for count simulation, corresponding bit
    will be set if the count is higher than the number provided for that spot
    /param foldedSize       size of the folded version of the fingerprint

   */
  AtomPairArguments(const bool countSimulation = true,
                    const bool includeChirality = false,
                    const bool use2D = true, const unsigned int minDistance = 1,
                    const unsigned int maxDistance = (maxPathLen - 1),
                    const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
                    const std::uint32_t foldedSize = 2048);
};

/*!
  /brief class that holds atom-environment data needed for atom-pair fingerprint
  generation

 */
template <typename OutputType>
class AtomPairAtomEnv : public AtomEnvironment<OutputType> {
  const unsigned int d_atomIdFirst;
  const unsigned int d_atomIdSecond;
  const unsigned int d_distance;

 public:
  std::uint32_t getBitId(FingerprintArguments<OutputType> *arguments,
                         const std::vector<std::uint32_t> *atomInvariants,
                         const std::vector<std::uint32_t> *bondInvariants,
                         const AdditionalOutput *additionalOutput) const;

  /*!
    /brief construct a new AtomPairAtomEnv object

    /param atomIdFirst      id of the first atom of the atom-pair
    /param atomIdSecond     id of the second atom of the atom-pair
    /param distance         distance between the atoms
   */
  AtomPairAtomEnv(const unsigned int atomIdFirst,
                  const unsigned int atomIdSecond, const unsigned int distance);
};

/*!
  /brief class that generates atom-environments for atom-pair fingerprint

 */
template <typename OutputType>
class AtomPairEnvGenerator : public AtomEnvironmentGenerator<OutputType> {
 public:
  std::vector<AtomEnvironment<OutputType> *> getEnvironments(
      const ROMol &mol, FingerprintArguments<OutputType> *arguments,
      const std::vector<std::uint32_t> *fromAtoms,
      const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
      const AdditionalOutput *additionalOutput,
      const std::vector<std::uint32_t> *atomInvariants,
      const std::vector<std::uint32_t> *bondInvariants) const;

  std::string infoString() const;
};

/*!
  /brief helper function that generates a /c FingerprintGenerator that generates
  atom-pair fingerprints


  /param minDistance            minimum distance between atoms to be considered
  in a pair, default is 1 bond
  /param maxDistance            maximum distance between atoms to be considered
  in a pair, default is maxPathLen-1 bonds
  /param includeChirality       if set, chirality will be used in the atom
  invariants, this is ignored if atomInvariantsGenerator is provided
  /param use2D                  if set, the 2D (topological) distance matrix
  will be used
  /param useCountSimulation         if set, use count simulation while
  generating the fingerprint
  /param atomInvariantsGenerator    atom invariants to be used during
  fingerprint generation
  /param bondInvariantsGenerator    bond invariants to be used during
  fingerprint generation

  /return FingerprintGenerator that generates atom-pair fingerprints
 */
FingerprintGenerator<std::uint32_t> *getAtomPairGenerator(
    const unsigned int minDistance = 1,
    const unsigned int maxDistance = maxPathLen - 1,
    const bool includeChirality = false, const bool use2D = true,
    const bool useCountSimulation = true,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    BondInvariantsGenerator *bondInvariantsGenerator = nullptr);
}  // namespace AtomPair
}  // namespace RDKit

#endif
