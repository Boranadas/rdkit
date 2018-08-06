#ifndef RD_MORGANGEN_H_2018_07
#define RD_MORGANGEN_H_2018_07

#include <GraphMol/Fingerprints/FingerprintGenerator.h>
#include <cstdint>

namespace RDKit {

namespace MorganFingerprint {

/**
 /brief Default atom invariants generator for Morgan fingerprint, generates
 ECFP-type invariants

 */
class MorganAtomInvGenerator : public AtomInvariantsGenerator {
  const bool df_includeRingMembership;

 public:
  /**
   /brief Construct a new MorganAtomInvGenerator object

   /param includeRingMembership : if set, whether or not the atom is in a ring
   will be used in the invariant list.
   */
  MorganAtomInvGenerator(const bool includeRingMembership = true);

  std::vector<std::uint32_t> *getAtomInvariants(const ROMol &mol) const;

  std::string infoString() const;
  MorganAtomInvGenerator *clone() const;
};

/**
 /brief Alternative atom invariants generator for Morgan fingerprint, generate
 FCFP-type invariants

 */
class MorganFeatureAtomInvGenerator : public AtomInvariantsGenerator {
  std::vector<const ROMol *> *dp_patterns;

 public:
  /**
   /brief Construct a new MorganFeatureAtomInvGenerator object

   /param patterns : if provided should contain the queries used to assign
   atom-types. if not provided, feature definitions adapted from reference:
   Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998) will be used for
   Donor, Acceptor, Aromatic, Halogen, Basic, Acidic.
   */
  MorganFeatureAtomInvGenerator(std::vector<const ROMol *> *patterns = nullptr);

  std::vector<std::uint32_t> *getAtomInvariants(const ROMol &mol) const;

  std::string infoString() const;
  MorganFeatureAtomInvGenerator *clone() const;
};

/**
 /brief Bond invariants generator for Morgan fingerprint

 */
class MorganBondInvGenerator : public BondInvariantsGenerator {
  const bool df_useBondTypes;
  const bool df_useChirality;

 public:
  /**
   /brief Construct a new MorganBondInvGenerator object

   /param useBondTypes : if set, bond types will be included as a part of the
   bond invariants
   /param useChirality : if set, chirality information will be included as a
   part of the bond invariants
   */
  MorganBondInvGenerator(const bool useBondTypes = true,
                         const bool useChirality = false);

  std::vector<std::uint32_t> *getBondInvariants(const ROMol &mol) const;

  std::string infoString() const;
  MorganBondInvGenerator *clone() const;
};

/**
 /brief Class for holding Morgan fingerprint specific arguments

 */
template <typename OutputType>
class MorganArguments : public FingerprintArguments<OutputType> {
 public:
  const bool df_includeChirality;
  const bool df_onlyNonzeroInvariants;
  const unsigned int d_radius;

  OutputType getResultSize() const;

  std::string infoString() const;

  /**
   /brief Construct a new MorganArguments object

   /param radius : the number of iterations to grow the fingerprint
   /param countSimulation : if set, use count simulation while generating the
    fingerprint
   /param includeChirality : if set, chirality information will be added to the
   generated bit id, independently from bond invariants
   /param onlyNonzeroInvariants : if set, bits will only be set from atoms that
   have a nonzero invariant
   /param countBounds : boundaries for count simulation, corresponding bit will
   be set if the count is higher than the number provided for that spot
   /param foldedSize : size of the folded version of the fingerprints
   */
  MorganArguments(const unsigned int radius, const bool countSimulation = true,
                  const bool includeChirality = false,
                  const bool onlyNonzeroInvariants = false,
                  const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
                  const std::uint32_t foldedSize = 2048);
};

/**
 /brief Class for holding the bit-id created from Morgan fingerprint
 environments and the additional data necessary extra outputs

 */
template <typename OutputType>
class MorganAtomEnv : public AtomEnvironment<OutputType> {
  const OutputType d_code;
  const unsigned int d_atomId;
  const unsigned int d_layer;

 public:
  OutputType getBitId(FingerprintArguments<OutputType> *arguments,
                      const std::vector<std::uint32_t> *atomInvariants,
                      const std::vector<std::uint32_t> *bondInvariants,
                      const AdditionalOutput *additionalOutput, const bool hashResults = false) const;

  /**
   /brief Construct a new MorganAtomEnv object

   /param code : bit id generated from this environment
   /param atomId : atom id of the atom at the center of this environment
   /param layer : radius of this environment
   */
  MorganAtomEnv(const std::uint32_t code, const unsigned int atomId,
                const unsigned int layer);
};

/**
 /brief Class that generates atom environments for Morgan fingerprint

 */
template <typename OutputType>
class MorganEnvGenerator : public AtomEnvironmentGenerator<OutputType> {
 public:
  std::vector<AtomEnvironment<OutputType> *> getEnvironments(
      const ROMol &mol, FingerprintArguments<OutputType> *arguments,
      const std::vector<std::uint32_t> *fromAtoms,
      const std::vector<std::uint32_t> *ignoreAtoms, const int confId,
      const AdditionalOutput *additionalOutput,
      const std::vector<std::uint32_t> *atomInvariants,
      const std::vector<std::uint32_t> *bondInvariants, const bool hashResults = false) const;

  std::string infoString() const;
};

/**
 /brief Get a fingerprint generator for Morgan fingerprint

 /tparam OutputType determines the size of the bitIds and the result, can be 32
 or 64 bit unsigned integer
 /param radius : the number of iterations to grow the fingerprint
 /param countSimulation : if set, use count simulation while generating the
 fingerprint
 /param includeChirality : if set, chirality information will be added to the
 generated bit id, independently from bond invariants
 /param onlyNonzeroInvariants : if set, bits will only be set from atoms that
 have a nonzero invariant
 /param countBounds : boundaries for count simulation, corresponding bit willbe
 set if the count is higher than the number provided for that spot
 /param foldedSize : size of the folded version of the fingerprints
 /param countSimulation : countSimulation for MorganArguments
 /param includeChirality : sets includeChirality flag for both MorganArguments
 and the default bond generator MorganBondInvGenerator
 /param useBondTypes : useBondTypes for MorganBondInvGenerator
 /param onlyNonzeroInvariants : onlyNonzeroInvariants for MorganArguments
 /param atomInvariantsGenerator : custom atom invariants generator to use
 /param bondInvariantsGenerator : custom bond invariants generator to use
 /param ownsAtomInvGen  if set atom invariants generator is destroyed with the
 fingerprint generator
 /param ownsBondInvGen  if set bond invariants generator is destroyed with the
 fingerprint generator

 /return FingerprintGenerator<OutputType>* that generates Morgan fingerprints
 */
template <typename OutputType>
FingerprintGenerator<OutputType> *getMorganGenerator(
    const unsigned int radius, const bool countSimulation = true,
    const bool includeChirality = false, const bool useBondTypes = true,
    const bool onlyNonzeroInvariants = false,
    AtomInvariantsGenerator *atomInvariantsGenerator = nullptr,
    BondInvariantsGenerator *bondInvariantsGenerator = nullptr,
    const std::uint32_t foldedSize = 2048,
    const std::vector<std::uint32_t> countBounds = {1, 2, 4, 8},
    const bool ownsAtomInvGen = false, const bool ownsBondInvGen = false);

}  // namespace MorganFingerprint
}  // namespace RDKit

#endif