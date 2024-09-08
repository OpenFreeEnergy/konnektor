"""Clustering compounds based on scaffolds

This clusterer attempts to cluster compounds based on their scaffolds.
It is built on rdkit's rdScaffoldNetwork module.
"""

from collections import defaultdict
import itertools

from rdkit import Chem
from rdkit.Chem import rdMolHash
from rdkit.Chem.Scaffolds import rdScaffoldNetwork

from gufe import Component

from ._abstract_clusterer import _AbstractClusterer


class ScaffoldClusterer(_AbstractClusterer):
    scaffold_looseness: int
    cid_scaffold: dict[int, str]
    cid_components: dict[int, list[Component]]

    def __init__(self, scaffold_looseness: int = 9):
        """
        Parameters
        ----------
        scaffold_looseness : int
          a heuristic to define to what extent alternate/smaller scaffolds can
          be used to match a certain molecule.
          this value decides how many heavy atoms from the *largest* scaffold
          other scaffolds may be permitted.
          a too high value may result in inappropriately generic scaffolds being
          used, while too low will result in too many scaffolds being identified
        """
        self.scaffold_looseness = scaffold_looseness

    @staticmethod
    def normalise_molecules(mols: list[Component]) -> dict[Component, Chem.Mol]:
        # Convert SMC to a normalised (reduced & cleaned up) version in rdkit
        # This is anonymous, so no bond orders, charges or elements
        # This makes comparing scaffolds more suited to RBFEs
        def normalised_rep(mol):
            smi = rdMolHash.MolHash(
                Chem.RemoveHs(mol.to_rdkit()), rdMolHash.HashFunction.AnonymousGraph
            )
            return Chem.MolFromSmiles(smi)

        # returns mapping of molecules to their normalised rep
        mols2anonymous = {mol: normalised_rep(mol) for mol in mols}

        return mols2anonymous

    @staticmethod
    def generate_scaffold_network(
        mols: list[Chem.Mol],
    ) -> rdScaffoldNetwork.ScaffoldNetwork:
        # generates the scaffold network from the rdkit mol objects
        params = rdScaffoldNetwork.ScaffoldNetworkParams()
        params.includeScaffoldsWithAttachments = False
        params.flattenChirality = True
        params.pruneBeforeFragmenting = True
        net = rdScaffoldNetwork.CreateScaffoldNetwork(mols, params)

        return net

    @staticmethod
    def match_scaffolds_to_source(
        network: rdScaffoldNetwork.ScaffoldNetwork,
        mols: list[Chem.Mol],
        hac_heuristic: int,
    ) -> dict[Chem.Mol, list[str]]:
        # match scaffolds in network back to normalised input molecules
        # i.e. for each molecule, which scaffolds can apply

        # will store for each molecule, potential scaffolds and their size
        mols2scaffolds = defaultdict(list)
        for scaff in network.nodes:
            # determine size of scaffold
            q = Chem.MolFromSmarts(scaff)
            natoms = q.GetNumAtoms()

            for m in mols:
                if m.HasSubstructMatch(q):
                    mols2scaffolds[m].append((scaff, natoms))

        # then filter these scaffolds to only allow those which are large enough
        mols2candidates = defaultdict(list)
        for m, scaffs in mols2scaffolds.items():
            # determine what the largest scaffold for this molecule was
            largest_scaff = max(scaffs, key=lambda x: x[1])
            # for each molecule, a size ordered list of scaffolds that they could be assigned to
            mols2candidates[m] = sorted(
                [s for s in scaffs if (largest_scaff[1] - s[1]) <= hac_heuristic],
                key=lambda x: x[1],
                reverse=True,
            )

        return mols2candidates

    @staticmethod
    def find_solution(
        mol_to_candidates: dict[Chem.Mol, list[str]]
    ) -> list[tuple[str, int]]:
        # returns the best scaffolds that cover all mols
        # returns a list of (scaffold smiles, n heavy atoms)

        # reverse mapping of scaffolds onto the mols they cater for
        scaffold2mols = defaultdict(list)
        anon_mols = set()
        for mol, scaffs in mol_to_candidates.items():
            anon_mols.add(mol)
            for scaff in scaffs:
                scaffold2mols[scaff].append(mol)

        def scaffold_coverage(scaffolds, scaff2mol, all_mols) -> bool:
            """Does this combination of scaffolds cover all ligands"""
            covered_mols = set()

            for scaff in scaffolds:
                covered_mols |= set(scaff2mol[scaff])

            return covered_mols == set(all_mols)

        candidate_scaffolds = set(
            itertools.chain.from_iterable(mol_to_candidates.values())
        )
        # try one scaffold to see if it catches all molecules
        # then try all combinations of two scaffolds to see if we cover
        # etc until we find a solution
        # then pick the solution with the largest scaffolds
        for i in range(1, len(candidate_scaffolds)):
            solutions = []

            for scaffolds in itertools.combinations(candidate_scaffolds, i):
                if not scaffold_coverage(scaffolds, scaffold2mols, anon_mols):
                    continue

                solutions.append(scaffolds)

            if solutions:
                # pick the best, based on HAC
                solution = max(solutions, key=lambda x: sum(s[1] for s in x))
                return solution

    @staticmethod
    def formulate_answer(
        solution: list[tuple[str, int]], mols_to_norm: dict[Component, Chem.Mol]
    ) -> dict[str, list[Component]]:
        # relate the solution scaffolds back to the input SMC

        # for each molecule, pick the largest scaffold that matches
        relationship = []
        for input_mol, anon_mol in mols_to_norm.items():
            best = -1
            best_scaff = None
            for scaff, natoms in solution:
                if natoms < best:
                    continue
                if anon_mol.HasSubstructMatch(Chem.MolFromSmarts(scaff)):
                    best_scaff = scaff
            relationship.append((input_mol, best_scaff))

        final_answer = defaultdict(list)
        for input_mol, scaff in relationship:
            final_answer[scaff].append(input_mol)

        return dict(final_answer)

    def cluster_compounds(self, components: list[Component]):
        # first normalise the molecules to a generic representation
        mols_to_norm = self.normalise_molecules(components)

        # then create a scaffold network from the normalised molecules
        network = self.generate_scaffold_network(list(mols_to_norm.values()))

        # reassign the scaffolds in the network back to the normalised reps
        mol_to_candidates = self.match_scaffolds_to_source(
            network,
            list(mols_to_norm.values()),
            self.scaffold_looseness,
        )

        # then try increasingly larger number of scaffolds
        # until we hit full coverage of molecules
        solution = self.find_solution(mol_to_candidates)

        # finally, relate this solution back to the input set and store res.
        scaffold_components = self.formulate_answer(solution, mols_to_norm)
        self.cid_scaffold = {}
        self.cid_components = {}

        for i, (scaff, components) in enumerate(scaffold_components.items()):
            self.cid_scaffold[i] = scaff
            self.cid_components[i] = components

        return self.cid_components
